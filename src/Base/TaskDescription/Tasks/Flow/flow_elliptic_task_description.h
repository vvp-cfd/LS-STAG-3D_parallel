/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_elliptic_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 18:18:32
* \brief     flow_elliptic_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ELLIPTIC_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_ELLIPTIC_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowEllipticTaskDescription

Данный файл содержит в себе реализацию класса FlowEllipticTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetElliptic, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, rad,
		const double, ANGLE, const double, coef, const double, coef2, const double, coef3, const double, coef4);
/*!
\brief Класс, реализующий работу с задачей об обтекании эллиптического профиля

Данный класс реализует работу с описанием решаемой задачи об обтекании эллиптического профиля
*/
	struct FlowEllipticTaskDescription : FlowTaskDescription
	{
	  public:
		FlowEllipticTaskDescription():
			FlowTaskDescription(AIRFOIL_ELLIPTIC, MOTION_NO, (LSSTAGUniformMesh*)&mesh, sizeof(FlowEllipticTaskDescription)),
			rad(0.5),
			ThRatio(DEFAULT_ELLIPTIC_TH_RATIO)//default
		{
			discontin = 1;
		}
		LS_STAG_KERNEL_2D(calcLevelSetElliptic, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double rad, const double ANGLE, const double coef, const double coef2,
			const double coef3, const double coef4)
		{
			int curNum = _NM * j + _N * i + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = curNum + hf; curNum < end; curNum++)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				_phi[curNum] = rad - ((fabs(ANGLE) < EQ_DOUBLE) ? sqrt(pToC.x * pToC.x + coef * pToC.y * pToC.y)
					: sqrt(coef2 * pToC.x * pToC.x + coef3 * pToC.x * pToC.y + coef4 * pToC.y * pToC.y));
				curD++;
				curPoint.x += curD->x;
			};
		}
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> ThRatio;
			in >> ANGLE;
		};
	  private:
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << ThRatio << "\n";
			out << ANGLE << "\n";
		};
		/*!
		Вычисляет значения функции уровня в углах ячеек основной сетки
		\param zs номер ряда по z, с которого рассматриваем ячейки
		\param zf номер ряда по z, до которого рассматриваем ячейки
		\param vs номер ряда по вертикали, с которого рассматриваем ячейки
		\param vf номер ряда по вертикали, до которого рассматриваем ячейки
		\param hs номер ряда по горизонтали, с которого рассматриваем ячейки
		\param hf номер ряда по горизонтали, до которого рассматриваем ячейки
		*/
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf)
		{
			LS_STAG_PROFILE("buildLS");
			LS_STAG_PAR_FOR_2D_RANGE(FlowEllipticTaskDescription, calcLevelSetElliptic, vs, vf, zs, zf, hs, hf, mesh.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad, ANGLE, coef, coef2, coef3, coef4)
		};
		void initMeshParams()
		{
			L_BASE = 1.0;// хорда профиля не зависит от угла атаки (cos(ANGLE)=0 при ANGLE=90 deg)
			fcoef /= (L_BASE * ThRatio);
			FlowTaskDescription::initMeshParams();
			coef = 1.0 / (ThRatio * ThRatio);
			double sinPhi = sin(ANGLE);
			double cosPhi = cos(ANGLE);
			coef2 = cosPhi * cosPhi + coef * sinPhi * sinPhi;
			coef3 = 2.0 * (coef - 1.0) * sinPhi * cosPhi;
			coef4 = coef * cosPhi * cosPhi + sinPhi * sinPhi;
		};
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
		double rad, coef, coef2, coef3, coef4;
	  public:
		double ThRatio; ///< для эллиптического цилиндра: отношение малой полуоси к большой
	}; 
}

#endif //_LS_STAG_FLOW_ELLIPTIC_TASK_DESCRIPTION_H_
