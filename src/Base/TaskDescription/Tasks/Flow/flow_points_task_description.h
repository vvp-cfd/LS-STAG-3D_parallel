/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_points_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 18:23:26
* \brief     flow_points_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_POINTS_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_POINTS_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_task_description.h"
#include "../../../../Bezier/airfoil_bezier_approximation.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowPointsTaskDescription

Данный файл содержит в себе реализацию класса FlowPointsTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetPoints, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, ANGLE,
		const double, len, const double, inv_len, const double, delta, const double, _LE, const int, cU_L2,
		LS_STAG_CR_PD, cU_coef, const int, cdU_L2, LS_STAG_CR_PD, cdU_coef, const int, cD_L2, LS_STAG_CR_PD, cD_coef,
		const int, cdD_L2, LS_STAG_CR_PD, cdD_coef, const double, _eps);
/*!
\brief Класс, реализующий работу с задачей об обтекании профиля, заданного по точкам

Данный класс реализует работу с описанием решаемой задачи об обтекании профиля, который задан массивом координат точек 
границы. Для аппрокимации погруженной границы используются кривые Безье
*/
	struct FlowPointsTaskDescription : FlowTaskDescription
	{
	  public:
		FlowPointsTaskDescription():
		FlowTaskDescription(AIRFOIL_POINTS, MOTION_NO, (LSSTAGUniformMesh*)&mesh, sizeof(FlowPointsTaskDescription))
		{
			discontin = 1;
		}
		LS_STAG_KERNEL_2D(calcLevelSetPoints, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double ANGLE, const double len, const double inv_len,
			const double delta, const double _LE, const int cU_L2, LS_STAG_CR_PD cU_coef, const int cdU_L2,
			LS_STAG_CR_PD cdU_coef, const int cD_L2, LS_STAG_CR_PD cD_coef, const int cdD_L2, LS_STAG_CR_PD cdD_coef,
			const double _eps)
		{
			int curNum = _NM * j + _N * i + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = curNum + hf; curNum < end; curNum++, curD++, curPoint.x += curD->x)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				const double cosa = cos(ANGLE);
				const double sina = sin(ANGLE);
				const double x_rot = P_C.x + pToC.x * cosa - pToC.y * sina;
				const double y_rot = P_C.y + pToC.x * sina + pToC.y * cosa;
_phi[curNum] = AirfoilBezierApproximation::getValue(x_rot, y_rot, len,
				inv_len, delta, _LE, cU_L2, cU_coef, cdU_L2, cdU_coef, cD_L2, cD_coef, cdD_L2, cdD_coef, _eps);
			};
		}
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> ANGLE;
			airfoil.load(in);
		};
	  private:
		void initMeshParams()
		{
			L_BASE = 1.0;// хорда профиля не зависит от угла атаки (cos(ANGLE)=0 при ANGLE=90 deg ломает построение тела)
			FlowTaskDescription::initMeshParams();
			airfoil.set(this);
		};
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << ANGLE << "\n";
			airfoil.print(out);
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
			LS_STAG_PAR_FOR_2D_RANGE(FlowPointsTaskDescription, calcLevelSetPoints, vs, vf, zs, zf, hs, hf, mesh.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, ANGLE, airfoil.param.length, airfoil.param.inv_length,
				airfoil.param.delta, airfoil.param.LE, airfoil.curveU.curve.L2, airfoil.curveU.curve.coef, 
				airfoil.curveU.curve_derivative.L2, airfoil.curveU.curve_derivative.coef, airfoil.curveD.curve.L2, 
				airfoil.curveD.curve.coef, airfoil.curveD.curve_derivative.L2, airfoil.curveD.curve_derivative.coef, 
				airfoil.param.eps)
			{
			}
		};
		AirfoilBezierApproximation airfoil;///< аппроксимация границы профиля кривыми Безье
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
	};
}

#endif //_LS_STAG_FLOW_POINTS_TASK_DESCRIPTION_H_
