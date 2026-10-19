/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_rotor_darrieus_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 20:38:15
* \brief     flow_rotor_darrieus_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ROTOR_DARRIEUS_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_ROTOR_DARRIEUS_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_rotor_task_description.h"
#include "../../../../../Bezier/airfoil_bezier_approximation.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowRotorDarrieusTaskDescription

Данный файл содержит в себе реализацию класса FlowRotorDarrieusTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetRotorDarrieus, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, BaseRadius,
		const int, numBlade, const double, radius, const double, step2, const double, step, const double, ANGLE,
		const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor, const double, len, const double, inv_len, const double, delta,
		const double, _LE, const int, cU_L2, LS_STAG_CR_PD, cU_coef, const int, cdU_L2, LS_STAG_CR_PD, cdU_coef,
		const int, cD_L2, LS_STAG_CR_PD, cD_coef, const int, cdD_L2, LS_STAG_CR_PD, cdD_coef, const double, _eps);
/*!
\brief Класс, реализующий работу с задачей об обтекании ротора Дарье

Данный класс реализует работу с описанием решаемой задачи об обтекании ротора Дарье. Ротор Дарье может 
авторотировать или совершать автоколебания
*/
	struct FlowRotorDarrieusTaskDescription : FlowRotorTaskDescription
	{
	  public:
		FlowRotorDarrieusTaskDescription():
		FlowRotorTaskDescription(ROTOR_DARRIEUS, sizeof(FlowRotorDarrieusTaskDescription))
		{
			radius = 0.5 * L_BASE;
		}
		LS_STAG_KERNEL_2D(calcLevelSetRotorDarrieus, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double BaseRadius, const int numBlade, const double radius,
			const double step2, const double step, const double ANGLE, const double level, const RotorKinds id_rotor,
			const double len, const double inv_len, const double delta, const double _LE, const int cU_L2,
			LS_STAG_CR_PD cU_coef, const int cdU_L2, LS_STAG_CR_PD cdU_coef, const int cD_L2,
			LS_STAG_CR_PD cD_coef, const int cdD_L2, LS_STAG_CR_PD cdD_coef, const double _eps)
		{
			int curNum = _NM * j + _N * i + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);// по z просто параллельный перенос
			for (int end = curNum + hf; curNum < end; curNum++, curD++, curPoint.x += curD->x)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				const double r = pToC.norm();
				if (BaseRadius < r)
				{
					double phi = 0.0;
					const double phiBlade = FlowRotorTaskDescription::CurBlade(pToC, phi, numBlade, step2, step, ANGLE, level, id_rotor);
					pToC.diff(cos(phiBlade), sin(phiBlade), 0.0, radius);
					const double attackBlade = -phiBlade + PI_2;// 11.06.13 ---  дополняющий к углу атаки //attackBlade = phiBlade + pi2;
					const double cosa = cos(attackBlade);
					const double sina = sin(attackBlade);
					const double x_rot = P_C.x + pToC.x * cosa - pToC.y * sina;
					const double y_rot = P_C.y + pToC.x * sina + pToC.y * cosa;
_phi[curNum] = AirfoilBezierApproximation::getValue(x_rot, y_rot, len,
					inv_len, delta, _LE, cU_L2, cU_coef, cdU_L2, cdU_coef, cD_L2, cD_coef, cdD_L2, cdD_coef, _eps);
				}
				else _phi[curNum] = BaseRadius - r;
			};
		}
		void load(std::ifstream &in)
		{
			FlowRotorTaskDescription::load(in);
			airfoil.load(in);
		};
	  private:
		void initMeshParams()
		{
			BaseRadius = 0.1 * L_BASE;// радиус опорной башни ВЭУ
			FlowRotorTaskDescription::initMeshParams();
			airfoil.set(this);// + здесь нужно считать координаты
		};
		void printExtParams(std::ofstream &out)
		{
			FlowRotorTaskDescription::printExtParams(out);
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
			LS_STAG_PAR_FOR_2D_RANGE(FlowRotorDarrieusTaskDescription, calcLevelSetRotorDarrieus, vs, vf, zs, zf, hs, hf, meshMov.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, BaseRadius, numBlade, radius, step2, step, ANGLE, level,
				id_rotor, airfoil.param.length, airfoil.param.inv_length, airfoil.param.delta, airfoil.param.LE, 
				airfoil.curveU.curve.L2, airfoil.curveU.curve.coef, airfoil.curveU.curve_derivative.L2, 
				airfoil.curveU.curve_derivative.coef, airfoil.curveD.curve.L2, airfoil.curveD.curve.coef, 
				airfoil.curveD.curve_derivative.L2, airfoil.curveD.curve_derivative.coef, airfoil.param.eps)
		};
		AirfoilBezierApproximation airfoil;///< аппроксимация профиля лопасти ротора кривой Безье 
		double radius;
		double BaseRadius;///< радиус опорной башни ВЭУ
	};
}

#endif //_LS_STAG_FLOW_ROTOR_DARRIEUS_TASK_DESCRIPTION_H_
