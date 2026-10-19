/*
* \copyright (C) 2013 Puzikova V.V.
* \file      airfoil_bezier_approximation.h
* \author    Puzikova Valeria (vvp) 
* \date      21.10.2013 10:27:05
* \brief     airfoil_bezier_approximation.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_AIRFOIL_BEZIER_APPROXIMATION_H_
#define _LS_STAG_AIRFOIL_BEZIER_APPROXIMATION_H_

// ================================================================== INCLUDES
#include "bezier_approximation.h"
#include "../Base/TaskDescription/Tasks/Flow/flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса AirfoilBezierApproximation

Данный файл содержит в себе реализацию класса AirfoilBezierApproximation
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий аппроксимацию границы профиля кривыми Безье

Данный класс реализует построение аппроксимации границы профиля кривыми Безье и численный расчет функции уровня.

В 3D достаточно вычислить это расстояние для первого слоя, выше просто используем параллельный перенос.
*/
	struct AirfoilBezierApproximation
	{// для всех поворачиваем пр-во, а не профиль, т.е. здесь работаем всегда с 0-м положением
		typedef FlowTaskDescription::AirfoilKind kind;
		AirfoilBezierApproximation()
		{
			;
		}
		/*!
		\return Значение функции уровня в точке с координатами (x,y)
		*/
		LS_STAG_KERNEL_INLINE double getValue(const double x, const double y,
			const double len, const double inv_len, const double delta, const double _LE, const int cU_L2,
			LS_STAG_CR_PD cU_coef, const int cdU_L2, LS_STAG_CR_PD cdU_coef, const int cD_L2,
			LS_STAG_CR_PD cD_coef, const int cdD_L2, LS_STAG_CR_PD cdD_coef, const double _eps)
		{
			double t = (x - _LE) * inv_len;
			if ((t < -DEFAULT_THETA_EPS) || (t > DEFAULT_THETA_ONE_EPS))  return -1.0;
			if (t < 0.0) t = 0.0;
			if (t > 1.0) t = 1.0;
			double resU = Polynomial::getValue(1.0 - t, cU_L2, cU_coef) - y;
			double resD = y - Polynomial::getValue(t, cD_L2, cD_coef);
			const bool up = fabs(resU) < fabs(resD);
			const double res = up ? resU : resD;
			if (fabs(res) > delta) return res;
			double sgn = (res < 0.0) ? -1.0 : 1.0;
			return sgn * ((up) ? PolynomialDistanceFunction::getValue(x, y, up, len, inv_len, delta, _LE, cU_L2, cU_coef,
				cdU_L2, cdU_coef, _eps) : PolynomialDistanceFunction::getValue(x, y, up, len, inv_len, delta, _LE, cD_L2,
					cD_coef, cdD_L2, cdD_coef, _eps));
		};
		void load(std::ifstream &in)
		{
			curveU.load(in);
			curveD.load(in);
		};
		void print(std::ofstream &out)
		{
			curveU.print(out);
			curveD.print(out);
		};
		void set(FlowTaskDescription* _task)
		{
			m_task = _task;
			param.delta = 0.25 * m_task->aa / m_task->nn;
			if(m_task->id_airfoil == kind::AIRFOIL_ROTOR) param.delta *= 0.25;// т.к. блок равномерной сетки для ротора измельчён в 4 раза
			param.length = m_task->L_BASE;
			param.LE = m_task->P_C.x - 0.5 * param.length;
			param.inv_length = 1.0 / param.length;
			param.eps = param.delta / 16.0;
			curveU.set(&basisBase, param.LE, m_task->P_C.y);
			curveD.set(&basisBase, param.LE, m_task->P_C.y);
		};
		BezierApproximation curveU;///< Аппроксимация верхней части границы профиля кривой Безье
		BezierApproximation curveD;///< Аппроксимация нижней части границы профиля кривой Безье
		BaseForBasisBezier basisBase;///< Базисная матрица Безье --- строится в конструкторе
		ForNumericalLSFunction param;///< Параметры численного расчета функции уровня
	  private:
		FlowTaskDescription* m_task;
	}; 
}

#endif //_LS_STAG_AIRFOIL_BEZIER_APPROXIMATION_H_
