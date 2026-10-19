/*
* \copyright (C) 2013 Puzikova V.V.
* \file      polynomial_distance_function.h
* \author    Puzikova Valeria (vvp) 
* \date      18.10.2013 0:44:32
* \brief     polynomial_distance_function.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_POLYNOMIAL_DISTANCE_FUNCTION_H_
#define _LS_STAG_POLYNOMIAL_DISTANCE_FUNCTION_H_

// ================================================================== INCLUDES
#include "polynomial_derivative.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием классов PolynomialDistanceFunction и ForNumericalLSFunction

Данный файл содержит в себе реализацию классов PolynomialDistanceFunction и ForNumericalLSFunction
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Структура параметров для численного расчета функции уровня

Вспомогательная структура параметров для численного расчета функции уровня
*/
	struct ForNumericalLSFunction
	{
		double delta;
		double length;
		double inv_length;// = 1 / length
		double LE;
		double eps; // delta/16.0
	};
/*!
\brief Класс, реализующий работу с полиномиальной функцией расстояния

Данный класс реализует расчет расстояния от точки расчетной области до погруженной границы,
аппроксимированной полиномиальной функцией, и хранение вспомогательных данных для данной процедуры.

В 3D достаточно вычислить это расстояние для первого слоя, выше просто используем параллельный перенос.
*/
	struct PolynomialDistanceFunction
	{
		PolynomialDistanceFunction():
		length(0)
		{
			;
		}
		LS_STAG_KERNEL_INLINE double getDistanceDerivative(const double m_x, const double m_y, const double ts, const double tus,
			const double len, const double inv_len, const double _LE, const int c_L2, LS_STAG_CR_PD c_coef, const int cd_L2,
			LS_STAG_CR_PD cd_coef)
		{
			return m_x - ts * len + _LE + (Polynomial::getValue(tus, c_L2, c_coef) - m_y) *
				Polynomial::getValue(tus, cd_L2, cd_coef) * inv_len;
		}
		/*!
		\return Расстояние от точки с координатами (_x,_y) до погруженной границы, аппроксимированной полиномом
		*/
		LS_STAG_KERNEL_INLINE double getValue(const double _x, const double _y, const bool up, const double len,
			const double inv_len, const double delta, const double _LE, const int c_L2, LS_STAG_CR_PD c_coef, const int cd_L2,
			LS_STAG_CR_PD cd_coef, const double _eps)
		{
			double x1 = _x - delta;
			double x2 = _x + delta;
			double t1 = (x1 - _LE) * inv_len;
			double tu1 = (up) ? (1.0 - t1) : t1; // setUpDown();
			double aux = PolynomialDistanceFunction::getDistanceDerivative(_x, _y, t1, tu1, len, inv_len, _LE, c_L2, c_coef,
				cd_L2, cd_coef);
			double t2 = (x2 - _LE) * inv_len;
			double tu2 = (up) ? (1.0 - t2) : t2; // setUpDown();
			double aux2 = PolynomialDistanceFunction::getDistanceDerivative(_x, _y, t2, tu2, len, inv_len, _LE, c_L2, c_coef,
				cd_L2, cd_coef);
			double ts = 0.5 * (t1 + t2); // setMiddle();
			double tus = (up) ? (1.0 - ts) : ts; //setUpDown(); 
			double auxs = PolynomialDistanceFunction::getDistanceDerivative(_x, _y, ts, tus, len, inv_len, _LE, c_L2, c_coef,
				cd_L2, cd_coef);
			while (fabs(x2 - x1) > _eps)
			{
				if (aux * auxs < 0.0)
				{
					t2 = ts; tu2 = tus; aux2 = auxs; // shiftLeft();
				}
				else
				{
					t1 = ts; tu1 = tus; aux = auxs; // shiftRight();
				}
				ts = 0.5 * (t1 + t2); tus = (up) ? (1.0 - ts) : ts; //setUpDown(); // setMiddle();
				x1 = t1 * len + _LE;
				x2 = t2 * len + _LE;
				auxs = PolynomialDistanceFunction::getDistanceDerivative(_x, _y, ts, tus, len, inv_len, _LE, c_L2, c_coef,
					cd_L2, cd_coef);
			}
			x1 = ts * len + _LE;
			x2 = Polynomial::getValue(tus, c_L2, c_coef);
			aux = (_x - x1) * (_x - x1) + (_y - x2) * (_y - x2);
			aux2 = sqrt(aux);
			return aux2;
		};
	  protected:
		virtual void build()
		{
			curve_derivative.set(&curve);
		};
		void init()
		{// length уже считано вместе с массивом координат
			curve.init(length);
		};
	  public:
		Polynomial curve;                      ///< Полином, аппроксимирующий погруженную границу                  
		PolynomialDerivative curve_derivative; ///< Производная полинома, аппроксимирующего погруженную границу
		int length;
	};
}

#endif //_LS_STAG_POLYNOMIAL_DISTANCE_FUNCTION_H_
