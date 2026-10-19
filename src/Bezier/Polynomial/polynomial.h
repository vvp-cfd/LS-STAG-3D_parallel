/*
* \copyright (C) 2013 Puzikova V.V.
* \file      polynomial.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 22:53:40
* \brief     polynomial.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_POLYNOMIAL_H_
#define _LS_STAG_POLYNOMIAL_H_

// ================================================================== INCLUDES
#include <string.h>
#include <fstream>
#include "../../Base/MeshFunctions/array.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса Polynomial

Данный файл содержит в себе реализацию класса Polynomial
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий работу с полиномами

Данный базовый класс реализует работу с полиномами
*/
	struct Polynomial
	{
		Polynomial():
		length(0)
		{
			;
		}
		/*!
		\param t точка, в которой вычисляется значение полинома
		\return Значение полинома в точке t
		*/
		LS_STAG_KERNEL_INLINE double getValue(const double t, const int _L2, LS_STAG_CR_PD _coef)
		{
			double res = 0.0;
			double aux = 1.0;
			for (int i = _L2; i > -1; i--)//0-й эл-т coef --- коэф-т при старшей степени 
			{
				res += aux * _coef[i];
				aux *= t;
			}
			return res;
		};
		virtual void init(const int _length)
		{
			length = _length;
			L2 = length - 1;
			m_coef.init(length, coef);
			m_coef.setCopy(pCoef);
		};
		void print(const std::string &name)
		{
			m_coef.print(name);
		};
		Array<double> m_coef;
		double* coef;///< Массив, в котором хранятся коэффициенты полинома
		double* pCoef;// CPU версия coef
		int L2;   ///< Степень полинома
		int length;///< Число коэффициентов полинома
	};
}

#endif //_LS_STAG_POLYNOMIAL_H_
