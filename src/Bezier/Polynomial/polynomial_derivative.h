/*
* \copyright (C) 2013 Puzikova V.V.
* \file      polynomial_derivative.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 23:05:55
* \brief     polynomial_derivative.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_POLYNOMIAL_DERIVATIVE_H_
#define _LS_STAG_POLYNOMIAL_DERIVATIVE_H_

// ================================================================== INCLUDES
#include "../../Base/ls_stag_parallel.h"
#include "polynomial.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса PolynomialDerivative

Данный файл содержит в себе реализацию класса PolynomialDerivative
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calc, LS_STAG_R_PD, _coef, LS_STAG_CR_PD, _antiderivative_coef, const int, _length);
/*!
\brief Класс, реализующий работу с производной полинома

Данный класс реализует построение и хранение полинома, являющегося производной для исходного полинома
*/
	struct PolynomialDerivative : Polynomial
	{
	  private:
		void init(const int _length)
		{
			Polynomial::init(_length - 1);
		};
	  public:
		  LS_STAG_KERNEL(calc, LS_STAG_R_PD _coef, LS_STAG_CR_PD _antiderivative_coef, const int _length)
		  {
			  _coef[i] = _antiderivative_coef[i] * (_length - i);
		  }
		  /*!
		  Расчет коэффициентов полинома, являющегося производной для исходного
		  */
		void build()
		{
			LS_STAG_PAR_FOR(PolynomialDerivative, calc, length, coef, m_antiderivative->coef, length)
			m_coef.setCopy(pCoef);
		};
		void set(Polynomial* _antiderivative)
		{
			m_antiderivative = _antiderivative;
			init(m_antiderivative->length);
			build();
		};
		Polynomial* m_antiderivative;///< Указатель на полином, для которого строится производная
	};
}

#endif //_LS_STAG_POLYNOMIAL_DERIVATIVE_H_
