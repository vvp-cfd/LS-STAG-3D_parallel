/*
* \copyright (C) 2013 Puzikova V.V.
* \file      gradient_discrete_operator_y.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 13:15:16
* \brief     gradient_discrete_operator_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_Y_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_Y_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorY

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calc, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell, LS_STAG_R_PC, _theta, 
		LS_STAG_R_PC, _delta, LS_STAG_CR_PCT, _type, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildGy, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, baseN, 
		const int, NM, const int, m_K);
/*!
\brief Класс, реализующий дискретный оператор градиента на у-сетке

Данный класс реализует работу с дискретным аналогом оператора градиента на у-сетке
*/
	struct GradientDiscreteOperatorY : GradientDiscreteOperatorBase
	{
		LS_STAG_SINGLE_KERNEL(portraitBuildGy, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int baseN, const int NM,
		const int m_K)
		{
		for(int k = 0, p = 0, n = baseN, numA2 = 0, num = 0; k < m_K; k++, p += baseN, n += baseN)
		{
		for (int i = 0; i < NM; i++, p++, n++, numA2++, num++)
		{
				if (m_type[num] == SolidCell)
					Portrait[numA2] = p;
				else
				{
					Portrait[numA2] = p;
					numA2++;
					Portrait[numA2] = n;
				}
			}
		}
		}
		LS_STAG_KERNEL(calc, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, LS_STAG_R_PC _theta, LS_STAG_R_PC _delta,
			LS_STAG_CR_PCT _type, const int _N, const int _NM)
		{
			const int numF = i * _NM; // на y-сетке
			const int numFbase = numF + i * _N;
			coordinate* curD = _delta + numFbase;
			coordinate* curT = _theta + numFbase;
			for (int num = numF, m_end = num + _NM; num < m_end; num++, curD++, curT++)
			{
			const int numA2 = _Number[num];
			if (_type[num] == SolidCell)
				_Cell[numA2] = 0.0;
			else
			{
				const double cell = curT->y * curD->x * curD->z; 
				_Cell[numA2] = cell;
				_Cell[numA2 + 1] = -cell;
			}
			}
		}
	  protected:
		  /*!
		  Строит портрет разреженной матрицы дискретного аналога оператора градиента на у-сетке
		  */
		void buildPortait()
		{
			LS_STAG_SINGLE(GradientDiscreteOperatorY, portraitBuildGy, Portrait, m_type, baseN, m_N * m_M, m_K)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора градиента на у-сетке
		  */
		virtual void calcElem()
		{
			LS_STAG_PAR_FOR(GradientDiscreteOperatorY, calc, m_K, Number, Cell, m_theta, m_delta, m_type, m_N, m_N * m_M)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_Y_H_
