/*
* \copyright (C) 2013 Puzikova V.V.
* \file      gradient_discrete_operator_x.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 13:13:07
* \brief     gradient_discrete_operator_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_X_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_X_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorX

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(calc, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
			LS_STAG_R_PC, _delta, LS_STAG_R_PC, _theta, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildGx, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, m_M, 
		const int, m_N, const int, m_K);
/*!
\brief Класс, реализующий дискретный оператор градиента на х-сетке

Данный класс реализует работу с дискретным аналогом оператора градиента на х-сетке
*/
	struct GradientDiscreteOperatorX : GradientDiscreteOperatorBase
	{
		LS_STAG_SINGLE_KERNEL(portraitBuildGx, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int m_M, 
			const int m_N, const int m_K)
		{
			for(int k = 0, p = 0, e = 1, num = 0, numA2 = 0; k < m_K; k++)
			{
			for (int j = 0; j < m_M; j++, p++, e++)
			{
				for (int i = 0; i < m_N; i++, p++, e++, num++, numA2++)
				{
					if (m_type[num] == SolidCell)
						Portrait[numA2] = p;
					else
					{
						Portrait[numA2] = p;
						numA2++;
						Portrait[numA2] = e;
					}
				}
			}
		}
		}
		LS_STAG_KERNEL_2D(calc, LS_STAG_CR_PCT _type, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell,
			LS_STAG_R_PC _delta, LS_STAG_R_PC _theta, const int _N, const int _M, const int _NM)
		{
			const int numF = _N * i + j * _NM; // на x-сетке
			const int numFbase = numF + i + j * _M;
			coordinate* curD = _delta + numFbase;
			coordinate* curT = _theta + numFbase;
			for (int num = numF, m_end = num + _N, numA2 = _Number[num]; num < m_end; num++, numA2++, curD++, curT++)
			{
				if (_type[num] == SolidCell)
					_Cell[numA2] = 0.0;
				else
				{
					const double cell = curT->x * curD->y * curD->z;
					_Cell[numA2] = cell;
					numA2++;
					_Cell[numA2] = -cell;
				}
			}
		}
	  private:
		void setBaseN(){ baseN = m_N + 1; };
	  protected:
		  /*!
		  Строит портрет разреженной матрицы дискретного аналога оператора градиента на х-сетке
		  */
		void buildPortait()
		{
			LS_STAG_SINGLE(GradientDiscreteOperatorX, portraitBuildGx, Portrait, m_type, m_M, m_N, m_K)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора градиента на х-сетке
		  */
		virtual void calcElem()
		{
			LS_STAG_PAR_FOR_2D(GradientDiscreteOperatorX, calc, m_M, m_K, m_type, Number, Cell, m_delta, m_theta, m_N, m_M, m_N * m_M)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_X_H_
