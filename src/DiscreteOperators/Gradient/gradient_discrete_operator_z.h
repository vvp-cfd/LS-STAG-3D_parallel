/*
* \copyright (C) 2026 Puzikova V.V.
* \file      gradient_discrete_operator_z.h
* \author    Puzikova Valeria (vvp) 
* \date      05.06.2026 14:54:23
* \brief     gradient_discrete_operator_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_Z_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_Z_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorZ

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calc, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _delta, LS_STAG_CR_PCT, _type);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildGz, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, baseNM, const int, NMK);
/*!
\brief Класс, реализующий дискретный оператор градиента на z-сетке

Данный класс реализует работу с дискретным аналогом оператора градиента на z-сетке
*/
	struct GradientDiscreteOperatorZ : GradientDiscreteOperatorBase
	{
		LS_STAG_SINGLE_KERNEL(portraitBuildGz, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int baseNM, const int NMK)
		{
			for (int p = 0, t = baseNM, numA2 = 0; p < NMK; p++, t++, numA2++)
			{
				if (m_type[p] == SolidCell)
					Portrait[numA2] = p;
				else
				{
					Portrait[numA2] = p;
					numA2++;
					Portrait[numA2] = t;
				}
			}
		}
		LS_STAG_KERNEL(calc, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, LS_STAG_CR_PC _theta, LS_STAG_CR_PC _delta,
			LS_STAG_CR_PCT _type)
		{
			const int numA2 = _Number[i];
			if (_type[i] == SolidCell)
				_Cell[numA2] = 0.0;
			else
			{
				const coordinate* curD = _delta + i;
				const double cell = (_theta + i)->z * curD->x * curD->y;
				_Cell[numA2] = cell;
				_Cell[numA2 + 1] = -cell;
			}
		}
	  protected:
		  /*!
		  Строит портрет разреженной матрицы дискретного аналога оператора градиента на z-сетке
		  */
		void buildPortait()
		{
			LS_STAG_SINGLE(GradientDiscreteOperatorZ, portraitBuildGz, Portrait, m_type, m_N * m_M, NMK)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора градиента на z-сетке
		  */
		virtual void calcElem()
		{
			LS_STAG_PAR_FOR(GradientDiscreteOperatorZ, calc, NMK, Number, Cell, m_theta, m_delta, m_type)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_Y_H_
