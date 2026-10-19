/*
* \copyright (C) 2026 Puzikova V.V.
* \file      gradient_discrete_operator_flux_x.h
* \author    Puzikova Valeria (vvp) 
* \date      04.06.2026 23:39:37
* \brief     gradient_discrete_operator_flux_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_X_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_X_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_flux_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorFluxX

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorFluxX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(computeFluxGx, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_inv_delta_ib, LS_STAG_CR_PC, m_theta,
		const int, m_N, const int, m_M, const int, NM);
/*!
\brief Класс, реализующий дискретный оператор градиента на х-сетке с основной сетки

Данный класс реализует работу с дискретным аналогом оператора градиента на х-сетке с основной сетки
*/
	struct GradientDiscreteOperatorFluxX : GradientDiscreteOperatorFluxBase
	{
		LS_STAG_KERNEL_2D(computeFluxGx, LS_STAG_R_PD flux, LS_STAG_CR_PC m_inv_delta_ib, LS_STAG_CR_PC m_theta,
			const int m_N, const int m_M, const int NM)
		{
			for(int numX = i * m_N + j * NM + 1, num = numX + i + j * m_M, m_end = numX + m_N - 1; numX < m_end; numX++, num++)
			{
				const coordinate* curT = m_theta + num;
				flux[numX] = (curT->x < EQ_DOUBLE) ? -m_inv_delta_ib[num].x : (((curT - 1)->x < EQ_DOUBLE) ?  m_inv_delta_ib[num].x : 0.0);
			}
		}
		/*!
		Вычисляет значения на погруженных границах
		*/
		void calcElem()
		{		
			LS_STAG_PAR_FOR_2D(GradientDiscreteOperatorFluxX, computeFluxGx, m_M, m_K, flux, m_inv_delta_ib, m_theta, m_N, m_M, NM)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_X_H_
