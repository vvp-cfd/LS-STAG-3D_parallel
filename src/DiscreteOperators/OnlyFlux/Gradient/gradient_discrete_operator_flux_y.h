/*
* \copyright (C) 2026 Puzikova V.V.
* \file      gradient_discrete_operator_flux_y.h
* \author    Puzikova Valeria (vvp) 
* \date      05.06.2026 00:01:05
* \brief     gradient_discrete_operator_flux_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Y_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Y_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_flux_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorFluxY

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorFluxY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(computeFluxGy, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_inv_delta_ib, LS_STAG_CR_PC, m_theta,
		const int, m_N, const int, NM);
/*!
\brief Класс, реализующий дискретный оператор градиента на y-сетке с основной сетки

Данный класс реализует работу с дискретным аналогом оператора градиента на y-сетке с основной сетки
*/
	struct GradientDiscreteOperatorFluxY : GradientDiscreteOperatorFluxBase
	{
		LS_STAG_KERNEL_2D(computeFluxGy, LS_STAG_R_PD flux, LS_STAG_CR_PC m_inv_delta_ib, LS_STAG_CR_PC m_theta,
			const int m_N, const int NM)
		{
			for(int numY = i * m_N + j * NM + 1, num = numY + j * m_N, m_end = numY + m_N - 1; numY < m_end; numY++, num++)
			{
				const coordinate* curT = m_theta + num;
				flux[numY] = (curT->y < EQ_DOUBLE) ? -m_inv_delta_ib[num].y : (((curT - m_N)->y < EQ_DOUBLE) ?  m_inv_delta_ib[num].y : 0.0);
			}
		}
		/*!
		Вычисляет значения на погруженных границах
		*/
		void calcElem()
		{		
			LS_STAG_PAR_FOR_2D_RANGE(GradientDiscreteOperatorFluxY, computeFluxGy, 1, m_M - 1, 0, m_K, flux,
				m_inv_delta_ib, m_theta, m_N, NM)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Y_H_
