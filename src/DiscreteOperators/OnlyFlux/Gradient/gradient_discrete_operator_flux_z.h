/*
* \copyright (C) 2026 Puzikova V.V.
* \file      gradient_discrete_operator_flux_z.h
* \author    Puzikova Valeria (vvp) 
* \date      05.06.2026 00:09:37
* \brief     gradient_discrete_operator_flux_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Z_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Z_H_

// ================================================================== INCLUDES
#include "gradient_discrete_operator_flux_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorFluxZ

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorFluxZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(computeFluxGz, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_inv_delta_ib, LS_STAG_CR_PC, m_theta,
		const int, m_N, const int, NM);
/*!
\brief Класс, реализующий дискретный оператор градиента на z-сетке с основной сетки

Данный класс реализует работу с дискретным аналогом оператора градиента на z-сетке с основной сетки
*/
	struct GradientDiscreteOperatorFluxZ : GradientDiscreteOperatorFluxBase
	{
		LS_STAG_KERNEL_2D(computeFluxGz, LS_STAG_R_PD flux, LS_STAG_CR_PC m_inv_delta_ib, LS_STAG_CR_PC m_theta,
			const int m_N, const int NM)
		{
			for(int num = i * m_N + j * NM + 1, bNM = (j) ? NM : 0, m_end = num + m_N - 1; num < m_end; num++)
			{
				const coordinate* curT = m_theta + num;
				flux[num] = (curT->z < EQ_DOUBLE) ? -m_inv_delta_ib[num].z : (((curT - bNM)->z < EQ_DOUBLE) ?  m_inv_delta_ib[num].z : 0.0);
			}
		}
		/*!
		Вычисляет значения на погруженных границах
		*/
		void calcElem()
		{		
			LS_STAG_PAR_FOR_2D_RANGE(GradientDiscreteOperatorFluxZ, computeFluxGz, 1, m_M - 1, 0, m_K, flux,
				m_inv_delta_ib, m_theta, m_N, NM)
		};
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_Z_H_
