/*
* \copyright (C) 2026 Puzikova V.V.
* \file      smoother_adlj_z.h
* \author    Puzikova Valeria (vvp) 
* \date      20.05.2026 15:15:15
* \brief     smoother_adlj_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SMOOTHER_ADLJ_Z_H_
#define _LS_STAG_SMOOTHER_ADLJ_Z_H_

// ================================================================== INCLUDES
#include "smoother_adlj_y.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SmootherADLJz

Данный файл содержит в себе реализацию класса SmootherADLJz
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(initOrderZ, LS_STAG_R_PI, _ho, LS_STAG_R_PI, _vo, const int, _N, const int, _K, const int, _NM, const int, _NK);
/*!
\brief Класс, реализующий релаклированную вертикальную линию Якоби

Данный класс реализует итерацию релаклированной вертикальной линии Якоби
*/
	struct SmootherADLJz : SmootherADLJy
	{
		LS_STAG_KERNEL_2D(initOrderZ, LS_STAG_R_PI _ho, LS_STAG_R_PI _vo, const int _N, const int _K, const int _NM, const int _NK)
		{
			for (int j2 = 0, number = j * _NM + i * _N, numS = j + i *_NK; j2 < _N; j2++, number++)
			{
				int num = numS + j2 * _K;
				_ho[num] = number;
				_vo[number] = num;
			}
		}
	protected:
	void initOrder()
	{
		M2 = m_K - 1;
		NK = m_N * m_K;
		LS_STAG_PAR_FOR_2D(SmootherADLJz, initOrderZ, m_M, m_K, ho, vo, m_N, m_K, NM, NK)
	}
	public:
		/*!
		Выполняет вертикальную прогонку
		*/
		void TDMA()
		{
			LS_STAG_PROFILE("TDMAz");
			LS_STAG_PAR_FOR_2D(SmootherADLJy, tdma, m_M, m_N, solution, aux_sol, d, c, e, right_member, ho, m_K, NK)
		};
		int NK;
	};
}

#endif //_LS_STAG_SMOOTHER_ADLJ_Y_H_
