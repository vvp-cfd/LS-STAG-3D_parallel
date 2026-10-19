/*
* \copyright (C) 2013 Puzikova V.V.
* \file      smoother_adlj_x.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 14:17:13
* \brief     smoother_adlj_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SMOOTHER_ADLJ_X_H_
#define _LS_STAG_SMOOTHER_ADLJ_X_H_

// ================================================================== INCLUDES
#include "smoother_adlj_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SmootherADLJx

Данный файл содержит в себе реализацию класса SmootherADLJx
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(tdma, LS_STAG_R_PD, _solution, LS_STAG_CR_PD, _d, LS_STAG_CR_PD, _c,
		LS_STAG_CR_PD, _e, LS_STAG_CR_PD, _right_member, const int, _N, const int, _NM);
/*!
\brief Класс, реализующий релаклированную горизонтальную линию Якоби

Данный класс реализует итерацию релаклированной горизонтальной линии Якоби
*/
	struct SmootherADLJx : SmootherADLJInterface
	{
		/*!
		Выполняет горизонтальную прогонку
		*/
		LS_STAG_KERNEL_2D(tdma, LS_STAG_R_PD _solution, LS_STAG_CR_PD _d, LS_STAG_CR_PD _c, LS_STAG_CR_PD _e, 
			LS_STAG_CR_PD _right_member, const int _N, const int _NM)
		{
			int num = j * _N + i * _NM;// 27/09/15 M;
			_solution[num] = _right_member[num] * _d[num];// 25/07/15 / d[num];
			for (int i2 = num, j2 = i2 + 1, number = i2 + _N; j2 < number; j2++, i2++)
				_solution[j2] = (_right_member[j2] - _c[i2] * _solution[i2]) * _d[j2];
			num--;
			for (int j2 = num + _N - 1; j2 > num; j2--)// 27/09/15 M - 1; i > num; i--)
				_solution[j2] -= _e[j2] * _solution[j2 + 1];
		}
		void TDMA()
		{
			LS_STAG_PROFILE("TDMAx");
			LS_STAG_PAR_FOR_2D(SmootherADLJx, tdma, m_K, m_M, solution, d, c, e, right_member, m_N, NM)
		};
	};
}

#endif //_LS_STAG_SMOOTHER_ADLJ_X_H_
