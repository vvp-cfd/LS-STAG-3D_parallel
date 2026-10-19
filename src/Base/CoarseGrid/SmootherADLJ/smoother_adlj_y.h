/*
* \copyright (C) 2013 Puzikova V.V.
* \file      smoother_adlj_y.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 14:17:38
* \brief     smoother_adlj_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SMOOTHER_ADLJ_Y_H_
#define _LS_STAG_SMOOTHER_ADLJ_Y_H_

// ================================================================== INCLUDES
#include "smoother_adlj_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SmootherADLJy

Данный файл содержит в себе реализацию класса SmootherADLJy
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(tdma, LS_STAG_R_PD, _solution, LS_STAG_R_PD, _aux_sol, LS_STAG_CR_PD, _d,
		LS_STAG_CR_PD, _c, LS_STAG_CR_PD, _e, LS_STAG_CR_PD, _right_member, LS_STAG_CR_PI, _ho, const int, _M2, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(initOrder, LS_STAG_R_PI, _ho, LS_STAG_R_PI, _vo, const int, _N, const int, _M, const int, _NM);
/*!
\brief Класс, реализующий релаклированную вертикальную линию Якоби

Данный класс реализует итерацию релаклированной вертикальной линии Якоби
*/
	struct SmootherADLJy : SmootherADLJInterface
	{
		LS_STAG_KERNEL_2D(tdma, LS_STAG_R_PD _solution, LS_STAG_R_PD _aux_sol, LS_STAG_CR_PD _d, LS_STAG_CR_PD _c, LS_STAG_CR_PD _e,
			LS_STAG_CR_PD _right_member, LS_STAG_CR_PI _ho, const int _M2, const int _NM)
		{
			int is = i * _NM;
			int iY = is + j;
			int num = is + j * _M2;
			const int number2 = num + _M2 - 1;
			_aux_sol[num] = _right_member[iY] * _d[num];// 25/07/15  / d[num];
			for (int i2 = num, j2 = i2 + 1, number = number2 + 1; j2 < number; j2++, i2++)
				_aux_sol[j2] = (_right_member[_ho[j2]] - _c[i2] * _aux_sol[i2]) * _d[j2];// 25/07/15  / d[i];
			_solution[_ho[number2]] = _aux_sol[number2];
			num--;
			for (int j2 = number2 - 1; j2 > num; j2--)
				_solution[_ho[j2]] = _aux_sol[j2] - _e[j2] * _solution[_ho[j2 + 1]];
		}
		LS_STAG_KERNEL_2D(initOrder, LS_STAG_R_PI _ho, LS_STAG_R_PI _vo, const int _N, const int _M, const int _NM)
		{
			for (int j2 = 0, numS = i * _NM, number = numS + j * _N, num2 = numS + j; j2 < _N; j2++, number++)
			{
				int num = num2 + j2 * _M;
				_ho[num] = number;
				_vo[number] = num;
			}
		}
	protected:
	virtual void initOrder()
	{
		M2 = m_M - 1;
		LS_STAG_PAR_FOR_2D(SmootherADLJy, initOrder, m_K, m_M, ho, vo, m_N, m_M, NM)
	}
	public:
		void setSolver(SparseMatrix* _AH, extDouble *_r, extDouble *_c, RealValuedMeshFunction *_q, extDouble *_ch, const int _n_min)
		{
			SmootherADLJInterface::setSolver(_AH, _r, _c, _q, _ch, _n_min);
			m_q->getCopy(aux_sol);
		}
		void init(const int _N, const int _M, const int _K)
		{
			TriDiagonalMatrix::init(_N, _M, _K);
			m_ho.init(NMK, ho); 
			m_vo.init(NMK, vo);
			initOrder();
		};
		/*!
		Выполняет вертикальную прогонку
		*/
		void TDMA()
		{
			LS_STAG_PROFILE("TDMAy");
			LS_STAG_PAR_FOR_2D(SmootherADLJy, tdma, m_K, m_N, solution, aux_sol, d, c, e, right_member, ho, m_M, NM)
		};
		double* aux_sol;
		Array<int> m_ho, m_vo;
		int *ho;///< Массив с горизонтально упорядоченной нумерация узлов на предыдущей подробной сетке
		int *vo;///< Массив с номерами на вертикально упорядоченной сетке (для ускорения перестроения)
		int M2;
	};
}

#endif //_LS_STAG_SMOOTHER_ADLJ_Y_H_
