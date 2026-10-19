/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_even_even_even.h
* \author    Puzikova Valeria (vvp) 
* \date      23.05.2026 16:01:44
* \brief     transfer_operators_even_even_even.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_EVEN_H_
#define _LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_EVEN_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsEvenEvenEven

Данный файл содержит в себе реализацию класса TransferOperatorsEvenEvenEven
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numP, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numR, LS_STAG_R_PI, _Number, const int, NMR);
	LS_STAG_DECLARATION_S_KERNEL(portraitP, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR, const int, NMR);
	LS_STAG_DECLARATION_S_KERNEL(portraitR, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR,
		const int, NP, const int, NR, const int, KR, const int, NMP);
/*!
\brief Класс, реализующий переход с сеточного уровня с четным числом узлов

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит четное число узлов как по вертикали, так и по горизонтали, z
*/
	struct TransferOperatorsEvenEvenEven : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsEvenEvenEven():
		TransferOperatorsInterface(EvenEvenEven, sizeof(TransferOperatorsEvenEvenEven))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numP, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			TransferOperatorsInterface::NumWidth4xy(1, 2, 4, num, _Number, N2, M2);
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::NumWidth4xy(2, 4, 8, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth4xy(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth4xy(1, 2, 4, num, _Number, N2, M2);
		}
		LS_STAG_SINGLE_KERNEL(numR, LS_STAG_R_PI _Number, const int NMR)
		{
			for (int j = 0; j < NMR; j++)
				_Number[j + 1] = _Number[j] + 8;
		}
		LS_STAG_SINGLE_KERNEL(portraitP, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR, 
			const int NMR)
		{
			int a = 0;
			int b = 1;
			int c = NR;
			int d = c + 1;
			int aa = NMR;
			int bb = aa + 1;
			int cc = aa + NR;
			int dd = cc + 1;
			int num = 0;
			TransferOperatorsInterface::PortraitWidth2D4D(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
			a = 0; b = 1; 
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::PortraitWidth4D8D(0.75, 0.25, 0.5625, 0.1875, 0.1875, 0.0625, 0.421875, 0.140625, 
				0.140625, 0.046875, 0.140625, 0.046875, 0.046875, 0.015625, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
			aa = a; bb = b; a -= NMR; b = a + 1; 
			TransferOperatorsInterface::PortraitWidth4D8D(0.25, 0.75, 0.1875, 0.0625, 0.5625, 0.1875, 0.140625, 0.046875, 
				0.046875, 0.015625, 0.421875, 0.140625, 0.140625, 0.046875, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
		}
			c = a + NR; d = c + 1;
			TransferOperatorsInterface::PortraitWidth2D4D(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
		};
		LS_STAG_SINGLE_KERNEL(portraitR, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int MR, const int NP, const int NR, const int KR, const int NMP)
		{
			for (int k = 0, num = 0, a  = 0, b = 1, c = NP, d = c + 1, aa = NMP, bb = aa + 1, cc = aa + NP, dd = cc + 1; 
				k < KR; k++, a = aa, b = bb, c = cc, d = dd, aa += NMP, bb = aa + 1, cc = aa + NP, dd = cc + 1)
			{
				TransferOperatorsInterface::PortraitWidth8_1xy(0.125, 0.125, 0.125, 0.125, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, NR, MR, NP);
		}
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenEven, numP, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenEven, portraitP, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR, m_NR * m_MR)
		};
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenEven, numR, m_R->Number, m_NMKR)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenEven, portraitR, m_R->Portrait, m_R->Cell, m_MR, m_NP, m_NR, m_KR, m_NP * m_MP)
		};
		void buildR() override
		{
			buildFullR();
		}
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_EVEN_H_
