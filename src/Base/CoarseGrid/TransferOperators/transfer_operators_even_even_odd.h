/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_even_even_odd.h
* \author    Puzikova Valeria (vvp) 
* \date      29.05.2026 22:23:44
* \brief     transfer_operators_even_even_odd.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_ODD_H_
#define _LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_ODD_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsEvenEvenOdd

Данный файл содержит в себе реализацию класса TransferOperatorsEvenEvenOdd
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPeeo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numReeo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitPeeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR, const int, NMR);
	LS_STAG_DECLARATION_S_KERNEL(portraitReeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR,
		const int, NP, const int, NR, const int, K2, const int, NMP);
/*!
\brief Класс, реализующий переход с сеточного уровня с четным числом узлов

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит четное число узлов как по вертикали, так и по горизонтали, но нечетное по z
*/
	struct TransferOperatorsEvenEvenOdd : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsEvenEvenOdd():
		TransferOperatorsInterface(EvenEvenOdd, sizeof(TransferOperatorsEvenEvenOdd))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPeeo, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::NumWidth4xy(1, 2, 4, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth4xy(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth4xy(1, 2, 4, num, _Number, N2, M2);
		}
	LS_STAG_KERNEL_INLINE void NumWidth_1xy(const int N2, const int M2, const int m_a, int& num, int* _Number)
    {
			for (int j = 0; j < M2; j++)
			SparseMatrix::NumWidth_1(N2, m_a, num, _Number);
    }	
		LS_STAG_SINGLE_KERNEL(numReeo, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			TransferOperatorsEvenEvenOdd::NumWidth_1xy(N2, M2, 8, num, _Number);
			for (int k = 1; k < K2; k++)
				TransferOperatorsEvenEvenOdd::NumWidth_1xy(N2, M2, 12, num, _Number);
			TransferOperatorsEvenEvenOdd::NumWidth_1xy(N2, M2, 8, num, _Number);
		}
		LS_STAG_SINGLE_KERNEL(portraitPeeo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR, const int NMR)
		{
			int num = 0, a = 0, b = 1, c = NR, d = c + 1;
			for (int k = 0, aa = NMR, bb = aa + 1, cc = aa + NR, dd = cc + 1; 
				k < K2; k++, c = a + NR, d = c + 1, cc = aa + NR, dd = cc + 1)
			{
			TransferOperatorsInterface::PortraitWidth2D4D(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
			a = b - NMR; b = a + 1; c = a + NR, d = c + 1;
			TransferOperatorsInterface::PortraitWidth4D8D(0.5, 0.5, 0.375, 0.125, 0.375, 0.125, 0.28125, 0.09375, 0.09375,
				0.03125, 0.28125, 0.09375, 0.09375, 0.03125, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
		}
			TransferOperatorsInterface::PortraitWidth2D4D(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
		};
		LS_STAG_SINGLE_KERNEL(portraitReeo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int MR, const int NP, const int NR, const int K2, const int NMP)
		{
			int aa = 0;
			int bb = 1;
			int cc = NP;
			int dd = cc + 1;
			int a = NMP;
			int b = a + 1;
			int c = a + NP;
			int d = c + 1;
			int aaa = NMP;
			int bbb = aaa + 1;
			int ccc = aaa + NP;
			int ddd = ccc + 1;
			int num = 0;
			TransferOperatorsInterface::PortraitWidth8_1xy(0.25, 0.25, 0.125, 0.125, aa, bb, cc, dd, aaa, bbb, ccc, ddd, num, _Cell, _Portrait, NR, MR, NP);
			aa = aaa; bb = bbb; cc = ccc; dd = ddd; aaa += NMP; bbb = aaa + 1; ccc = aaa + NP; ddd = ccc + 1;
			for (int k = 1; k < K2; k++, a = aa, b = bb, c = cc, d = dd, aa = aaa, bb = bbb, cc = ccc, dd = ddd, aaa += NMP, bbb = aaa + 1, ccc = aaa + NP, ddd = ccc + 1)
			{
			for (int i = 0; i < MR; i++, a = c, aa = cc, aaa = ccc, b = a + 1, bb = aa + 1, bbb = aaa + 1, c = a + NP, cc = aa + NP, ccc = aaa + NP, d = c + 1, dd = cc + 1, ddd = ccc + 1)
			{
				TransferOperatorsInterface::PortraitWidth12_1(0.125, 0.125, 0.25, 0.25, 0.125, 0.125, a, b, c, d, aa, bb, cc, dd, aaa, bbb, ccc, ddd, num, _Cell, _Portrait, NR);
			}//for j
		}
		TransferOperatorsInterface::PortraitWidth8_1xy(0.125, 0.125, 0.25, 0.25, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, NR, MR, NP);	
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenOdd, numPeeo, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenOdd, portraitPeeo, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR, m_NR * m_MR)
		};
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenOdd, numReeo, m_R->Number, m_NR, m_MR, K2)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenEvenOdd, portraitReeo, m_R->Portrait, m_R->Cell, m_MR, m_NP, m_NR, K2, m_NP * m_MP)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_EVEN_EVEN_ODD_H_
