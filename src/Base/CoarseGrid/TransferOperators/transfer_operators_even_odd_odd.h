/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_even_odd_odd.h
* \author    Puzikova Valeria (vvp) 
* \date      23.05.2026 19:21:20
* \brief     transfer_operators_even_odd_odd.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_ODD_H_
#define _LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_ODD_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsEvenOddOdd

Данный файл содержит в себе реализацию класса TransferOperatorsEvenOddOdd
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPeoo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numReoo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitPeoo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR);
	LS_STAG_DECLARATION_S_KERNEL(portraitReoo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2, const int, K2,
		const int, NP, const int, NR, const int, NMP);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным M, K и четным N

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов по вертикали и z и четное --- по горизонтали
*/
	struct TransferOperatorsEvenOddOdd : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsEvenOddOdd():
		TransferOperatorsInterface(EvenOddOdd, sizeof(TransferOperatorsEvenOddOdd))
		{
			;
		}
	LS_STAG_KERNEL_INLINE void NumWidth1_o(const int m_a, const int m_b, int& num, int* _Number, const int _NR, const int _M2)
    {
			SparseMatrix::NumWidth_1(_NR, m_a, num, _Number);
			for (int j = 1; j < _M2; j++)
				SparseMatrix::NumWidth_1(_NR, m_b, num, _Number);
			SparseMatrix::NumWidth_1(_NR, m_a, num, _Number);
    }	
		LS_STAG_SINGLE_KERNEL(numPeoo, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			for (int k = 0; k < K2; k++)
			{
				TransferOperatorsInterface::NumWidth4_o(1, 2, 4, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth4_o(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth4_o(1, 2, 4, num, _Number, N2, M2);
		}
		LS_STAG_SINGLE_KERNEL(numReoo, LS_STAG_R_PI _Number, const int NR, const int M2, const int K2)
		{
			int num = 0;
			TransferOperatorsEvenOddOdd::NumWidth1_o(8, 12, num, _Number, NR, M2);
			for (int k = 1; k < K2; k++)
			   TransferOperatorsEvenOddOdd::NumWidth1_o(12, 18, num, _Number, NR, M2);
			TransferOperatorsEvenOddOdd::NumWidth1_o(8, 12, num, _Number, NR, M2);
		}
		LS_STAG_SINGLE_KERNEL(portraitPeoo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR)
		{
			int a = 0;
			int b = 1;
			int c = NR;
			int d = c + 1;
			int aa = 0;
			int bb = 1;
			int cc = NR;
			int dd = cc + 1;
			int num = 0;
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::PortraitWidth2D4D_o(aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			TransferOperatorsInterface::PortraitWidth4D8D_o(0.5, 0.5, 0.375, 0.125, 0.375, 0.125, 0.25, 0.25, 
				0.1875, 0.0625, 0.1875, 0.0625, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
				aa = a; bb = b; cc = c; dd = d;
		}
			TransferOperatorsInterface::PortraitWidth2D4D_o(a, b, c, d, num, _Cell, _Portrait, N2, M2);
		};
		LS_STAG_SINGLE_KERNEL(portraitReoo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int M2, const int K2, const int NP, const int NR, const int NMP)
		{
			int a = NMP + NP;
			int b = a + 1;
			int c = NMP;
			int d = c + 1;
			int e = a;
			int f = e + 1;
			int aa = NP;
			int bb = aa + 1;
			int cc = 0;
			int dd = cc + 1;
			int ee = NP;
			int ff = ee + 1;
			int aaa = a;
			int bbb = aaa + 1;
			int ccc = NMP;
			int ddd = ccc + 1;
			int eee = aaa;
			int fff = eee + 1;
			int num = 0;
			TransferOperatorsInterface::PortraitWidth8_1(0.5, 0.25, 0.25, 0.125, cc, dd, ee, ff, ccc, ddd, eee, fff, num, _Cell, _Portrait, NR);
			cc += NP; dd = cc + 1; ee += NP; ff = ee + 1; ccc += NP; ddd = ccc + 1; eee += NP; fff = eee + 1;
			for (int i = 1; i < M2; i++, aa += NP, bb = aa + 1, aaa += NP, bbb = aaa + 1, ee += NP, ff = ee + 1, eee += NP, fff = eee + 1,
				cc += NP, dd = cc + 1, ccc += NP, ddd = ccc + 1)
			    TransferOperatorsInterface::PortraitWidth12_1(0.25, 0.5, 0.25, 0.125, 0.25, 0.125, aa, bb, cc, dd, ee, ff, aaa, bbb, ccc, ddd, eee, fff,
			num, _Cell, _Portrait, NR);
			TransferOperatorsInterface::PortraitWidth8_1(0.25, 0.5, 0.125, 0.25, aa, bb, cc, dd, aaa, bbb, ccc, ddd, num, _Cell, _Portrait, NR);
			cc += NMP, aa = cc + NP, bb = aa + 1, dd = cc + 1, ee = aa, ff = ee + 1;
			ccc += NMP, aaa = ccc + NP, bbb = aaa + 1, ddd = ccc + 1, eee = aaa, fff = eee + 1;
			for (int k = 1; k < K2; k++, c = cc, a = c + NP, b = a + 1, d = c + 1, e = a, f = e + 1, 
				cc = c + NMP, aa = cc + NP, bb = aa + 1, dd = cc + 1, ee = a + NMP, ff = ee + 1,
				ccc = cc + NMP, aaa = ccc + NP, bbb = aaa + 1, ddd = ccc + 1, eee = aa + NMP, fff = eee + 1)
			{
			TransferOperatorsInterface::PortraitWidth12_1(0.25, 0.125, 0.5, 0.25, 0.25, 0.125, c, d, e, f, cc, dd, ee, ff, ccc, ddd, eee, fff,
			num, _Cell, _Portrait, NR);
			c += NP; d = c + 1; e += NP; f = e + 1; cc += NP; dd = cc + 1; ee += NP; ff = ee + 1; ccc += NP; ddd = ccc + 1; eee += NP; fff = eee + 1;
			for (int j = 1; j < M2; j++, a += NP, b = a + 1, e += NP, f = e + 1, c += NP, d = c + 1, aa += NP, bb = aa + 1, aaa += NP, bbb = aaa + 1, 
				ee += NP, ff = ee + 1, eee += NP, fff = eee + 1, cc += NP, dd = cc + 1, ccc += NP, ddd = ccc + 1)
			{
				for (int i = 0; i < NR; i++, a += 2,  b += 2,  c += 2,  d += 2,  e += 2,  f += 2, 
				aa += 2,  bb += 2,  cc += 2,  dd += 2,  ee += 2,  ff += 2, 
				aaa += 2,  bbb += 2,  ccc += 2,  ddd += 2,  eee += 2,  fff += 2)
				{
					SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, c, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, d, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, e, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, f, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, aa, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, bb, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, cc, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, dd, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, ee, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, ff, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, aaa, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, bbb, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, ccc, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, ddd, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, eee, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, fff, 0.125);
				}
			}
			TransferOperatorsInterface::PortraitWidth12_1(0.125, 0.25, 0.25, 0.5, 0.125, 0.25, a, b, c, d, aa, bb, cc, dd, aaa, bbb, ccc, ddd,
			num, _Cell, _Portrait, NR);
		}
		TransferOperatorsInterface::PortraitWidth8_1(0.25, 0.125, 0.5, 0.25, c, d, e, f, cc, dd, ee, ff, num, _Cell, _Portrait, NR);
		c += NP; d = c + 1; e += NP; f = e + 1; cc += NP; dd = cc + 1; ee += NP; ff = ee + 1;
			for (int i = 1; i < M2; i++, a += NP, b = a + 1, e += NP, f = e + 1, c += NP, d = c + 1, aa += NP, bb = aa + 1, ee += NP, ff = ee + 1, cc += NP, dd = cc + 1)
			    TransferOperatorsInterface::PortraitWidth12_1(0.125, 0.25, 0.125, 0.25, 0.5, 0.25, a, b, c, d, e, f, aa, bb, cc, dd, ee, ff,
			num, _Cell, _Portrait, NR);
			TransferOperatorsInterface::PortraitWidth8_1(0.125, 0.25, 0.25, 0.5, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, NR);
			
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddOdd, numPeoo, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddOdd, portraitPeoo, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR)
		}
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddOdd, numReoo, m_R->Number, m_NR, M2, K2)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddOdd, portraitReoo, m_R->Portrait, m_R->Cell, M2, K2, m_NP, m_NR, m_NP * m_MP)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_ODD_H_
