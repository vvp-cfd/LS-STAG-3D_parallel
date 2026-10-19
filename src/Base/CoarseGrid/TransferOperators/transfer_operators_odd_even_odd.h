/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_odd_even_odd.h
* \author    Puzikova Valeria (vvp) 
* \date      24.05.2026 11:42:51
* \brief     transfer_operators_odd_even_odd.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_ODD_H_
#define _LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_ODD_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsOddEvenOdd

Данный файл содержит в себе реализацию класса TransferOperatorsOddEvenOdd
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPoeo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numRoeo, LS_STAG_R_PI, _Number, const int, N2, const int, MR, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitPoeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR);
	LS_STAG_DECLARATION_S_KERNEL(portraitRoeo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR,
		const int, NP, const int, K2, const int, nn);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным N, K и четным M

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов как горизонтали, z и четное --- по вертикали
*/
	struct TransferOperatorsOddEvenOdd : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsOddEvenOdd():
		TransferOperatorsInterface(OddEvenOdd, sizeof(TransferOperatorsOddEvenOdd))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPoeo, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::NumWidth_2o(1, 2, 4, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth_2o(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth_2o(1, 2, 4, num, _Number, N2, M2);
		}
		LS_STAG_SINGLE_KERNEL(numRoeo, LS_STAG_R_PI _Number, const int N2, const int MR, const int K2)
		{
			int num = 0;
			TransferOperatorsInterface::NumWidth1_oe(8, 12, num, _Number, N2, MR);
			for (int k = 1; k < K2; k++)
			   TransferOperatorsInterface::NumWidth1_oe(12, 18, num, _Number, N2, MR);
			TransferOperatorsInterface::NumWidth1_oe(8, 12, num, _Number, N2, MR);
		}
		LS_STAG_SINGLE_KERNEL(portraitPoeo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR)
		{
			int a = 0;
			int b = 1;
			int c = 0;
			int d = c + 1;
			int aa = 0;
			int bb = 1;
			int cc = 0;
			int dd = cc + 1;
			int num = 0;
			for (int k = 0; k < K2; k++, cc = c, dd = d, aa = a, bb = b)
			{
				TransferOperatorsInterface::PortraitWidth24(aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
			TransferOperatorsInterface::PortraitWidth4(0.5, 0.5, 0.25, 0.25, c, d, cc, dd, num, _Cell, _Portrait, N2);
			for (int j = 0; j < M2; j++)
			{
			TransferOperatorsInterface::PortraitWidth8(0.375, 0.125, 0.1875, 0.0625, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2);
				a -= NR;
				b = a + 1;
				c -= NR;
				d = c + 1;
				aa -= NR;
				bb = aa + 1;
				cc -= NR;
				dd = cc + 1;
			TransferOperatorsInterface::PortraitWidth8(0.125, 0.375, 0.0625, 0.1875, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2);
			}
			TransferOperatorsInterface::PortraitWidth4(0.5, 0.5, 0.25, 0.25, a, b, aa, bb, num, _Cell, _Portrait, N2);
		}
			TransferOperatorsInterface::PortraitWidth24(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
		}; 
		LS_STAG_SINGLE_KERNEL(portraitRoeo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int MR, const int NP, const int K2, const int nn)
		{
			int NMP = NP * 2 * MR;
			int c = NMP;
			int e = c + 1;
			int a = e;
			int b = a + NP;
			int d = c + NP;
			int f = d + 1;
			int aa = 1;
			int cc = 0;
			int ee = 1;
			int bb = aa + NP;
			int dd = NP;
			int ff = dd + 1;
			int ccc = NMP;
			int eee = ccc + 1;
			int aaa = eee;
			int bbb = aaa + NP;
			int ddd = ccc + NP;
			int fff = ddd + 1;
			int num = 0;
			for (int j = 0; j < MR; j++)
			{
				TransferOperatorsInterface::PortraitWidth8_1o(0.5, 0.25, 0.25, 0.125, cc, ee, dd, ff, ccc, eee, ddd, fff, 
					num, _Cell, _Portrait, 1);
			    TransferOperatorsInterface::PortraitWidth12_tb1o(0.25, 0.5, 0.125, 0.25, aa, cc, ee, bb, dd, ff, aaa, ccc, eee, bbb, ddd, fff,
			num, _Cell, _Portrait, nn);
			TransferOperatorsInterface::PortraitWidth8_1o(0.25, 0.5, 0.125, 0.25, aa, cc, bb, dd, aaa, ccc, bbb, ddd, 
				num, _Cell, _Portrait, 1);
			cc = bb; ee = bb + 1; aa = ee;  bb = aa + NP; dd = cc + NP; ff = dd + 1;
			ccc = bbb; eee = bbb + 1; aaa = eee;  bbb = aaa + NP; ddd = ccc + NP; fff = ddd + 1;
			}
			cc = ccc; ee = eee; aa = ee; bb = aa + NP; dd = ddd; ff = fff; ccc += NMP; eee += NMP; aaa = eee; bbb = aaa + NP; ddd += NMP; fff += NMP;
			for (int k = 1; k < K2; k++, c = cc, e = ee, a = e, b = a + NP, d = dd, f = ff, cc = ccc, ee = eee, aa = ee, bb = aa + NP, dd = ddd, ff = fff,
				ccc += NMP, eee += NMP, aaa = eee, bbb = aaa + NP, ddd += NMP, fff += NMP)
			{
			for (int j = 0; j < MR; j++)
			{
			TransferOperatorsInterface::PortraitWidth12_1o(0.25, 0.125, 0.5, 0.25, 0.25, 0.125, c, e, d, f, cc, ee, dd, ff, 
				ccc, eee, ddd, fff, num, _Cell, _Portrait, 1);
				for (int i = 0; i < nn; i++, a += 2,  b += 2,  c += 2,  d += 2,  e += 2,  f += 2, 
				aa += 2,  bb += 2,  cc += 2,  dd += 2,  ee += 2,  ff += 2, 
				aaa += 2,  bbb += 2,  ccc += 2,  ddd += 2,  eee += 2,  fff += 2)
				{
					SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, c, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, e, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, d, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, f, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, aa, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, cc, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, ee, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, bb, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, dd, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, ff, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, aaa, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, ccc, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, eee, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, bbb, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, ddd, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, fff, 0.125);
				}
			TransferOperatorsInterface::PortraitWidth12_1o(0.125, 0.25, 0.25, 0.5, 0.125, 0.25, a, c, b, d, aa, cc, bb, dd, 
				aaa, ccc, bbb, ddd, num, _Cell, _Portrait, 1);
			ccc = bbb; eee = bbb + 1; aaa = eee;  bbb = aaa + NP; ddd = ccc + NP; fff = ddd + 1;
				cc = bb; ee = bb + 1; aa = ee;  bb = aa + NP; dd = cc + NP; ff = dd + 1;
				c = b; e = b + 1; a = e;  b = a + NP; d = c + NP; f = d + 1;
			}
		}
			for (int j = 0; j < MR; j++)
			{
			TransferOperatorsInterface::PortraitWidth8_1o(0.25, 0.125, 0.5, 0.25, c, e, d, f, cc, ee, dd, ff, 
					num, _Cell, _Portrait, 1);
			    TransferOperatorsInterface::PortraitWidth12_tb1o(0.125, 0.25, 0.25, 0.5, a, c, e, b, d, f, aa, cc, ee, bb, dd, ff,
			num, _Cell, _Portrait, nn);
			TransferOperatorsInterface::PortraitWidth8_1o(0.125, 0.25, 0.25, 0.5, a, c, b, d, aa, cc, bb, dd, 
				num, _Cell, _Portrait, 1);
			cc = bb; ee = bb + 1; aa = ee;  bb = aa + NP; dd = cc + NP; ff = dd + 1;
			c = b; e = b + 1; a = e;  b = a + NP; d = c + NP; f = d + 1;
			}
			
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenOdd, numPoeo, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenOdd, portraitPoeo, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR)
		};
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenOdd, numRoeo, m_R->Number, N2, m_MR, K2)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenOdd, portraitRoeo, m_R->Portrait, m_R->Cell, m_MR, m_NP, K2, N2 - 1)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_ODD_H_
