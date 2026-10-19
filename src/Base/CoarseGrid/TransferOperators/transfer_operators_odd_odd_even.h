/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_odd_odd_even.h
* \author    Puzikova Valeria (vvp) 
* \date      29.05.2026 10:23:41
* \brief     transfer_operators_odd_odd_even.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_ODD_ODD_EVEN_H_
#define _LS_STAG_TRANSFER_OPERATORS_ODD_ODD_EVEN_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsOddOddEven

Данный файл содержит в себе реализацию класса TransferOperatorsOddOddEven
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPooe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitPooe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numRooe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitRooe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2,
		const int, NP, const int, NMP, const int, KR, const int, nn);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным числом узлов

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов как по вертикали, так и по горизонтали, но четное - по z
*/
	struct TransferOperatorsOddOddEven : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsOddOddEven():
		TransferOperatorsInterface(OddOddEven, sizeof(TransferOperatorsOddOddEven))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPooe, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			TransferOperatorsInterface::NumWidth24(1, 2, 4, num, _Number, N2, M2);
			for (int k = 0; k < K2; k++)
			{
				TransferOperatorsInterface::NumWidth24(2, 4, 8, num, _Number, N2, M2);
				TransferOperatorsInterface::NumWidth24(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth24(1, 2, 4, num, _Number, N2, M2);
		}
		LS_STAG_SINGLE_KERNEL(numRooe, LS_STAG_R_PI _Number, const int N2, const int MR, const int K2)
		{
			for (int k = 0, num = 0; k < K2; k++)// KR
			{
			TransferOperatorsInterface::NumWidth1_oe(8, 12, num, _Number, N2, 1);
			   TransferOperatorsInterface::NumWidth1_oe(12, 18, num, _Number, N2, MR);
			TransferOperatorsInterface::NumWidth1_oe(8, 12, num, _Number, N2, 1);
			}
		}
		LS_STAG_KERNEL_INLINE void PortraitWidth24_base(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, 
			const int _N2, const int _M2)
		{
			for (int j = 0; j < _M2; j++)
			{
			TransferOperatorsInterface::PortraitWidth2(c, d, num, _Cell, _Portrait, _N2);
			TransferOperatorsInterface::PortraitWidth4(0.5, 0.5, 0.25, 0.25, a, b, c, d, num, _Cell, _Portrait, _N2);
			c = a; d = b;
			}
			TransferOperatorsInterface::PortraitWidth2(c, d, num, _Cell, _Portrait, _N2);
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8e(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellC);// ooo 
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
				a++;
				b++;
				c++;
				d++;
				aa++;
				bb++;
				cc++;
				dd++;
			}
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellB);
			SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellB);
			a++;
			b++;
			c++;
			d++;
			aa++;
			bb++;
			cc++;
			dd++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth48_base(const double cellA, const double cellB, const double cellC, 
			const double cellD, const double cellC2, const double cellD2, int& a, int& b, int& c, int& d, int& aa, 
			int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2, const int _M2)
		{
			for (int j = 0, _NR = _N2 + 1; j < _M2; j++)
			{
				TransferOperatorsOddOddEven::PortraitWidth4(cellA, cellB, cellC, cellD, c, d, cc, dd, num, _Cell, _Portrait, _N2);
			aa = cc - _NR; bb = aa + 1;
			TransferOperatorsOddOddEven::PortraitWidth8e(cellC, cellD, cellC2, cellD2, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
			c = a; d = b; cc = aa; dd = bb;
			}
			TransferOperatorsOddOddEven::PortraitWidth4(cellA, cellB, cellC, cellD, c, d, cc, dd, num, _Cell, _Portrait, _N2);
		}
		LS_STAG_SINGLE_KERNEL(portraitPooe, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2)
		{	
			int a = 0;
			int b = 1;
			int c = 0;
			int d = 1;
			int aa = 0;
			int bb = 1;
			int cc = 0;
			int dd = 1;
			int num = 0;
			TransferOperatorsOddOddEven::PortraitWidth24_base(aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			for (int k = 0, _NMR = (N2 + 1) * (M2 + 1); k < K2; k++)
			{
			TransferOperatorsOddOddEven::PortraitWidth48_base(0.75, 0.25, 0.375, 0.125, 0.1875, 0.0625, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			aa = a; bb = b; cc = c; dd = d;  a = c - _NMR; b = a + 1; c = a; d = c + 1;
			TransferOperatorsOddOddEven::PortraitWidth48_base(0.25, 0.75, 0.125, 0.375, 0.0625, 0.1875, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			a = c; b = a + 1;
		}
			TransferOperatorsOddOddEven::PortraitWidth24_base(a, b, c, d, num, _Cell, _Portrait, N2, M2);
		};
		LS_STAG_SINGLE_KERNEL(portraitRooe, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int M2, const int NP, const int NMP, const int KR, 
			const int nn)
		{
			int a = NP + 1;
			int c = NP;
			int e = c + 1;
			int aa = 1;
			int cc = 0;
			int ee = 1;
			int aaa = a;
			int ccc = c;
			int eee = ccc + 1;
			int b = a + NMP;
			int d = c + NMP;
			int f = d + 1;
			int bb = NMP + 1;
			int dd = NMP;
			int ff = dd + 1;
			int bbb = bb + NP;
			int ddd = dd + NP;
			int fff = ddd + 1;
			int num = 0;
			for (int k = 0; k < KR; k++)
			{
				TransferOperatorsInterface::PortraitWidth8_1oo(0.5, 0.25, 0.25, 0.125, cc, ee, ccc, eee, dd, ff, ddd, fff, 
					num, _Cell, _Portrait, 1);
			    TransferOperatorsInterface::PortraitWidth12_1oo(0.25, 0.5, 0.25, 0.125, 0.25, 0.125, aa, cc, ee, aaa, ccc, eee, bb, dd, ff, bbb, ddd, fff,
			num, _Cell, _Portrait, nn);
			TransferOperatorsInterface::PortraitWidth8_1oo(0.25, 0.5, 0.125, 0.25, aa, cc, aaa, ccc, bb, dd, bbb, ddd, 
				num, _Cell, _Portrait, 1);
				cc = aaa; ee = ccc; aa = ccc; bb = ddd; ff = ddd; dd = fff;
				aaa = aa + NP; bbb = bb + NP; ccc = cc + NP; ddd = dd + NP; eee = ee + NP; fff = ff + NP;
			for (int j = 1; j < M2; j++)
			{
			TransferOperatorsInterface::PortraitWidth12_1oo(0.25, 0.125, 0.5, 0.25, 0.25, 0.125, c, e, cc, ee, ccc, eee, 
				d, f, dd, ff, ddd, fff, num, _Cell, _Portrait, 1);
				for (int i = 0; i < nn; i++, a += 2,  b += 2,  c += 2,  d += 2,  e += 2,  f += 2, 
				aa += 2,  bb += 2,  cc += 2,  dd += 2,  ee += 2,  ff += 2, 
				aaa += 2,  bbb += 2,  ccc += 2,  ddd += 2,  eee += 2,  fff += 2)
				{
					SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, c, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, e, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, aa, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, cc, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, ee, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, aaa, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, ccc, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, eee, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, d, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, f, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, bb, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, dd, 0.5);
					SparseMatrix::setPC(num, _Cell, _Portrait, ff, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, bbb, 0.125);
					SparseMatrix::setPC(num, _Cell, _Portrait, ddd, 0.25);
					SparseMatrix::setPC(num, _Cell, _Portrait, fff, 0.125);
				}
			TransferOperatorsInterface::PortraitWidth12_1oo(0.125, 0.25, 0.25, 0.5, 0.125, 0.25, a, c, aa, cc, aaa, ccc, 
				b, d, bb, dd, bbb, ddd, num, _Cell, _Portrait, 1);
				c = aa; e = cc; a = cc; b = dd; f = dd; d = ff;
				cc = aaa; ee = ccc; aa = ccc; bb = ddd; ff = ddd; dd = fff;
				aaa = aa + NP; bbb = bb + NP; ccc = cc + NP; ddd = dd + NP; eee = ee + NP; fff = ff + NP; 
			}
				TransferOperatorsInterface::PortraitWidth8_1oo(0.25, 0.125, 0.5, 0.25, c, e, cc, ee, d, f, dd, ff, 
					num, _Cell, _Portrait, 1);
			    TransferOperatorsInterface::PortraitWidth12_1oo(0.125, 0.25, 0.125, 0.25, 0.5, 0.25, a, c, e, aa, cc, ee, b, d, f, bb, dd, ff,
			num, _Cell, _Portrait, nn);
			TransferOperatorsInterface::PortraitWidth8_1oo(0.125, 0.25, 0.25, 0.5, a, c, aa, cc, b, d, bb, dd, 
				num, _Cell, _Portrait, 1);
			cc = ff; ee = cc + 1; aa = ee; a = aa + NP; c = cc + NP; e = c + 1; aaa = a; ccc = c; eee = e; 
			b = a + NMP; d = c + NMP;  f = d + 1;  dd = cc + NMP; ff = dd + 1; bb = ff;
				 bbb = bb + NP; ddd = dd + NP; fff = ddd + 1;
		}
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddEven, numPooe, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddEven, portraitPooe, m_P->Portrait, m_P->Cell, N2, M2, K2)
		};
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddEven, numRooe, m_R->Number, N2, M2 - 1, m_KR)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddEven, portraitRooe, m_R->Portrait, m_R->Cell, M2, m_NP, m_NP * m_MP, m_KR, N2 - 1)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_ODD_ODD_EVEN_H_
