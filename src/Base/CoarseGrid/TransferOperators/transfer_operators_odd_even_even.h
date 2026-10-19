/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_odd_even_even.h
* \author    Puzikova Valeria (vvp) 
* \date      30.05.2026 16:03:51
* \brief     transfer_operators_odd_even_even.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_EVEN_H_
#define _LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_EVEN_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsOddEvenEven

Данный файл содержит в себе реализацию класса TransferOperatorsOddEvenEven
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPoee, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numRoee, LS_STAG_R_PI, _Number, const int, nn, const int, MR, const int, KR);
	LS_STAG_DECLARATION_S_KERNEL(portraitPoee, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR, const int, NMR);
	LS_STAG_DECLARATION_S_KERNEL(portraitRoee, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, MR, const int, KR,
		const int, NP, const int, NR, const int, NMP);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным N и четным M, K

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов как по горизонтали и четное --- по вертикали, z
*/
	struct TransferOperatorsOddEvenEven : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsOddEvenEven():
		TransferOperatorsInterface(OddEvenEven, sizeof(TransferOperatorsOddEvenEven))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPoee, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
    {	
        int num = 0;
			TransferOperatorsInterface::NumWidth_2o(1, 2, 4, num, _Number, N2, M2);
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::NumWidth_2o(2, 4, 8, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth_2o(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth_2o(1, 2, 4, num, _Number, N2, M2);
    }
		LS_STAG_SINGLE_KERNEL(numRoee, LS_STAG_R_PI _Number, const int nn, const int MR, const int KR)
		{
			for (int k = 0, num = 0; k < KR; k++)
			{
			for (int j = 0; j < MR; j++)
			{
				_Number[num + 1] = _Number[num] + 8;
				num++;
				SparseMatrix::NumWidth_1(nn, 12, num, _Number);
				_Number[num + 1] = _Number[num] + 8;
				num++;
			}
		}
		}
		LS_STAG_KERNEL_INLINE void PortraitWidth8_baseO(const double cellA, const double cellB, const double cellC, const double cellD, 
			const double cellA2, const double cellB2, const double cellC2, const double cellD2, int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA2);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellA2);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC2);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellC2);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD2);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD2);
				a++;
				b++;
				c++;
				d++;
				aa++;
				bb++;
				cc++;
				dd++;
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth24_baseO(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int N2, const int M2, const int NR)
		{
			TransferOperatorsInterface::PortraitWidth2(c, d, num, _Cell, _Portrait, N2);
			for (int j = 0; j < M2; j++)
			{
			TransferOperatorsInterface::PortraitWidth4(0.75, 0.25,  0.375, 0.125, a, b, c, d, num, _Cell, _Portrait, N2);
				c = a; d = c + 1; a = c - NR; b = a + 1;
			TransferOperatorsInterface::PortraitWidth4(0.25, 0.75,  0.125, 0.375, a, b, c, d, num, _Cell, _Portrait, N2);
			}
			TransferOperatorsInterface::PortraitWidth2(a, b, num, _Cell, _Portrait, N2);
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8e(const double cellA, const double cellB, const double cellC, const double cellD,
			const double cellA2, const double cellB2, const double cellC2, const double cellD2,  int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA2);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellA2);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB2);// ooo 
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC2);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellC2);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD2);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD2);
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
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
			SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
			SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
			a++;
			b++;
			c++;
			d++;
			aa++;
			bb++;
			cc++;
			dd++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth48_baseO(const double cellA, const double cellB, const double cellC, const double cellD, 
			const double cellA2, const double cellB2, const double cellC2, const double cellD2, 
			const double cellA3, const double cellB3, const double cellC3, const double cellD3, int& a, int& b, int& c, int& d, int& aa, 
			int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2, const int _M2, const int NR)
		{
			TransferOperatorsInterface::PortraitWidth4(cellA, cellB, cellC, cellD, c, d, cc, dd, num, _Cell, _Portrait, _N2);
			for (int j = 0; j < _M2; j++)
			{
			TransferOperatorsOddEvenEven::PortraitWidth8e(cellA2, cellB2, cellC2, cellD2, cellA3, cellB3, cellC3, cellD3, a, b, c, d,
				aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
				c = a; d = c + 1; a = c - NR; b = a + 1; cc = aa; dd = cc + 1; aa = cc - NR; bb = aa + 1;
				TransferOperatorsOddEvenEven::PortraitWidth8e(cellB2, cellA2, cellD2, cellC2, cellB3, cellA3, cellD3, cellC3, a, b, c, d,
				aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
			}
			TransferOperatorsInterface::PortraitWidth4(cellA, cellB, cellC, cellD, a, b, aa, bb, num, _Cell, _Portrait, _N2);
		}
		LS_STAG_SINGLE_KERNEL(portraitPoee, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR, 
			const int NMR)
		{
			 int a = 0, b = 1, c = 0, d = 1;
        int aa = 0, bb = 1, cc = 0, dd = 1;
        int num = 0;
        TransferOperatorsOddEvenEven::PortraitWidth24_baseO(aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
        for (int k = 0; k < K2; ++k) {
			TransferOperatorsOddEvenEven::PortraitWidth48_baseO(0.75, 0.25, 0.375, 0.125, 0.5625, 0.1875, 0.1875, 0.0625,
				 0.28125, 0.09375, 0.09375, 0.03125, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
				 cc = c; dd = d; aa = c; bb = d; c -= NMR; a = c; b = a + 1; d = b;
				 TransferOperatorsOddEvenEven::PortraitWidth48_baseO(0.25, 0.75, 0.125, 0.375, 0.1875, 0.0625, 0.5625, 0.1875,
                 0.09375, 0.03125, 0.28125, 0.09375, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2, NR);
        }
        TransferOperatorsOddEvenEven::PortraitWidth24_baseO(a, b, c, d, num, _Cell, _Portrait, N2, M2, NR);
		}; 
		LS_STAG_SINGLE_KERNEL(portraitRoee, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int MR, const int KR, const int NP, const int NR, const int NMP)
		{
			for (int k = 0, num = 0, a = 1, aa = 0, aaa = 1, bb = NP, bbb = bb + 1, b = bbb, cc = NMP, ccc = cc + 1, c = ccc, dd = cc + NP, ddd = dd + 1, d = ddd; 
				k < KR; k++, b = ddd, aa = cc, aaa = ccc, a = aaa, bb = dd, bbb = ddd, cc += NMP, ccc = cc + 1, dd += NMP, ddd = dd + 1, c = ccc, d = ddd)
			{
			for (int i = 0; i < MR; i++, aa = bb - 1, aaa = aa + 1, bb = aa + NP, bbb = bb + 1, cc = aa + NMP, ccc = cc + 1, dd = cc + NP, ddd = dd + 1,
				a = aaa, b = bbb, c = ccc, d = ddd)
			{
				TransferOperatorsInterface::PortraitWidth8_1o(0.25, 0.125, 0.25, 0.125, aa, aaa, bb, bbb, cc, ccc, dd, ddd, num, _Cell, _Portrait, 1);
				TransferOperatorsInterface::PortraitWidth12_1oo(0.125, 0.25, 0.125, 0.125, 0.25, 0.125, a, aa, aaa, b, bb, bbb, c, cc, ccc, d, dd, ddd, num, _Cell, _Portrait, NR);
			TransferOperatorsInterface::PortraitWidth8_1o(0.125, 0.25, 0.125, 0.25, a, aa, b, bb, c, cc, d, dd, num, _Cell, _Portrait, 1);	
			}//for j
		}
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenEven, numPoee, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenEven, portraitPoee, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR, m_NR * m_MR)
		};
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenEven, numRoee, m_R->Number, N2 - 1, m_MR, m_KR)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsOddEvenEven, portraitRoee, m_R->Portrait, m_R->Cell, m_MR, m_KR, m_NP, N2 - 1, m_NP * m_MP)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_ODD_EVEN_EVEN_H_
