/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_even_odd_even.h
* \author    Puzikova Valeria (vvp) 
* \date      30.05.2026 19:57:20
* \brief     transfer_operators_even_odd_even.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_EVEN_H_
#define _LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_EVEN_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsEvenOddEven

Данный файл содержит в себе реализацию класса TransferOperatorsEvenOddEven
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPeoe, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(numReoe, LS_STAG_R_PI, _Number, const int, NR, const int, M2, const int, KR);
	LS_STAG_DECLARATION_S_KERNEL(portraitPeoe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2, const int, NR);
	LS_STAG_DECLARATION_S_KERNEL(portraitReoe, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, M2, const int, K2,
		const int, NP, const int, NR, const int, MR, const int, KR, const int, NMP);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным M и четным N, K

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов как по вертикали и четное --- по горизонтали, z
*/
	struct TransferOperatorsEvenOddEven : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsEvenOddEven():
		TransferOperatorsInterface(EvenOddEven, sizeof(TransferOperatorsEvenOddEven))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPeoe, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
		TransferOperatorsInterface::NumWidth4_o(1, 2, 4, num, _Number, N2, M2);
        for (int k = 0; k < K2; ++k) {
			TransferOperatorsInterface::NumWidth4_o(2, 4, 8, num, _Number, N2, M2);
			TransferOperatorsInterface::NumWidth4_o(2, 4, 8, num, _Number, N2, M2);
        }
		TransferOperatorsInterface::NumWidth4_o(1, 2, 4, num, _Number, N2, M2);// NumWidth4xy_o
		}
		LS_STAG_SINGLE_KERNEL(numReoe, LS_STAG_R_PI _Number, const int NR, const int M2, const int KR)
		{	
			for (int k = 0, num = 0; k < KR; k++)
			{
				SparseMatrix::NumWidth_1(NR, 8, num, _Number);
			for (int j = 1; j < M2; j++) SparseMatrix::NumWidth_1(NR, 12, num, _Number);
				SparseMatrix::NumWidth_1(NR, 8, num, _Number);
			}
		}
		LS_STAG_SINGLE_KERNEL(portraitPeoe, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2, const int NR)
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
			int NMR = (N2 + 1) * (M2 + 1);
			TransferOperatorsInterface::PortraitWidth2D4D_o(aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			for (int k = 0; k < K2; k++)
			{
			TransferOperatorsInterface::PortraitWidth4D8D_o(0.75, 0.25, 0.5625, 0.1875, 0.1875, 0.0625, 0.375, 0.125, 
				0.28125, 0.09375, 0.09375, 0.03125, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
				aa = a; bb = b; cc = c; dd = d; a -= NMR; b = a + 1; c = a + NR; d = c + 1;
			TransferOperatorsInterface::PortraitWidth4D8D_o(0.25, 0.75, 0.1875, 0.0625, 0.5625, 0.1875, 0.125, 0.375, 
				0.09375, 0.03125, 0.28125, 0.09375, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2, M2);
			}
			TransferOperatorsInterface::PortraitWidth2D4D_o(a, b, c, d, num, _Cell, _Portrait, N2, M2);
		};
		LS_STAG_SINGLE_KERNEL(portraitReoe, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int M2, const int K2, const int NP, const int NR, 
			const int MR, const int KR, const int NMP)
		{
			for (int k = 0, num = 0, aa = 0, bb = 1, aaa = NP, bbb = aaa + 1, cc = NMP, dd = cc + 1, ccc = cc + NP, 
				ddd = ccc + 1, a = aaa, b = bbb, c = ccc, d = ddd; k < KR; k++, aa = cc, bb = dd, aaa = aa + NP, bbb = aaa + 1, 
				cc += NMP, dd = cc + 1, ccc = cc + NP, ddd = ccc + 1, a = aaa, b = bbb, c = ccc, d = ddd)
			{
				TransferOperatorsInterface::PortraitWidth8_1(0.25, 0.125, 0.25, 0.125, aa, bb, aaa, bbb, cc, dd, ccc, ddd, num, _Cell, _Portrait, NR);
				aa = aaa; bb = bbb; cc = ccc; dd = ddd; aaa += NP; bbb = aaa + 1; ccc += NP; ddd = ccc + 1;
				for (int j = 1; j < M2; j++, a = aa, b = bb, c = cc, d = dd, aa = aaa, bb = bbb, cc = ccc, dd = ddd, 
					aaa += NP, bbb = aaa + 1, ccc += NP, ddd = ccc + 1) 
					TransferOperatorsInterface::PortraitWidth12_1(0.125, 0.25, 0.125, 0.125, 
					0.25, 0.125, a, b, aa, bb, aaa, bbb, c, d, cc, dd, ccc, ddd, num, _Cell, _Portrait, NR);
				TransferOperatorsInterface::PortraitWidth8_1(0.125, 0.25, 0.125, 0.25, a, b, aa, bb, c, d, cc, dd, num, _Cell, _Portrait, NR);
			}
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddEven, numPeoe, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddEven, portraitPeoe, m_P->Portrait, m_P->Cell, N2, M2, K2, m_NR)
		}
		void buildNumR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddEven, numReoe, m_R->Number, m_NR, M2, m_KR)
		};
		void buildPortraitR()
		{
			LS_STAG_SINGLE(TransferOperatorsEvenOddEven, portraitReoe, m_R->Portrait, m_R->Cell, M2, K2, m_NP, m_NR, m_MR, m_KR, m_NP * m_MP)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_EVEN_ODD_EVEN_H_
