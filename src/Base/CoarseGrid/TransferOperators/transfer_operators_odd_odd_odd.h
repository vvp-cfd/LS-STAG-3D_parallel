/*
* \copyright (C) 2026 Puzikova V.V.
* \file      transfer_operators_odd_odd_odd.h
* \author    Puzikova Valeria (vvp) 
* \date      22.05.2026 12:17:41
* \brief     transfer_operators_odd_odd_odd.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_ODD_ODD_ODD_H_
#define _LS_STAG_TRANSFER_OPERATORS_ODD_ODD_ODD_H_

// ================================================================== INCLUDES
#include "transfer_operators_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsOddOddOdd

Данный файл содержит в себе реализацию класса TransferOperatorsOddOddOdd
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numPooo, LS_STAG_R_PI, _Number, const int, N2, const int, M2, const int, K2);
	LS_STAG_DECLARATION_S_KERNEL(portraitPooo, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell, const int, N2,
		const int, M2, const int, K2);
/*!
\brief Класс, реализующий переход с сеточного уровня с нечетным числом узлов

Данный класс реализует построение дискретных аналогов операторов перехода с сеточного уровня многосеточного метода, 
который содержит нечетное число узлов как по вертикали, так и по горизонтали, z
*/
	struct TransferOperatorsOddOddOdd : TransferOperatorsInterface
	{
	  public:
		TransferOperatorsOddOddOdd():
		TransferOperatorsInterface(OddOddOdd, sizeof(TransferOperatorsOddOddOdd))
		{
			;
		}
		LS_STAG_KERNEL_INLINE void PortraitWidth24(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2, 
			const int _M2, const int _N)
		{
			for (int j = 0; j < _M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth2(a, b, num, _Cell, _Portrait, _N2);
				a -= _N;
				b = a + 1;
				TransferOperatorsInterface::PortraitWidth4(a, b, c, d, num, _Cell, _Portrait, _N2);
			}
			TransferOperatorsInterface::PortraitWidth2(a, b, num, _Cell, _Portrait, _N2);
		};
		LS_STAG_SINGLE_KERNEL(numPooo, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			int num = 0;
			TransferOperatorsInterface::numPoo(num, _Number, N2, M2, K2);
		}
		LS_STAG_SINGLE_KERNEL(portraitPooo, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell, const int N2, const int M2, const int K2)
		{
			const int N = N2 + 1;
			const int NM = N * (M2 + 1);
			int a = 0;
			int b = 1;
			int c = N;
			int d = c + 1;
			int aa = NM;
			int bb = aa + 1;
			int cc = c + NM;
			int dd = cc + 1;
			int num = 0;
			for (int k = 0; k < K2; k++)
			{
				TransferOperatorsOddOddOdd::PortraitWidth24(a, b, c, d, num, _Cell, _Portrait, N2, M2, N);
			a -= NM;
			b = a + 1;
			c = a + N;
			d = c + 1;
			for (int j = 0; j < M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth4(a, b, aa, bb, num, _Cell, _Portrait, N2);
				a -= N;
				b = a + 1;
				aa -= N;
				bb = aa + 1;
				TransferOperatorsInterface::PortraitWidth8(a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, N2);
			}
				TransferOperatorsInterface::PortraitWidth4(a, b, aa, bb, num, _Cell, _Portrait, N2);
			c += N;
			d = c + 1;
			cc += N;
			dd = cc + 1;
		}
				TransferOperatorsOddOddOdd::PortraitWidth24(a, b, c, d, num, _Cell, _Portrait, N2, M2, N);
		};
	  private:
		void buildNumP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddOdd, numPooo, m_P->Number, N2, M2, K2)
		};
		void buildPortraitP()
		{
			LS_STAG_SINGLE(TransferOperatorsOddOddOdd, portraitPooo, m_P->Portrait, m_P->Cell, N2, M2, K2)
		};
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_ODD_ODD_ODD_H_
