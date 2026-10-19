/*
* \copyright (C) 2013 Puzikova V.V.
* \file      transfer_operators_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      05.10.2013 1:41:04
* \brief     transfer_operators_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSFER_OPERATORS_INTERFACE_H_
#define _LS_STAG_TRANSFER_OPERATORS_INTERFACE_H_

// ================================================================== INCLUDES
#include "../../SparseMatrix/transposed_sparse_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransferOperatorsInterface

Данный файл содержит в себе реализацию класса TransferOperatorsInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий переход между сеточными уровнями

Данный базовый класс реализует интерфейс для построения дискретных аналогов операторов для перехода между 
соседними сеточными уровнями в многосеточном методе, а также для работы с ними
*/
	struct TransferOperatorsInterface
	{
	  public:	
		  /// Тип сетки
		  enum OddEvenGrid3D
		{
			UnknownOddEven = -1,///< Неизвестен
			OddOddOdd,   ///<  Нечетное число узлов и по горизонтали, и по вертикали, и по z      
			OddOddEven,  ///<  Нечетное число узлов по горизонтали и по вертикали, четное по z     
			OddEvenEven, ///<  Нечетное число узлов по горизонтали, четное по вертикали и по z         
			OddEvenOdd,  ///<  Нечетное число узлов по горизонтали и z, четное --- по вертикали
			EvenOddEven, ///<  Четное число узлов по горизонтали и по z, нечетное --- по вертикали
			EvenOddOdd,  ///<  Четное число узлов по горизонтали, нечетное --- по вертикали и по z
			EvenEvenOdd, ///<  Четное число узлов по горизонтали и по вертикали, нечетное по z
			EvenEvenEven ///<  Четное число узлов и по горизонтали, и по вертикали, и по z
		};
		explicit TransferOperatorsInterface(OddEvenGrid3D _id = UnknownOddEven, int _size = 0):
		id(_id),
		size(_size)
		{
			;
		}
		LS_STAG_KERNEL_INLINE void PortraitWidth2_base(int& a, int& b, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 1.0);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.5);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.5);
				a++;
				b++;
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth2(int& a, int& b, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			TransferOperatorsInterface::PortraitWidth2_base(a, b, num, _Cell, _Portrait, _N2);
			SparseMatrix::setPC(num, _Cell, _Portrait, a, 1.0);
			a++;
			b++;
		};
	LS_STAG_KERNEL_INLINE void NumWidth_2o(const int m_a, const int m_b, const int m_c, int& num, int* _Number, const int _N2, const int _MR)
    {
		SparseMatrix::NumWidth2(m_a, m_b, num, _Number, _N2);//_2(_N2, m_a, m_b, num, _Number);
			for (int j = 0; j < _MR; j++)
			{
				SparseMatrix::NumWidth2(m_b, m_c, num, _Number, _N2);//_2(_N2, m_b, m_c, num, _Number);
				SparseMatrix::NumWidth2(m_b, m_c, num, _Number, _N2);//_2(_N2, m_b, m_c, num, _Number);
			}
			SparseMatrix::NumWidth2(m_a, m_b, num, _Number, _N2);//_2(_N2, m_a, m_b, num, _Number);
    }		
    LS_STAG_KERNEL_INLINE void NumWidth4xy_o(const int m_b, const int m_c, int& num, int* _Number, const int _N2, const int _M2)
    {
			for (int j = 0; j < _M2; j++)
			{
				SparseMatrix::NumWidth4(m_b, m_c, num, _Number, _N2);
				SparseMatrix::NumWidth4(m_b, m_c, num, _Number, _N2);
			}
    }
	LS_STAG_KERNEL_INLINE void NumWidth4_o(const int m_a, const int m_b, const int m_c, int& num, int* _Number, const int _N2, const int _M2)
    {
			for (int j = 0; j < _M2; j++)
			{
				SparseMatrix::NumWidth4(m_a, m_b, num, _Number, _N2);
				SparseMatrix::NumWidth4(m_b, m_c, num, _Number, _N2);
			}
			SparseMatrix::NumWidth4(m_a, m_b, num, _Number, _N2);// eoo
    }	
    LS_STAG_KERNEL_INLINE void NumWidth4xy(const int m_a, const int m_b, const int m_c, int& num, int* _Number, const int _N2, const int _M2)
    {
        SparseMatrix::NumWidth4(m_a, m_b, num, _Number, _N2);
		TransferOperatorsInterface::NumWidth4xy_o(m_b, m_c, num, _Number, _N2, _M2);
		SparseMatrix::NumWidth4(m_a, m_b, num, _Number, _N2);
    }
	LS_STAG_KERNEL_INLINE void NumWidth1_oe(const int m_a, const int m_b, int& num, int* _Number, const int _N2, const int _MR)
    {
		for (int j = 0; j < _MR; j++) {
			_Number[num + 1] = _Number[num] + m_a;
            num++;
			for (int i = 1; i < _N2; i++)
			{
				_Number[num + 1] = _Number[num] + m_b;
            num++;
			}
			_Number[num + 1] = _Number[num] + m_a;
            num++;
        }
    }
    LS_STAG_KERNEL_INLINE void NumWidth24(const int m_a, const int m_b, const int m_c, int& num, int* _Number, const int _N2, const int _M2)
    {
        for (int j = 0; j < _M2; j++)
			{
				SparseMatrix::NumWidth2(m_a, m_b, num, _Number, _N2);
				SparseMatrix::NumWidth2(m_b, m_c, num, _Number, _N2);
			}
			SparseMatrix::NumWidth2(m_a, m_b, num, _Number, _N2);
    }
		LS_STAG_KERNEL_INLINE void numPoo(int& num, LS_STAG_R_PI _Number, const int N2, const int M2, const int K2)
		{
			for (int k = 0; k < K2; k++)
			{
				TransferOperatorsInterface::NumWidth24(1, 2, 4, num, _Number, N2, M2);
				TransferOperatorsInterface::NumWidth24(2, 4, 8, num, _Number, N2, M2);
			}
			TransferOperatorsInterface::NumWidth24(1, 2, 4, num, _Number, N2, M2);
		}
		LS_STAG_KERNEL_INLINE void PortraitWidth4_base(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++, a++, b++, c++, d++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellD);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth4(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			TransferOperatorsInterface::PortraitWidth4_base(0.5, 0.5, 0.25, 0.25, a, b, c, d, num, _Cell, _Portrait, _N2);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.5);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, 0.5);
			a++;
			b++;
				c++;
				d++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth4(const double cellA, const double cellB, const double cellC, const double cellD, 
			int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			TransferOperatorsInterface::PortraitWidth4_base(cellA, cellB, cellC, cellD, a, b, c, d, num, _Cell, _Portrait, _N2);
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
			a++;
			b++;
			c++;
			d++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth24(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int N2,
			const int M2, const int NR)
		{
			TransferOperatorsInterface::PortraitWidth2(c, d, num, _Cell, _Portrait, N2);
			for (int j = 0; j < M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth4(0.75, 0.25, 0.375, 0.125, a, b, c, d, num, _Cell, _Portrait, N2);
				a -= NR;
				b = a + 1;
				c -= NR;
				d = c + 1;
				TransferOperatorsInterface::PortraitWidth4(0.25, 0.75, 0.125, 0.375, a, b, c, d, num, _Cell, _Portrait, N2);
			}
			TransferOperatorsInterface::PortraitWidth2(a, b, num, _Cell, _Portrait, N2);
		}; 
		LS_STAG_KERNEL_INLINE void PortraitWidth8_base(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellD);// ooo 
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellC);
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
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8(int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2)
		{
			TransferOperatorsInterface::PortraitWidth8_base(0.25, 0.25, 0.125, 0.125, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.25);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, 0.25);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, 0.25);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, 0.25);
			a++;
			b++;
				c++;
				d++;
				aa++;
				bb++;
				cc++;
				dd++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8(const double cellA, const double cellB, const double cellC, const double cellD, 
			int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell,
			LS_STAG_R_PI _Portrait, const int _N2)
		{
			TransferOperatorsInterface::PortraitWidth8_base(cellA, cellB, cellC, cellD, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
			SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellA);
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
		LS_STAG_KERNEL_INLINE void PortraitWidth2D(int& a, int& b, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _N2)
		{
			SparseMatrix::setPC(num, _Cell, _Portrait, a, 1.0);
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.75);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.25);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, 0.25);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, 0.75);
				a++;
				b++;
			}
			SparseMatrix::setPC(num, _Cell, _Portrait, a, 1.0);
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth4D(const double cellA, const double cellA2, const double cellB, const double cellB2, 
			const double cellB3, const double cellB4, int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2)
		{
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA2);
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB3);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB4);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB4);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB3);
				a++;
				b++;
				c++;
				d++;
			}
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA2);
			a++;
			b++;
			c++;
			d++;
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8D(const double cellA, const double cellA2, const double cellA3, const double cellA4, 
			const double cellB, const double cellB2, const double cellB3, const double cellB4, 
			const double cellB5, const double cellB6, const double cellB7, const double cellB8, 
			int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2)
		{// for eee
			SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA2);
			SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellA3);
			SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellA4);
			for (int i = 0; i < _N2; i++)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB3);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB4);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellB5);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellB6);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellB7);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellB8);
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellB2);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB4);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB3);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellB6);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellB5);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellB8);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellB7);
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
			SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA2);
			SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellA3);
			SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellA4);
			a++;
			b++;
			c++;
			d++;
				aa++;
				bb++;
				cc++;
				dd++;
		};
    LS_STAG_KERNEL_INLINE void PortraitWidth2D4D(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2, const int _M2, const int NR)
    {
        TransferOperatorsInterface::PortraitWidth2D(a, b, num, _Cell, _Portrait, _N2);
			a -= _N2;
			b = a + 1;
			for (int j = 0; j < _M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth4D(0.75, 0.25, 0.5625, 0.1875, 0.1875, 0.0625, a, b, c, d, num,
					_Cell, _Portrait, _N2);
				a -= NR;
				b = a + 1;
				c -= NR;
				d = c + 1;
				TransferOperatorsInterface::PortraitWidth4D(0.25, 0.75, 0.1875, 0.0625, 0.5625, 0.1875, a, b, c, d, num,
					_Cell, _Portrait, _N2);
			}
			TransferOperatorsInterface::PortraitWidth2D(a, b, num, _Cell, _Portrait, _N2);
    }
    LS_STAG_KERNEL_INLINE void PortraitWidth2D4D_o(int& a, int& b, int& c, int& d, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2, const int _M2)
    {
			for (int j = 0; j < _M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth2D(a, b, num, _Cell, _Portrait, _N2);
				a -= _N2;
				b = a + 1;
				TransferOperatorsInterface::PortraitWidth4D(0.5, 0.5, 0.375, 0.125, 0.375, 0.125, a, b, c, d, num,
					_Cell, _Portrait, _N2);
				b = a + 1;
				d = c + 1;
			}
			TransferOperatorsInterface::PortraitWidth2D(a, b, num, _Cell, _Portrait, _N2);
			a++; b++; c += 2; d = c + 1;
    }
    LS_STAG_KERNEL_INLINE void PortraitWidth4D8D_o(const double cellA, const double cellA2, const double cellB, const double cellB2, 
		 const double cellB3, const double cellB4, const double cellC, const double cellC2, const double cellD, const double cellD2, 
		 const double cellD3, const double cellD4, int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2, const int _M2)
    {
			for (int j = 0; j < _M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth4D(cellA, cellA2, cellB, cellB2, cellB3, cellB4, a, b, aa, bb, num,
					_Cell, _Portrait, _N2);
				b = a - _N2;
				a = b - 1;
				bb = aa - _N2;
				aa = bb - 1; 
				TransferOperatorsInterface::PortraitWidth8D(cellC, cellC, cellC2, cellC2, 
					cellD, cellD2, cellD, cellD2, cellD3, cellD4, cellD3, cellD4, 
					a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
				b = a + 1;
				d = c + 1;
			}
			TransferOperatorsInterface::PortraitWidth4D(cellA, cellA2, cellB, cellB2, cellB3, cellB4, a, b, aa, bb, num,
				_Cell, _Portrait, _N2);
			c += 2; d = c + 1;
    }
    LS_STAG_KERNEL_INLINE void PortraitWidth4D8D(const double cellA, const double cellA2, const double cellB, const double cellB2, 
		 const double cellB3, const double cellB4, const double cellC, const double cellC2, const double cellC3, 
		 const double cellC4, const double cellD, const double cellD2, const double cellD3, const double cellD4, 
		 int& a, int& b, int& c, int& d, int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, 
			LS_STAG_R_PI _Portrait, const int _N2, const int _M2, const int NR)
			{
				c = a + NR;
			d = c + 1;
			cc = aa + NR;
			dd = cc + 1;
				TransferOperatorsInterface::PortraitWidth4D(cellA, cellA2, cellB, cellB2, cellB3, cellB4, a, b, aa, bb, num,
					_Cell, _Portrait, _N2);
			a -= NR;
			aa -= NR;
			b = a + 1;
			bb = aa + 1;
			for (int j = 0; j < _M2; j++)
			{
				TransferOperatorsInterface::PortraitWidth8D(cellB, cellB2, cellB3, cellB4, cellC, cellC2, cellC3, cellC4,
					cellD, cellD2, cellD3, cellD4, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
					c = a; cc = aa;
				a -= NR;// CHECK
				aa -= NR;
				b = a + 1;
				bb = aa + 1;
			d = c + 1;
			dd = cc + 1;
				TransferOperatorsInterface::PortraitWidth8D(cellB2, cellB, cellB4, cellB3, cellC3, cellC4, cellC, cellC2, 
					cellD3, cellD4, cellD, cellD2, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _N2);
			}
				TransferOperatorsInterface::PortraitWidth4D(cellA, cellA2, cellB, cellB2, cellB3, cellB4, a, b, aa, bb, num,// eee order
					_Cell, _Portrait, _N2);
			}
		LS_STAG_KERNEL_INLINE void PortraitWidth8_1o(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8_1oo(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8_1(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth8_1xy(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR, const int _MR, const int NP)
		{
			for (int i = 0; i < _MR; i++, a = c, aa = cc, b = a + 1, bb = aa + 1, c = a + NP, d = c + 1, cc = aa + NP, dd = cc + 1)
			{
				TransferOperatorsInterface::PortraitWidth8_1(cellA, cellB, cellC, cellD, a, b, c, d, aa, bb, cc, dd, num, _Cell, _Portrait, _NR);
			}//for j
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth12_tb1o(const double cellA, const double cellB, const double cellC, const double cellD, int& a, int& b, int& c, int& d,
			int& aa, int& bb, int& cc, int& dd, int& aaa, int& bbb, int& ccc, int& ddd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2, aaa += 2, bbb += 2, ccc += 2, ddd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, aaa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bbb, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, ccc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, ddd, cellC);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth12_1o(const double cellA, const double cellB, const double cellC, const double cellD, const double cellE, const double cellF, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& aaa, int& bbb, int& ccc, int& ddd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2, aaa += 2, bbb += 2, ccc += 2, ddd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, aaa, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, bbb, cellF);
				SparseMatrix::setPC(num, _Cell, _Portrait, ccc, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, ddd, cellF);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth12_1oo(const double cellA, const double cellB, const double cellC, const double cellD, const double cellE, const double cellF, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& aaa, int& bbb, int& ccc, int& ddd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2, aaa += 2, bbb += 2, ccc += 2, ddd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellF);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aaa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bbb, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, ccc, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, ddd, cellF);
			}
		};
		LS_STAG_KERNEL_INLINE void PortraitWidth12_1(const double cellA, const double cellB, const double cellC, const double cellD, const double cellE, const double cellF, int& a, int& b, int& c, int& d, 
			int& aa, int& bb, int& cc, int& dd, int& aaa, int& bbb, int& ccc, int& ddd, int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait, const int _NR)
		{
			for (int j = 0; j < _NR; j++, a += 2, b += 2, c += 2, d += 2, aa += 2, bb += 2, cc += 2, dd += 2, aaa += 2, bbb += 2, ccc += 2, ddd += 2)
			{
				SparseMatrix::setPC(num, _Cell, _Portrait, a, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, b, cellA);
				SparseMatrix::setPC(num, _Cell, _Portrait, c, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, d, cellB);
				SparseMatrix::setPC(num, _Cell, _Portrait, aa, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, bb, cellC);
				SparseMatrix::setPC(num, _Cell, _Portrait, cc, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, dd, cellD);
				SparseMatrix::setPC(num, _Cell, _Portrait, aaa, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, bbb, cellE);
				SparseMatrix::setPC(num, _Cell, _Portrait, ccc, cellF);
				SparseMatrix::setPC(num, _Cell, _Portrait, ddd, cellF);
			}
		};
	  protected:
		virtual void buildNumP(){};
		/*!
		Cтроит и заполняет портрет разреженной матрицы дискретного аналога оператора продолжения
		*/
		virtual void buildPortraitP(){};
		virtual void buildNumR(){};
		/*!
		Cтроит и заполняет портрет разреженной матрицы дискретного аналога оператора проекции
		*/
		virtual void buildPortraitR(){};
		/*!
		Cтроит и заполняет разреженную матрицу дискретного аналога оператора проекции
		*/
		virtual void buildR()
		{
			((TransposedSparseMatrix*)m_R)->set(m_P);
			m_R->m_N = m_NR;//17.01.14
			m_R->m_M = m_MR;//17.01.14
			m_R->m_K = m_KR;
		};
		void buildFullR()
		{// для сочетаний четностей, где R != P^T
			m_R->set(m_NR, m_MR, m_KR, m_NMKR);
			buildNumR();
			m_R->resize();
			m_R->maxElemCount = m_R->elemCount;
			m_R->reserve();
			buildPortraitR();
		}
		/*!
		Cтроит и заполняет разреженную матрицу дискретного аналога оператора продолжения
		*/
		void buildP()
		{
			m_P->set(m_NP, m_MP, m_KP, m_NMKP);
			buildNumP();
			m_P->resize();
			m_P->maxElemCount = m_P->elemCount;
			m_P->reserve();
			buildPortraitP();
		};
	  public:
	   virtual void build(SparseMatrix* _P, TransposedSparseMatrix* _R,
                       const int _NP, const int _MP, const int _KP,
                       const int _NR, const int _MR, const int _KR)
    {
        m_R = (SparseMatrix*)_R;
        m_P = _P;
        m_NP = _NP; m_MP = _MP; m_KP = _KP;
        m_NR = _NR; m_MR = _MR; m_KR = _KR;
        m_NMKP = m_NP * m_MP * m_KP;
        m_NMKR = m_NR * m_MR * m_KR;
        N2 = m_NR - 1; M2 = m_MR - 1; K2 = m_KR - 1;
			buildP();
			buildR();
		};
	  protected:
		SparseMatrix* m_P;///< Указатель на разреженную матрицу оператора продолжения (prolongation operator)
		SparseMatrix* m_R;///< Указатель на разреженную матрицу оператора проекции (restriction operator), нет при id=OddOdd
		int m_NP, m_MP, m_KP, m_NR, m_MR, m_KR, m_NMKP, m_NMKR, N2, M2, K2;
	  public:
		OddEvenGrid3D id;///< Идентификатор структуры операторов перехода между сетками 
		int size;      ///< Размер структуры операторов перехода между сетками
	};
}

#endif //_LS_STAG_TRANSFER_OPERATORS_INTERFACE_H_
