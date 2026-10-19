/*
* \copyright (C) 2013 Puzikova V.V.
* \file      product_sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2013 0:03:09
* \brief     product_sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_PRODUCT_SPARSE_MATRIX_H_
#define _LS_STAG_PRODUCT_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include "../profiler.h"
#include "sparse_matrix.h"
#include "../MeshFunctions/mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ProductSparseMatrix

Данный файл содержит в себе реализацию класса ProductSparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(init, LS_STAG_R_PI, _Portrait, LS_STAG_R_PD, _Cell);
	LS_STAG_DECLARATION_S_KERNEL(numSpMM, LS_STAG_R_PI, Number, LS_STAG_R_PI, mark,
		LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Number2, LS_STAG_CR_PI, m_Portrait, LS_STAG_CR_PI, m_Portrait2, const int, NMK);
	LS_STAG_DECLARATION_KERNEL(spMM, LS_STAG_R_PI, Number, LS_STAG_R_PI, Portrait, LS_STAG_R_PD, Cell,
		LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Number2, LS_STAG_CR_PI, m_Portrait, LS_STAG_CR_PI, m_Portrait2,
		LS_STAG_CR_PD, m_Cell, LS_STAG_CR_PD, m_Cell2);
/*!
\brief Класс, реализующий умножение разреженных матриц

Данный класс реализует расчет и построение разреженной матрицы, которая является произведением двух разреженных матриц,
а также ее хранение
*/
	struct ProductSparseMatrix : SparseMatrix
	{
		/*!
		При перестроении матрицы также вычисляется произведение разреженных матриц
		*/
		LS_STAG_KERNEL(init, LS_STAG_R_PI _Portrait, LS_STAG_R_PD _Cell)
		{
			_Portrait[i] = -1;
			_Cell[i] = 0.0;
		}
		LS_STAG_SINGLE_KERNEL(numSpMM, LS_STAG_R_PI Number, LS_STAG_R_PI mark, LS_STAG_CR_PI m_Number, LS_STAG_CR_PI m_Number2, 
			LS_STAG_CR_PI m_Portrait, LS_STAG_CR_PI m_Portrait2, const int NMK)
		{
			int total_nnz = 0;
			Number[0] = 0;
			for (int i = 0; i < NMK; ++i) 
			{
				for (int h = m_Number[i], hEnd = m_Number[i + 1]; h < hEnd; ++h)
				{
					int nL = m_Portrait[h];          // столбец j в A
					for (int k = m_Number2[nL], kEnd = m_Number2[nL + 1]; k < kEnd; ++k)
					{
						int UP = m_Portrait2[k];     // столбец в B
						if (mark[UP] != i) 
						{
							mark[UP] = i;
							++total_nnz;
						}
					}
				}
				Number[i + 1] = total_nnz;
			}
		}
		LS_STAG_KERNEL(spMM, LS_STAG_R_PI Number, LS_STAG_R_PI Portrait, LS_STAG_R_PD Cell,
			LS_STAG_CR_PI m_Number, LS_STAG_CR_PI m_Number2, LS_STAG_CR_PI m_Portrait, LS_STAG_CR_PI m_Portrait2,
			LS_STAG_CR_PD m_Cell, LS_STAG_CR_PD m_Cell2)
		{
			int start = Number[i];                     // начало строки в Portrait/Cell
			int count = 0;                             // сколько элементов уже вставлено
			for (int h = m_Number[i], numL = m_Number[i + 1]; h < numL; h++)
			{
				double a_val = m_Cell[h];
				for (int nL = m_Portrait[h], k = m_Number2[nL], numU = m_Number2[nL + 1]; k < numU; k++)
				{
					int col = m_Portrait2[k];
					double val = a_val * m_Cell2[k];
					// линейный поиск позиции вставки (столбцы в [start, start+count) уже отсортированы)
					int pos = start, shift = start + count;
					for (; (pos < shift) && (Portrait[pos] < col); ++pos) {}
					if ((pos < shift) && (Portrait[pos] == col))
						Cell[pos] += val;       // столбец уже существует — суммируем
					else
					{// сдвигаем хвост вправо, освобождая место
						for (int q = shift; q > pos; --q)
						{
							Portrait[q] = Portrait[q - 1];
							Cell[q] = Cell[q - 1];
						}
						Portrait[pos] = col;
						Cell[pos] = val;
						++count;
					}
				};//for(k=0;k<U.Number[nL];k++)
			};//for(h=0;h<L.Number[i];h++)
		}
		void resize()// (m_Matrix)^T
		{// buildNum
			LS_STAG_PROFILE("numSpMM");
			LS_STAG_SINGLE(ProductSparseMatrix, numSpMM, Number, mark, m_Number, m_Number2, m_Portrait, m_Portrait2, NMK)
			SparseMatrix::resize();
			calc();
		};
		void calc()
		{
			LS_STAG_PROFILE("SpMM");
			LS_STAG_PAR_FOR(ProductSparseMatrix, init, maxElemCount, Portrait, Cell)
			LS_STAG_PAR_FOR(ProductSparseMatrix, spMM, NMK, Number, Portrait, Cell, m_Number, m_Number2, m_Portrait,
					m_Portrait2, m_Cell, m_Cell2)
		}
		using SparseMatrix::set;
		void set(SparseMatrix* _Matrix, SparseMatrix* _Matrix2)
		{
			m_Cell = _Matrix->Cell;
			m_Number = _Matrix->Number;
			m_Portrait = _Matrix->Portrait;
			m_Cell2 = _Matrix2->Cell;
			m_Number2 = _Matrix2->Number;
			m_Portrait2 = _Matrix2->Portrait;
			SparseMatrix::setCopy(_Matrix);
			maxElemCount = _Matrix->maxElemCount + _Matrix2->maxElemCount;// как в 2D: запас на произведение (сумма резервов сомножителей)
			reserve();
			portrait_capacity = maxElemCount * sizeof(int);
			N_cols = _Matrix2->m_Portrait.getElem(_Matrix2->elemCount - 1) + 1;
			m_mark.init(N_cols, mark);
			m_mark.set(-1);
			resize();
		};
	  private:
		  MeshFunction<int> m_mark;
		// SparseMatrix *m_Matrix, *m_Matrix2;// указатели на матрицы-сомножители
		double *m_Cell, *m_Cell2;
		int *m_Number, *m_Number2, *m_Portrait, *m_Portrait2, *mark;
		size_t portrait_capacity;
		int N_cols;// число столбцов второй матрицы
	};
}

#endif //_LS_STAG_PRODUCT_SPARSE_MATRIX_H_
