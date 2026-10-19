/*
* \copyright (C) 2013 Puzikova V.V.
* \file      transposed_sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      07.10.2013 22:21:03
* \brief     transposed_sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRANSPOSED_SPARSE_MATRIX_H_
#define _LS_STAG_TRANSPOSED_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include "../profiler.h"
#include "sparse_matrix.h"
#include "../MeshFunctions/mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TransposedSparseMatrix

Данный файл содержит в себе реализацию класса TransposedSparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(transpose, LS_STAG_R_PI, CNT, LS_STAG_R_PI, Number,
		LS_STAG_R_PI, Portrait, LS_STAG_R_PD, Cell, LS_STAG_CR_PI, m_Number, LS_STAG_CR_PI, m_Portrait,
		LS_STAG_CR_PD, m_Cell, const int, elemCount, const int, NM, const int, CNM);
/*!
\brief Класс, реализующий транспонирование разреженной матрицы

Данный класс реализует транспонирование разреженной матрицы и ее хранение
*/
	struct TransposedSparseMatrix : SparseMatrix
	{
		TransposedSparseMatrix():
			SparseMatrix(),
		tmp_capacity(0)
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(transpose, LS_STAG_R_PI CNT, LS_STAG_R_PI Number, LS_STAG_R_PI Portrait,
			LS_STAG_R_PD Cell, LS_STAG_CR_PI m_Number, LS_STAG_CR_PI m_Portrait, LS_STAG_CR_PD m_Cell,
			const int elemCount, const int NM, const int CNM)
		{
			for (int i = 0; i < elemCount; i++)
				CNT[m_Portrait[i]]++;
			for (int i = 0; i < NM; i++)
			{
				Number[i + 1] = CNT[i] + Number[i];
				CNT[i] = 0;
			} //memset((void*)CNT, 0, tmp_capacity);
			for (int i = 0; i < CNM; i++)
			{
				for (int j = m_Number[i], num2 = m_Number[i + 1]; j < num2; j++)
				{
					const int LP = m_Portrait[j];
					const int k = Number[LP] + CNT[LP];
					Portrait[k] = i;
					Cell[k] = m_Cell[j];
					CNT[LP]++;
				}
			}
		}
		/*!
		При перестроении матрицы также выполняет транспонирование
		*/
		void resize()// (m_Matrix)^T
		{
			LS_STAG_PROFILE("SpT");
			if(tmp_capacity > 0)
			{// else работали как с обычной SparseMatrix
				LS_STAG_PAR_FOR(MeshFunction<int>, set, NMK, 0, CNT)
				LS_STAG_SINGLE(TransposedSparseMatrix, transpose, CNT, Number, Portrait, Cell, m_Number, 
					m_Portrait, m_Cell, m_Matrix->elemCount, NMK, CNM)
			}
			else SparseMatrix::resize();
		};
		using SparseMatrix::set;
		void set(SparseMatrix* _Matrix)
		{
			m_Matrix = _Matrix;
			m_N = m_Matrix->m_N;// вообще, они здесь не нужны
			m_M = m_Matrix->m_M;// вообще, они здесь не нужны
			m_K = m_Matrix->m_K;// вообще, они здесь не нужны
			m_Cell = m_Matrix->Cell;
			m_Number = m_Matrix->Number;
			m_Portrait = m_Matrix->Portrait;
			NMK = m_Matrix->m_Portrait.getElem(m_Matrix->elemCount-1) + 1; 
			CNM  = m_Matrix->NMK;
			init();
			elemCount = m_Matrix->elemCount;
			maxElemCount = m_Matrix->maxElemCount;
			reserve();
			m_CNT.init(NMK, CNT);
			tmp_capacity = NMK * sizeof(int);
			resize();
		};
	  private:
		SparseMatrix* m_Matrix;///< Указатель на транспонируемую матрицу
		double* m_Cell;
		int *m_Number, *m_Portrait;
		int CNM;// м-цы м.б. прямоугольными
		Array<int> m_CNT;
		int* CNT;///< Вспомогательный массив для определения числа элементов в строках транспонированной матрицы
		size_t tmp_capacity;
	}; 
}

#endif //_LS_STAG_TRANSPOSED_SPARSE_MATRIX_H_
