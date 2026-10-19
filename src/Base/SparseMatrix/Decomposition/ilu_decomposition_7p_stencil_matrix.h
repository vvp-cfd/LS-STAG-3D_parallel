/*
* \copyright (C) 2026 Puzikova V.V.
* \file      ilu_decomposition_7p_stencil_matrix.h
* \author    Puzikova Valeria (vvp)
* \date      11.05.2026 22:15:11
* \brief     ilu_decomposition_7p_stencil_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ILU_DECOMPOSITION_7P_STENCIL_MATRIX_H_
#define _LS_STAG_ILU_DECOMPOSITION_7P_STENCIL_MATRIX_H_

// ================================================================== INCLUDES
#include "ilu_decomposition_sparse_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ILUDecomposition7pStencilMatrix

Данный файл содержит в себе реализацию класса ILUDecomposition7pStencilMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcInverseGauss7p, LS_STAG_R_PD, m_x, LS_STAG_CR_PD, Linv_b, LS_STAG_CR_PD, U_Cell,
		LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, U_Portrait, const int, levelNum, const int, _N, const int, _K, const int, _NM);
	LS_STAG_DECLARATION_R_KERNEL_2D(calcDirectGauss7p, LS_STAG_R_PD, Linv_b, LS_STAG_CR_PD, L_Cell, LS_STAG_CR_PI, L_Number,
		LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PD, _b, const int, levelNum, const int, _N, const int, _K, const int, _NM);
	/*!
	\brief Класс для работы с ILU-разложением разреженной матрицы, полученной при использовании 7-точечного шаблона

	В данном классе реализована работа с ILU-разложением (неполным LU-разложением) разреженной матрицы.
	*/
	struct ILUDecomposition7pStencilMatrix : ILUDecompositionSparseMatrix
	{
		ILUDecomposition7pStencilMatrix():
			levelNum(0),
			levelNum2(0),
			N(0),
			N2(0),
			K(0), 
			M(0), 
			K2(0), 
			M2(0), 
			KM2(0), 
			NM(0), 
			i_start(0), 
			i_end(0), 
			j_start(0), 
			j_end(0)
		{
			;
		}
		LS_STAG_KERNEL_2D(calcInverseGauss7p, LS_STAG_R_PD m_x, LS_STAG_CR_PD Linv_b, LS_STAG_CR_PD U_Cell,
			LS_STAG_CR_PI U_Number, LS_STAG_CR_PI U_Portrait, const int levelNum, const int _N, const int _K, const int _NM)
		{
			const int k = levelNum - i - j;
			if((k > -1) && (k < _K))
			{
			const int num = k * _NM + j * _N + i;
			DecompositionSparseMatrixInterface::_k_inverseGaussRow(num, m_x, Linv_b, U_Cell, U_Number, U_Portrait);
			}
		}
		LS_STAG_KERNEL_2D(calcDirectGauss7p, LS_STAG_R_PD Linv_b, LS_STAG_CR_PD L_Cell, LS_STAG_CR_PI L_Number,
			LS_STAG_CR_PI L_Portrait, LS_STAG_CR_PD _b, const int levelNum, const int _N, const int _K, const int _NM)
		{
			const int k = levelNum - i - j;
			if((k > -1) && (k < _K))
			{
			const int num = k * _NM + j * _N + i;
			DecompositionSparseMatrixInterface::_k_directGaussRow(num, Linv_b, L_Cell, L_Number, L_Portrait, _b);
			}
		}
		void set(SparseMatrix* _LU)
		{
			N = _LU->m_N;
			K = _LU->m_K;
			M = _LU->m_M;
			NM = N * M;
			M2 = M - 1;
			N2 = N - 1;
			K2 = K - 1; 
			KM2 = M2 + K2;
			levelNum = N + KM2;   // количество диагональных плоскостей
			levelNum2 = levelNum - 1;
			DecompositionSparseMatrixInterface::set(_LU);
		}
	private:
		inline void levelRange(const int level) 
		{ 
				i_start = (level < KM2) ? 0 : (level - KM2);
				i_end = ((N2 < level) ? N2 : level) + 1;
				j_start = level - i_end - K;
				if (j_start < 0) j_start = 0;
				j_end = level - i_start + 1;
				if (j_end > M) j_end = M;
		}
		/*!
		Выполняет прямой ход метода Гаусса
		*/
		void DirectGauss()
		{
			LS_STAG_PROFILE("DirectGauss");
			for(int level = 0; level < levelNum; level++)
			{
				levelRange(level);
				LS_STAG_PAR_FOR_2D_RANGE(ILUDecomposition7pStencilMatrix, calcDirectGauss7p, i_start, i_end, j_start, j_end,
					Linv_b, L_Cell, L_Number, L_Portrait, m_b, level, N, K, NM)
			}
		};
		void DirectGauss(const double* _b)
		{
			LS_STAG_PROFILE("DirectGauss");
			for (int level = 0; level < levelNum; level++)
			{
				levelRange(level);
				LS_STAG_PAR_FOR_2D_RANGE(ILUDecomposition7pStencilMatrix, calcDirectGauss7p, i_start, i_end, j_start, j_end,
					Linv_b, L_Cell, L_Number, L_Portrait, _b, level, N, K, NM)
			}
		};
		/*!
		Выполняет обратный ход метода Гаусса
		*/
		void InverseGauss()
		{
			LS_STAG_PROFILE("InverseGauss");
			for (int level = levelNum2; level > -1; level--)
			{
				levelRange(level);
				LS_STAG_PAR_FOR_2D_RANGE(ILUDecomposition7pStencilMatrix, calcInverseGauss7p, i_start, i_end, j_start, j_end,
					m_x, Linv_b, U_Cell, U_Number, U_Portrait, level, N, K, NM)
			}
		};
		int levelNum, levelNum2, N, K, M, N2, K2, M2, KM2, NM, i_start, i_end, j_start, j_end;
	};
}

#endif //_LS_STAG_ILU_DECOMPOSITION_7P_STENCIL_MATRIX_H_
