/*
* \copyright (C) 2025 Puzikova V.V.
* \file      ilu_decomposition_5p_stencil_matrix.h
* \author    Puzikova Valeria (vvp)
* \date      01.04.2025 16:23:29
* \brief     ilu_decomposition_5p_stencil_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ILU_DECOMPOSITION_5P_STENCIL_MATRIX_H_
#define _LS_STAG_ILU_DECOMPOSITION_5P_STENCIL_MATRIX_H_

// ================================================================== INCLUDES
#include "ilu_decomposition_sparse_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ILUDecomposition5pStencilMatrix

Данный файл содержит в себе реализацию класса ILUDecomposition5pStencilMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calcInverseGauss5p, LS_STAG_R_PD, m_x, LS_STAG_CR_PD, Linv_b, LS_STAG_CR_PD, U_Cell,
		LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, U_Portrait, const int, levelNum, const int, _N, const int, _N2);
	LS_STAG_DECLARATION_KERNEL(calcDirectGauss5p, LS_STAG_R_PD, Linv_b, LS_STAG_CR_PD, L_Cell, LS_STAG_CR_PI, L_Number,
		LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PD, _b, const int, levelNum, const int, _N, const int, _N2);
	/*!
	\brief Класс для работы с ILU-разложением разреженной матрицы, полученной при использовании 5-точечного шаблона

	В данном классе реализована работа с ILU-разложением (неполным LU-разложением) разреженной матрицы.
	*/
	struct ILUDecomposition5pStencilMatrix : ILUDecompositionSparseMatrix
	{
		ILUDecomposition5pStencilMatrix() :
			levelNum(0),
			levelNum2(0),
			levelConstNum(0),
			levelDecreaseNum(0),
			N(0),
			N2(0)
		{
			;
		}
		LS_STAG_KERNEL(calcInverseGauss5p, LS_STAG_R_PD m_x, LS_STAG_CR_PD Linv_b, LS_STAG_CR_PD U_Cell,
			LS_STAG_CR_PI U_Number, LS_STAG_CR_PI U_Portrait, const int levelNum, const int _N, const int _N2)
		{
			const int k = levelNum * _N - (((_N2 < levelNum) ? _N2 : levelNum) - i) * _N2;
			DecompositionSparseMatrixInterface::_k_inverseGaussRow(k, m_x, Linv_b, U_Cell, U_Number, U_Portrait);
		}
		LS_STAG_KERNEL(calcDirectGauss5p, LS_STAG_R_PD Linv_b, LS_STAG_CR_PD L_Cell, LS_STAG_CR_PI L_Number,
			LS_STAG_CR_PI L_Portrait, LS_STAG_CR_PD _b, const int levelNum, const int _N, const int _N2)
		{
			const int k = levelNum * _N - (((_N2 < levelNum) ? _N2 : levelNum) - i) * _N2;
			DecompositionSparseMatrixInterface::_k_directGaussRow(k, Linv_b, L_Cell, L_Number, L_Portrait, _b);
		}
		void set(SparseMatrix* _LU)
		{
			N = _LU->m_N;
			N2 = N - 1;
			levelNum = _LU->m_M + N2;
			levelNum2 = levelNum - 1;
			if (_LU->m_M > N)
			{
				levelConstNum = N;
				levelDecreaseNum = _LU->m_M;
			}
			else
			{
				levelConstNum = _LU->m_M;
				levelDecreaseNum = N;
			}
			DecompositionSparseMatrixInterface::set(_LU);
		}
	private:
		inline int levelSize(const int curLevelNum) 
		{ 
			return (curLevelNum < levelConstNum) ? (curLevelNum + 1) :
				((curLevelNum < levelDecreaseNum) ? levelConstNum : (levelNum - curLevelNum));
		}
		/*!
		Выполняет прямой ход метода Гаусса
		*/
		void DirectGauss()
		{
			LS_STAG_PROFILE("DirectGauss");
			for(int level = 0; level < levelNum; level++)
				LS_STAG_PAR_FOR(ILUDecomposition5pStencilMatrix, calcDirectGauss5p, levelSize(level), Linv_b, L_Cell, 
					L_Number, L_Portrait, m_b, level, N, N2)
		};
		void DirectGauss(const double* _b)
		{
			LS_STAG_PROFILE("DirectGauss");
			for (int level = 0; level < levelNum; level++)
				LS_STAG_PAR_FOR(ILUDecomposition5pStencilMatrix, calcDirectGauss5p, levelSize(level), Linv_b, L_Cell, 
					L_Number, L_Portrait, _b, level, N, N2)
		};
		/*!
		Выполняет обратный ход метода Гаусса
		*/
		void InverseGauss()
		{
			LS_STAG_PROFILE("InverseGauss");
			for (int level = levelNum2; level > -1; level--)
				LS_STAG_PAR_FOR(ILUDecomposition5pStencilMatrix, calcInverseGauss5p, levelSize(level), m_x, Linv_b, 
					U_Cell, U_Number, U_Portrait, level, N, N2)
		};
		int levelNum, levelNum2, levelConstNum, levelDecreaseNum, N, N2;
	};
}

#endif //_LS_STAG_ILU_DECOMPOSITION_5P_STENCIL_MATRIX_H_
