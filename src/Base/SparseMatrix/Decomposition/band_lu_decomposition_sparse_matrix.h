/*
* \copyright (C) 2013 Puzikova V.V.
* \file      band_lu_decomposition_sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2013 0:25:50
* \brief     band_lu_decomposition_sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BAND_LU_DECOMPOSITION_SPARSE_MATRIX_H_
#define _LS_STAG_BAND_LU_DECOMPOSITION_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include "math.h"
#include "decomposition_sparse_matrix_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BandLUDecompositionSparseMatrix

Данный файл содержит в себе реализацию класса BandLUDecompositionSparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(portraitLU, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, L_Number,
		LS_STAG_R_PI, U_Portrait, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PD, L_Cell, const int, L_m_N);
	LS_STAG_DECLARATION_S_KERNEL(numBuildBLU, LS_STAG_R_PI, L_Number, LS_STAG_R_PI, U_Number,
		const int, L_m_N, const int, width_band, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(calcLU, LS_STAG_R_PD, L_Cell, LS_STAG_R_PD, U_Cell,
		LS_STAG_CR_PD, LU_Cell, LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, LU_Number,
		LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PI, LU_Portrait, const int, width_band, const int, NM, const int, L_m_N);
/*!
\brief Класс для работы с ленточным LU-разложением разреженной матрицы

В данном классе реализована работа с LU-разложением разреженной матрицы, которая рассматривается как ленточная.
*/
	struct BandLUDecompositionSparseMatrix : DecompositionSparseMatrixInterface
	{
		LS_STAG_KERNEL(portraitLU, LS_STAG_CR_PI U_Number, LS_STAG_CR_PI L_Number, LS_STAG_R_PI U_Portrait,
			LS_STAG_R_PI L_Portrait, LS_STAG_R_PD L_Cell, const int L_m_N)
		{// для LU_band портрет можно построить сразу
			const int numL = L_Number[i + 1];
			for (int numL2 = L_Number[i], j = numL2, num2 = i - numL2 - L_m_N; j < numL; j++)
				L_Portrait[j] = (i < L_m_N) ? (j - numL2) : (num2 + j);
			for (int j = U_Number[i], num = i - j, numU = U_Number[i + 1]; j < numU; j++)
				U_Portrait[j] = j + num;//портрет U
			L_Cell[numL - 1] = 1.0;//диагональ
		}
		LS_STAG_SINGLE_KERNEL(numBuildBLU, LS_STAG_R_PI L_Number, LS_STAG_R_PI U_Number, const int L_m_N,
			const int width_band, const int NM)
		{
			for (int i = 0, i2 = 1; i < NM; i++, i2++)
			{
				L_Number[i2] = L_Number[i] + ((i < L_m_N) ? i2 : width_band);
				U_Number[i2] = U_Number[i] + (((i + L_m_N) < NM) ? width_band : (NM - i));
			}
		}
		LS_STAG_SINGLE_KERNEL(calcLU, LS_STAG_R_PD L_Cell, LS_STAG_R_PD U_Cell, LS_STAG_CR_PD LU_Cell,
			LS_STAG_CR_PI L_Number, LS_STAG_CR_PI U_Number, LS_STAG_CR_PI LU_Number, LS_STAG_CR_PI L_Portrait,
			LS_STAG_CR_PI LU_Portrait, const int width_band, const int NM, const int L_m_N)
		{
			for (int j = 0, num2 = LU_Number[0], num = LU_Number[1]; j < width_band; j++)// 1-я строка U
			{
				double val = 0.0;
				for (int i = num2; i < num; i++)
				{
					if (LU_Portrait[i] == j)
					{
						val = LU_Cell[i];
						i = num;
					}
				}
				U_Cell[j] = val;
			};
			double aux = 1.0 / U_Cell[0];
			for (int j = 1; j < width_band; j++)
			{// в этих строках есть 1-й столбец в L (кроме 1-й строки, там только диагональный элемент)
				const int num = LU_Number[j];
				L_Cell[L_Number[j]] = (LU_Portrait[num]) ? 0.0 : (LU_Cell[num] * aux);
			};
			for (int i = 1, i3 = 0, numL2 = 1 + width_band; i < NM; i++, i3++, numL2++)// заполняем остальные строки
			{
				int numL = L_Number[i];
				int numU = U_Number[i];
				int numU2 = i + U_Number[i3] - numU;
				if (numL2 > NM) numL2 = NM;
				//вычисление элементов U
				for (int j = i, num2 = LU_Number[i], num = LU_Number[i3]; j < numU2; j++)
				{
					double aux = 0.0;
					int k_min = 0;
					if (j > L_m_N) k_min = j - L_m_N;
					int j2 = j - i;
					for (int k = k_min, k2 = numL - k_min, i2 = j2 + 1; k < i; k++)
						aux -= U_Cell[U_Number[k] + i2] * L_Cell[k2 + k];
					for (int h = num2; h < num; h++)
					{
						if (LU_Portrait[h] == j)
						{
							aux += LU_Cell[h];
							h = num;
						}
					}
					U_Cell[numU + j2] = aux;
				};
				double coef = 1.0 / U_Cell[numU];
				//вычисление элементов L
				for (int j = i3; j < numL2; j++)
				{
					double aux = 0.0;
					int numL = L_Number[j];
					int k_min = 0;
					if (j > L_m_N) k_min = j - L_m_N;
					int k2 = numL - k_min;
					for (int k = k_min; k < i; k++)
						aux -= U_Cell[U_Number[k]] * L_Cell[k2 + k];
					for (int h = LU_Number[j], num = LU_Number[j + 1]; h < num; h++)
					{
						if (LU_Portrait[h] == i)
						{
							aux += LU_Cell[h];
							h = num;
						}
					}
					L_Cell[numL + i - L_Portrait[numL]] = aux * coef;
				};//вычисление элементов 
			}
		}
	  private:
		/*!
		Исходя из того, что ленточное LU-разложение приводит к полному заполнению портрета матриц 
		внутри ленты выбранной ширины, 
		оценивает максимальное число ненулевых элементов матриц, составляющих LU-разложение разреженной матрицы
		*/
		void estimateNeedMemory()
		{
			L_m_N = L.m_N;
			width_band = L_m_N + 1;
			LS_STAG_SINGLE(BandLUDecompositionSparseMatrix, numBuildBLU, L.Number, U.Number, L_m_N, width_band, NMK)
			reserveLU();
			LS_STAG_PAR_FOR(BandLUDecompositionSparseMatrix, portraitLU, NMK, U.Number, L.Number, U.Portrait, L.Portrait, L.Cell, L_m_N)
		};
	  public:
		/*!
		Производит расчет ленточного LU-разложения разреженной матрицы
		*/
		void decompose()
		{
			LS_STAG_PROFILE("BandLU");
			LS_STAG_SINGLE(BandLUDecompositionSparseMatrix, calcLU, L_Cell, U_Cell, LU_Cell, L_Number, U_Number, 
				LU_Number, L_Portrait, LU_Portrait, width_band, NMK, L_m_N)
		};
	  private:
		int width_band;///< ширина ленты
		int L_m_N;
	}; 
}


#endif //_LS_STAG_BAND_LU_DECOMPOSITION_SPARSE_MATRIX_H_
