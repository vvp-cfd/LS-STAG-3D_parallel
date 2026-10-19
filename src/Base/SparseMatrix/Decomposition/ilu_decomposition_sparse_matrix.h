/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ilu_decomposition_sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2013 0:24:04
* \brief     ilu_decomposition_sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ILU_DECOMPOSITION_SPARSE_MATRIX_H_
#define _LS_STAG_ILU_DECOMPOSITION_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include "math.h"
#include "../../ls_stag_real_constants.h"
#include "decomposition_sparse_matrix_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ILUDecompositionSparseMatrix

Данный файл содержит в себе реализацию класса ILUDecompositionSparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(portraitILU, LS_STAG_R_PI, L_Number, LS_STAG_R_PI, U_Number,
		LS_STAG_CR_PI, LU_Number, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PI, U_Portrait, LS_STAG_CR_PI, LU_Portrait, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(factorization, LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, U_Number,
		LS_STAG_CR_PI, LU_Number, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PI, U_Portrait, LS_STAG_R_PD, L_Cell,
		LS_STAG_R_PD, U_Cell, LS_STAG_CR_PD, LU_Cell, const int, NM);
/*!
\brief Класс для работы с ILU-разложением разреженной матрицы

В данном классе реализована работа с ILU-разложением (неполным LU-разложением) разреженной матрицы.
*/
	struct ILUDecompositionSparseMatrix : DecompositionSparseMatrixInterface
	{
		LS_STAG_SINGLE_KERNEL(portraitILU, LS_STAG_R_PI L_Number, LS_STAG_R_PI U_Number, LS_STAG_CR_PI LU_Number,
			LS_STAG_R_PI L_Portrait, LS_STAG_R_PI U_Portrait, LS_STAG_CR_PI LU_Portrait, const int NM)
		{
			for (int i = 0, numL = 0, numU = 0, i2 = -1, i3 = 1; i < NM; i++, i2++, i3++)
			{
				for (int j = LU_Number[i], num2 = LU_Number[i3]; j < num2; j++)
				{
					const int k = LU_Portrait[j];
					if (k < i3)
					{
						L_Portrait[numL] = k;//портрет L
						numL++;    //количество элементов нижней  треугольной части А
					}
					if (k > i2)
					{
						U_Portrait[numU] = k;//портрет U
						numU++;    //количество элементов вехней  треугольной части А
					}
				}
				L_Number[i3] = numL;
				U_Number[i3] = numU;
			} //for i
		}
		LS_STAG_SINGLE_KERNEL(factorization, LS_STAG_CR_PI L_Number, LS_STAG_CR_PI U_Number, LS_STAG_CR_PI LU_Number,
			LS_STAG_R_PI L_Portrait, LS_STAG_R_PI U_Portrait, LS_STAG_R_PD L_Cell, LS_STAG_R_PD U_Cell, LS_STAG_CR_PD LU_Cell,
			const int NM)
		{
			for (int k = 0; k < NM; k++)
			{
				const int num2 = L_Number[k];
				int num = L_Number[k + 1] - 1;
				const int i2 = LU_Number[k];
				//вычисление элементов L
				for (int j = num2, i3 = i2 - num2; j < num; j++)
				{
					const int number = L_Portrait[j];
					double aux = 0.0;
					for (int i = num2; i < j; i++)
					{
						double coef = 0.0;
						for (int numL = L_Portrait[i], h = U_Number[numL], numU = U_Number[numL + 1]; h < numU; h++)
						{
							if (U_Portrait[h] == number)
							{
								coef = U_Cell[h];
								h = numU;
							}
						}
						aux += L_Cell[i] * coef;
					}//for i
					double curU = U_Cell[U_Number[number]];
					if (fabs(curU) < EQ_TINY_DOUBLE) curU = (curU < 0.0) ? -EQ_TINY_DOUBLE : EQ_TINY_DOUBLE;
					L_Cell[j] = (LU_Cell[i3 + j] - aux) / curU;
				}// for j
				L_Cell[num] = 1.0;//вычисление элементов L
				//вычисление элементов U
				num++;//L.Number[k+1]
				for (int numU2 = U_Number[k], j = numU2, i3 = i2 - 1 + num - num2 - numU2, numU = U_Number[k + 1]; j < numU; j++)
				{
					double aux = 0.0;
					const int number = U_Portrait[j];
					for (int i = num2; i < num; i++)
					{
						double coef = 0.0;
						for (int numL = L_Portrait[i], h = U_Number[numL], j2 = U_Number[numL + 1]; h < j2; h++)
						{
							if (U_Portrait[h] == number)
							{
								coef = U_Cell[h];
								h = j2;
							}
						}
						aux += L_Cell[i] * coef;
					}
					double curU = LU_Cell[i3 + j] - aux;
					if ((number == k) && (fabs(curU) < PLUS_TINY)) curU = (curU < 0.0) ? -PLUS_TINY : PLUS_TINY;
					U_Cell[j] = curU;
				}//for i
			}//for k
		}
	  private:
		/*!
		Исходя из того, что ILU-разложение не приводит к заполнению портрета матриц, 
		оценивает максимальное число ненулевых элементов разреженных матриц, составляющих ILU-разложение
		*/
		void estimateNeedMemory()
		{
			int num = NMK + (int)floor(0.5 * LU->maxElemCount);
			L.maxElemCount = num;
			L.reserve();
			U.maxElemCount = num;
			U.reserve();
		};
	  public:
		/*!
		Производит расчет ILU-разложения разреженной матрицы
		*/
		void decompose()
		{
			LS_STAG_PROFILE("ILUDecomp");
			LS_STAG_SINGLE(ILUDecompositionSparseMatrix, portraitILU, L_Number, U_Number, LU_Number, L_Portrait, U_Portrait, 
				LU_Portrait, NMK)//создание портретов
		  L.resize();
		  U.resize();
		  //портреты готовы
		  LS_STAG_PAR_FOR(MeshFunction<double>, set, U.elemCount, 0.0, U_Cell)
			LS_STAG_SINGLE(ILUDecompositionSparseMatrix, factorization, L_Number, U_Number, LU_Number, L_Portrait, U_Portrait, 
				L_Cell, U_Cell, LU_Cell,  NMK)// факторизация
		};
	};
}

#endif //_LS_STAG_ILU_DECOMPOSITION_SPARSE_MATRIX_H_
