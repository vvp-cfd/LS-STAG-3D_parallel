/*
* \copyright (C) 2013 Puzikova V.V.
* \file      lu_decomposition_sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2013 0:22:15
* \brief     lu_decomposition_sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LU_DECOMPOSITION_SPARSE_MATRIX_H_
#define _LS_STAG_LU_DECOMPOSITION_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include "math.h"
#include "decomposition_sparse_matrix_interface.h"
#include "../../ls_stag_real_constants.h"
#include "../../MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LUDecompositionSparseMatrix

Данный файл содержит в себе реализацию класса LUDecompositionSparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(portraitLU, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, L_Number,
		LS_STAG_R_PI, U_Portrait, LS_STAG_R_PI, L_Portrait, LS_STAG_R_PD, L_Cell);
	LS_STAG_DECLARATION_KERNEL(initScalePerm, LS_STAG_CR_PI, LU_Number, LS_STAG_CR_PI, LU_Portrait, LS_STAG_CR_PD, LU_Cell, LS_STAG_R_PD, scale, LS_STAG_R_PI, perm);
	LS_STAG_DECLARATION_S_KERNEL(numBuildLU, LS_STAG_R_PI, L_Number, LS_STAG_R_PI, U_Number, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(calcLU, LS_STAG_CR_PD, scale, LS_STAG_R_PI, perm, LS_STAG_R_PD, col, LS_STAG_R_PD, L_Cell, LS_STAG_R_PD, U_Cell, LS_STAG_CR_PD, LU_Cell,
		LS_STAG_CR_PI, L_Number, LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, LU_Number, LS_STAG_CR_PI, LU_Portrait, const int, NM);
	/*!
	\brief Класс для работы с LU-разложением разреженной матрицы

	В данном классе реализована работа с LU-разложением разреженной матрицы.
	*/
	struct LUDecompositionSparseMatrix : DecompositionSparseMatrixInterface
	{
		static const int MAX_IT_REFINEMENT = 5;
		LS_STAG_KERNEL(portraitLU, LS_STAG_CR_PI U_Number, LS_STAG_CR_PI L_Number, LS_STAG_R_PI U_Portrait,
			LS_STAG_R_PI L_Portrait, LS_STAG_R_PD L_Cell)
		{// для LU портрет можно построить сразу
			const int i3 = i + 1;
			const int numU2 = U_Number[i];
			const int numL2 = L_Number[i];
			const int numU = U_Number[i3];
			const int numL = L_Number[i3];
			for (int j = numL2; j < numL; j++)
				L_Portrait[j] = j - numL2;//портрет L
			for (int j = numU2, num = i - numU2; j < numU; j++)
				U_Portrait[j] = j + num;//портрет U
			L_Cell[numL - 1] = 1.0;//диагональ
		}
		LS_STAG_SINGLE_KERNEL(numBuildLU, LS_STAG_R_PI L_Number, LS_STAG_R_PI U_Number, const int NM)
		{
			for (int i = 0; i < NM; i++)
			{
				const int numL = i + 1;
				L_Number[numL] = L_Number[i] + numL;
				U_Number[numL] = U_Number[i] + NM - i;
			}
		}
		LS_STAG_KERNEL(initScalePerm, LS_STAG_CR_PI LU_Number, LS_STAG_CR_PI LU_Portrait, LS_STAG_CR_PD LU_Cell, LS_STAG_R_PD scale, LS_STAG_R_PI perm)
		{// Вычисляем диагональные масштабы d[i] = 1 / sqrt(|A[i][i]|)
			double diag = 0.0;
			for (int k = LU_Number[i], kEnd = LU_Number[i + 1]; k < kEnd; ++k)
			{
				if (LU_Portrait[k] == i)
				{
					diag = fabs(LU_Cell[k]);
					k = kEnd;
				}
			}
			scale[i] = (diag > EQ_DOUBLE) ? (1.0 / sqrt(diag)) : 1.0;// EQ_TINY3_DOUBLE
			perm[i] = i;
		}
		LS_STAG_SINGLE_KERNEL(calcLU, LS_STAG_CR_PD scale, LS_STAG_R_PI perm, LS_STAG_R_PD col, LS_STAG_R_PD L_Cell, LS_STAG_R_PD U_Cell, LS_STAG_CR_PD LU_Cell,
			LS_STAG_CR_PI L_Number, LS_STAG_CR_PI U_Number, LS_STAG_CR_PI LU_Number, LS_STAG_CR_PI LU_Portrait,
			const int NM)
		{
			// Вычисляем maxAbsA по исходной матрице(немасштабированной) + i=0
			double maxAbsA_local = 0.0;
			double maxVal = 0.0;
			int pivot = 0;
			double si = scale[0];
			for (int i = 0; i < NM; ++i)
			{
				int k = LU_Number[i];
				double curA = LU_Cell[k];
				if (LU_Portrait[k] == 0)
				{
					double aik = curA * scale[i] * si;
					double absVal = fabs(aik);
					if (absVal > maxVal) { maxVal = absVal; pivot = i; }
				}
				double absVal = fabs(curA);
				if (absVal > maxAbsA_local) maxAbsA_local = absVal;
				int kEnd = LU_Number[i + 1];
				k++;
				if (k < kEnd)
				{
					for (; k < kEnd; ++k)
					{
						double absVal = fabs(LU_Cell[k]);
						if (absVal > maxAbsA_local) maxAbsA_local = absVal;
					}
				}
			}
			const double tiny = EQ_TINY3_DOUBLE * maxAbsA_local;
			if (pivot != 0)
			{
				int tmp = perm[0];
				perm[0] = perm[pivot];
				perm[pivot] = tmp;
			}
			double* U0 = U_Cell + U_Number[0];
			int permi = perm[0];
			si = scale[permi];
			int hStart = LU_Number[permi], hEnd = LU_Number[permi + 1];
			for (int j = 0; j < NM; ++j)
			{
				double A_0j = 0.0;
				for (int h = hStart; h < hEnd; ++h)
				{
					if (LU_Portrait[h] == j) { A_0j = LU_Cell[h]; h = hEnd; }
				}
				U0[j] = A_0j * si * scale[j];
			}
			double U00 = U0[0];
			if (fabs(U00) < tiny) U00 = (U00 < 0) ? -tiny : tiny;
			U0[0] = U00;
			double invU00 = scale[0] / U00;
			for (int k = 1; k < NM; ++k)
			{
				double A_k0 = 0.0;
				int  permk = perm[k];
				for (int h = LU_Number[permk], hEnd = LU_Number[permk + 1]; h < hEnd; ++h)
				{
					if (LU_Portrait[h] == 0) { A_k0 = LU_Cell[h]; h = hEnd; }
				}
				L_Cell[L_Number[k]] = A_k0 * scale[permk] * invU00;
			}
			L_Cell[L_Number[0]] = 1.0;
			for (int i = 1; i < NM; ++i)
			{// 1. Вычисляем столбец i для строк k >= i
				maxVal = 0.0;
				pivot = i;
				si = scale[i];
				for (int k = i; k < NM; ++k)
				{
					double sum = 0.0;
					double* Lk = L_Cell + L_Number[k];
					for (int t = 0; t < i; ++t)
					{
						double U_ti = U_Cell[U_Number[t] + (i - t)];
						sum += Lk[t] * U_ti;
					}
					double A_ki = 0.0;
					int permk = perm[k];
					for (int h = LU_Number[permk], hEnd = LU_Number[permk + 1]; h < hEnd; ++h)
					{
						if (LU_Portrait[h] == i) { A_ki = LU_Cell[h]; h = hEnd; }
					}
					double val = A_ki * scale[permk] * si - sum;
					col[k] = val;
					double absVal = fabs(val);
					if (absVal > maxVal) { maxVal = absVal; pivot = k; }
				}
				// 2. Перестановка
				if (pivot != i)
				{
					int tmp = perm[i]; perm[i] = perm[pivot]; perm[pivot] = tmp;
					double* Li = L_Cell + L_Number[i];
					double* Lp = L_Cell + L_Number[pivot];
					for (int t = 0; t < i; ++t)
					{
						double dtmp = Li[t]; Li[t] = Lp[t]; Lp[t] = dtmp;
					}
					double vtmp = col[i]; col[i] = col[pivot]; col[pivot] = vtmp;
				}
				// 3. Строка U(i, j)
				double* Ui = U_Cell + U_Number[i];
				double* Li = L_Cell + L_Number[i];
				permi = perm[i];
				si = scale[permi];
				for (int j = i; j < NM; ++j)
				{
					double sum = 0.0;
					for (int t = 0; t < i; ++t)
					{
						double U_tj = U_Cell[U_Number[t] + (j - t)];
						sum += Li[t] * U_tj;
					}
					double A_ij = 0.0;
					for (int h = LU_Number[permi], hEnd = LU_Number[permi + 1]; h < hEnd; ++h)
					{
						if (LU_Portrait[h] == j) { A_ij = LU_Cell[h]; h = hEnd; }
					}
					Ui[j - i] = A_ij * si * scale[j] - sum;
				}
				double Uii = Ui[0];
				if (fabs(Uii) < tiny) Uii = (Uii < 0) ? -tiny : tiny;
				Ui[0] = Uii;
				double invUii = 1.0 / Uii;
				// 4. Столбец L(k,i)
				for (int k = i + 1; k < NM; ++k)
				{
					L_Cell[L_Number[k] + i] = col[k] * invUii;
				}
				L_Cell[L_Number[i] + i] = 1.0;
			}
		}
	private:
		/*!
		Исходя из того, что LU-разложение приводит к полному заполнению портрета матриц,
		оценивает максимальное число ненулевых элементов матриц, составляющих LU-разложение разреженной матрицы
		*/
		void estimateNeedMemory()
		{
			LS_STAG_SINGLE(LUDecompositionSparseMatrix, numBuildLU, L.Number, U.Number, NMK)
				reserveLU();
			LS_STAG_PAR_FOR(LUDecompositionSparseMatrix, portraitLU, NMK, U.Number, L.Number, U.Portrait, L.Portrait, L.Cell)
		};
	public:
		/*!
		Производит расчет LU-разложения разреженной матрицы
		*/
		void decompose()
		{
			LS_STAG_PROFILE("LUDecomp");
			LS_STAG_PAR_FOR(LUDecompositionSparseMatrix, initScalePerm, NMK, LU_Number, LU_Portrait, LU_Cell, scale, perm)
				LS_STAG_SINGLE(LUDecompositionSparseMatrix, calcLU, scale, perm, Linv_b, L_Cell, U_Cell, LU_Cell, L_Number, U_Number,
					LU_Number, LU_Portrait, NMK)
		};
		void set(SparseMatrix * _LU) override
		{
			m_scale.init(_LU->NMK, scale);
			m_bScaled.init(_LU->NMK, bScaled);
			m_perm.init(_LU->NMK, perm);
			m_bScaled.setComponentwiseMultiplier(scale);
			DecompositionSparseMatrixInterface::set(_LU);
			m_bScaled.setMatr(_LU);
		}
	private:
		void InverseGauss() override
		{
			LS_STAG_PROFILE("InverseGauss");
			LS_STAG_SINGLE(DecompositionSparseMatrixInterface, calcInverseGauss, bScaled, Linv_b, U_Cell, U_Number,
				U_Portrait, NM2)
			m_bScaled.copyComponentwiseProductTo(m_x);//double norm_r = m_bScaled.setProductWithResidual(m_x, m_b);
		};
		void DirectGauss() override
		{
			m_bScaled.copyComponentwiseProductPerm(m_b, perm);
			DecompositionSparseMatrixInterface::DirectGauss(bScaled);
		}
		RealValuedMeshFunction m_scale, m_bScaled;
		Array<int> m_perm;
		double* scale, * bScaled;
		int* perm;
	};
}

#endif //_LS_STAG_LU_DECOMPOSITION_SPARSE_MATRIX_H_
