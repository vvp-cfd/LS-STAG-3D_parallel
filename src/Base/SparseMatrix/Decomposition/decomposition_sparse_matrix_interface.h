/*
* \copyright (C) 2013 Puzikova V.V.
* \file      decomposition_sparse_matrix_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2013 0:13:16
* \brief     decomposition_sparse_matrix_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DECOMPOSITION_SPARSE_MATRIX_INTERFACE_H_
#define _LS_STAG_DECOMPOSITION_SPARSE_MATRIX_INTERFACE_H_

// ================================================================== INCLUDES
#include "../transposed_sparse_matrix.h"
#include "../../MeshFunctions/mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DecompositionSparseMatrixInterface

Данный файл содержит в себе реализацию класса DecompositionSparseMatrixInterface
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(calcDirectGauss, LS_STAG_R_PD, Linv_b, LS_STAG_CR_PD, L_Cell, LS_STAG_CR_PI, L_Number,
		LS_STAG_CR_PI, L_Portrait, LS_STAG_CR_PD, _b, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(calcInverseGauss, LS_STAG_R_PD, m_x, LS_STAG_CR_PD, Linv_b, LS_STAG_CR_PD, U_Cell,
		LS_STAG_CR_PI, U_Number, LS_STAG_CR_PI, U_Portrait, const int, NM2);
	/*!
	\brief Базовый класс-интерфейс для работы с разложением разреженной матрицы

	В данном базовом классе определен интерфейс для работы с разложением разреженной матрицы
	*/
	struct DecompositionSparseMatrixInterface
	{
		DecompositionSparseMatrixInterface() :
			NMK(0)
		{
			;
		}
		LS_STAG_KERNEL(inverseGaussRow, LS_STAG_R_PD m_x, LS_STAG_CR_PD Linv_b, LS_STAG_CR_PD U_Cell,
			LS_STAG_CR_PI U_Number, LS_STAG_CR_PI U_Portrait)
		{
			double aux = Linv_b[i];
			int num = U_Number[i];
			for (int j = num + 1, num2 = U_Number[i + 1]; j < num2; j++)
				aux -= U_Cell[j] * m_x[U_Portrait[j]];
			m_x[i] = aux / U_Cell[num];
		}
		LS_STAG_KERNEL(directGaussRow, LS_STAG_R_PD Linv_b, LS_STAG_CR_PD L_Cell, LS_STAG_CR_PI L_Number,
			LS_STAG_CR_PI L_Portrait, LS_STAG_CR_PD _b)
		{
			double aux = _b[i];
			for (int j = L_Number[i], num2 = L_Number[i + 1] - 1; j < num2; j++)
				aux -= L_Cell[j] * Linv_b[L_Portrait[j]];
			Linv_b[i] = aux;
		}
		LS_STAG_SINGLE_KERNEL(calcInverseGauss, LS_STAG_R_PD m_x, LS_STAG_CR_PD Linv_b, LS_STAG_CR_PD U_Cell,
			LS_STAG_CR_PI U_Number, LS_STAG_CR_PI U_Portrait, const int NM2)
		{
			for (int i = NM2; i > -1; i--)
				DecompositionSparseMatrixInterface::_k_inverseGaussRow(i, m_x, Linv_b, U_Cell, U_Number, U_Portrait);
		}
		LS_STAG_SINGLE_KERNEL(calcDirectGauss, LS_STAG_R_PD Linv_b, LS_STAG_CR_PD L_Cell, LS_STAG_CR_PI L_Number,
			LS_STAG_CR_PI L_Portrait, LS_STAG_CR_PD _b, const int NM)
		{
			for (int i = 0; i < NM; i++)
				DecompositionSparseMatrixInterface::_k_directGaussRow(i, Linv_b, L_Cell, L_Number, L_Portrait, _b);
		}
		/*!
		Производит расчет элементов разреженных матриц, составляющих разложение
		*/
		virtual void decompose() {};
	protected:
		/*!
		Оценивает максимальное число ненулевых элементов разреженных матриц, составляющих разложение
		*/
		virtual void estimateNeedMemory() {};
		/*!
		Выделяет память под разреженные матрицы из разложения.
		Таким образом, при изменении портрета матриц из-за движения погруженной границы не требуется каждый раз
		освобождать и выделять память для ее разложения.
		*/
		void reserveLU()
		{// вызывается из LU и bandLU
			L.finResize();
			U.finResize();
		};
	protected:
		/*!
		Выполняет прямой ход метода Гаусса
		*/
		virtual void DirectGauss()
		{
			LS_STAG_PROFILE("DirectGauss");
			LS_STAG_SINGLE(DecompositionSparseMatrixInterface, calcDirectGauss, Linv_b, L_Cell, L_Number,
				L_Portrait, m_b, NMK)
		};
		virtual void DirectGauss(const double* _b)
		{
			LS_STAG_PROFILE("DirectGauss");
			LS_STAG_SINGLE(DecompositionSparseMatrixInterface, calcDirectGauss, Linv_b, L_Cell, L_Number,
				L_Portrait, _b, NMK)
		};
		/*!
		Выполняет обратный ход метода Гаусса
		*/
		virtual void InverseGauss()
		{
			LS_STAG_PROFILE("InverseGauss");
			LS_STAG_SINGLE(DecompositionSparseMatrixInterface, calcInverseGauss, m_x, Linv_b, U_Cell, U_Number,
				U_Portrait, NM2)
		};
	private:
		/*!
		Выполняет метода Гаусса
		*/
		void GaussMethod()
		{
			LS_STAG_PROFILE("GaussMeth");
			DirectGauss();
			InverseGauss();
		};
		void GaussMethod(const double* _b)
		{
			LS_STAG_PROFILE("GaussMeth");
			DirectGauss(_b);
			InverseGauss();
		};
	public:
		/*!
		Производит решение СЛАУ (системы линейных алгебраических уравнений),
		для матрицы которой было построено данное разложение, методом Гаусса с использованием хранимого разложения.
		\param _b указатель на массив, в котором хранится правая часть системы
		\param _x указатель на массив, в котором хранится решение системы
		*/
		void solve(const double* _b, double* _x)
		{
			LS_STAG_PROFILE("GaussSolve");
			m_x = _x;
			GaussMethod(_b);
		};
		void setSolAndRight(double* _b, double* _x)
		{
			m_b = _b;
			m_x = _x;
		};
		void setSolAndRight(extDouble * _b, extDouble * _x)
		{
			_b->getCopy(m_b);
			_x->getCopy(m_x);
		};
		void solve()
		{
			LS_STAG_PROFILE("GaussSolve");
			GaussMethod();
		};
		/*!
		Вызывает пересчет декомпозиции при перестроении разреженной матрицы
		*/
		void rebuild()//24.11.14
		{
			decompose();
		};
		virtual void set(SparseMatrix * _LU)
		{
			LU = _LU;
			NMK = LU->NMK;
			NM2 = NMK - 1;
			m_Linv_b.init(NMK, Linv_b);
			solution_capacity = NMK * sizeof(double);
			L.setCopy(LU);
			U.setCopy(LU);
			estimateNeedMemory();
			L_Number = L.Number;
			L_Portrait = L.Portrait;
			L_Cell = L.Cell;
			U_Number = U.Number;
			U_Portrait = U.Portrait;
			U_Cell = U.Cell;
			LU_Number = LU->Number;
			LU_Portrait = LU->Portrait;
			LU_Cell = LU->Cell;
			decompose();
		};
		void print(const std::string & name)
		{
			L.print(name + "_L");
			U.print(name + "_U");
		}
		void printNotSparse(const std::string & name)
		{
			LU->printNotSparse(name);
			L.printNotSparse(name + "_L");
			U.printNotSparse(name + "_U");
		}
	protected:
		SparseMatrix L, U;///< Разреженные матрицы из разложения (декомпозиции) разреженной матрицы
		SparseMatrix * LU;///< Указатель на разреженную матрицу, для которой строится декомпозиция
		double* L_Cell, *U_Cell, *LU_Cell;
		int* L_Number, *U_Number, *LU_Number, *L_Portrait, *U_Portrait, *LU_Portrait;
		Array<double> m_Linv_b;
		double* Linv_b, *m_x, *m_b;
		int NMK, NM2;
	private:
		size_t solution_capacity;
	};
}

#endif //_LS_STAG_DECOMPOSITION_SPARSE_MATRIX_INTERFACE_H_
