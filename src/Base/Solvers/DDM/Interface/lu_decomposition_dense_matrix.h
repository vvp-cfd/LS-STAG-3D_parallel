/*
* \copyright (C) 2019 Puzikova V.V.
* \file      dense_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      03.08.2019 06:31:15
* \brief     dense_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LU_DECOMPOSITION_DENSE_MATRIX_H_
#define _LS_STAG_LU_DECOMPOSITION_DENSE_MATRIX_H_

// ================================================================== INCLUDES
#include "dense_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LUDecompositionDenseMatrix

Данный файл содержит в себе реализацию класса LUDecompositionDenseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий LU-разложение

Данный класс реализует построение LU-разложения заполненной матрицы и выполнение метода Гаусса при помощи него
*/
	struct LUDecompositionDenseMatrix
	{
	  LUDecompositionDenseMatrix():
		m_N(0)
		{
			;
		}
		~LUDecompositionDenseMatrix()
		{
			if(m_N > 0) delete[] Linv_b;
		}
		/*!
		Вычисляет LU-разложение заполненной матрицы и сохраняет его
		*/
	    void decompose()
		{
			LS_STAG_PROFILE("interfaceLUDecomp");
			double aux = 1.0 / LU[0];
			for(int j = m_N; j < elemCount; j += m_N)
			{
				LU[j] *= aux;
			};
			for(int start = m_N, i = N3, i2 = 1, end = start + m_N; i2 < m_N; i += N3, start = end, end += m_N, i2++)
			{	  
				//вычисление элементов U
				for(int j = i; j < end; j++)
				{
					double aux = 0.0;
					for(int k = start, k2 = j - start; k < i; k++, k2 += m_N)
						aux += LU[k2] * LU[k];
					LU[j] -= aux;
				};
				double coef = 1.0 / LU[i];
				//вычисление элементов L
				for(int j = i + m_N; j < elemCount; j += m_N)
				{
					double aux = 0.0;
					for(int k = j - i2, k2 = i2; k < j; k++, k2 += m_N)
						aux += LU[k2] * LU[k];
					LU[j] -= aux;
					LU[j] *= coef;
				};//вычисление элементов 
			}
		};
	private:
		/*!
		Выполняет прямой ход метода Гаусса
		*/
		void DirectGauss()
		{
			Linv_b[0] = m_b[0];
			for(int i = 1, start = m_N; i < m_N; i++, start += m_N)
			{
				  double aux = 0.0;
				  for(int j = 0, k = start; j < i; j++, k++)
					  aux += LU[k] * Linv_b[j]; 
				  Linv_b[i] = m_b[i] - aux;
			}
		}; 
		void DirectGauss(const double* _b)
		{
			Linv_b[0] = _b[0];
			for(int i = 1, start = m_N; i < m_N; i++, start += m_N)
			{
				  double aux = 0.0;
				  for(int j = 0, k = start; j < i; j++, k++)
					  aux += LU[k] * Linv_b[j]; 
				  Linv_b[i] = _b[i] - aux;
			}
		}; 
		/*!
		Выполняет обратный ход метода Гаусса
		*/
		void InverseGauss()
		{
			m_x[N2] = Linv_b[N2] / LU[NM2];
			for(int i = N4, p = NM3, p2 = p + 1; i > -1; i--, p -= N3, p2 = p + 1)
			{
				double aux = 0.0;
				for(int j = i + 1, k = p2; j < m_N; j++, k++)
					aux += LU[k] * m_x[j]; 
				m_x[i] = (Linv_b[i] - aux) / LU[p];
			}
		};
		/*!
		Выполняет решение СЛАУ методом Гаусса
		*/
		void GaussMethod()
		{
			LS_STAG_PROFILE("interfaceGauss");
			DirectGauss();
			InverseGauss();
		};
		void GaussMethod(const double* _b)
		{
			LS_STAG_PROFILE("interfaceGauss");
			DirectGauss(_b);
			InverseGauss();
		};
	public:
		void solve(const double* _b, double* _x)
		{
			m_x = _x;
			GaussMethod(_b);
		};
		void solve(const double* _b)
		{
			GaussMethod(_b);
		};
		void setSolAndRight(double* _b, double* _x)
		{
			m_b = _b;
			m_x = _x;
		};
		/*!
		Выполняет решение СЛАУ методом Гаусса
		*/
		void solve()
		{
			GaussMethod();
		};
		void rebuild()
		{
			decompose();
		};
		void set(DenseMatrix* _A)
		{
			LU = _A->Cell;
			elemCount = _A->elemCount;
			m_N = _A->m_N;
			N2 = m_N - 1;
			N3 = m_N + 1;
			N4 = N2 - 1;
			NM2 = elemCount - 1;
			NM3 = NM2 - N3;
			Linv_b = new double[m_N];
		};
	private:
		double *Linv_b, *m_x, *m_b;
	    double* LU; ///< элементы матрицы и в ней же записываем разложение
		int m_N, N2, elemCount, N3, NM2, N4, NM3;
	};
}
#endif //_LS_STAG_LU_DECOMPOSITION_DENSE_MATRIX_H_
