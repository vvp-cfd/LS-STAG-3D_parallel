/*
* \copyright (C) 2019 Puzikova V.V.
* \file      interface_block_solver.h
* \author    Puzikova Valeria (vvp) 
* \date      04.08.2019 22:44:05
* \brief     interface_block_solver.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_INTERFACE_BLOCK_SOLVER_H_
#define _LS_STAG_INTERFACE_BLOCK_SOLVER_H_

// ================================================================== INCLUDES
#include "lu_decomposition_dense_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса InterfaceBlockSolver

Данный файл содержит в себе реализацию класса InterfaceBlockSolver
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий решение задачи на первой интерфейсной границе

Данный класс реализует построение LU-разложения для для блоков интерфейсной матрицы, соответствующих 
 первой интерфейсной границе, а также прямой и обратный ход метода Гаусса для них. Кроме того, данный 
 класс является базовым для классов, реализующих решение интерфейсной задачи на последующих интерфейсных границах
*/
	struct InterfaceBlockSolver
	{
		virtual void set(DenseMatrix* _L, DenseMatrix* _D, DenseMatrix* _U, 
			             double* _xb, double* _bb, double* _yb, double* _xpr, LUDecompositionDenseMatrix* _prevLU)
		{
			LU.set(_D);
			LU.setSolAndRight(_yb, _xb);
			m_xb = _xb;
			m_bb = _bb;
			m_yb = _yb;
			m_N = _D->m_N;
			capasity = m_N * sizeof(double);
		}
		/*!
		Выполняет LU-разложение для блоков интерфейсной матрицы, соответствующей данной интерфейсной границе
		*/
		virtual void buildLU()
		{// для первой подобласти
			LU.rebuild();
		}
		/*!
		Выполняет прямой ход метода Гаусса для блоков интерфейсной матрицы, соответствующей данной интерфейсной границе
		*/
		virtual void directGauss()
		{// для первой подобласти
			memcpy((void*)m_yb, (void*)m_bb, capasity);
			LU.solve();
		}
		/*!
		Выполняет обратный ход метода Гаусса для блоков интерфейсной матрицы, соответствующей данной интерфейсной границе
		*/
		virtual void inverseGauss()
		{// для первой подобласти
			prepareIG();
			LU.solve();
		}
	 protected:
		void prepareIG()
		{
			for(int i = 0; i < m_N; i++)
				m_yb[i] -= m_xb[i];
		}
	public:
		LUDecompositionDenseMatrix LU;///< LU-разложение заполненной матрицы, соответствующей блоку с главной диагонали
	 protected:
	  double *m_xb, *m_bb, *m_yb;
	  int m_N, capasity;
	};
}
#endif //_LS_STAG_INTERFACE_BLOCK_SOLVER_H_