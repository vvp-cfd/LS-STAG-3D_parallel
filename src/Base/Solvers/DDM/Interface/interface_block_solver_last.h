/*
* \copyright (C) 2019 Puzikova V.V.
* \file      interface_block_solver_last.h
* \author    Puzikova Valeria (vvp) 
* \date      09.08.2019 22:10:15
* \brief     interface_block_solver_last.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_INTERFACE_BLOCK_SOLVER_LAST_H_
#define _LS_STAG_INTERFACE_BLOCK_SOLVER_LAST_H_

// ================================================================== INCLUDES
#include "interface_block_solver.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса InterfaceBlockSolverLast

Данный файл содержит в себе реализацию класса InterfaceBlockSolverLast
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий решение задачи на последующих интерфейсных границах

Данный класс реализует построение LU-разложения для для блоков интерфейсной матрицы, соответствующих 
 интерфейсным границам начиная со второй, а также прямой и обратный ход метода Гаусса для блоков интерфейсной 
 матрицы, соответствующих последней интерфейсной границе
*/
	struct InterfaceBlockSolverLast : InterfaceBlockSolver
	{
		void set(DenseMatrix* _L, DenseMatrix* _D, DenseMatrix* _U, 
			             double* _xb, double* _bb, double* _yb, double* _xpr, LUDecompositionDenseMatrix* _prevLU)
		{
			InterfaceBlockSolver::set(_L, _D, _U, _xb, _bb, _yb, _xpr, _prevLU);
			m_L = _L;
			m_U = _U;
			m_xpr = _xpr;
			m_prevLU = _prevLU;
			m_CellL = m_L->Cell;
			m_CellU = m_U->Cell;
			m_CellD = _D->Cell;
			m_L->setMultProd(m_xpr, m_xb);
			m_U->setMultProd(m_xb, m_xpr);
		}
		void buildLU()
		{// для интерфейсов после первого
			for(int i = 0; i < m_N; i++)
			{
				for(int j = 0, num = i; j < m_N; j++, num += m_N)
					m_xb[j] = m_CellU[num];// правая часть - i-й столбец U
				m_prevLU->solve(m_xb);
				for(int j = 0, num = i, num2 = 0; j < m_N; j++, num += m_N)
				{
					double aux = 0.0;
					for(int k = 0; k < m_N; k++, num2++)
						aux += m_xpr[k] * m_CellL[num2];
					m_CellD[num] -= aux;
				}
			}
			LU.rebuild();
		}
		virtual void directGauss()
		{// для последней подобласти
			m_L->product();
			for(int i = 0; i < m_N; i++)
				m_yb[i] = m_bb[i] - m_xb[i];
		}
		virtual void inverseGauss()
		{// для последней подобласти
			LU.solve();
			m_U->product();
		}
	private:
	  DenseMatrix *m_L, *m_U;
	  LUDecompositionDenseMatrix* m_prevLU;      ///< Указатель на LU-разложение диагонального блока с предыдущей границы
	  double *m_xpr, *m_CellL, *m_CellU, *m_CellD;
	};
}
#endif //_LS_STAG_INTERFACE_BLOCK_SOLVER_LAST_H_