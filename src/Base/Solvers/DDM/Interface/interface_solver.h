/*
* \copyright (C) 2019 Puzikova V.V.
* \file      interfac_solver.h
* \author    Puzikova Valeria (vvp) 
* \date      10.08.2019 20:58:35
* \brief     interface_solver.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_INTERFACE_SOLVER_H_
#define _LS_STAG_INTERFACE_SOLVER_H_

// ================================================================== INCLUDES
#include "interface_block_solver_mid.h"
#include "interface_block_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса InterfaceSolver

Данный файл содержит в себе реализацию класса InterfaceSolver
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий решение интерфейсной задачи

Данный класс реализует построение LU-разложение блоков с главной диагонали интерфейсной матрицы и решение интерфейсной 
задачи с учетом ее блочной трехдиагональной структуры (т.е. прямым методом)
*/
	struct InterfaceSolver
	{
		InterfaceSolver():
		interfaceCount(0)
		{
			;
		}
		~InterfaceSolver()
		{
			if(interfaceCount)
			{
				for(int i = 0; i < interfaceCount; i++)
					delete interfaces[i];
				delete[] interfaces;
				delete[] y;
				interfaceCount = 0;
			}
		};
		void set(InterfaceBlockMatrix* _A, double* _b, double* _x)
		{
			interfaceCount = _A->iCount;
			iCount2 = interfaceCount - 1;
			if(interfaceCount)
			{
				interfaces = new InterfaceBlockSolver*[interfaceCount];
				interfaces[0] = new InterfaceBlockSolver;
				if(iCount2)
				{
					for(int i = 1; i < iCount2; i++)
						interfaces[i] = (InterfaceBlockSolver*)(new InterfaceBlockSolverMid);
					interfaces[iCount2] = (InterfaceBlockSolver*)(new InterfaceBlockSolverLast);
				}
				DenseMatrix *m_L = _A->L;
				DenseMatrix *m_D = _A->D;
				DenseMatrix *m_U = _A->U;
				int m_N = m_D->m_N;
				int NM = m_N * interfaceCount;
				y = new double[NM];
				double* curX = _x;
				double* curB = _b;
				double* curY = y;
				interfaces[0]->set((DenseMatrix*)NULL, m_D, (DenseMatrix*)NULL, curX, curB, curY, NULL, NULL);
				for(int i = 1; i < interfaceCount; i++, m_L++, m_U++)
				{
					double* prevX = curX;
					curX += m_N;
					curB += m_N;
					curY += m_N;
					m_D++;
					interfaces[i]->set(m_L, m_D, m_U, curX, curB, curY, prevX, &(interfaces[i-1]->LU));
				}
				rebuild();
			}
		};
		/*!
		Строит LU-декомпозицию блоков с главной диагонали интерфейсной матрицы
		*/
		void rebuild()
		{
			LS_STAG_PROFILE("rebuildInterface");
			for(int i = 0; i < interfaceCount; i++)
				interfaces[i]->buildLU();
		};
		/*!
		Выполняет решение интерфейсной системы
		*/
		void solve()
		{
			LS_STAG_PROFILE("InterfaceSolver");
			for(int i = 0; i < interfaceCount; i++)
				interfaces[i]->directGauss();
			if(iCount2)
				for(int i = iCount2; i > -1; i--)
					interfaces[i]->inverseGauss();
		};
		InterfaceBlockSolver** interfaces;///< Массив решателей для блоков, соответствующих интерфейсным границам
		double* y;
		int interfaceCount;///< Число интерфейсных границ между подобластями
		int iCount2;
	};
}
#endif //_LS_STAG_INTERFACE_SOLVER_H_