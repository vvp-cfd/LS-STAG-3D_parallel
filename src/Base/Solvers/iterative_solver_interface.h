/*
* \copyright (C) 2015 Puzikova V.V.
* \file      iterative_solver_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      08.08.2015 0:58:13
* \brief     iterative_solver_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ITERATIVE_SOLVER_INTERFACE_H_
#define _LS_STAG_ITERATIVE_SOLVER_INTERFACE_H_

// ================================================================== INCLUDES
#include "../TaskDescription/base_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса IterativeSolverInterface

Данный файл содержит в себе реализацию класса IterativeSolverInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий интерфейс итерационного метода решения СЛАУ

Данный базовый класс реализует интерфейс итерационного метода решения СЛАУ (в т.ч. с использованием предобуславливания)
*/
	struct IterativeSolverInterface
	{
		IterativeSolverInterface():
		it(1),
		NM(0)
		{
			;
		}
		virtual ~IterativeSolverInterface(){};
		virtual void rebuildPreconditioner(){};
		/*!
		Решает СЛАУ итерационным методом
		\return it число итераций
		*/
		virtual int solve(){ return it; };
		virtual void set(SparseMatrix* _A, double* _rightM, double *_sol, BaseTaskDescription* _task)
		{
			m_A = _A;
			m_rightM = _rightM;
			m_sol = _sol;
			eps = _task->eps;
			max_it = _task->MAX_ITER;
			setPreconditioner(m_A, _task);
			NM = m_A->NMK;
		}
		virtual void print(const std::string &name){};
		virtual void printPreconditioner(const std::string &name){};
	  protected:
		virtual void setPreconditioner(SparseMatrix* _A, BaseTaskDescription* _task){};
		virtual void preconditioner(const double* _b, double* _x)
		{// если предобуславливателя нет, просто копируем
			LS_STAG_PAR_FOR(MeshFunction<double>, copy, NM, _b, _x)
		};
	  public:
		double eps;              ///< Точность решения 
		double norm_r;           ///< Норма невязки
		int it;                  ///< Число итераций
	  protected:
		SparseMatrix *m_A;       ///< Указатель на разреженную матрицу СЛАУ
		double *m_sol, *m_rightM;///< Указатели на векторы решения и правой части
		int max_it;              ///< Максимальное число итераций
		int NM;                  ///< Размерность СЛАУ
	};
};

#endif //_LS_STAG_ITERATIVE_SOLVER_INTERFACE_H_
