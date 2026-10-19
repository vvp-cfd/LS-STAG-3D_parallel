/*
* \copyright (C) 2019 Puzikova V.V.
* \file      ddm_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2019 21:07:29
* \brief     ddm_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DDM_INTERFACE_H_
#define _LS_STAG_DDM_INTERFACE_H_

// ================================================================== INCLUDES
#include "../iterative_solver_interface.h"
#include "SubDomains/domain_decompositor.h"
#include "Interface/interface_solver.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DDMInterface

Данный файл содержит в себе реализацию класса DDMInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий метод декомпозиции области

Данный базовый класс реализует интерфейс метода декомпозиции области
*/
	struct DDMInterface : IterativeSolverInterface
	{
		DDMInterface():
		IterativeSolverInterface(),
		domainCount(0),
		coef(1.0)
		{
			;
		}
		~DDMInterface()
		{
			if(domainCount)
			{
				for(int i = 0; i < domainCount; i++) delete domainSolvers[i];
				delete [] domainSolvers;
				domainCount = 0;
			}
		};
		void set(SparseMatrix* _A, double* _rightM, double *_sol, BaseTaskDescription* _task)
		{
			IterativeSolverInterface::set(_A, _rightM, _sol, _task);
			domainCount = _task->NUM_SUBDOMAINS;
			domainCount2 = domainCount - 1;
			m_N = _A->m_N * _A->m_M;
			if (domainCount2)
			{
				decompositor.set(_A, _sol, _rightM, domainCount);//build();
				w.init(m_A->NMK);
				w.setOperator(m_A);
				subTasks = decompositor.subTasks;
				interTask = &(decompositor.interTask);
			}
			createDomainSolvers();
			if (domainCount2)
			{
				{
					LS_STAG_PROFILE("domainSolversSet");
#pragma omp parallel for
					for (int i = 0; i < domainCount; i++)
						domainSolvers[i]->set(subTasks[i]->App, subTasks[i]->right, subTasks[i]->sol, _task);
				}
				buildInterface();
				interfaceSolver.set(&(interTask->Ass), interTask->right, interTask->sol);
			}
			else domainSolvers[0]->set(_A, _rightM, _sol, _task);
		};
		/*!
		Решает СЛАУ методом декомпозиции области
		\return it среднее число итераций по всем подобластям
		*/
		int solve()
		{
			if (domainCount < 2) 
			{
			    it = domainSolvers[0]->solve();
				norm_r = domainSolvers[0]->norm_r;
				return it;
			}
			it = 0;
			norm_r = w.setProductWithResidual(m_sol, m_rightM);
			if (norm_r < eps) return it;
			interTask->initInterface();// подготовка bp, bs
			int m_it = 0;
			{
				LS_STAG_PROFILE("subdomainsSolve");
#pragma omp parallel for reduction(+:m_it)
				for (int i = 0; i < domainCount; i++)
				{// расчет \tilde{bs}
					subTasks[i]->initSubDomain();
					m_it += domainSolvers[i]->solve();
					subTasks[i]->recalcInterfaceBHS();
				}
			}
			interfaceSolver.solve();// 10/08/19 it += interfaceSolver->solve();// решение интерфейсной задачи
			{
				LS_STAG_PROFILE("subdomainsSolve");
#pragma omp parallel for reduction(+:m_it)
				for (int i = 0; i < domainCount; i++)
				{
					subTasks[i]->recalcRHS();// подготовка правой части для каждой подобласти
					m_it += domainSolvers[i]->solve();
					subTasks[i]->sendSubDomain();
				}
			}
			interTask->sendInterface();// собираем решение по подобластям	
			// вообще, норму надо считать хотя бы как среднее по подобластям
			norm_r = w.setProductWithResidual(m_sol, m_rightM);// пока так, но если матрица хранится распределенно, изменить!
			it += ((int)((m_it + 1) * coef));// 10/08/19 it = ((int)(it * coef));
			return it;// среднее число итераций по всем подобластям
		}
		void rebuildPreconditioner()
		{
			if (domainCount2)
			{
				{
					LS_STAG_PROFILE("subdomainsRebuild");
#pragma omp parallel for
					for (int i = 0; i < domainCount; i++)
					{
						subTasks[i]->rebuild();
						domainSolvers[i]->rebuildPreconditioner();
					}
				}
				interTask->rebuild();
				buildInterface();
				interfaceSolver.rebuild();//10/08/19 interfaceSolver->rebuildPreconditioner();
			}
			else domainSolvers[0]->rebuildPreconditioner();
		};
	  protected:
		virtual IterativeSolverInterface* getNewSolver(){ return ((IterativeSolverInterface*)NULL); } ;
		/*!
		Достраивает интерфейсную матрицу при помощи решателей из подобластей
		*/
		void buildInterface()
		{// достроение интерфейсной матрицы при помощи решателей из подобластей
			LS_STAG_PROFILE("buildInterface");
#pragma omp parallel for
			for(int i = 1; i < domainCount; i++) 
			{
				for(int j = 0; j < m_N; j++)
				{
					subTasks[i]->setSouthRHS(j);
					domainSolvers[i]->solve();
					subTasks[i]->getSouthInterface();//10/08/19
				}
			}
#pragma omp parallel for
			for(int i = 0; i < domainCount2; i++) 
			{
				for(int j = 0; j < m_N; j++)
				{
					subTasks[i]->setNorthRHS(j);
					domainSolvers[i]->solve();
					subTasks[i]->getNorthInterface();//10/08/19
				}
			}
		}
		/*!
		Создает решатели для задач в подобластях
		*/
		void createDomainSolvers()
		{
			if(domainCount) 
			{
				domainSolvers = new IterativeSolverInterface*[domainCount];
				coef = 1.0 / (1 + 2 * domainCount);
			}
			for(int i = 0; i < domainCount; i++) domainSolvers[i] = getNewSolver();
		};
		RealValuedMeshFunction w;
		DomainDecompositor decompositor;///< Задачи для подобластей и интерфейса пока храним тут же, потом перенести
		double coef;
		IterativeSolverInterface** domainSolvers;///< Массив решателей для задач в подобластях
		InterfaceSolver interfaceSolver;///< Решатель для интерфейсной задачи
		ShurSubdomain **subTasks;///< Массив задач в подобластях
		ShurInterfaceTask *interTask;///< Указатель на интерфейсную задачу
		int domainCount;       ///< Число подобластей
		int domainCount2, m_N;
	};
};

#endif //_LS_STAG_DDM_INTERFACE_H_
