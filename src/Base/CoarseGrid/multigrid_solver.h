/*
* \copyright (C) 2013 Puzikova V.V.
* \file      multigrid_solver.h
* \author    Puzikova Valeria (vvp) 
* \date      30.10.2013 18:25:57
* \brief     multigrid_solver.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_MULTIGRID_SOLVER_H_
#define _LS_STAG_MULTIGRID_SOLVER_H_

// ================================================================== INCLUDES
#include "coarse_grids_hierarchy.h"
#include "../SparseMatrix/Decomposition/lu_decomposition_sparse_matrix.h"
#include "../TaskDescription/base_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса MultiGridSolver

Данный файл содержит в себе реализацию класса MultiGridSolver
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий многосеточный метод 

Данный класс реализует решение СЛАУ многосеточным методом, а также хранение и перестроение всех вспомогательных параметров
*/
	struct MultiGridSolver : CoarseGridsHierarchy
	{// n_post пока не берём из BaseTaskDescription
		~MultiGridSolver()
		{
			if (numLevMG > 0)
			{
				delete[] r;
				numLevMG = 0;
			}
		};
		/*!
		Выполняет решение СЛАУ многосеточным методом 
		\param _b указатель на вектор правой части
		\param _x указатель на вектор решения
		*/
		void solve(const double* _b, double* _x)
		{// n_pre пока нет
			{
				LS_STAG_PROFILE("copy");
				LS_STAG_PAR_FOR(MeshFunction<double>, copy, NMK, _b, m_r)
			}
			RealValuedMeshFunction* curR = r + 1;
			for (int k = 0; k < numLevMG; k++, curR++)
			{	
				curR->setProduct(); 
			}
			LU.solve();
			CoarseGrid* curGrid = grids + lastLevel;
			for(int k = lastLevel; k > -1; k--, curGrid--)
			{
				curNMK = curGrid->m_A->NMK;
				resize();
				ch.setProduct(&curGrid->P);
				c.copy(m_ch);
				for (int i = 0; i < n_post; i++)//постсглаживание
					curGrid->ADLJ.solve();
				cHH.copy(m_ch);
			}
			{
				LS_STAG_PROFILE("copy");
				LS_STAG_PAR_FOR(MeshFunction<double>, copy, NMK, m_ch, _x)
			}
		};
		void set(SparseMatrix* _A, BaseTaskDescription* _task)
		{
			CoarseGridsHierarchy::set(_A);
			LU.set((SparseMatrix*)&grids[lastLevel].AH);
			NMK = m_A->NMK;
			m_task = _task;
			n_pre = m_task->n_pre;
			n_post = m_task->n_post;
			n_min = m_task->N_MIN;
			c.init(NMK);
			ch.init(NMK, m_ch);
			cHH.init(NMK);
			q.init(NMK);
			r = new RealValuedMeshFunction[numLevMG + 1];
			r[0].init(NMK, m_r);
			for(int i = 0, i2 = 1; i < numLevMG; i++, i2++)
			{
				r[i2].init(grids[i].NMK);
				r[i2].setMultipliers((SparseMatrix*)&grids[i].R, &(r[i]));
			}
			q.setMultiplier(&c);
			ch.setMultiplier(&cHH);
			setSolverADLJ();
			LU.setSolAndRight(&(r[numLevMG]), &cHH);
			capasity = NMK * sizeof(double); 
		};
		void rebuild()
		{// 02.12.14
			CoarseGridsHierarchy::rebuild();
			LU.rebuild();//02.12.14
		};
		void print(const std::string &name) 
		{ 
			LU.print(name); //LU.printNotSparse(name); 
			for(int i = 0; i < numLevMG; i++)
				grids[i].ADLJ.print(name + std::to_string(i));
		}
	  private:
		void setSolverADLJ()
		{
			for (int i = 0; i < numLevMG; i++) 
				grids[i].ADLJ.setSolver(&r[i], &c, &q, &ch, n_min);
		};
		void resize()
		{
			LS_STAG_PROFILE("resize");
			c.resize(curNMK);
			ch.resize(curNMK);
			cHH.resize(curNMK);
			q.resize(curNMK);
		};
		LUDecompositionSparseMatrix LU;///< LU-декомпозиция разреженной матрицы системы на самом грубом сеточном уровне
	    RealValuedMeshFunction ch, q;
		extDouble c, cHH;
		RealValuedMeshFunction *r;//невязки и правые части
		BaseTaskDescription* m_task;
		double *m_ch, *m_r;
		int n_post, n_pre, NMK, n_min, curNMK;
		size_t capasity;
	};
}

#endif //_LS_STAG_MULTIGRID_SOLVER_H_
