/*
* \copyright (C) 2013 Puzikova V.V.
* \file      discrete_equation_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      31.10.2013 20:05:08
* \brief     discrete_equation_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DISCRETE_EQUATION_INTERFACE_H_
#define _LS_STAG_DISCRETE_EQUATION_INTERFACE_H_

// ================================================================== INCLUDES
#include "../DiscreteOperators/WithFlux/discrete_operator_with_flux_interface.h"
#include "../Base/Solvers/iterative_solver_interface.h"
#include "../Base/MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DiscreteEquationInterface

Данный файл содержит в себе реализацию класса DiscreteEquationInterface
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(sym, LS_STAG_R_PD, _sol, const int, _N, const int, _NM, const int, M2);
/*!
\brief Базовый класс для работы с разностными аналогами уравнений

Данный базовый класс-интерфейс реализует работу с разностными аналогами уравнений
*/
	struct DiscreteEquationInterface
	{
		explicit DiscreteEquationInterface(DiscreteOperatorWithFluxInterface* _A = (DiscreteOperatorWithFluxInterface*)0, IterativeSolverInterface* _solver = (IterativeSolverInterface*)0):
		m_A(_A),
		m_solver(_solver),
		needRebuild(false),
		isMov(false),
		needSym(false)
		{
			;
		}
		LS_STAG_KERNEL_2D(sym, LS_STAG_R_PD _sol, const int _N, const int _NM, const int M2)
		{
			for (int numS = j * _NM, num_sym = (M2 - i) * _N + numS, num = i * _N + numS, m_end = num + _N; num < m_end; num++, num_sym++)
				_sol[num_sym] = _sol[num];
		}
	  protected:
		  /*!
		  Производит пересчет правой части
		  */
		virtual void recalcRightMember(){};
		  /*!
		  Проверяет, присутствует ли движение погруженных границ и надо ли перестроить матрицу СЛАУ  
		  */
		virtual void checkNeedRebuild()
		{
			LS_STAG_PROFILE("baseEqRebuild");
			if(!needRebuild) needRebuild = m_meshes->movBody();
			if(needRebuild)
				m_A->build();
			else 
				ifNotRebuild();
		};
		virtual void ifNotRebuild() { m_A->calcFlux(); };
		virtual bool getNeedRebuildSolver() { return needRebuild; };
		virtual void initCurMesh(){};
		virtual void symmetric()
		{
			LS_STAG_PAR_FOR_2D(DiscreteEquationInterface, sym, center_row, m_cur_mesh->K, m_sol, m_cur_mesh->N, m_cur_mesh->M * m_cur_mesh->N, m_cur_mesh->M - 1)
		};
	  public:
		  /*!
		  Решает разностное уравнение на текущем шаге по времени
		  */
		virtual void solve()
		{
			LS_STAG_PROFILE("solve");
			recalcRightMember();
			if(getNeedRebuildSolver()) m_solver->rebuildPreconditioner();// 24.11.14
			m_solver->solve();
			if (needSym) symmetric();
		};
		virtual void set(BaseTaskDescription* _task, RealValuedMeshFunction* _sol)
		{
			_sol->getCopy(m_sol);
			m_meshes = _task->m_mesh;
			initCurMesh();
			m_A->init(_task, m_cur_mesh);// + build   // 0419 --- по-хорошему здесь задачи по подобластям
			NMK = m_A->NMK;
			rightMember.init(NMK, m_right);
			m_solver->set(m_A, m_right, m_sol, _task);
			aux.init(NMK, m_aux);
			m_flux = m_A->flux;
			m_type = m_cur_mesh->typeOfCell;
			if((_task->id_motion == BaseTaskDescription::MotionKinds::MOTION_MOV) && 
				(_task->id == BaseTaskDescription::TaskKinds::TASK_FLOW))
			{
				FlowTaskDescription* f_task = (FlowTaskDescription*)_task;
				if((f_task->id_airfoil == FlowTaskDescription::AirfoilKind::AIRFOIL_ROTOR) ||
					(f_task->id_airfoil == FlowTaskDescription::AirfoilKind::AIRFOIL_POINTS))
					needRebuild = true;
			}
			if(_task->id_motion == BaseTaskDescription::MotionKinds::MOTION_MOV) isMov = true;
			if((_task->id == BaseTaskDescription::TaskKinds::TASK_FLOW) && (!(_task->discontin)) && 
				(_task->id_motion == BaseTaskDescription::MotionKinds::MOTION_NO) && (_task->Re > RE_CRIT) && (!(_task->useTM)))
			{// включить только если нужно смоделировать стационарное "базовое" неустойчивое течение
				// needSym = true;// нужна ли симметризация решения
				center_row = (int)floor(0.5 * m_cur_mesh->M);
			}
		};
		RealValuedMeshFunction rightMember;///< Сеточная функция --- правая часть уравнения
		RealValuedMeshFunction aux;
		DiscreteOperatorWithFluxInterface* m_A;///< Указатель на дискретный оператор, дающий матрицу СЛАУ
		IterativeSolverInterface* m_solver;///< Указатель на базовую структуру итерационного метода решения СЛАУ
		LSSTAGUniformMesh* m_meshes;///< Указатель на структуру, в которой хранятся данные о всех разнесенных сетках LS-STAG-сетки
		LSSTAGBaseMeshParams* m_cur_mesh;///< Указатель на структуру с параметрами сетки, на которой строится разностный аналог уравнения
		CellTypes *m_type;
		double* m_sol;
		double* m_right;
	  protected:
		double* m_flux;
		double* m_aux;
		int NMK, center_row;
		bool needRebuild, isMov, needSym;
	};
}

#endif //_LS_STAG_DISCRETE_EQUATION_INTERFACE_H_
