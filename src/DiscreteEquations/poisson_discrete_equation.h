/*
* \copyright (C) 2013 Puzikova V.V.
* \file      laplace_discrete_equation.h
* \author    Puzikova Valeria (vvp) 
* \date      04.11.2013 14:15:01
* \brief     poisson_discrete_equation.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_POISSON_DISCRETE_EQUATION_H_
#define _LS_STAG_POISSON_DISCRETE_EQUATION_H_

// ================================================================== INCLUDES
#include "discrete_equation_interface.h"
#include "../DiscreteOperators/WithFlux/Sym/laplace_discrete_operator.h"
#include "../DiscreteOperators/Divergence/divergence_discrete_operator_x.h"
#include "../DiscreteOperators/Divergence/divergence_discrete_operator_y.h"
#include "../DiscreteOperators/Divergence/divergence_discrete_operator_z.h"
#ifdef LS_STAG_USE_CUDA
	#ifdef  LS_STAG_USE_BICGSTAB
        #ifdef LS_STAG_FULL_ILU
#include "../Base/Solvers/BiCGSTab/ilu_bicgstab.h"
#else
		#include "../Base/Solvers/BiCGSTab/mg_bicgstab.h"
#endif
	#else
#ifdef LS_STAG_FULL_ILU
#include "../Base/Solvers/FGMRES/ilu_fgmres.h"
#else
		#include "../Base/Solvers/FGMRES/mg_fgmres.h"
#endif
	#endif
#elif  LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
#ifdef LS_STAG_FULL_ILU 
#include "../Base/Solvers/DDM/ddm_ilu_bicgstab.h"
#else
	#include "../Base/Solvers/DDM/ddm_mg_bicgstab.h" 
#endif
#else
#ifdef LS_STAG_FULL_ILU 
#include "../Base/Solvers/BiCGSTab/ilu_bicgstab.h"
#else
#include "../Base/Solvers/BiCGSTab/mg_bicgstab.h"
#endif
#endif
#else
#ifdef LS_STAG_USE_DDM
#ifdef LS_STAG_FULL_ILU 
#include "../Base/Solvers/DDM/ddm_ilu_fgmres.h"
#else
	#include "../Base/Solvers/DDM/ddm_mg_fgmres.h" 
#endif
#else
#ifdef LS_STAG_FULL_ILU
#include "../Base/Solvers/FGMRES/ilu_fgmres.h"
#else
#include "../Base/Solvers/FGMRES/mg_fgmres.h"
#endif
#endif
#endif
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса PoissonDiscreteEquation

Данный файл содержит в себе реализацию класса PoissonDiscreteEquation
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(fixPoint, LS_STAG_R_PD, m_right, const int, num, const bool, isFlow, const int, num2);
/*!
\brief Класс для работы с разностным аналогом уравнения Пуассона

Данный класс реализует работу с разностным аналогом уравнения Пуассона для поправки давления
*/
// ===================== PoissonDiscreteEquation (3D) =====================
struct PoissonDiscreteEquation : DiscreteEquationInterface
{
    typedef BaseTaskDescription::TaskKinds task;
    PoissonDiscreteEquation()
        : DiscreteEquationInterface((DiscreteOperatorWithFluxInterface*)&A, (IterativeSolverInterface*)&solver),
        existOutFlow(false) {}

    LS_STAG_SINGLE_KERNEL(fixPoint, LS_STAG_R_PD m_right, const int num, const bool isFlow, const int num2)
    {
        m_right[num] = 0.0;
        if (isFlow) m_right[num2] = 0.0;  // пины на углах выхода (E-S-B, E-N-B)
    }

private:
    void checkNeedRebuild() override
    {
        LS_STAG_PROFILE("rebuild");
        if (!needRebuild) needRebuild = m_meshes->movBody();
        if (needRebuild) {
            m_A->build();
            Dx.build();
            Dy.build();
            Dz.build();                // <-- новый дивергентный оператор по z
        }
    }

    void initCurMesh() override
    {
        m_cur_mesh = &m_meshes->baseMesh;  // основная сетка (центры ячеек)
    }
public:
    void set(BaseTaskDescription* _task, RealValuedMeshFunction* _sol) override
    {
        // p-решатель получает собственное окно рестарта FGMRES (FGMRES_RESTART_NUM_P):
        // на finest-сетках малое окно (13) вызывает стагнацию --- p не сходится за MAX_ITER
        // и расчёт расходится; u/v/w-решатели остаются с обычным окном
        const int savedRestart = _task->FGMRES_RESTART_NUM;
        _task->FGMRES_RESTART_NUM = _task->FGMRES_RESTART_NUM_P;
        DiscreteEquationInterface::set(_task, _sol);
        _task->FGMRES_RESTART_NUM = savedRestart;
        Dx.initMeshParams(m_cur_mesh, m_meshes);
        Dy.initMeshParams(m_cur_mesh, m_meshes);
        Dz.initMeshParams(m_cur_mesh, m_meshes);   // инициализация Dz

        aux.getCopy(m_aux);
        m_predUx = m_meshes->m_predUx;
        m_predUy = m_meshes->m_predUy;
        m_predUz = m_meshes->m_predUz;             // новая z-компонента прогнозной скорости

        aux.setMultipliers((SparseMatrix*)&Dx, m_predUx);
        rightMember.setMultipliers((SparseMatrix*)&Dy, m_predUy);

        // Добавляем член от w
        auxZ.init(NMK, m_auxZ);
        auxZ.setMultipliers((SparseMatrix*)&Dz, m_predUz);
        m_task = _task;
        yM = m_meshes->yMesh.M;
        zK = m_meshes->zMesh.K;
        bN = m_meshes->baseMesh.N;
        bNM = bN * m_meshes->baseMesh.M;
        if (m_task->id == task::TASK_FLOW)
            existOutFlow = true;
    }

private:
    void recalcRightMember() override
    {
        LS_STAG_PROFILE("recalcRM");
        checkNeedRebuild();
        // Dx * u_pred + Dy * v_pred
        aux.setProduct();
        rightMember.setProductWithSum(m_aux, m_flux);
        // добавляем Dz * w_pred
        auxZ.setProduct();
        rightMember += m_auxZ;
        // фиксация
        LS_STAG_SINGLE(PoissonDiscreteEquation, fixPoint, m_right, A.fixPointNum, A.isFlow, A.fixPointFlowNum)
        if (isMov)
            LS_STAG_PAR_FOR(LSSTAGUniformMesh, corr, NMK, m_right, m_type)
    }

public:
    LaplaceDiscreteOperator A;
    DivergenceDiscreteOperatorX Dx;
    DivergenceDiscreteOperatorY Dy;
    DivergenceDiscreteOperatorZ Dz;   // новый оператор
public:
#ifdef LS_STAG_USE_CUDA
    #ifdef LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_FULL_ILU
        ILUBiCGStab solver;
#else
		MGBiCGStab solver;
#endif
    #else
#ifdef LS_STAG_FULL_ILU
        ILUFGMRES solver;
#else
	MGFGMRES solver;
    #endif
#endif
#elif LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
#ifdef LS_STAG_FULL_ILU
        DDMILUBiCGStab solver;
#else
	DDMMGBiCGStab solver;
#endif
#else
#ifdef LS_STAG_FULL_ILU
	ILUBiCGStab solver;
#else
	MGBiCGStab solver;
#endif
#endif
#else
#ifdef LS_STAG_USE_DDM
#ifdef LS_STAG_FULL_ILU
        DDMILUFGMRES solver;
#else
	DDMMGFGMRES solver;
#endif
#else
#ifdef LS_STAG_FULL_ILU
	ILUFGMRES solver;
#else
	MGFGMRES solver;
#endif
#endif
#endif
private:
    RealValuedMeshFunction auxZ;
    double *m_aux, *m_predUx, *m_predUy, *m_predUz, *m_auxZ;
    int yM, zK, bN, bNM, num2;
    BaseTaskDescription* m_task;
    bool existOutFlow;
};
}

#endif //_LS_STAG_POISSON_DISCRETE_EQUATION_H_
