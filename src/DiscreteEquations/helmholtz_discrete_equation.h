/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_equation.h
* \author    Puzikova Valeria (vvp) 
* \date      04.11.2013 15:42:49
* \brief     helmholtz_discrete_equation.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_H_

// ================================================================== INCLUDES
#include "discrete_equation_interface.h"
#include "../DiscreteOperators/WithFlux/Sym/helmholtz_discrete_operator_base.h"
#ifdef LS_STAG_FULL_MG
#ifdef LS_STAG_USE_CUDA
#ifdef  LS_STAG_USE_BICGSTAB
#include "../Base/Solvers/BiCGSTab/mg_bicgstab.h"
#else
#include "../Base/Solvers/FGMRES/mg_fgmres.h"
#endif
#elif  LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
#include "../Base/Solvers/DDM/ddm_mg_bicgstab.h"
#else
#include "../Base/Solvers/BiCGSTab/mg_bicgstab.h"
#endif
#else
#ifdef LS_STAG_USE_DDM
#include "../Base/Solvers/DDM/ddm_mg_fgmres.h" 
#else
#include "../Base/Solvers/FGMRES/mg_fgmres.h"
#endif
#endif
#else
#ifdef LS_STAG_USE_CUDA
	#ifdef  LS_STAG_USE_BICGSTAB
		#include "../Base/Solvers/BiCGSTab/ilu_bicgstab.h"
	#else
		#include "../Base/Solvers/FGMRES/ilu_fgmres.h"
	#endif
#elif  LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
	#include "../Base/Solvers/DDM/ddm_ilu_bicgstab.h"
#else
#include "../Base/Solvers/BiCGSTab/ilu_bicgstab.h"
#endif
#else
#ifdef LS_STAG_USE_DDM
	#include "../Base/Solvers/DDM/ddm_ilu_fgmres.h" // #include "../Base/Solvers/DDM/ddm_mg_fgmres.h"
#else
#include "../Base/Solvers/FGMRES/ilu_fgmres.h"
#endif
#endif
#endif
#include "../DiscreteOperators/Gradient/gradient_discrete_operator_base.h"
#include "../DiscreteOperators/WithFlux/SkewSym/skew_sym_discrete_operator.h"
#include "../Base/LS_STAG_Mesh/ls_stag_mesh_mov.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteEquation

Данный файл содержит в себе реализацию класса HelmholtzDiscreteEquation
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(initRM, LS_STAG_R_PD, _right, LS_STAG_CR_PD, _lastSol, LS_STAG_CR_PD, _prev_sol,
		LS_STAG_CR_PD, _VOF, const double, _coef2);
/*!
\brief Класс для работы с разностным аналогом уравнения Гельмгольца

Данный класс реализует работу с разностным аналогом уравнения Гельмгольца для прогноза скорости
*/
// ===================== HelmholtzDiscreteEquation (3D) =====================
struct HelmholtzDiscreteEquation : DiscreteEquationInterface
{
    explicit HelmholtzDiscreteEquation(
        DiscreteOperatorWithFluxInterface* _A = nullptr,
        GradientDiscreteOperatorBase* _grad = nullptr,
        SkewSymDiscreteOperator* _conv = nullptr,
        DiscreteOperatorInterface* _div1 = nullptr,   // первый shear (τ_xy или τ_xz)
        DiscreteOperatorInterface* _div2 = nullptr)   // второй shear (τ_yz)
        : DiscreteEquationInterface(_A, (IterativeSolverInterface*)&solver),
          m_grad(_grad), m_conv(_conv), m_div1(_div1), m_div2(_div2),
          useTM(false), useNNM(false), useFOS(false), useVEF(false), useGNF(false), useTDF(false) {}

    LS_STAG_KERNEL(initRM, LS_STAG_R_PD _right, LS_STAG_CR_PD _lastSol, LS_STAG_CR_PD _prev_sol,
        LS_STAG_CR_PD _VOF, const double _coef2)
    {
        _right[i] = (-2.0 * _lastSol[i] + 0.5 * _prev_sol[i]) * _VOF[i] * _coef2;
    }

    void solve() override
    {
        LS_STAG_PROFILE("helmSolve");
        DiscreteEquationInterface::solve();
        initRightMemberIb();
    }
    virtual void setDiffusionGNF(double* _fd)
    {
        m_fd = _fd;
        useNNM = true;
        useGNF  = true;
    }
    void setReynoldsStresses(double* _normalTS, RealValuedMeshFunction* _shearTS1, RealValuedMeshFunction* _shearTS2)
    {
        m_normalTS = _normalTS;
        _shearTS1->getCopy(m_shearTS1);
         _shearTS2->getCopy(m_shearTS2);
        useTM = true;
    }
    void setNonNewtonianStresses(double* _normalNNS, double* _shearNNS1, double* _shearNNS2)
    {
        m_normalNNS = _normalNNS;
        m_shearNNS1 = _shearNNS1;
        m_shearNNS2 = _shearNNS2;
        useNNM = true;
        useVEF = true;
    }
    void set(BaseTaskDescription* _task, RealValuedMeshFunction* _sol) override
    {
        m_task = _task;
        DiscreteEquationInterface::set(_task, _sol);
        m_VOF = m_cur_mesh->VOF;
        m_grad->initMeshParams(m_cur_mesh, m_meshes);
        m_conv->init(_task, m_cur_mesh);
        if (useTM || useVEF) {
            m_div1->initMeshParams(m_cur_mesh, m_meshes);
            m_div2->initMeshParams(m_cur_mesh, m_meshes);
        }
        m_flux2 = m_conv->flux;
        m_P = m_meshes->m_P;
        aux.setMultipliers((SparseMatrix*)m_conv, m_lastSol);
        aux.set(0.0);
        aux2.init(NMK, m_aux2);
        aux2.setMultipliers((SparseMatrix*)m_grad, m_P);
        useFOS = (isMov || (useNNM && (!useTDF)));
        if (m_meshes->id == LSSTAGMeshMov::MeshKinds::MESH_MOVING) {
            LSSTAGMeshMov* aux_mov = (LSSTAGMeshMov*)m_meshes;
            m_hs = aux_mov->hsMov;  m_hf = aux_mov->hfMov;
            m_vs = aux_mov->vsMov;  m_vf = aux_mov->vfMov;
            m_zs = aux_mov->zsMov;  m_zf = aux_mov->zfMov;
            m_theta   = aux_mov->theta;
            m_theta_N = aux_mov->theta_N;
        }
        coef2 = 1.0 / _task->t;
        if (useFOS) {
            coef = coef2;
            rightMember.setComponentwiseMultiplier(m_lastSol);
            if (useTM) {
                prev_aux.init(NMK, m_prev_aux);
                aux3.init(NMK, m_aux3);
                aux5.init(NMK, m_aux5);
                prev_aux.setMultipliers((SparseMatrix*)m_grad, m_normalTS);
                aux3.setMultipliers((SparseMatrix*)m_div1, m_shearTS1);
                aux5.setMultipliers((SparseMatrix*)m_div2, m_shearTS2);
            }
            if (useVEF) {
                prev_aux.init(NMK, m_prev_aux);
                aux3.init(NMK, m_aux3);
                aux5.init(NMK, m_aux5);
                prev_aux.setMultipliers((SparseMatrix*)m_grad, m_normalNNS);
                aux3.setMultipliers((SparseMatrix*)m_div1, m_shearNNS1);
                aux5.setMultipliers((SparseMatrix*)m_div2, m_shearNNS2);
            }
        } else {
            coef = 1.5 * coef2;
            prev_aux.init(NMK, m_prev_aux);
            prev_aux.copy(m_lastSol);
            LS_STAG_PAR_FOR(MeshFunction<double>, copy, NMK, m_prev_sol, m_lastSol)
            m_conv->recalc();
            aux.setProductWithSum(m_flux2);
            LS_STAG_PAR_FOR(MeshFunction<double>, copy, NMK, m_prev_aux, m_lastSol)
            prev_aux.copy(m_aux);
            if (useTM) {
                aux4.init(NMK, m_aux4);
                aux3.init(NMK, m_aux3);
				aux5.init(NMK, m_aux5);
                aux4.setMultipliers((SparseMatrix*)m_grad, m_normalTS);
                aux3.setMultipliers((SparseMatrix*)m_div1, m_shearTS1);;
                aux5.setMultipliers((SparseMatrix*)m_div2, m_shearTS2);
            }
            if (useVEF) {
                aux4.init(NMK, m_aux4);
                aux3.init(NMK, m_aux3);
                aux5.init(NMK, m_aux5);
                aux4.setMultipliers((SparseMatrix*)m_grad, m_normalNNS);
                aux3.setMultipliers((SparseMatrix*)m_div1, m_shearNNS1);
                aux5.setMultipliers((SparseMatrix*)m_div2, m_shearNNS2);
            }
        }
        initRightMemberIb();
    }
protected:
    LS_STAG_KERNEL(correctionOnImB, LS_STAG_R_PD _Cell, LS_STAG_CR_PI _Number)
    {
        const int numA2 = _Number[i];
        _Cell[numA2] = 0.0;
        _Cell[numA2 + 1] = 0.0;
        _Cell[numA2 + 2] = 0.0;
        _Cell[numA2 + 3] = 1.0;
        _Cell[numA2 + 4] = 0.0;
        _Cell[numA2 + 5] = 0.0;
        _Cell[numA2 + 6] = 0.0;
    }
    virtual void correctionRightMemberImB() {}

private:
    void recalcRightMember() override
    {
        LS_STAG_PROFILE("recalcRM");
        aux2.setProduct();
        if (useTM || useVEF) {
            if (useFOS) prev_aux.setProduct();
            else aux4.setProduct();
            aux3.setProduct();
            aux5.setProduct();
        }
        checkNeedRebuild();
        aux.setProductWithSum(m_flux2);
        if (useFOS) {
            if (useGNF) {
                LS_STAG_PAR_FOR(LSSTAGUniformMesh, corr, NMK, m_fd, m_type)
                rightMember.componentwiseProductAddAndSubtract(m_aux, m_flux, m_aux2, m_fd);
            } else
                rightMember.componentwiseProductAddAndSubtract(m_aux, m_flux, m_aux2);
            if (useTM || useVEF) {
                rightMember.addAndSubtract(m_prev_aux, m_aux3);
                rightMember -= m_aux5;   // второе касательное напряжение
            }
            if (isMov) {
                LS_STAG_PROFILE("correction");
                LS_STAG_PAR_FOR(LSSTAGUniformMesh, corr, NMK, m_right, m_type)
                correctionRightMemberImB();
            }
        } else {
            {
                LS_STAG_PROFILE("initRM");
                LS_STAG_PAR_FOR(HelmholtzDiscreteEquation, initRM, NMK, m_right, m_lastSol, m_prev_sol, m_VOF, coef2)
            }
            rightMember.subtract(m_aux2, m_flux, m_prev_aux);
            prev_aux.copy(m_aux);
            aux *= 2.0;
            rightMember += m_aux;
            if (useTM || useVEF) {
                rightMember.addAndSubtract(m_aux4, m_aux3);
                rightMember -= m_aux5;
            }
        }
    }

    void checkNeedRebuild() override
    {
        LS_STAG_PROFILE("rebuildHelm");
        DiscreteEquationInterface::checkNeedRebuild();
        if (needRebuild) {
            LS_STAG_PROFILE("rebuild");
            m_conv->build();
            m_grad->build();
            if (useTM || useVEF) {
                m_div1->build();
                m_div2->build();
            }
        } else
            m_conv->recalc();
    }

    bool getNeedRebuildSolver() override { return (needRebuild || useGNF || useTDF); }
    void ifNotRebuild() override
    {
        if (useGNF) m_A->recalc();
        else if (useTDF) m_A->recalcTDF();
        else DiscreteEquationInterface::ifNotRebuild();
    }
    void initRightMemberIb()
    {
        if (useFOS) LS_STAG_PAR_FOR(extDouble, copy, NMK, m_VOF, m_right, -coef)
    }

public:
#ifdef LS_STAG_FULL_MG
#ifdef LS_STAG_USE_CUDA
#ifdef LS_STAG_USE_BICGSTAB
	MGBiCGStab solver;
#else
	MGFGMRES solver;
#endif
#elif LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
	DDMMGBiCGStab solver;
#else
	MGBiCGStab solver;
#endif
#else
#ifdef LS_STAG_USE_DDM
	DDMMGFGMRES solver;
#else
	MGFGMRES solver;
#endif
#endif
#else
#ifdef LS_STAG_USE_CUDA
    #ifdef LS_STAG_USE_BICGSTAB
        ILUBiCGStab solver;
    #else
        ILUFGMRES solver;
    #endif
#elif LS_STAG_USE_BICGSTAB
#ifdef LS_STAG_USE_DDM
        DDMILUBiCGStab solver;
#else
	ILUBiCGStab solver;
#endif
#else
#ifdef LS_STAG_USE_DDM
        DDMILUFGMRES solver;
#else
	ILUFGMRES solver;
#endif
#endif
#endif
    RealValuedMeshFunction aux2, prev_aux, aux3, aux5, aux4;
    GradientDiscreteOperatorBase* m_grad;
    SkewSymDiscreteOperator* m_conv;
    DiscreteOperatorInterface* m_div1;// shear1
    DiscreteOperatorInterface* m_div2;// shear2   
    double *m_normalTS, *m_shearTS1, *m_shearTS2;
    double *m_normalNNS, *m_shearNNS1, *m_shearNNS2;
    double coef;
protected:
    double *m_prev_sol, *m_lastSol;
    coordinate *m_theta_N, *m_theta;
    BaseTaskDescription* m_task;
    int m_hs, m_hf, m_vs, m_vf, m_zs, m_zf;
    bool useTM, useNNM, useFOS, useVEF, useGNF, useTDF;
private:
    double *m_fd, *m_aux2, *m_flux2, *m_prev_aux, *m_aux3, *m_aux5, *m_aux4;
    double *m_VOF, *m_P;
    double coef2;
};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_H_
