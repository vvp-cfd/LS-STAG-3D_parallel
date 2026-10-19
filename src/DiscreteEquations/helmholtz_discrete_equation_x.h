/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_equation_x.h
* \author    Puzikova Valeria (vvp) 
* \date      06.11.2013 14:46:18
* \brief     helmholtz_discrete_equation_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_X_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_X_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_equation.h"
#include "../DiscreteOperators/WithFlux/Sym/helmholtz_discrete_operator_x.h"
#include "../DiscreteOperators/Gradient/gradient_discrete_operator_x.h"
#include "../DiscreteOperators/WithFlux/SkewSym/convective_discrete_operator_x.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteEquationX

Данный файл содержит в себе реализацию класса HelmholtzDiscreteEquationX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(correctionX, LS_STAG_R_PD, _right, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _Number, LS_STAG_CR_PCT, _type,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _theta_N, LS_STAG_CR_PC, _curV_C, LS_STAG_CR_PI, _CellOnAirfoill, const int, _hs,
		const int, _hf, const int, _N, const int, _M, const int, _NM);
/*!
\brief Класс для работы с разностным аналогом уравнения Гельмгольца на х-й сетке

Данный класс реализует работу с разностным аналогом уравнения Гельмгольца для прогноза проекции скорости на горизонталь
*/
// ===================== HelmholtzDiscreteEquationX (3D) =====================
struct HelmholtzDiscreteEquationX : HelmholtzDiscreteEquation
{
    HelmholtzDiscreteEquationX()
        : HelmholtzDiscreteEquation((DiscreteOperatorWithFluxInterface*)&Kx,
              (GradientDiscreteOperatorBase*)&DxT,
              (SkewSymDiscreteOperator*)&Cx,
              nullptr, nullptr) {}

    LS_STAG_KERNEL_2D(correctionX, LS_STAG_R_PD _right, LS_STAG_R_PD _Cell, LS_STAG_CR_PI _Number,
        LS_STAG_CR_PCT _type, LS_STAG_CR_PC _theta, LS_STAG_CR_PC _theta_N,
        LS_STAG_CR_PC _curV_C, LS_STAG_CR_PI _CellOnAirfoill, const int _hs, const int _hf, const int _N, const int _M, const int _NM)
    {
        for (int k = i * _N + j * _NM, sym_num = k + _hs, m_end = k + _hf, num = sym_num + i + j * _M; sym_num < m_end; num++, sym_num++)
        {
            if (_type[sym_num] > SolidCell)
            {
                if ((_theta_N[num].y < EQ_DOUBLE) && (_theta[num].y > EQ_DOUBLE))
                {
                    const int airNum = LS_STAG_VIRTUAL_INLINE_CALL(CLASS_NAME_AIRFOIL_NUM, getAirfoilNum, num, _CellOnAirfoill);
                    _right[sym_num] = _curV_C[airNum].x;
                    HelmholtzDiscreteEquation::_k_correctionOnImB(sym_num, _Cell, _Number);
                }
            }
        }
    }

private:
    void initCurMesh() override
    {
        m_cur_mesh = &m_meshes->xMesh;
        m_prev_sol = m_meshes->m_prevUx;
        m_lastSol  = m_meshes->m_Ux;
    }

    void correctionRightMemberImB() override
    {
        LS_STAG_PROFILE("correctionX");
        LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteEquationX, correctionX, m_vs, m_vf, m_zs, m_zf,
            m_right, m_A->Cell, m_A->Number, m_type,
            m_theta, m_theta_N, m_meshes->curV_C, m_task->getPCellOnAirfoill(),
            m_hs, m_hf, m_cur_mesh->N, m_cur_mesh->M, m_cur_mesh->M * m_cur_mesh->N)
    }

    HelmholtzDiscreteOperatorX Kx;
    GradientDiscreteOperatorX DxT;
    ConvectiveDiscreteOperatorX Cx;
}; 
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_X_H_
