/*
* \copyright (C) 2026 Puzikova V.V.
* \file      helmholtz_discrete_equation_z.h
* \author    Puzikova Valeria (vvp) 
* \date      11.05.2026 19:04:35
* \brief     helmholtz_discrete_equation_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Z_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Z_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_equation.h"
#include "../DiscreteOperators/WithFlux/Sym/helmholtz_discrete_operator_z.h"
#include "../DiscreteOperators/Gradient/gradient_discrete_operator_z.h"
#include "../DiscreteOperators/WithFlux/SkewSym/convective_discrete_operator_z.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteEquationZ

Данный файл содержит в себе реализацию класса HelmholtzDiscreteEquationZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(correctionZ, LS_STAG_R_PD, _right, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _Number, LS_STAG_CR_PCT, _type,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _theta_N, LS_STAG_CR_PC, _curV_C, LS_STAG_CR_PI, _CellOnAirfoill, const int, _hs,
		const int, _hf, const int, _N, const int, _NM);
/*!
\brief Класс для работы с разностным аналогом уравнения Гельмгольца на z-й сетке

Данный класс реализует работу с разностным аналогом уравнения Гельмгольца для прогноза проекции скорости на Oz
*/
// ===================== HelmholtzDiscreteEquationZ (3D) =====================
struct HelmholtzDiscreteEquationZ : HelmholtzDiscreteEquation
{
    HelmholtzDiscreteEquationZ()
        : HelmholtzDiscreteEquation((DiscreteOperatorWithFluxInterface*)&Kz,
              (GradientDiscreteOperatorBase*)&DzT,
              (SkewSymDiscreteOperator*)&Cz,
              nullptr, nullptr) {}

    LS_STAG_KERNEL_2D(correctionZ, LS_STAG_R_PD _right, LS_STAG_R_PD _Cell, LS_STAG_CR_PI _Number, LS_STAG_CR_PCT _type, LS_STAG_CR_PC _theta, 
		LS_STAG_CR_PC _theta_N, LS_STAG_CR_PC _curV_C, LS_STAG_CR_PI _CellOnAirfoill, const int _hs, const int _hf, const int _N, const int _NM)
    {
        for (int k = i * _N + j * _NM, num = k + _hs, m_end = k + _hf; num < m_end; ++num) {
            if (_type[num] > SolidCell) {
                if ((_theta_N[num].z < EQ_DOUBLE) && (_theta[num].z > EQ_DOUBLE)) {
                    const int airNum = LS_STAG_VIRTUAL_INLINE_CALL(CLASS_NAME_AIRFOIL_NUM, getAirfoilNum, num, _CellOnAirfoill);
                    _right[num] = _curV_C[airNum].z;
                    HelmholtzDiscreteEquation::_k_correctionOnImB(num, _Cell, _Number);
                }
            }
        }
    }
private:
    void initCurMesh() override
    {
        m_cur_mesh = &m_meshes->zMesh;
        m_prev_sol = m_meshes->m_prevUz;
        m_lastSol  = m_meshes->m_Uz;
    }
    void correctionRightMemberImB() override
    {
        LS_STAG_PROFILE("correctionZ");
        LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteEquationZ, correctionZ, m_vs, m_vf, m_zs, m_zf-1, m_right, m_A->Cell, m_A->Number, m_type,
            m_theta, m_theta_N, m_meshes->curV_C, m_task->getPCellOnAirfoill(), m_hs, m_hf, m_cur_mesh->N, m_cur_mesh->N * m_cur_mesh->M)
    }

    HelmholtzDiscreteOperatorZ Kz;
    GradientDiscreteOperatorZ DzT;
    ConvectiveDiscreteOperatorZ Cz;
};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Z_H_
