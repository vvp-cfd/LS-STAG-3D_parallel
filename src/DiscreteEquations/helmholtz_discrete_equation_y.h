/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_equation_y.h
* \author    Puzikova Valeria (vvp) 
* \date      06.11.2013 14:53:41
* \brief     helmholtz_discrete_equation_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Y_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Y_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_equation.h"
#include "../DiscreteOperators/WithFlux/Sym/helmholtz_discrete_operator_y.h"
#include "../DiscreteOperators/Gradient/gradient_discrete_operator_y.h"
#include "../DiscreteOperators/WithFlux/SkewSym/convective_discrete_operator_y.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteEquationY

Данный файл содержит в себе реализацию класса HelmholtzDiscreteEquationY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(symY, LS_STAG_R_PD, _sol, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_R_KERNEL_2D(correctionY, LS_STAG_R_PD, _right, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _Number, LS_STAG_CR_PCT, _type,
		LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _theta_N, LS_STAG_CR_PC, _curV_C, LS_STAG_CR_PI, _CellOnAirfoill, const int, _hs,
		const int, _hf, const int, _N, const int, _NM);
/*!
\brief Класс для работы с разностным аналогом уравнения Гельмгольца на у-й сетке

Данный класс реализует работу с разностным аналогом уравнения Гельмгольца для прогноза проекции скорости на вертикаль
*/
// ===================== HelmholtzDiscreteEquationY (3D) =====================
struct HelmholtzDiscreteEquationY : HelmholtzDiscreteEquation
{
    HelmholtzDiscreteEquationY(): 
		HelmholtzDiscreteEquation((DiscreteOperatorWithFluxInterface*)&Ky,
              (GradientDiscreteOperatorBase*)&DyT,
              (SkewSymDiscreteOperator*)&Cy,
              nullptr, nullptr) 
	{
		;
	}

    LS_STAG_KERNEL_2D(symY, LS_STAG_R_PD _sol, const int _N, const int _NM)
    {
		for (int k = 0, last = (j + 1) * _NM - 1, num = i * _N + j * _NM; k < _N; k++, num++)
			_sol[last - num] = -_sol[num];
    }

    LS_STAG_KERNEL_2D(correctionY, LS_STAG_R_PD _right, LS_STAG_R_PD _Cell, LS_STAG_CR_PI _Number,
        LS_STAG_CR_PCT _type, LS_STAG_CR_PC _theta, LS_STAG_CR_PC _theta_N,
        LS_STAG_CR_PC _curV_C, LS_STAG_CR_PI _CellOnAirfoill, const int _hs, const int _hf, const int _N, const int _NM)
    {
        for (int k = i * _N + j * _NM, num = k + _hs, numB = num + j * _N, m_end = k + _hf; num < m_end; num++, numB++)
        {
            if (_type[num] > SolidCell)
            {
                if ((_theta_N[numB].x < EQ_DOUBLE) && (_theta[numB].x > EQ_DOUBLE))
                {
                    const int airNum = LS_STAG_VIRTUAL_INLINE_CALL(CLASS_NAME_AIRFOIL_NUM, getAirfoilNum, numB, _CellOnAirfoill);
                    _right[num] = _curV_C[airNum].y;
                    HelmholtzDiscreteEquation::_k_correctionOnImB(num, _Cell, _Number);
                }
            }
        }
    }

private:
    void initCurMesh() override
    {
        m_cur_mesh = &m_meshes->yMesh;
        m_prev_sol = m_meshes->m_prevUy;
        m_lastSol  = m_meshes->m_Uy;
    }

    void symmetric()
    {
        LS_STAG_PAR_FOR_2D(HelmholtzDiscreteEquationY, symY, center_row, m_cur_mesh->K, m_sol, m_cur_mesh->N, m_cur_mesh->M * m_cur_mesh->N)
    }

    void correctionRightMemberImB() override
    {
        LS_STAG_PROFILE("correctionY");
        LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteEquationY, correctionY, m_vs, m_vf, m_zs, m_zf,
            m_right, m_A->Cell, m_A->Number, m_type,
            m_theta, m_theta_N, m_meshes->curV_C, m_task->getPCellOnAirfoill(),
            m_hs, m_hf, m_cur_mesh->N, m_cur_mesh->N * m_cur_mesh->M)
    }

    HelmholtzDiscreteOperatorY Ky;
    GradientDiscreteOperatorY DyT;
    ConvectiveDiscreteOperatorY Cy;
};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_EQUATION_Y_H_
