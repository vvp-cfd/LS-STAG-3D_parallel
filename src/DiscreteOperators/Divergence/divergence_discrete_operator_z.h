/*
* \copyright (C) 2026 Puzikova V.V.
* \file      divergence_discrete_operator_z.h
* \author    Puzikova Valeria (vvp) 
* \date      04.06.2026 11:07:35
* \brief     divergence_discrete_operator_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Z_H_
#define _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Z_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DivergenceDiscreteOperatorZ

Данный файл содержит в себе реализацию класса DivergenceDiscreteOperatorZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calcDz, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, const coordinate*, _delta, const coordinate*, _theta, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(numBuildDz, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildDz, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, NM, const int, NMK);
/*!
\brief Класс, реализующий дискретный оператор дивергенции на z-сетке

Данный класс реализует работу с дискретным аналогом оператора дивергенции на z-сетке
*/
struct DivergenceDiscreteOperatorZ : DiscreteOperatorInterface
{
    using DiscreteOperatorInterface::initMeshParams;
    virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes) override
    {
        DiscreteOperatorInterface::initMeshParams(_cur_mesh, _meshes);
        build();
    }
    LS_STAG_SINGLE_KERNEL(numBuildDz, LS_STAG_R_PI Number, LS_STAG_CR_PCT m_type, const int NMK)
    {
        for (int i = 0, num = 0; i < NMK; ++i, Number[i] = num) {
            switch (m_type[i]) {
            case SolidCell:
            case B_Border: case T_Border:
            case BN_Edge: case BS_Edge: case BE_Edge: case BW_Edge:
            case TN_Edge: case TS_Edge: case TE_Edge: case TW_Edge:
            case BNE_Corner: case BNW_Corner: case BSE_Corner: case BSW_Corner:
            case TNE_Corner: case TNW_Corner: case TSE_Corner: case TSW_Corner:
                ++num;
                break;
            default:
                num += 2;   // B и P
            }
        }
    }
    LS_STAG_SINGLE_KERNEL(portraitBuildDz, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type,
                          const int NM, const int NMK)
    {
		for (int p = 0; p < NM; p++) Portrait[p] = p; // только P
        for (int b = 0, p = NM, numA2 = p; p < NMK; ++numA2, p++, b++) 
        {
                    switch (m_type[p]) 
                    {
					case SolidCell:
                    case T_Border:
                    case TN_Edge: case TS_Edge: case TE_Edge: case TW_Edge:
                    case TNE_Corner: case TNW_Corner: case TSE_Corner: case TSW_Corner:
                        Portrait[numA2] = b;         // только B
                        break;
					case B_Border:
					case BN_Edge: case BS_Edge: case BE_Edge: case BW_Edge:
					case BNE_Corner: case BNW_Corner: case BSE_Corner: case BSW_Corner:
						Portrait[numA2] = p;       // только P
						break;
                    default:
                        Portrait[numA2] = b;
                        ++numA2;
                        Portrait[numA2] = p;       // B и P
                    }
        }
    }
    LS_STAG_KERNEL(calcDz, LS_STAG_CR_PCT _type, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell,
                   const coordinate* _delta, const coordinate* _theta, const int NM)
    {
        const coordinate* curD = _delta + i;
        const coordinate* curT = _theta + i;
        const int numA2 = _Number[i];
        switch (_type[i]) {
        case SolidCell:
            _Cell[numA2] = 0.0;
            break;
        case B_Border:
        case BN_Edge: case BS_Edge: case BE_Edge: case BW_Edge:
        case BNE_Corner: case BNW_Corner: case BSE_Corner: case BSW_Corner:
            _Cell[numA2] = curT->z * curD->x * curD->y;
            break;
        case T_Border:
        case TN_Edge: case TS_Edge: case TE_Edge: case TW_Edge:
        case TNE_Corner: case TNW_Corner: case TSE_Corner: case TSW_Corner:
            _Cell[numA2] = -(curT - NM)->z * curD->x * curD->y;
            break;
        default:
            double dxy = curD->x * curD->y;
            _Cell[numA2]     = -(curT - NM)->z * dxy;   // B
            _Cell[numA2 + 1] =  curT->z * dxy;         // P
        }
    }
protected:
    void estimateNeedMemory() override
    {
        maxElemCount = 2 * (NMK - m_N * m_M);
    }
    void buildNum() override
    {
        LS_STAG_SINGLE(DivergenceDiscreteOperatorZ, numBuildDz, Number, m_type, NMK)
    }
    void buildPortait() override
    {
        LS_STAG_SINGLE(DivergenceDiscreteOperatorZ, portraitBuildDz, Portrait, m_type, m_M * m_N, NMK)
    }
    void calcElem() override
    {
        LS_STAG_PAR_FOR(DivergenceDiscreteOperatorZ, calcDz, NMK, m_type, Number, Cell, m_delta, m_theta, m_M * m_N)
    }
};
}

#endif //_LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Z_H_
