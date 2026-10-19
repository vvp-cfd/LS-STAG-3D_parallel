/*
* \copyright (C) 2013 Puzikova V.V.
* \file      divergence_discrete_operator_y.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 11:45:45
* \brief     divergence_discrete_operator_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Y_H_
#define _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Y_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DivergenceDiscreteOperatorY

Данный файл содержит в себе реализацию класса DivergenceDiscreteOperatorY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calcDy, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, const coordinate*, _delta, const coordinate*, _theta, const int, _N);
	LS_STAG_DECLARATION_S_KERNEL(numBuildDy, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildDy, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, NM, const int, m_N, const int, m_K);
/*!
\brief Класс, реализующий дискретный оператор дивергенции на у-сетке

Данный класс реализует работу с дискретным аналогом оператора дивергенции на у-сетке
*/
	struct DivergenceDiscreteOperatorY : DiscreteOperatorInterface
	{
		using DiscreteOperatorInterface::initMeshParams;
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorInterface::initMeshParams(_cur_mesh, _meshes);
			build();
		};
		LS_STAG_SINGLE_KERNEL(numBuildDy, LS_STAG_R_PI Number, LS_STAG_CR_PCT m_type, const int NMK)
    {
        for (int i = 0, num = 0; i < NMK; ++i, Number[i] = num) {
            switch (m_type[i]) {
            case SolidCell:
            case N_Border: case S_Border:
            case NE_Edge: case NW_Edge: case SE_Edge: case SW_Edge:
            case TN_Edge: case TS_Edge: case BN_Edge: case BS_Edge:
            case TNE_Corner: case TNW_Corner: case TSE_Corner: case TSW_Corner:
            case BNE_Corner: case BNW_Corner: case BSE_Corner: case BSW_Corner:
                ++num;
                break;
            default:
                num += 2;   // S и P //твёрдая или с твёрдой N/S стенкой
            }
        }
    }	
		LS_STAG_SINGLE_KERNEL(portraitBuildDy, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int NM, const int m_N, const int m_K)
    {
		for(int k = 0, p = 0, s = -m_N, numA2 = 0, num = 0; k < m_K; k++, p-= m_N, s -=m_N)
		{
			for(int i = 0; i < m_N; i++, p++, s++, numA2++, num++) Portrait[numA2] = p;
		for (int i = m_N; i < NM; i++, p++, s++, numA2++, num++)
		{
                    switch (m_type[num]) 
					{
					case SolidCell:
                    case N_Border:
                    case NE_Edge: case NW_Edge:
                    case TN_Edge: case BN_Edge:
                    case TNE_Corner: case TNW_Corner:
                    case BNE_Corner: case BNW_Corner:
                        Portrait[numA2] = s;         // только S
                        break;
					case S_Border:
					case SE_Edge: case SW_Edge:
					case TS_Edge: case BS_Edge:
					case TSE_Corner: case TSW_Corner:
					case BSE_Corner: case BSW_Corner:
						Portrait[numA2] = p;       // только P
						break;
                    default:
                        Portrait[numA2] = s;
                        ++numA2;
                        Portrait[numA2] = p;
                    }
        }
	}
    }
    LS_STAG_KERNEL(calcDy, LS_STAG_CR_PCT _type, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell,
                   const coordinate* _delta, const coordinate* _theta, const int _N)
    {
        const coordinate* curD = _delta + i;
        const coordinate* curT = _theta + i;
        const int numA2 = _Number[i];
        switch (_type[i]) {
        case SolidCell:
            _Cell[numA2] = 0.0;
            break;
        case N_Border:
        case NE_Edge: case NW_Edge:
        case TN_Edge: case BN_Edge:
        case TNE_Corner: case TNW_Corner:
        case BNE_Corner: case BNW_Corner:
            _Cell[numA2] = -(curT - _N)->y * curD->x * curD->z;
            break;
        case S_Border:
        case SE_Edge: case SW_Edge:
        case TS_Edge: case BS_Edge:
        case TSE_Corner: case TSW_Corner:
        case BSE_Corner: case BSW_Corner:
            _Cell[numA2] = curT->y * curD->x * curD->z;
            break;
        default:
            double dxz = curD->x * curD->z;
            _Cell[numA2]     = -(curT - _N)->y * dxz;
            _Cell[numA2 + 1] =  curT->y * dxz;
        }
    }
	  protected:
		  /*!
		  Оценивает максимальное количество ненулевых элементов матрицы дискретного аналога оператора дивергенции на у-сетке
		  */
		void estimateNeedMemory()
		{
			maxElemCount = 2 * (NMK - m_N * m_K);
		};
		void buildNum()
		{
			LS_STAG_SINGLE(DivergenceDiscreteOperatorY, numBuildDy, Number, m_type, NMK)
		};
		  /*!
		  Строит портрет разреженной матрицы дискретного аналога оператора дивергенции на у-сетке
		  */
		void buildPortait()
		{
			LS_STAG_SINGLE(DivergenceDiscreteOperatorY, portraitBuildDy, Portrait, m_type, m_N *  m_M, m_N, m_K)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора дивергенции на у-сетке
		  */
		virtual void calcElem()
		{
			LS_STAG_PAR_FOR(DivergenceDiscreteOperatorY, calcDy, NMK, m_type, Number, Cell, m_delta, m_theta, m_N)
		};
	};
}

#endif //_LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_Y_H_
