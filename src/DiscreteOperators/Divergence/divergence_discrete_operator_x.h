/*
* \copyright (C) 2013 Puzikova V.V.
* \file      divergence_discrete_operator_x.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 11:08:32
* \brief     divergence_discrete_operator_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_X_H_
#define _LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_X_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DivergenceDiscreteOperatorX

Данный файл содержит в себе реализацию класса DivergenceDiscreteOperatorX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calcDx, LS_STAG_CR_PCT, _type, LS_STAG_CR_PI, _Number,
		LS_STAG_R_PD, _Cell, const coordinate*, _delta, const coordinate*, _theta);
	LS_STAG_DECLARATION_S_KERNEL(numBuildDx, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuildDx, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, m_M, const int, m_N, const int, m_K);
/*!
\brief Класс, реализующий дискретный оператор дивергенции на х-сетке

Данный класс реализует работу с дискретным аналогом оператора дивергенции на х-сетке
*/
	struct DivergenceDiscreteOperatorX : DiscreteOperatorInterface
	{
		using DiscreteOperatorInterface::initMeshParams;
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorInterface::initMeshParams(_cur_mesh, _meshes);
			build();
		};
		LS_STAG_SINGLE_KERNEL(numBuildDx, LS_STAG_R_PI Number, LS_STAG_CR_PCT m_type, const int NMK)
		{
			for (int i = 0, num = 0; i < NMK; i++, Number[i] = num)
			{
				switch (m_type[i])
				{
				case SolidCell:
			case W_Border:
			case NW_Edge:
			case SW_Edge:
			case TW_Edge:
			case BW_Edge:
			case TNW_Corner:
			case TSW_Corner:
			case BNW_Corner:
			case BSW_Corner:
			case E_Border:
			case NE_Edge:
			case SE_Edge:
			case TE_Edge:
			case BE_Edge:
			case TNE_Corner:
			case TSE_Corner:
			case BNE_Corner:
			case BSE_Corner:
					num++;
					break;
				default:
					num += 2;
				}
			}
		}
		LS_STAG_SINGLE_KERNEL(portraitBuildDx, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int m_M, 
			const int m_N, const int m_K)
		{
			for(int k = 0, p = 0, w = -1, num = 0, numA2 = 0; k < m_K; k++)
			{
			for (int j = 0; j < m_M; j++, p--, w--)
			{
				Portrait[numA2] = p; p++, w++, num++, numA2++;
				for (int i = 1; i < m_N; i++, p++, w++, num++, numA2++)
				{
					switch (m_type[num])
					{
					case SolidCell:
					case E_Border:
					case NE_Edge:
					case SE_Edge:
					case TE_Edge:
					case BE_Edge:
					case TNE_Corner:
					case TSE_Corner:
					case BNE_Corner:
					case BSE_Corner:
						Portrait[numA2] = w;
						break;
			case W_Border:
			case NW_Edge:
			case SW_Edge:
			case TW_Edge:
			case BW_Edge:
			case TNW_Corner:
			case TSW_Corner:
			case BNW_Corner:
			case BSW_Corner:
						Portrait[numA2] = p;
						break;
					default:
						Portrait[numA2] = w;
						numA2++;
						Portrait[numA2] = p;
					}//switch
				}//for i
			}//for j
			}
		}
		LS_STAG_KERNEL(calcDx, LS_STAG_CR_PCT _type, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell,
			const coordinate* _delta, const coordinate* _theta)
		{
			const coordinate* curD = _delta + i;
			const coordinate* curT = _theta + i;
			const int numA2 = _Number[i];
			switch (_type[i])
			{
			case SolidCell:
				_Cell[numA2] = 0.0;
				break;
			case E_Border:
			case NE_Edge:
			case SE_Edge:
			case TE_Edge:
			case BE_Edge:
			case TNE_Corner:
			case TSE_Corner:
			case BNE_Corner:
			case BSE_Corner:
				_Cell[numA2] = -(curT - 1)->x * curD->y * curD->z;
				break;
			case W_Border:
			case NW_Edge:
			case SW_Edge:
			case TW_Edge:
			case BW_Edge:
			case TNW_Corner:
			case TSW_Corner:
			case BNW_Corner:
			case BSW_Corner:
				_Cell[numA2] = curT->x * curD->y * curD->z;
				break;
			default:
			{
				double dyz = curD->y * curD->z;
				_Cell[numA2] = -(curT - 1)->x * dyz;
				_Cell[numA2 + 1] = curT->x * dyz;
			}
			}//switch
		}
	  protected:
		  /*!
		  Оценивает максимальное количество ненулевых элементов матрицы дискретного аналога оператора дивергенции на х-сетке
		  */
		void estimateNeedMemory()
		{
			maxElemCount = 2 * (NMK - m_M * m_K);
		};
		void buildNum()
		{
			LS_STAG_SINGLE(DivergenceDiscreteOperatorX, numBuildDx, Number, m_type, NMK)
		};
		  /*!
		  Строит портрет разреженной матрицы дискретного аналога оператора дивергенции на х-сетке
		  */
		void buildPortait()
		{
			LS_STAG_SINGLE(DivergenceDiscreteOperatorX, portraitBuildDx, Portrait, m_type, m_M, m_N, m_K)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора дивергенции на х-сетке
		  */
		virtual void calcElem()
		{
			LS_STAG_PAR_FOR(DivergenceDiscreteOperatorX, calcDx, NMK, m_type, Number, Cell, m_delta, m_theta)
		};
	};
}

#endif //_LS_STAG_DIVERGENCE_DISCRETE_OPERATOR_X_H_
