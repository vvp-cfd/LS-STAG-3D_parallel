/*
* \copyright (C) 2013 Puzikova V.V.
* \file      sym_discrete_operator.h
* \author    Puzikova Valeria (vvp) 
* \date      01.10.2013 12:08:39
* \brief     sym_discrete_operator.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SYM_DISCRETE_OPERATOR_H_
#define _LS_STAG_SYM_DISCRETE_OPERATOR_H_

// ================================================================== INCLUDES
#include "../discrete_operator_with_flux_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SymDiscreteOperator

Данный файл содержит в себе реализацию класса SymDiscreteOperator
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numBuild, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(portraitBuild, LS_STAG_R_PI, Portrait, LS_STAG_CR_PCT, m_type, const int, m_N, 
		const int, m_NM, const int, NMK);
/*!
\brief Класс, реализующий симметрический дискретный оператор

Данный класс реализует интерфейс симметрического дискретного оператора с источниковым членом
*/
	struct SymDiscreteOperator : DiscreteOperatorWithFluxInterface
	{
		LS_STAG_SINGLE_KERNEL(numBuild, LS_STAG_R_PI Number, LS_STAG_CR_PCT m_type, const int NM)
		{
			for (int i = 0, num = 0; i < NM; i++, Number[i] = num)
			{
				switch (m_type[i])
				{
			case SolidCell:
                ++num;
                break;
            // Грани (отсутствует одна сторона) – остаётся 5 соседей
            case N_Border: case S_Border:
            case E_Border: case W_Border:
            case T_Border: case B_Border:
                num += 6;
                break;
            // Рёбра (отсутствуют две стороны) – остаётся 4 соседа
            case TN_Edge: case TS_Edge: case TE_Edge: case TW_Edge:
            case BN_Edge: case BS_Edge: case BE_Edge: case BW_Edge:
            case NE_Edge: case NW_Edge: case SE_Edge: case SW_Edge:
                num += 5;
                break;
            // Углы (отсутствуют три стороны) – остаётся 3 соседа
            case TNE_Corner: case TNW_Corner: case TSE_Corner: case TSW_Corner:
            case BNE_Corner: case BNW_Corner: case BSE_Corner: case BSW_Corner:
                num += 4;
                break;
				default:
					num += 7;
				}
			}
		}
		LS_STAG_SINGLE_KERNEL(portraitBuild, LS_STAG_R_PI Portrait, LS_STAG_CR_PCT m_type, const int m_N, const int m_NM, const int NMK)
		{
        for (int p = 0, n = m_N, s = -n, e = 1, w = -1, t = m_NM, b = -t, numA2 = 0; p < NMK; ++p, ++b, ++s, ++w, ++e, ++n, ++t) 
		{
            switch (m_type[p]) 
			{
            case SolidCell:
                Portrait[numA2++] = p;   // как в 2D – ставим e, нужно только одно значение
                break;
            // ---------- грани (5 соседей) ----------
            case N_Border:   // нет n
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = t;
                break;
            case S_Border:   // нет s
                Portrait[numA2++] = b;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case E_Border:   // нет e
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case W_Border:   // нет w
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case T_Border:   // нет t
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                break;
            case B_Border:   // нет b
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            // ---------- рёбра (4 соседа) ----------
            case TN_Edge:    // нет t,n
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                break;
            case TS_Edge:    // нет t,s
                Portrait[numA2++] = b;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                break;
            case TE_Edge:    // нет t,e
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                break;
            case TW_Edge:    // нет t,w
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                break;
            case BN_Edge:    // нет b,n
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = t;
                break;
            case BS_Edge:    // нет b,s
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case BE_Edge:    // нет b,e
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case BW_Edge:    // нет b,w
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case NE_Edge:    // нет n,e
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = t;
                break;
            case NW_Edge:    // нет n,w
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = t;
                break;
            case SE_Edge:    // нет s,e
                Portrait[numA2++] = b;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case SW_Edge:    // нет s,w
                Portrait[numA2++] = b;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            // ---------- углы (3 соседа) ----------
            case TNE_Corner: // нет t,n,e
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                break;
            case TNW_Corner: // нет t,n,w
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                break;
            case TSE_Corner: // нет t,s,e
                Portrait[numA2++] = b;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                break;
            case TSW_Corner: // нет t,s,w
                Portrait[numA2++] = b;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                break;
            case BNE_Corner: // нет b,n,e
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = t;
                break;
            case BNW_Corner: // нет b,n,w
                Portrait[numA2++] = s;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = t;
                break;
            case BSE_Corner: // нет b,s,e
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            case BSW_Corner: // нет b,s,w
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
                break;
            default:   // FluidCell и разрезанные ячейки – 7 соседей
                Portrait[numA2++] = b;
                Portrait[numA2++] = s;
                Portrait[numA2++] = w;
                Portrait[numA2++] = p;
                Portrait[numA2++] = e;
                Portrait[numA2++] = n;
                Portrait[numA2++] = t;
            }
		}
		}
		using DiscreteOperatorWithFluxInterface::set;
		virtual void set(const int _N, const int _M, const int _K)
		{//  к этому моменту уже д.б. theta и delta  (т.е. они построены) и V_IB для TA
			DiscreteOperatorWithFluxInterface::set(_N, _M, _K);//stLR = m_N * (m_M - 1);
		};
	  protected:
		virtual void initFlux(){};
		virtual void setRecalcFluxMarkers()
		{
			if((m_task->id_motion == motion::MOTION_OSC) || (m_task->id_motion == motion::MOTION_MOV))
				ibRecalcFlux = true;
			else// для MOV проверяем m_meshes->MovBody()
			{
				if(m_task->id == kind::TASK_FLOW)
					if((flow_task->id_airfoil == airfoil::AIRFOIL_CIRCLES) || 
						(flow_task->id_airfoil == airfoil::AIRFOIL_ROTOR)) 
						ibRecalcFlux = true;
			}
		};
	  public:
		virtual void calcFlux()
		{
			if(ibRecalcFlux) ibFlux();
		};
	  private:
		/*!
		Оценивает максимальное число ненулевых элементов разреженной симметрической матрицы дискретного оператора
		*/
		void estimateNeedMemory()
		{
			maxElemCount = 7 * NMK - 2 * ((m_N + m_M) * m_K + m_N * m_M);// max кол-во эл-в м-цы (если нет твёрдых ячеек)
		};
	  protected:
		void buildNum()
		{
			LS_STAG_SINGLE(SymDiscreteOperator, numBuild, Number, m_type, NMK)
		};
		/*!
		Строит портрет разреженной симметрической матрицы дискретного оператора
		*/
		void buildPortait()
		{
			LS_STAG_SINGLE(SymDiscreteOperator, portraitBuild, Portrait, m_type, m_N, m_N * m_M, NMK)
		};
	};
}
#endif //_LS_STAG_SYM_DISCRETE_OPERATOR_H_
