/*
* \copyright (C) 2026 Puzikova V.V.
* \file      conditions_wnstb.h
* \author    Puzikova Valeria (vvp) 
* \date      08.06.2026 22:22:25
* \brief     conditions_wnstb.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_CONDITIONS_WNSTB_H_
#define _LS_STAG_CONDITIONS_WNSTB_H_

// ================================================================== INCLUDES
#include "MeshFunctions/mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ConditionWNSTB

Данный файл содержит в себе реализацию класса ConditionWNSTB
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с г.у. на границах расчетной области на основной сетке

Данный класс реализует хранение г.у. на границах расчетной области на основной сетке и некоторые базовые операции с ними
*/
    struct ConditionWNSTB
		{
			void init(const int _N, const int M, const int K)
			{
				W.init(M * K);
				N.init(_N * K);
				S.init(_N * K);
				B.init(_N * M);
				T.init(_N * M);
				W.set(0.0);
				N.set(0.0);
				S.set(0.0);
				B.set(0.0);
				T.set(0.0);
			}
            void copy(const double* m_W, const double* m_N, const double* m_S, const double* m_B, const double* m_T)
            {
                W.copy(m_W);
                N.copy(m_N);
                S.copy(m_S);
                B.copy(m_B);
                T.copy(m_T);
            }
            void copy(const ConditionWNSTB* _wnsbt)
            {
                W.copy(&(_wnsbt->W));
                N.copy(&(_wnsbt->N));
                S.copy(&(_wnsbt->S));
                B.copy(&(_wnsbt->B));
                T.copy(&(_wnsbt->T));
            }
            void copyTo(extDouble* _W, extDouble* _N, extDouble* _S, extDouble* _T, extDouble* _B)
            {
				_W->copy(&W);
				_N->copy(&N);
				_S->copy(&S);
				_B->copy(&B);
				_T->copy(&T);
            }
			void operator *=(const double factor)
			{
				W *= factor;
				N *= factor;
				S *= factor;
				B *= factor;
				T *= factor;
			}
			void operator +=(const ConditionWNSTB* _add)
			{
				W += (extDouble*)(&(_add->W));
				N += (extDouble*)(&(_add->N));
				S += (extDouble*)(&(_add->S));
				B += (extDouble*)(&(_add->B));
				T += (extDouble*)(&(_add->T));
			}
			void operator -=(const ConditionWNSTB* _add)
			{
				W -= (extDouble*)(&(_add->W));
				N -= (extDouble*)(&(_add->N));
				S -= (extDouble*)(&(_add->S));
				B -= (extDouble*)(&(_add->B));
				T -= (extDouble*)(&(_add->T));
			}
			void copyWithScale(const ConditionWNSTB* values, const double factor)
			{
				W.copyWithScale((extDouble*)(&(values->W)), factor);
				N.copyWithScale((extDouble*)(&(values->N)), factor);
				S.copyWithScale((extDouble*)(&(values->S)), factor);
				B.copyWithScale((extDouble*)(&(values->B)), factor);
				T.copyWithScale((extDouble*)(&(values->T)), factor);
			}
			void addWithScale(const ConditionWNSTB* values, const double factor)
			{
				W.addWithScale((extDouble*)(&(values->W)), factor);
				N.addWithScale((extDouble*)(&(values->N)), factor);
				S.addWithScale((extDouble*)(&(values->S)), factor);
				B.addWithScale((extDouble*)(&(values->B)), factor);
				T.addWithScale((extDouble*)(&(values->T)), factor);
			}
			void getCopy(double*& m_W, double*& m_N, double*& m_S, double*& m_B, double*& m_T)
			{
				W.getCopy(m_W);
				N.getCopy(m_N);
				S.getCopy(m_S);
				B.getCopy(m_B);
				T.getCopy(m_T);
			}
			void setCopy(double*& m_W, double*& m_N, double*& m_S, double*& m_B, double*& m_T)
			{
				W.setCopy(m_W);
				N.setCopy(m_N);
				S.setCopy(m_S);
				B.setCopy(m_B);
				T.setCopy(m_T);
			}
			void loadState(std::ifstream &in)
			{
				W.loadState(in);
				N.loadState(in);
				S.loadState(in);
				B.loadState(in);
				T.loadState(in);
			}
			void saveState(std::ofstream &out)
			{
				W.saveState(out);
				N.saveState(out);
				S.saveState(out);
				B.saveState(out);
				T.saveState(out);
			}
			extDouble W, N, S, B, T;
		};

} // namespace LS_STAG

#endif //_LS_STAG_CONDITIONS_WNSTB_H_