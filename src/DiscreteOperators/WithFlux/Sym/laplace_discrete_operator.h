/*
* \copyright (C) 2013 Puzikova V.V.
* \file      laplace_discrete_operator.h
* \author    Puzikova Valeria (vvp) 
* \date      02.10.2013 0:21:19
* \brief     laplace_discrete_operator.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LAPLACE_DISCRETE_OPERATOR_H_
#define _LS_STAG_LAPLACE_DISCRETE_OPERATOR_H_

// ================================================================== INCLUDES
#include "sym_discrete_operator.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LaplaceDiscreteOperator

Данный файл содержит в себе реализацию класса LaplaceDiscreteOperator
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(calcLap, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell, LS_STAG_CR_PC, _theta, LS_STAG_CR_PC, _delta, 
		LS_STAG_R_PD, _inv_Mx, const double*, _inv_My, const double*, _inv_Mz, LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_R_KERNEL_2D(ibFlux, LS_STAG_R_PD, _flux, const coordinate*, _theta,
		LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(setBaseFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PC, _delta, LS_STAG_CR_PD, _Ux, 
		const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(initFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PC, _delta, 
		LS_STAG_CR_PD, _U_W, const double, _U_E, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxNSp, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, V_S, const double, V_N,const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxTBp, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_PC, m_theta, LS_STAG_R_PC, m_delta, 
		 const double, W_B, const double, W_T, const int, _N, const int, _NMK2);
	LS_STAG_DECLARATION_S_KERNEL(fixPointFlow, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, LS_STAG_CR_PI, Portrait, const int, numA);
	LS_STAG_DECLARATION_S_KERNEL(fixPointFlowBase, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(fixPointDef, LS_STAG_CR_PI, Number, LS_STAG_R_PD, Cell, const int, N2);
/*!
\brief Класс, реализующий дискретный оператор Лапласа

Данный класс реализует работу с дискретным аналогом оператора Лапласа на основной сетке
*/
	struct LaplaceDiscreteOperator : SymDiscreteOperator
	{
		LS_STAG_KERNEL_2D(initFlux, LS_STAG_R_PD _flux, LS_STAG_CR_PC _delta, LS_STAG_CR_PD _U_W, const double _U_E,
			const int _N, const int _M, const int _NM)
		{
			int numF = _N * i + _NM * j;
			const coordinate* curD = _delta + numF;
			double dyz = curD->y * curD->z;
			_flux[numF] = -_U_W[i + j * _M] * dyz;
			_flux[numF + _N - 1] = _U_E * dyz;
		}
	LS_STAG_KERNEL(addFluxNSp, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double V_S, const double V_N, 
		const int _N, const int _NM)
		{
			int numS = i * _NM;
			coordinate* curD = m_delta + numS;
			double dz = curD->z;
			for(int ii = 0, numN = numS + _NM - _N; ii < _N; ii++, numS++, numN++, curD++)
			{
				double coef = curD->x * dz;
				flux[numS] -= coef * V_S;
				flux[numN] += coef * V_N;
			}
		}
	LS_STAG_KERNEL(addFluxTBp, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_CR_PC m_theta, LS_STAG_R_PC m_delta, const double W_B, 
		const double W_T, const int _N, const int _NMK2)
		{
			int numB = i * _N;
			coordinate* curD = m_delta + numB;
			double dy = curD->y;
			for(int ii = 0, numT = _NMK2 + numB; ii < _N; ii++, numB++, numT++, curD++)
			{
				double coef = dy * curD->x;
				// theta.z базовых ячеек B/T: вклад стенки идёт только через жидкую (неусечённую) часть грани
				flux[numB] -= coef * (m_theta + numB)->z * W_B;
				flux[numT] += coef * (m_theta + numT)->z * W_T;
			}
		}	
		LS_STAG_SINGLE_KERNEL(fixPointFlow, LS_STAG_CR_PI Number, LS_STAG_R_PD Cell, LS_STAG_CR_PI Portrait, const int numA)
		{// одна точка: строка приводится к тождественной по фактическому портрету (E_Border для flow, BNE_Corner для artery)
			const int numA2 = Number[numA];
			for (int c = numA2, m_end = Number[numA + 1]; c < m_end; c++)
				Cell[c] = (Portrait[c] == numA) ? 1.0 : 0.0;
		}
        
		LS_STAG_SINGLE_KERNEL(fixPointFlowBase, LS_STAG_CR_PI Number, LS_STAG_R_PD Cell, const int NM)
		{// NM -- это NM - 1 (BNE)
			const int numA2 = Number[NM];
			Cell[numA2] = 0.0;
			Cell[numA2 + 1] = 0.0;
			Cell[numA2 + 2] = 1.0;
			Cell[numA2 + 3] = 0.0;
		}
		LS_STAG_SINGLE_KERNEL(fixPointDef, LS_STAG_CR_PI Number, LS_STAG_R_PD Cell, const int N2)
		{// (BSE)
			const int numA2 = Number[N2];
			Cell[numA2] = 0.0;
			Cell[numA2 + 1] = 1.0;
			Cell[numA2 + 2] = 0.0;
			Cell[numA2 + 3] = 0.0;
		}
		LS_STAG_KERNEL_2D(calcLap, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, LS_STAG_CR_PC _theta, LS_STAG_CR_PC _delta,
			LS_STAG_R_PD _inv_Mx, const double* _inv_My, const double* _inv_Mz, LS_STAG_CR_PCT _type, const int _N, 
			const int _M, const int _NM)
		{
			for (int num = _N * i + j * _NM, m_end = num + _N, numX = num - i - j * _M, numY = num - j * _N; num < m_end; 
				num++, numX++, numY++)
			{
				int numA2 = _Number[num];
				switch (_type[num])
				{
				case FluidCell:
				{
					const coordinate* curD = _delta + num;
					double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->y * curD->x;
					double aux = dyz * dyz;
					double a_w = aux * _inv_Mx[numX - 1];
					double a_e = aux * _inv_Mx[numX];
					aux = dxz * dxz;
					double a_s = aux * _inv_My[numY - _N];
					double a_n = aux * _inv_My[numY];
					aux = dxy * dxy;
					double a_b = aux * _inv_Mz[num - _NM];
					double a_t = aux * _inv_Mz[num];
					_Cell[numA2] = a_b;
					_Cell[numA2 + 1] = a_s;
					_Cell[numA2 + 2] = a_w;
					_Cell[numA2 + 3] = -a_b - a_s - a_w - a_e - a_n - a_t;
					_Cell[numA2 + 4] = a_e;
					_Cell[numA2 + 5] = a_n;
					_Cell[numA2 + 6] = a_t;

				}
				break;
				case SolidCell:
					_Cell[numA2] = PLUS_TINY;
					break;
				case N_Border:
				{
					const coordinate* curD = _delta + num;
					const coordinate* curT = _theta + num;
					double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->y * curD->x;
					double aux = (curT - 1)->x * dyz;
					double a_w = aux * aux * _inv_Mx[numX - 1];
					aux = curT->x * dyz;
					double a_e = aux * aux * _inv_Mx[numX];
					aux = (curT - _N)->y * dxz;
					double a_s = aux * aux * _inv_My[numY - _N];
					aux = (curT - _NM)->z * dxy;
					double a_b = aux * aux * _inv_Mz[num - _NM];
					aux = curT->z * dxy;
					double a_t = aux * aux * _inv_Mz[num];
					_Cell[numA2] = a_b;
					_Cell[numA2 + 1] = a_s;
					_Cell[numA2 + 2] = a_w;
					_Cell[numA2 + 3] = -a_b - a_s - a_w - a_e - a_t;
					_Cell[numA2 + 4] = a_e;
					_Cell[numA2 + 5] = a_t;
				}
				break;
				case S_Border:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_b - a_w - a_e - a_n - a_t;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_n;
    _Cell[numA2 + 5] = a_t;
}
break;
case E_Border:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w - a_n - a_t;
    _Cell[numA2 + 4] = a_n;
    _Cell[numA2 + 5] = a_t;
}
break;
case W_Border:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = -a_b - a_s - a_e - a_n - a_t;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_n;
    _Cell[numA2 + 5] = a_t;
}
break;
case T_Border:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w - a_e - a_n;
    _Cell[numA2 + 4] = a_e;
    _Cell[numA2 + 5] = a_n;
}
break;
case B_Border:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_s - a_w - a_e - a_n - a_t;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_n;
    _Cell[numA2 + 5] = a_t;
}
break;
// ---------- Рёбра (4 элемента) ----------
case TN_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w - a_e;
    _Cell[numA2 + 4] = a_e;
}
break;
case TS_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_b - a_w - a_e - a_n;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_n;
}
break;
case TE_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w - a_n;
    _Cell[numA2 + 4] = a_n;
}
break;
case TW_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = -a_b - a_s - a_e - a_n;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_n;
}
break;
case BN_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_s - a_w - a_e - a_t;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_t;
}
break;
case BS_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_w;
    _Cell[numA2 + 1] = -a_w - a_e - a_n - a_t;
    _Cell[numA2 + 2] = a_e;
    _Cell[numA2 + 3] = a_n;
    _Cell[numA2 + 4] = a_t;
}
break;
case BE_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_s - a_w - a_n - a_t;
    _Cell[numA2 + 3] = a_n;
    _Cell[numA2 + 4] = a_t;
}
break;
case BW_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = -a_s - a_e - a_n - a_t;
    _Cell[numA2 + 2] = a_e;
    _Cell[numA2 + 3] = a_n;
    _Cell[numA2 + 4] = a_t;
}
break;
case NE_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w - a_t;
    _Cell[numA2 + 4] = a_t;
}
break;
case NW_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = -a_b - a_s - a_e - a_t;
    _Cell[numA2 + 3] = a_e;
    _Cell[numA2 + 4] = a_t;
}
break;
case SE_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_b - a_w - a_n - a_t;
    _Cell[numA2 + 3] = a_n;
    _Cell[numA2 + 4] = a_t;
}
break;
case SW_Edge:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = -a_b - a_e - a_n - a_t;
    _Cell[numA2 + 2] = a_e;
    _Cell[numA2 + 3] = a_n;
    _Cell[numA2 + 4] = a_t;
}
break;
// ---------- Углы (3 элемента) ----------
case TNE_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = a_w;
    _Cell[numA2 + 3] = -a_b - a_s - a_w;
}
break;
case TNW_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_s;
    _Cell[numA2 + 2] = -a_b - a_s - a_e;
    _Cell[numA2 + 3] = a_e;
}
break;
case TSE_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_b - a_w - a_n;
    _Cell[numA2 + 3] = a_n;
}
break;
case TSW_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = (curT - _NM)->z * dxy;
    double a_b = aux * aux * _inv_Mz[num - _NM];
    _Cell[numA2]   = a_b;
    _Cell[numA2 + 1] = -a_b - a_e - a_n;
    _Cell[numA2 + 2] = a_e;
    _Cell[numA2 + 3] = a_n;
}
break;
case BNE_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = a_w;
    _Cell[numA2 + 2] = -a_s - a_w - a_t;
    _Cell[numA2 + 3] = a_t;
}
break;
case BNW_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = (curT - _N)->y * dxz;
    double a_s = aux * aux * _inv_My[numY - _N];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_s;
    _Cell[numA2 + 1] = -a_s - a_e - a_t;
    _Cell[numA2 + 2] = a_e;
    _Cell[numA2 + 3] = a_t;
}
break;
case BSE_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = (curT - 1)->x * dyz;
    double a_w = aux * aux * _inv_Mx[numX - 1];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = a_w;
    _Cell[numA2 + 1] = -a_w - a_n - a_t;
    _Cell[numA2 + 2] = a_n;
    _Cell[numA2 + 3] = a_t;
}
break;
case BSW_Corner:
{
    const coordinate* curD = _delta + num;
    const coordinate* curT = _theta + num;
    double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->x * curD->y;
    double aux = curT->x * dyz;
    double a_e = aux * aux * _inv_Mx[numX];
    aux = curT->y * dxz;
    double a_n = aux * aux * _inv_My[numY];
    aux = curT->z * dxy;
    double a_t = aux * aux * _inv_Mz[num];
    _Cell[numA2]   = -a_e - a_n - a_t;
    _Cell[numA2 + 1] = a_e;
    _Cell[numA2 + 2] = a_n;
    _Cell[numA2 + 3] = a_t;
}
break;
				default:
				{
					const coordinate* curD = _delta + num;
					const coordinate* curT = _theta + num;
					double dyz = curD->y * curD->z, dxz = curD->x * curD->z, dxy = curD->y * curD->x;
					double aux = (curT - 1)->x * dyz;
					double a_w = aux * aux * _inv_Mx[numX - 1];
					aux = curT->x * dyz;
					double a_e = aux * aux * _inv_Mx[numX];
					aux = (curT - _N)->y * dxz;
					double a_s = aux * aux * _inv_My[numY - _N];
					aux = curT->y * dxz;
					double a_n = aux * aux * _inv_My[numY];
					aux = (curT - _NM)->z * dxy;
					double a_b = aux * aux * _inv_Mz[num - _NM];
					aux = curT->z * dxy;
					double a_t = aux * aux * _inv_Mz[num];
					_Cell[numA2] = a_b;
					_Cell[numA2 + 1] = a_s;
					_Cell[numA2 + 2] = a_w;
					_Cell[numA2 + 3] = -a_b - a_s - a_w - a_e - a_n - a_t;
					_Cell[numA2 + 4] = a_e;
					_Cell[numA2 + 5] = a_n;
					_Cell[numA2 + 6] = a_t;
				}
				}//switch
			}//for i
		}
		LS_STAG_KERNEL_2D(ibFlux, LS_STAG_R_PD _flux, const coordinate* _theta, LS_STAG_CR_PC _delta,
			LS_STAG_CR_PC _Vib, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _baseSideType, const int _N, const int _NM)
		{
			for (int k = _N * i + j * _NM, num = k + 1, m_end = k + _N - 1; num < m_end; num++)// W и E границы не войдут
			{
				const CellTypes curType = _type[num];
				switch(curType)
				{// границы расчетной области сюда не отправляются
					case SolidCell:
					case FluidCell:
					_flux[num] = 0.0;
					break;
					case B_Border:
					{
					if (_baseSideType[num].inRange(SolidSide, FluidSide))
					{// усечённая нижняя грань — поток IB через усечение
						const coordinate* curD = _delta + num;
						double dx = curD->x, dy = curD->y, dz = curD->z;
						const coordinate* curVib = _Vib + num;
						const coordinate* curT = _theta + num;
						_flux[num] = (curVib->x * ((curT - 1)->x - curT->x) * dy +  curVib->y * ((curT - _N)->y - curT->y) * 
						              dx) * dz;
					}
					else _flux[num] = 0.0;
					}
					break;
					case T_Border:
					{
					if (_baseSideType[num].inRange(SolidSide, FluidSide))
					{// усечённая верхняя грань — поток IB через усечение
						const coordinate* curD = _delta + num;
						double dx = curD->x, dy = curD->y, dz = curD->z;
						const coordinate* curVib = _Vib + num;
						const coordinate* curT = _theta + num;
						_flux[num] = (curVib->x * ((curT - 1)->x - curT->x) * dy +  curVib->y * ((curT - _N)->y - curT->y) * 
						              dx) * dz;
					}
					else _flux[num] = 0.0;
					}
					break;// N, S грани в эту функцию не попадают
                    case TW_Edge:
                    case TE_Edge:
                    case BW_Edge:
                    case BE_Edge:
                    break;
					default:
					{
					const coordinate* curD = _delta + num;
					double dx = curD->x, dy = curD->y, dz = curD->z;
					const coordinate* curVib = _Vib + num;
					const coordinate* curT = _theta + num;
					_flux[num] = (curVib->x * ((curT - 1)->x - curT->x) * dy +  curVib->y * ((curT - _N)->y - curT->y) * 
					              dx) * dz + curVib->z * ((curT - _NM)->z - curT->z) * dx * dy;
					}
				}
			}
		}
		LS_STAG_KERNEL_2D(setBaseFlux, LS_STAG_R_PD _flux, LS_STAG_CR_PC _delta, LS_STAG_CR_PD _Ux, const int _N,
			const int _M, const int _NM)
		{
			const int num = (i + 1) * _N - 1 + j * _NM;
			const coordinate* curD = _delta + num;
			_flux[num] = _Ux[num - 1 - i - j * _M] * curD->y * curD->z;//за \delta t прошли \delta x
		}
		using DiscreteOperatorWithFluxInterface::initMeshParams;
		void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorWithFluxInterface::initMeshParams(_cur_mesh, _meshes);
			m_inv_Mx = m_meshes->xMesh.inv_VOF;
			m_inv_My = m_meshes->yMesh.inv_VOF;
			m_inv_Mz = m_meshes->zMesh.inv_VOF;
			outEndZ = m_K;
			m_NM = m_N * m_M;
			N2 = m_N - 1;
			fixPointNum = N2;
			fixPointFlowNum = m_NM - 1;
			isFlow = true;
			initFlux();// вклад от стенок РО, на которых г.у. I рода (соот-но строим один раз)
			if (baseRecalcFlux) m_Ux = m_meshes->m_Ux;// вызыв-м после setRecalcFluxMarkers()
			build();
			calcFlux();// 08.12.14
		};
		void calcFlux()
		{
			if(baseRecalcFlux) baseFlux();
			SymDiscreteOperator::calcFlux();
		};
		void build()
		{//пересчет потока нужно делать перед решением уравнения Гельмгольца для правильного расчета конвективного члена на ImB
			DiscreteOperatorInterface::build();
		};
	  private:
		void setRecalcFluxMarkers()
		{
			outEnd = m_task->m_mesh->baseMesh.M;
			if(m_task->id == kind::TASK_FLOW)
				baseRecalcFlux = true;
			SymDiscreteOperator::setRecalcFluxMarkers();
		};
		void baseFlux()//outflow
		{
			if(outEnd)
			LS_STAG_PAR_FOR_2D(LaplaceDiscreteOperator, setBaseFlux, outEnd, outEndZ, flux, m_delta, m_Ux, m_N, m_M, m_NM)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора Пуассона на основной сетке на погруженных границах 
		  */
		void ibFlux()
		{
			LS_STAG_PAR_FOR_2D_RANGE(LaplaceDiscreteOperator, ibFlux, 1, M2, 0, m_K, flux, m_theta, m_delta, m_Vib, 
				m_type, m_baseSideTypes, m_N, m_NM)
		}; 
		void initFlux()
		{// уже обнулили => присваиваем только ненулевые значения
			ibFlux();
			LS_STAG_PAR_FOR_2D(LaplaceDiscreteOperator, initFlux, m_M, m_K, flux, m_delta, 
				m_task->pU_W, m_task->U_E, m_N, m_M, m_NM)
			LS_STAG_PAR_FOR(LaplaceDiscreteOperator, addFluxNSp, m_K, flux, m_delta, m_task->V_S, m_task->V_N, m_N,
				 m_NM)
			LS_STAG_PAR_FOR(LaplaceDiscreteOperator, addFluxTBp, m_M, flux, m_type, m_theta, m_delta, m_task->W_B, m_task->W_T, m_N, 
				NMK - m_NM)
		};
		void fixPoint()
		{
			switch(m_task->id)
			{
			case kind::TASK_FLOW:
					if (isFlow) LS_STAG_SINGLE(LaplaceDiscreteOperator, fixPointFlow, Number, Cell, Portrait, fixPointFlowNum)
                    LS_STAG_SINGLE(LaplaceDiscreteOperator, fixPointFlow, Number, Cell, Portrait, fixPointNum)
					break;
				default:
					if (isFlow) LS_STAG_SINGLE(LaplaceDiscreteOperator, fixPointFlowBase, Number, Cell, fixPointFlowNum)
					LS_STAG_SINGLE(LaplaceDiscreteOperator, fixPointDef, Number, Cell, fixPointNum)
			}
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора Пуассона на основной сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(LaplaceDiscreteOperator, calcLap, m_M, m_K, Number, Cell, m_theta, m_delta, m_inv_Mx, 
				m_inv_My, m_inv_Mz, m_type, m_N, m_M, m_NM)
			fixPoint();
		};
		private:
			double *m_inv_Mx, *m_inv_My, *m_inv_Mz, *m_Ux;
			int outEnd, outEndZ, m_NM;
		public:
			int fixPointNum, fixPointFlowNum;
			bool isFlow;
	};
}
#endif //_LS_STAG_LAPLACE_DISCRETE_OPERATOR_H_
