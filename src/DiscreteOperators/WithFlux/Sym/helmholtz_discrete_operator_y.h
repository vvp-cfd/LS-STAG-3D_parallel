/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_operator_y.h
* \author    Puzikova Valeria (vvp) 
* \date      19.10.2013 17:34:59
* \brief     helmholtz_discrete_operator_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Y_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Y_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_operator_base.h"
// ================================================================== 
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteOperatorY

Данный файл содержит в себе реализацию класса HelmholtzDiscreteOperatorY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(calcHelmY, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, const coordinate*, _alpha, const coordinate*, _delta, 
		const coordinate*, _inv_delta_ib_ne, const coordinate*, _inv_delta_ib_tn, const coordinate*, _delta_ib_en,
		const coordinate*, _delta_ib_ws, const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, 
		LS_STAG_CR_PD, _inv_V, LS_STAG_CR_PD, _My, LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, 
		const double, _nu, const double, _coefIntegr, const bool, _isFlow);
	LS_STAG_DECLARATION_R_KERNEL_2D(ibFluxHelmY, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType,
		const coordinate*, _theta, const coordinate*, _alpha, LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _inv_delta_ib_ne, 
		LS_STAG_CR_PC, _inv_delta_ib_tn, const coordinate*, _delta_ib_en, const coordinate*, _delta_ib_ws, 
		const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PC, _VibZY, 
		 LS_STAG_CR_PC, _VibYX, LS_STAG_CR_PD, _inv_V, const int, _N, const int, _NM, const int, baseNM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(setPlusOutFlowY, LS_STAG_R_PD, _PlusOutFlow, LS_STAG_CR_PC, _delta,
		const double, _nu, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxNS, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, V_S, const double, V_N, const int, _N, const int, _NM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL(addFluxTB, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, V_B, const double, V_T, const int, _N, const int, _K2, const int, _NMK2, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(inFlowFluxHelmY, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, _V_W, const double, _V_E, const bool, isCavity, const int, _N, const int, _NM, const double, _nu);
/*!
\brief Класс, реализующий дискретный оператор Гельмгольца на у-сетке

Данный класс реализует работу с дискретным аналогом оператора Гельмгольца на у-сетке
*/
	struct HelmholtzDiscreteOperatorY : HelmholtzDiscreteOperatorBase
	{//m_param_x = My, m_param_y = 1/V
	LS_STAG_KERNEL(addFluxNS, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double V_S, const double V_N, 
		const int _N, const int _NM, const double _nu)
		{
			int numS = i * _NM, numN = numS + _NM - _N, im = i * _N;
			coordinate* curD = m_delta + (numS + im);
			double coefS = _nu * V_S * curD->z / curD->y;
			double coefN = _nu * V_N * curD->z / (m_delta + (numN + im))->y;
			for(int ii = 0; ii < _N; ii++, numS++, numN++, curD++)
			{
				flux[numS] += curD->x * coefS;
				flux[numN] += curD->x * coefN;
			}
		}
	LS_STAG_KERNEL(addFluxTB, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta,  const double V_B, const double V_T,
		const int _N, const int _K2, const int _NMK2, const double _nu)
		{
			int numB = i * _N, numT = _NMK2 + numB;
			coordinate* curD = m_delta + numB; 
			double dyy = (curD->y + (curD + _N)->y);
			double coefB = _nu * V_B * dyy / curD->z;
			double coefT = _nu * V_T * dyy / (m_delta + (numT + i + _K2 * _N))->z;
			for(int ii = 0; ii < _N; ii++, numB++, numT++, curD++)
			{
				if (_type[numB] != SolidCell) flux[numB] += curD->x * coefB;
				if (_type[numT] != SolidCell) flux[numT] += curD->x * coefT;
			}
		}
		LS_STAG_KERNEL_2D(inFlowFluxHelmY, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, const double _V_W, const double _V_E,
			const bool isCavity, const int _N, const int _NM, const double _nu)
		{
			int numF = _N * i + _NM * j, _N2 = _N - 1, iM = j * _N;
			coordinate* curD = _delta + numF + iM;
			double dyz = _nu * curD->z * (curD->y + (curD + _N)->y);
			_flux[numF] = _V_W * dyz / curD->x;
			if (isCavity) 
			{
			numF += _N2;
			_flux[numF] = _V_E * dyz / (curD + _N2)->x;
			}
		}	
		LS_STAG_KERNEL_2D(calcHelmY, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate* _theta,
			const coordinate* _alpha, const coordinate* _delta, const coordinate* _inv_delta_ib_ne, 
			const coordinate* _inv_delta_ib_tn, const coordinate* _delta_ib_en, const coordinate* _delta_ib_ws, 
			const coordinate* _delta_ib_nb, const coordinate* _delta_ib_st, LS_STAG_CR_PD _inv_V, LS_STAG_CR_PD _My, 
			LS_STAG_CR_PCT _type, const int _N, const int _M, const int _NM, const double _nu, const double _coefIntegr,
			const bool _isFlow)
		{
			for (int numY = _N * i + _NM * j, m_end = numY + _N, num = numY + j * _N, baseNM = _N * (_M + 1); 
			numY < m_end; numY++, num++)
			{
				int numA2 = _Number[numY];
				switch (_type[numY])
				{
				case FluidCell:
				{
					const coordinate* curD = _delta + num;
					const coordinate* curT = _theta + num;
					double dz = curD->z, dx = curD->x, dyy = _nu * ((curD + _N)->y + curD->y);
					double aux = dz * dyy;
					double k_w = aux / (dx + (curT - 1)->y * (curD - 1)->x);
					double k_e = aux / (dx + (curT + 1)->y * (curD + 1)->x);
					aux = dx * dyy;
					double k_b = aux / (dz + (curT - baseNM)->y * (curD - baseNM)->z);
					double k_t = aux / (dz + (curT + baseNM)->y * (curD + baseNM)->z);
					aux = dx * dz;
					dyy = _nu * aux * aux;
					double k_n = dyy * _inv_V[num+_N];
					double k_s = dyy * _inv_V[num];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s * (curT - _N)->y;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_s - k_w - k_e - k_n - k_t;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_n * (curT + _N)->y;
					_Cell[numA2 + 6] = k_t;
				}
				break;
				case SolidCell:
					_Cell[numA2] = PLUS_TINY;
					break;
				case N_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num; 
					const coordinate* nT = curT + _N; const coordinate* wT = curT - 1; const coordinate* bT = curT - baseNM;
					double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty  = curT->y;
					double aux = dz * _nu, aux2 = ty * dx;
					double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
					double k_e = aux * (curT->x * dy + nT->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
					aux = dx * _nu, aux2 = ty * dz;
					double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
					double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
					dx *= dz;
					double k_s =  _nu * ty * dx * dx * _inv_V[num];
					double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s * (curT - _N)->y;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_w - k_s * ty - k_e - k_t - k_wallN;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_t;
				}
				break;
				case S_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const  coordinate* wT = curT - 1; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx;
    double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
    double k_e = aux * (curT->x * dy + nT->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
    aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_w - k_e - k_n * ty - k_t - k_wallS;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_n * nT->y;
    _Cell[numA2 + 5] = k_t;
}
break;

case E_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* wT = curT - 1; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - (k_s + k_n) * ty - k_w - k_e - k_t;
    _Cell[numA2 + 4] = k_n * (curT+_N)->y;
    _Cell[numA2 + 5] = k_t;
}
break;

case W_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - (k_s + k_n) * ty - k_e - k_t;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_n * nT->y;
    _Cell[numA2 + 5] = k_t;
}
break;

case T_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const   coordinate* wT = curT - 1; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx, az = (curA - 1)->z;
    bool existIb = (az < EQ_DOUBLE);
    double k_w = (existIb) ? 0.0 : (aux * (wT->x * dy + (wT + _N)->x * dyy) * az / (aux2 + wT->y * (curD - 1)->x));
    az = curA->z; bool existIb2 = (az < EQ_DOUBLE);
    double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dy + (curT + _N)->x * dyy) * az / (aux2 + (curT + 1)->y * (curD + 1)->x));
    double ib = (existIb || existIb2) ? (_nu * (_delta_ib_en[num].z + _delta_ib_ws[num+_N].z) * _inv_delta_ib_ne[num].x) : 0.0;
    aux = dx * _nu, aux2 = ty * dz;
    az = (curA - baseNM)->x; existIb = (az < EQ_DOUBLE);
    double k_b = (existIb) ? 0.0 : (aux * (bT->z * dy + (bT + _N)->z * dyy) * az / (aux2 + bT->y * (curD - baseNM)->z));
    if(existIb) ib += _nu * (_delta_ib_nb[num].x + _delta_ib_st[num+_N].x) * _inv_delta_ib_tn[num].z;
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - (k_s + k_n) * ty - k_w - k_e - k_wallT - ib;
    _Cell[numA2 + 4] = k_e;
    _Cell[numA2 + 5] = k_n * (curT+_N)->y;
}
break;

case B_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const    coordinate* nT = curT + _N; const  coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx, az = (curA - 1)->z;
    bool existIb = (az < EQ_DOUBLE);
    double k_w = (existIb) ? 0.0 : (aux * (wT->x * dy + (wT + _N)->x * dyy) * az / (aux2 + wT->y * (curD - 1)->x));
    az = curA->z; bool existIb2 = (az < EQ_DOUBLE);
    double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dy + nT->x * dyy) * az / (aux2 + (curT + 1)->y * (curD + 1)->x));
    double ib = (existIb || existIb2) ? (_nu * (_delta_ib_en[num].z + _delta_ib_ws[num+_N].z) * _inv_delta_ib_ne[num].x) : 0.0;
    aux = dx * _nu, aux2 = ty * dz;
    az = curA->x; existIb2 = (az < EQ_DOUBLE);
    double k_t = (existIb2) ? 0.0 : (aux * (curT->z * dy + nT->z * dyy) * az / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z));
    if(existIb2) ib += _nu * (_delta_ib_nb[num].x + _delta_ib_st[num+_N].x) * _inv_delta_ib_tn[num].z;
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz)) : 0.0;
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - (k_s + k_n) * ty - k_w - k_e - k_t - k_wallB - ib;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_n * nT->y;
    _Cell[numA2 + 5] = k_t;
}
break;
// ---------- Рёбра (4 элемента) ----------
case TN_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* wT = curT - 1; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx;
    double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
    double k_e = aux * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_s * ty - k_w - k_e - k_wallT - k_wallN;
    _Cell[numA2 + 4] = k_e;
}
break;

case TS_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* wT = curT - 1; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx;
    double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
    double k_e = aux * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_w - k_e - k_n * ty - k_wallT - k_wallS;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_n * (curT+_N)->y;
}
break;

case TE_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* wT = curT - 1; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - (k_s + k_n) * ty - k_w - k_e - k_wallT;
    _Cell[numA2 + 4] = k_n * (curT+_N)->y;
}
break;

case TW_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - (k_s + k_n) * ty - k_e - k_wallT;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_n * (curT+_N)->y;
}
break;

case BN_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx;
    double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
    double k_e = aux * (curT->x * dy + nT->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_s * ty - k_w - k_e - k_t - k_wallB - k_wallN;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_t;
}
break;

case BS_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N; const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double aux = dz * _nu, aux2 = ty * dx;
    double k_w = aux * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (aux2 + wT->y * (curD - 1)->x);
    double k_e = aux * (curT->x * dy + nT->x * dyy) * curA->z / (aux2 + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_w;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - k_w - k_e - k_n * ty - k_t - k_wallB - k_wallS;
    _Cell[numA2 + 2] = k_e;
    _Cell[numA2 + 3] = k_n * nT->y;
    _Cell[numA2 + 4] = k_t;
}
break;
case BE_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - (k_s + k_n) * ty - k_w - k_e - k_t - k_wallB;
    _Cell[numA2 + 3] = k_n * (curT+_N)->y;
    _Cell[numA2 + 4] = k_t;
}
break;

case BW_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz; dyy = _nu * ty * dx * dx;
    double k_n = dyy * _inv_V[num + _N];
    double k_s = dyy * _inv_V[num];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - (k_s + k_n) * ty - k_e - k_t - k_wallB;
    _Cell[numA2 + 2] = k_e;
    _Cell[numA2 + 3] = k_n * nT->y;
    _Cell[numA2 + 4] = k_t;
}
break;

case NE_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM; const  coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_s * ty - k_w - k_e - k_t - k_wallN;
    _Cell[numA2 + 4] = k_t;
}
break;

case NW_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_s * ty - k_e - k_t - k_wallN;
    _Cell[numA2 + 3] = k_e;
    _Cell[numA2 + 4] = k_t;
}
break;
case SE_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const coordinate* bT = curT - baseNM; const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_w - k_e - k_n * ty - k_t - k_wallS;
    _Cell[numA2 + 3] = k_n * nT->y;
    _Cell[numA2 + 4] = k_t;
}
break;
case SW_Edge:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double aux = dx * _nu, aux2 = ty * dz;
    double k_b = aux * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (aux2 + bT->y * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dy + nT->z * dyy) * curA->x / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - k_b - k_e - k_n * ty - k_t - k_wallS;
    _Cell[numA2 + 2] = k_e;
    _Cell[numA2 + 3] = k_n * nT->y;
    _Cell[numA2 + 4] = k_t;
}
break;
// ---------- Углы (3 элемента) ----------
case TNE_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* bT = curT - baseNM; const  coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = k_w;
    _Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_s * ty - k_w - k_e - k_wallT - k_wallN;
}
break;

case TNW_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s * (curT - _N)->y;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_s * ty - k_e - k_wallT - k_wallN;
    _Cell[numA2 + 3] = k_e;
}
break;

case TSE_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM; const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_b - k_w - k_e - k_n * ty - k_wallT - k_wallS;
    _Cell[numA2 + 3] = k_n * nT->y;
}
break;

case TSW_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_b = dx * _nu * (bT->z * dy + (bT + _N)->z * dyy) * (curA - baseNM)->x / (ty * dz + bT->y * (curD - baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallT = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - k_b - k_e - k_n * ty - k_wallT - k_wallS;
    _Cell[numA2 + 2] = k_e;
    _Cell[numA2 + 3] = k_n * nT->y;
}
break;

case BNE_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = k_w;
    _Cell[numA2 + 2] = -(_coefIntegr * _My[numY]) - k_s * ty - k_w - k_e - k_t - k_wallB - k_wallN;
    _Cell[numA2 + 3] = k_t;
}
break;

case BNW_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_s = _nu * ty * dx * dx * _inv_V[num];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallN = _nu * ty * dx * dx * _inv_V[num + _N];
    _Cell[numA2]     = k_s * (curT - _N)->y;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - k_s * ty - k_e - k_t - k_wallB - k_wallN;
    _Cell[numA2 + 2] = k_e;
    _Cell[numA2 + 3] = k_t;
}
break;

case BSE_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* wT = curT - 1;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_w = dz * _nu * (wT->x * dy + (wT + _N)->x * dyy) * (curA - 1)->z / (ty * dx + wT->y * (curD - 1)->x);
    double k_e = dz * _nu * (curT->x * dy + (curT + _N)->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = k_w;
    _Cell[numA2 + 1] = -(_coefIntegr * _My[numY]) - k_w - k_e - k_n * ty - k_t - k_wallB - k_wallS;
    _Cell[numA2 + 2] = k_n * (curT+_N)->y;
    _Cell[numA2 + 3] = k_t;
}
break;
case BSW_Corner:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* nT = curT + _N;
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty = curT->y;
    double k_e = dz * _nu * (curT->x * dy + nT->x * dyy) * curA->z / (ty * dx + (curT + 1)->y * (curD + 1)->x);
    double k_t = dx * _nu * (curT->z * dy + nT->z * dyy) * curA->x / (ty * dz + (curT + baseNM)->y * (curD + baseNM)->z);
    dx *= dz;
    double k_n = _nu * ty * dx * dx * _inv_V[num + _N];
double k_wallB = _isFlow ? 0.0 : (ty > EQ_DOUBLE) ? (dx * _nu * (curT->z * dy + (curT + _N)->z * dyy) * curA->x / (ty * dz)) : 0.0;
    double k_wallS = _nu * ty * dx * dx * _inv_V[num];
    _Cell[numA2]     = -(_coefIntegr * _My[numY]) - k_e - k_n * ty - k_t - k_wallB - k_wallS;
    _Cell[numA2 + 1] = k_e;
    _Cell[numA2 + 2] = k_n * nT->y;
    _Cell[numA2 + 3] = k_t;
}
break;
				default:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num; 
					const coordinate* nT = curT + _N; const coordinate* wT = curT - 1; const coordinate* bT = curT - baseNM;
					double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD + _N)->y, ty  = curT->y;
					double aux = dz * _nu, aux2 = ty * dx, az = (curA - 1)->z;
					bool existIb = (az < EQ_DOUBLE);
					double k_w = (existIb) ? 0.0 : (aux * (wT->x * dy + (wT + _N)->x * dyy) * az / (aux2 + wT->y * (curD - 1)->x));
					az = curA->z; bool existIb2 = (az < EQ_DOUBLE);
					double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dy + nT->x * dyy) * az / (aux2 + (curT + 1)->y * (curD + 1)->x));
double ib = (existIb || existIb2) ? (_nu * (_delta_ib_en[num].z + _delta_ib_ws[num+_N].z) * _inv_delta_ib_ne[num].x) : 0.0;
					aux = dx * _nu, aux2 = ty * dz;
					az = (curA - baseNM)->x; existIb = (az < EQ_DOUBLE);
					double k_b = (existIb) ? 0.0 : (aux * (bT->z * dy + (bT + _N)->z * dyy) * az / (aux2 + bT->y * (curD - baseNM)->z));
					az = curA->x; existIb2 = (az < EQ_DOUBLE);
					double k_t = (existIb2) ? 0.0 : (aux * (curT->z * dy + nT->z * dyy) * az / (aux2 + (curT + baseNM)->y * (curD + baseNM)->z));
					if(existIb || existIb2) ib += _nu * (_delta_ib_nb[num].x + _delta_ib_st[num+_N].x) * _inv_delta_ib_tn[num].z;
					dx *= dz;
					dyy = _nu * ty * dx * dx;
					double k_n = dyy * _inv_V[num+_N];
					double k_s = dyy * _inv_V[num];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s * (curT - _N)->y;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _My[numY]) - k_b - k_w - (k_s + k_n) * ty - k_e - k_t - ib;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_n * nT->y;
					_Cell[numA2 + 6] = k_t;
				}
				}//switch
			}//for i
		}
		LS_STAG_KERNEL_2D(ibFluxHelmY, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _SideType, const coordinate* _theta, const coordinate* _alpha,
			LS_STAG_CR_PC _delta, LS_STAG_CR_PC _inv_delta_ib_ne, LS_STAG_CR_PC _inv_delta_ib_tn, const coordinate* _delta_ib_en, 
			const coordinate* _delta_ib_ws, const coordinate* _delta_ib_nb, const coordinate* _delta_ib_st, 
			LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibZY, LS_STAG_CR_PC _VibYX, const double* _inv_V, const int _N, const int _NM, const int baseNM,
			const double _nu)
		{
			for (int numY = i * _N + j * _NM, num = numY + j * _N, bNum = (j) ? baseNM : 0, m_end = numY + _N; numY < m_end; numY++, num++)//for(i = 0; i < m_N; i++)
			{
				const CellTypes curType = _type[numY];
				bool isCutCell = false;
				if((curType > SolidCell) && (curType < FluidCell))
				isCutCell = true;
				else
				{
					switch(curType)
					{
						case B_Border:
						case T_Border:
                        if ((_SideType[num].inRange(SolidSide, FluidSide)) || (_SideType[num + _N].inRange(SolidSide, FluidSide))) isCutCell = true;
					}
				}
				if (isCutCell)
				{
					const coordinate* curT = _theta + num; double ty = curT->y;
					if (ty < EQ_DOUBLE)
						_flux[numY] = 0.0;
					else
					{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curVib = _Vib + num;
					double dxz = curD->x * curD->z;
					double ib = ((curA->z < EQ_DOUBLE) || ((curA - 1)->z < EQ_DOUBLE)) ? 
					     ((_delta_ib_en[num].z + _delta_ib_ws[num+_N].z) * _inv_delta_ib_ne[num].x * _VibYX[num].y) : 0.0;
					if((curA->x < EQ_DOUBLE) || ((curA - bNum)->x < EQ_DOUBLE)) 
					     ib += (_delta_ib_nb[num].x + _delta_ib_st[num+_N].x) * _inv_delta_ib_tn[num].z * _VibZY[num].y;
					_flux[numY] = _nu * (ib + ty * dxz * dxz * (((ty - (curT + _N)->y) * (curVib + _N)->y * _inv_V[num + _N]) - 
					                  (((curT - _N)->y - ty) * curVib->y * _inv_V[num])));
					}
				}// Grid.TypeOfCellX[numX] > 9
				else _flux[numY] = 0.0;
			}//for i
		}
		LS_STAG_KERNEL_2D(setPlusOutFlowY, LS_STAG_R_PD _PlusOutFlow, LS_STAG_CR_PC _delta, const double _nu, const int _N,
			const int _M, const int _NM)
		{
			const int numE = i + j * _M;
			const coordinate* curD = _delta + ((i + j + 1) * _N + j * _NM -1);
			_PlusOutFlow[numE] = _nu * (curD->y + (curD + _N)->y) * curD->z / curD->x;
		}
		using HelmholtzDiscreteOperatorBase::initMeshParams;
		void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			HelmholtzDiscreteOperatorBase::initMeshParams(_cur_mesh, _meshes);
			m_My = m_meshes->yMesh.VOF;
			build();
		};
	  private:
		void calcPlusOutFlow()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorY, setPlusOutFlowY, m_M, m_K, PlusOutFlow, m_delta, m_nu, m_N, m_M, m_N * m_M)
		}
		void initInFlowFlux()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorY, inFlowFluxHelmY, m_M, m_K, flux, m_delta, m_task->V_W, 
				m_task->V_E, false, m_N, m_N * m_M, m_nu)
			LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorY, addFluxNS, m_K, flux, m_delta, m_task->V_S, m_task->V_N, 
				m_N, m_N * m_M, m_nu)
			if(!existOutFlow) LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorY, addFluxTB, m_M, flux, m_type, m_delta, m_task->V_B, m_task->V_T,
				m_N, m_K - 1, NMK - m_N * m_M, m_nu)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора Гельмгольца на у-сетке на погруженных границах 
		  */
		void ibFlux()
		{// испр. 31.05.14
			LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteOperatorY, ibFluxHelmY, 1, M2, 0, m_K, flux, m_type, m_baseSideTypes,
				m_theta, m_alpha, m_delta, m_inv_delta_ib_ne, m_inv_delta_ib_tn, m_delta_ib_en, m_delta_ib_ws, 
				m_delta_ib_nb, m_delta_ib_st, m_Vib, m_VibZY, m_VibYX, m_inv_V, m_N, m_N * m_M, m_N * M3, m_nu)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора Гельмгольца на у-сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorY, calcHelmY, m_M, m_K, Number, Cell, m_theta, m_alpha, m_delta, 
				m_inv_delta_ib_ne, m_inv_delta_ib_tn, m_delta_ib_en, m_delta_ib_ws, m_delta_ib_nb, m_delta_ib_st, m_inv_V, 
				m_My, m_type, m_N, m_M, m_N * m_M, m_nu, coefIntegr, (m_task->id == kind::TASK_FLOW))
		   addOutFlow();
		};
	  private:
		double *m_My;///< массив значений объемов ячеек с у-сетки
	};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Y_H_
