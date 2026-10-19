/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_operator_x.h
* \author    Puzikova Valeria (vvp) 
* \date      19.10.2013 16:54:36
* \brief     helmholtz_discrete_operator_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_X_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_X_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_operator_base.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteOperatorX

Данный файл содержит в себе реализацию класса HelmholtzDiscreteOperatorX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(calcHelmX, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, const coordinate*, _alpha, const coordinate*, _delta, const coordinate*, _inv_delta_ib, const coordinate*, _delta_ib_en,
		const coordinate*, _delta_ib_ws, LS_STAG_CR_PD, _inv_V, LS_STAG_CR_PD, _Mx, LS_STAG_CR_PCT, _type,
		const int, _N, const int, _M, const int, _NM, const double, _nu, const double, _coefIntegr, const bool, _isFlow);
	LS_STAG_DECLARATION_R_KERNEL_2D(ibFluxHelmX, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType,
		const coordinate*, _theta, const coordinate*, _alpha, LS_STAG_CR_PC, _delta, const coordinate*, _inv_delta_ib, 
		const coordinate*, _delta_ib_en, const coordinate*, _delta_ib_ws, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PC, _VibZY, 
		 LS_STAG_CR_PC, _VibYX, LS_STAG_CR_PD, _inv_V, const int, _N, const int, _M, const int, _NM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(setPlusOutFlow, LS_STAG_R_PD, _PlusOutFlow, LS_STAG_CR_PC, _delta,
		const double, _nu, const int, _M, const int, _N3);
	LS_STAG_DECLARATION_KERNEL(addFluxNS, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, U_S, const double, U_N, const int, _N, const int, _M, const int, _NM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL(addFluxTB, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, U_B, const double, U_T, const int, _N, const int, _M, const int, _K2, const int, _NMK2, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(inFlowFluxHelmX, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		LS_STAG_CR_PD, _U_W, const double, _U_E, const bool, isArtery, const bool, isCavity, const int, _N, const int, _M, 
		const int, _NM, const int, _vs, const int, _vf, const int, _zs, const int, _zf, const double, _nu);
/*!
\brief Класс, реализующий дискретный оператор Гельмгольца на х-сетке

Данный класс реализует работу с дискретным аналогом оператора Гельмгольца на х-сетке
*/
	struct HelmholtzDiscreteOperatorX : HelmholtzDiscreteOperatorBase
	{//m_param_x = Mx, m_param_y = 1/V	
		
	LS_STAG_KERNEL(addFluxNS, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double U_S, const double U_N, 
		const int _N, const int _M, const int _NM, const double _nu)
		{
			int numS = i * _NM, numN = numS + _NM - _N, im = i * _M;
			coordinate* curD = m_delta + (numS + im);
			double coefS = _nu * U_S * curD->z / curD->y;
			double coefN = _nu * U_N * curD->z / (m_delta + (numN + im))->y;
			double dx = curD->x; curD++;
			for(int ii = 0; ii < _N; ii++, numS++, numN++, curD++)
			{
				double dxx = curD->x, coef = dx + dxx;
				flux[numS] += coef * coefS;
				flux[numN] += coef * coefN;
				dx = dxx;
			}
		}
	LS_STAG_KERNEL(addFluxTB, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta,  const double U_B, const double U_T,
		const int _N, const int _M, const int _K2, const int _NMK2, const double _nu)
		{
			int numB = i * _N, numT = _NMK2 + numB;
			coordinate* curD = m_delta + (numB + i);
			double coefB = _nu * U_B * curD->y / curD->z;
			double coefT = _nu * U_T * curD->y / (m_delta + (numT + i + _K2 * _M))->z;
			double dx = curD->x; curD++;
			for(int ii = 0; ii < _N; ii++, numB++, numT++, curD++)
			{
				double dxx = curD->x, coef = dx + dxx;
				if (_type[numB] != SolidCell) flux[numB] += coef * coefB;
				if (_type[numT] != SolidCell) flux[numT] += coef * coefT;
				dx = dxx;
			}
		}
		LS_STAG_KERNEL_2D(inFlowFluxHelmX, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, LS_STAG_CR_PD _U_W, const double _U_E,
			const bool isArtery, const bool isCavity, const int _N, const int _M, const int _NM, const int _vs, 
			const int _vf, const int _zs, const int _zf, const double _nu)
		{
			int numF = _N * i + _NM * j, _N2 = _N - 1, iM = j * _M;
			coordinate* curD = _delta + (numF + iM + i);
			double dyz = _nu * curD->z * curD->y;
			_flux[numF] = ((isArtery && ((i <= _vs) || (i >= _vf) || ((j < _zs) || (j > _zf))))) ? 0.0 :
			             (_U_W[i + iM] * dyz / curD->x);
			if (isCavity) 
			{
			numF += _N2;
			_flux[numF] = _U_E * dyz / (curD + _N)->x;
			}
		}	
		LS_STAG_KERNEL_2D(calcHelmX, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate* _theta, const coordinate* _alpha,
			const coordinate* _delta, const coordinate* _inv_delta_ib, const coordinate* _delta_ib_en, 
			const coordinate* _delta_ib_ws, LS_STAG_CR_PD _inv_V, LS_STAG_CR_PD _Mx, LS_STAG_CR_PCT _type, const int _N,
			 const int _M, const int _NM, const double _nu, const double _coefIntegr, const bool _isFlow)
		{
			for (int numX = _N * i + _NM * j, m_end = numX + _N, num = numX + i + j * _M, baseN = _N + 1, 
				baseNM = baseN * _M; numX < m_end; numX++, num++)
			{
				int numA2 = _Number[numX];
				switch (_type[numX])
				{
				case FluidCell:
				{
					const coordinate* curD = _delta + num;
					const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dxx = _nu * ((curD + 1)->x + curD->x);
					double aux = dz * dxx;
					double k_s = aux / (dy + (curT - baseN)->x * (curD - baseN)->y);
					double k_n = aux / (dy + (curT + baseN)->x * (curD + baseN)->y);
					aux = dy * dxx;
					double k_b = aux / (dz + (curT - baseNM)->x * (curD - baseNM)->z);
					double k_t = aux / (dz + (curT + baseNM)->x * (curD + baseNM)->z);
					aux = dy * dz;
					dxx = _nu * aux * aux;
					double k_e = dxx * _inv_V[num+1];
					double k_w = dxx * _inv_V[num];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w * (curT - 1)->x;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - k_w - k_e - k_n - k_t;
					_Cell[numA2 + 4] = k_e * (curT + 1)->x;
					_Cell[numA2 + 5] = k_n;
					_Cell[numA2 + 6] = k_t;
				}
				break;
				case SolidCell:
					_Cell[numA2] = PLUS_TINY;
					break;
				case N_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num; 
					const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
					double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx  = curT->x;
					double k_s = _nu * dz * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
					double aux = dy * _nu, aux2 = tx * dz;
					double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
					double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
					double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + eT->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
					aux = dy * dz;
					dxx = _nu * tx * aux * aux;
					double k_e = dxx * _inv_V[num+1];
					double k_w =  dxx * _inv_V[num];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w * (curT - 1)->x;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_t - k_wallN;
					_Cell[numA2 + 4] = k_e * eT->x;
					_Cell[numA2 + 5] = k_t;
				}
				break;
				case S_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
    double aux = dy * _nu, aux2 = tx * dz;
    double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
    double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * dxx) * curA->z / (tx * dy)) : 0.0;
    aux = dy * dz;
    dxx = _nu * tx * aux * aux;
    double k_e = dxx * _inv_V[num + 1];
    double k_w = dxx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_w * (curT - 1)->x;
    _Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - (k_w + k_e) * tx - k_n - k_t - k_wallS;
    _Cell[numA2 + 3] = k_e * eT->x;
    _Cell[numA2 + 4] = k_n;
    _Cell[numA2 + 5] = k_t;
}
break;
case E_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const  coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
    double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
    double aux = dy * _nu, aux2 = tx * dz;
    double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
    dy *= dz;
    double k_w = _nu * tx * dy * dy * _inv_V[num];
    double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s;
    _Cell[numA2 + 2] = k_w * (curT - 1)->x;
    _Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_n - k_t;
    _Cell[numA2 + 4] = k_n;
    _Cell[numA2 + 5] = k_t;
}
break;
case W_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const   coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const  coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
    double aux = dz * _nu, aux2 = tx * dy;
    double k_s = aux * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (aux2 + sT->x * (curD - baseN)->y);
    double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
    aux = dy * _nu, aux2 = tx * dz;
    double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
    double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
    double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s;
    _Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - k_e * tx - k_n - k_t;
    _Cell[numA2 + 3] = k_e * eT->x;
    _Cell[numA2 + 4] = k_n;
    _Cell[numA2 + 5] = k_t;
}
break;
case T_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const   coordinate* eT = curT + 1; const  coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
    double aux = dz * _nu, aux2 = tx * dy, az = (curA - baseN)->z;
    bool existIb = (az < EQ_DOUBLE);
    double k_s = (existIb) ? 0.0 : (aux * (sT->y * dx + (sT + 1)->y * dxx) * az / (aux2 + sT->x * (curD - baseN)->y));
    az = curA->z; bool existIb2 = (az < EQ_DOUBLE);
    double k_n = (existIb2) ? 0.0 : (aux * (curT->y * dx + eT->y * dxx) * az / (aux2 + (curT + baseN)->x * (curD + baseN)->y));
// tex: s^{ib,x,y} при alpha^z_{i,j-1,k}=0 ИЛИ alpha^z_{i,j,k}=0 --- дельта своей и E-ячейки
    double ib = (existIb || existIb2) ? (_nu * (_delta_ib_en[num].y + _delta_ib_ws[num+1].y) * _inv_delta_ib[num].y) : 0.0;
    aux = dy * _nu, aux2 = tx * dz;
    az = (curA - baseNM)->y; existIb = (az < EQ_DOUBLE);
    double k_b = (existIb) ? 0.0 : (aux * (bT->z * dx + (bT + 1)->z * dxx) * az / (aux2 + bT->x * (curD - baseNM)->z));
    if(existIb || (curA->y < EQ_DOUBLE)) ib += _nu * (_delta_ib_en[num].x + _delta_ib_ws[num+1].x) * _inv_delta_ib[num].z;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + eT->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
    aux = dy * dz;
    dxx = _nu * tx * aux * aux;
    double k_e = dxx * _inv_V[num + 1];
    double k_w = dxx * _inv_V[num];
    _Cell[numA2]     = k_b;
    _Cell[numA2 + 1] = k_s;
    _Cell[numA2 + 2] = k_w * (curT - 1)->x;
    _Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_n - k_wallT - ib;
    _Cell[numA2 + 4] = k_e * eT->x;
    _Cell[numA2 + 5] = k_n;
}
break;
case B_Border:
{
    const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
	const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN;
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
    double aux = dz * _nu, aux2 = tx * dy, az = (curA - baseN)->z;
    bool existIb = (az < EQ_DOUBLE);
    double k_s = (existIb) ? 0.0 : (aux * (sT->y * dx + (sT + 1)->y * dxx) * az / (aux2 + sT->x * (curD - baseN)->y));
    az = curA->z; bool existIb2 = (az < EQ_DOUBLE);
    double k_n = (existIb2) ? 0.0 : (aux * (curT->y * dx + eT->y * dxx) * az / (aux2 + (curT + baseN)->x * (curD + baseN)->y));
// tex: s^{ib,x,y} при alpha^z_{i,j-1,k}=0 ИЛИ alpha^z_{i,j,k}=0 --- дельта своей и E-ячейки
    double ib = (existIb || existIb2) ? (_nu * (_delta_ib_en[num].y + _delta_ib_ws[num+1].y) * _inv_delta_ib[num].y) : 0.0;
    aux = dy * _nu, aux2 = tx * dz;
    az = curA->y; existIb2 = (az < EQ_DOUBLE);
    double k_t = (existIb2) ? 0.0 : (aux * (curT->z * dx + eT->z * dxx) * az / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z));
    if(existIb2) ib += _nu * (_delta_ib_en[num].x + _delta_ib_ws[num+1].x) * _inv_delta_ib[num].z;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + eT->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
    aux = dy * dz;
    dxx = _nu * tx * aux * aux;
    double k_e = dxx * _inv_V[num + 1];
    double k_w = dxx * _inv_V[num];
    _Cell[numA2]     = k_s;
    _Cell[numA2 + 1] = k_w * (curT - 1)->x;
    _Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_s - (k_w + k_e) * tx - k_n - k_t - k_wallB - ib;
    _Cell[numA2 + 3] = k_e * eT->x;
    _Cell[numA2 + 4] = k_n;
    _Cell[numA2 + 5] = k_t;
}
break;
		case TN_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
dy *= dz;
dxx = _nu * tx * dy * dy;
double k_e = dxx * _inv_V[num + 1];
double k_w = dxx * _inv_V[num];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
_Cell[numA2 + 2] = k_w * (curT - 1)->x;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? ((dy / dz) * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * (dy / dz))) : 0.0;
_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_wallT - k_wallN;
_Cell[numA2 + 4] = k_e * eT->x;
}
break;
case TS_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
double aux = dz * _nu, aux2 = tx * dy;
double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
aux = dy * dz;
dxx = _nu * tx * aux * aux;
double k_e = dxx * _inv_V[num + 1];
double k_w = dxx * _inv_V[num];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - (k_w + k_e) * tx - k_n - k_wallT - k_wallS;
_Cell[numA2 + 3] = k_e * eT->x;
_Cell[numA2 + 4] = k_n;
}
break;
case TE_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double aux = dz * _nu, aux2 = tx * dy;
double k_s = aux * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (aux2 + sT->x * (curD - baseN)->y);
double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
_Cell[numA2 + 2] = k_w * (curT - 1)->x;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_n - k_wallT;
_Cell[numA2 + 4] = k_n;
}
break;
case TW_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double aux = dz * _nu, aux2 = tx * dy;
double k_s = aux * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (aux2 + sT->x * (curD - baseN)->y);
double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - k_e * tx - k_n - k_wallT;
_Cell[numA2 + 3] = k_e * eT->x;
_Cell[numA2 + 4] = k_n;
}
break;
case BN_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
dy *= dz;
dxx = _nu * tx * dy * dy;
double k_e = dxx * _inv_V[num + 1];
double k_w = dxx * _inv_V[num];
_Cell[numA2] = k_s;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_s - (k_w + k_e) * tx - k_t - k_wallB - k_wallN;
_Cell[numA2 + 3] = k_e * eT->x;
_Cell[numA2 + 4] = k_t;
}
break;
case BS_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
dy *= dz;
dxx = _nu * tx * dy * dy;
double k_e = dxx * _inv_V[num + 1];
double k_w = dxx * _inv_V[num];
_Cell[numA2] = k_w * (curT - 1)->x;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - (k_w + k_e) * tx - k_n - k_t - k_wallB - k_wallS;
_Cell[numA2 + 2] = k_e * eT->x;
_Cell[numA2 + 3] = k_n;
_Cell[numA2 + 4] = k_t;
}
break;
case BE_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double aux = dz * _nu, aux2 = tx * dy;
double k_s = aux * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (aux2 + sT->x * (curD - baseN)->y);
double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_s;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_s - (k_w + k_e) * tx - k_n - k_t - k_wallB;
_Cell[numA2 + 3] = k_n;
_Cell[numA2 + 4] = k_t;
}
break;
case BW_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double aux = dz * _nu, aux2 = tx * dy;
double k_s = aux * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (aux2 + sT->x * (curD - baseN)->y);
double k_n = aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_s;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - k_s - k_e * tx - k_n - k_t - k_wallB;
_Cell[numA2 + 2] = k_e * eT->x;
_Cell[numA2 + 3] = k_n;
_Cell[numA2 + 4] = k_t;
}
break;
case NE_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double aux = dy * _nu, aux2 = tx * dz;
double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
double k_t = aux * (curT->z * dx + (curT + 1)->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
_Cell[numA2 + 2] = k_w * (curT - 1)->x;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_t - k_wallN;
_Cell[numA2 + 4] = k_t;
}
break;
case NW_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double aux = dy * _nu, aux2 = tx * dz;
double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - k_e * tx - k_t - k_wallN;
_Cell[numA2 + 3] = k_e * eT->x;
_Cell[numA2 + 4] = k_t;
}
break;
case SE_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double aux = dy * _nu, aux2 = tx * dz;
double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - (k_w + k_e) * tx - k_n - k_t - k_wallS;
_Cell[numA2 + 3] = k_n;
_Cell[numA2 + 4] = k_t;
}
break;
case SW_Edge:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double aux = dy * _nu, aux2 = tx * dz;
double k_b = aux * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (aux2 + bT->x * (curD - baseNM)->z);
double k_t = aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - k_b - k_e * tx - k_n - k_t - k_wallS;
_Cell[numA2 + 2] = k_e * eT->x;
_Cell[numA2 + 3] = k_n;
_Cell[numA2 + 4] = k_t;
}
break;
// ---------- Углы (3 элемента) ----------
case TNE_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
_Cell[numA2 + 2] = k_w * (curT - 1)->x;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_wallT - k_wallN;
}
break;
case TNW_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_s;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - k_e * tx - k_wallT - k_wallN;
_Cell[numA2 + 3] = k_e * eT->x;
}
break;
case TSE_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_b - (k_w + k_e) * tx - k_n - k_wallT - k_wallS;
_Cell[numA2 + 3] = k_n;
}
break;
case TSW_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* bT = curT - baseNM;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double k_b = dy * _nu * (bT->z * dx + (bT + 1)->z * dxx) * (curA - baseNM)->y / (tx * dz + bT->x * (curD - baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_b;
double k_wallT = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - k_b - k_e * tx - k_n - k_wallT - k_wallS;
_Cell[numA2 + 2] = k_e * eT->x;
_Cell[numA2 + 3] = k_n;
}
break;
case BNE_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* sT = curT - baseN;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_t = dy * _nu * (curT->z * dx + (curT+1)->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_s;
_Cell[numA2 + 1] = k_w * (curT - 1)->x;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 2] = -(_coefIntegr * _Mx[numX]) - k_s - (k_w + k_e) * tx - k_t - k_wallB - k_wallN;
_Cell[numA2 + 3] = k_t;
}
break;
case BNW_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_s = dz * _nu * (sT->y * dx + (sT + 1)->y * dxx) * (curA - baseN)->z / (tx * dy + sT->x * (curD - baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
    dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_s;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallN = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - k_s - k_e * tx - k_t - k_wallB - k_wallN;
_Cell[numA2 + 2] = k_e * eT->x;
_Cell[numA2 + 3] = k_t;
}
break;
case BSE_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
dy *= dz;
double k_w = _nu * tx * dy * dy * _inv_V[num];
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
_Cell[numA2] = k_w * (curT - 1)->x;
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2 + 1] = -(_coefIntegr * _Mx[numX]) - (k_w + k_e) * tx - k_n - k_t - k_wallB - k_wallS;
_Cell[numA2 + 2] = k_n;
_Cell[numA2 + 3] = k_t;
}
break;
case BSW_Corner:
{
const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
const coordinate* eT = curT + 1;
double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx = curT->x;
double k_n = dz * _nu * (curT->y * dx + eT->y * dxx) * curA->z / (tx * dy + (curT + baseN)->x * (curD + baseN)->y);
double k_t = dy * _nu * (curT->z * dx + eT->z * dxx) * curA->y / (tx * dz + (curT + baseNM)->x * (curD + baseNM)->z);
dy *= dz;
double k_e = _nu * tx * dy * dy * _inv_V[num + 1];
double k_wallB = _isFlow ? 0.0 : (tx > EQ_DOUBLE) ? (curD->y * _nu * (curT->z * dx + (curT + 1)->z * (curD + 1)->x) * curA->y / (tx * dz)) : 0.0;
double k_wallS = (tx > EQ_DOUBLE) ? (_nu * dz * (curT->y * dx + (curT + 1)->y * (curD + 1)->x) * curA->z / (tx * curD->y)) : 0.0;
_Cell[numA2] = -(_coefIntegr * _Mx[numX]) - k_e * tx - k_n - k_t - k_wallB - k_wallS;
_Cell[numA2 + 1] = k_e * eT->x;
_Cell[numA2 + 2] = k_n;
_Cell[numA2 + 3] = k_t;
}
break;
				default:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num; 
					const coordinate* eT = curT + 1; const coordinate* sT = curT - baseN; const coordinate* bT = curT - baseNM;
					double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD + 1)->x, tx  = curT->x;
					double aux = dz * _nu, aux2 = tx * dy, az = (curA - baseN)->z;
					const bool existIbS = (az < EQ_DOUBLE);
					double k_s = (existIbS) ? 0.0 : (aux * (sT->y * dx + (sT + 1)->y * dxx) * az / (aux2 + sT->x * (curD - baseN)->y));
					az = curA->z; const bool existIbN = (az < EQ_DOUBLE);
					double k_n = (existIbN) ? 0.0 : (aux * (curT->y * dx + eT->y * dxx) * curA->z / (aux2 + (curT + baseN)->x * (curD + baseN)->y));
// tex: s^{ib,x,y}_{i,j,k}=(dS^{ib,e}_{xz,i,j,k}+dS^{ib,w}_{xz,i+1,j,k})/dy^{ib,e}_{i,j,k} при alpha^z_{i,j-1,k}=0 ИЛИ alpha^z_{i,j,k}=0
					double ib = (existIbS || existIbN) ? (_nu * (_delta_ib_en[num].y + _delta_ib_ws[num+1].y) * _inv_delta_ib[num].y) : 0.0;
					aux = dy * _nu, aux2 = tx * dz;
					az = (curA - baseNM)->y; const bool existIbB = (az < EQ_DOUBLE);
					double k_b = (existIbB) ? 0.0 : (aux * (bT->z * dx + (bT + 1)->z * dxx) * az / (aux2 + bT->x * (curD - baseNM)->z));
					az = curA->y; const bool existIbT = (az < EQ_DOUBLE);
					double k_t = (existIbT) ? 0.0 : (aux * (curT->z * dx + eT->z * dxx) * curA->y / (aux2 + (curT + baseNM)->x * (curD + baseNM)->z));
// tex: s^{ib,x,z}_{i,j,k}=(dS^{ib,e}_{xy,i,j,k}+dS^{ib,w}_{xy,i+1,j,k})/dz^{ib,e}_{i,j,k} при alpha^y_{i,j,k-1}=0 ИЛИ alpha^y_{i,j,k}=0
					if (existIbB || existIbT) ib += _nu * (_delta_ib_en[num].x + _delta_ib_ws[num+1].x) * _inv_delta_ib[num].z;
					aux = dy * dz;
					dxx = _nu * tx * aux * aux;
					double k_e = dxx * _inv_V[num+1];
					double k_w =  dxx * _inv_V[num];
					_Cell[numA2] = k_b;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w * (curT - 1)->x;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mx[numX]) - k_b - k_s - (k_w + k_e) * tx - k_n - k_t - ib;
					_Cell[numA2 + 4] = k_e * eT->x;
					_Cell[numA2 + 5] = k_n;
					_Cell[numA2 + 6] = k_t;
				}
				}//switch
			}//for i
		}
		LS_STAG_KERNEL_2D(ibFluxHelmX, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _SideType, const coordinate* _theta, const coordinate* _alpha,
			LS_STAG_CR_PC _delta, LS_STAG_CR_PC _inv_delta_ib, const coordinate* _delta_ib_en, const coordinate* _delta_ib_ws, 
			LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibZY, LS_STAG_CR_PC _VibYX, LS_STAG_CR_PD _inv_V, const int _N, const int _M, 
			const int _NM, const double _nu)
		{
			for (int numX = i * _N + j * _NM, num = numX + i + j * _M, baseN = _N + 1, baseNM = _NM + _M, bNum = (j) ? baseNM : 0, m_end = numX + _N; numX < m_end; numX++, num++)//for(i = 0; i < m_N; i++)
			{
				const CellTypes curType = _type[numX];
				bool isCutCell = false;
				if((curType > SolidCell) && (curType < FluidCell))
				isCutCell = true;
				else
				{
					switch(curType)
					{
						case B_Border:
						case T_Border:
                        if ((_SideType[num].inRange(SolidSide, FluidSide)) || (_SideType[num + 1].inRange(SolidSide, FluidSide))) isCutCell = true;
					}
				}
				if (isCutCell)
				{
					const coordinate* curT = _theta + num; double tx = curT->x;
					if (tx < EQ_DOUBLE)
						_flux[numX] = 0.0;
					else
					{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curVib = _Vib + num;
					double dyz = curD->y * curD->z;
					double ib = ((curA->z < EQ_DOUBLE) || ((curA - baseN)->z < EQ_DOUBLE)) ? 
					     ((_delta_ib_en[num].y + _delta_ib_ws[num+1].y) * _inv_delta_ib[num].y * _VibYX[num].x) : 0.0;
					if((curA->y < EQ_DOUBLE) || ((curA - bNum)->y < EQ_DOUBLE)) 
					     ib += (_delta_ib_en[num].x + _delta_ib_ws[num+1].x) * _inv_delta_ib[num].z * _VibZY[num].x;
					_flux[numX] = _nu * (ib + tx * dyz * dyz * (((tx - (curT + 1)->x) * (curVib + 1)->x * _inv_V[num + 1]) - 
					                  (((curT - 1)->x - tx) * curVib->x * _inv_V[num])));
					}
				}// Grid.TypeOfCellX[numX] > 9
				else _flux[numX] = 0.0;
			}//for i
		}
		LS_STAG_KERNEL_2D(setPlusOutFlow, LS_STAG_R_PD _PlusOutFlow, LS_STAG_CR_PC _delta, const double _nu, const int _M, const int _N3)
		{
			const int numE = i + j * _M;
			const coordinate* curD = _delta - 1 + (numE + 1) * _N3;
			_PlusOutFlow[numE] = _nu * curD->y * curD->z / curD->x;
		}
		using HelmholtzDiscreteOperatorBase::initMeshParams;
		void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			HelmholtzDiscreteOperatorBase::initMeshParams(_cur_mesh, _meshes);
			m_Mx = m_meshes->xMesh.VOF;
			build();
		};
	protected:
		void calcPlusOutFlow()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorX, setPlusOutFlow, m_M, m_K, PlusOutFlow, m_delta, m_nu, m_M, N3)
		}
		void initInFlowFlux()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorX, inFlowFluxHelmX, m_M, m_K, flux, m_delta, m_task->pU_W, 
				m_task->U_E, false, false, 
				m_N, m_M, m_N * m_M, vs, vf, zs, zf, m_nu)
			LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorX, addFluxNS, m_K, flux, m_delta, m_task->U_S, m_task->U_N, 
				m_N, m_M, m_N * m_M, m_nu)
			if(!existOutFlow) LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorX, addFluxTB, m_M, flux, m_type, m_delta, m_task->U_B, m_task->U_T,
				m_N, m_M, m_K - 1, NMK - m_N * m_M, m_nu)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора Гельмгольца на х-сетке на погруженных границах 
		  */
		void ibFlux()
		{
			LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteOperatorX, ibFluxHelmX, 1, M2, 0, m_K, flux, m_type, m_baseSideTypes,
				m_theta, m_alpha, m_delta, m_inv_delta_ib_ne, m_delta_ib_en, m_delta_ib_ws, m_Vib, m_VibZY, m_VibYX, 
				m_inv_V, m_N, m_M, m_N * m_M, m_nu)
		}
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора Гельмгольца на х-сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorX, calcHelmX, m_M, m_K, Number, Cell, m_theta, m_alpha, m_delta, 
				m_inv_delta_ib_ne, m_delta_ib_en, m_delta_ib_ws, m_inv_V, m_Mx, m_type, m_N, m_M, m_N * m_M, m_nu, coefIntegr,
				(m_task->id == kind::TASK_FLOW))
			addOutFlow();
		};
		double *m_Mx;///< массив значений объемов ячеек с х-сетки
	};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_X_H_
