/*
* \copyright (C) 2026 Puzikova V.V.
* \file      helmholtz_discrete_operator_z.h
* \author    Puzikova Valeria (vvp) 
* \date      15.06.2026 11:13:23
* \brief     helmholtz_discrete_operator_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Z_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Z_H_

// ================================================================== INCLUDES
#include "helmholtz_discrete_operator_base.h"
// ================================================================== 
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteOperatorZ

Данный файл содержит в себе реализацию класса HelmholtzDiscreteOperatorZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(calcHelmZ, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, const coordinate*, _alpha, const coordinate*, _delta, 
		const coordinate*, _inv_delta_ib_tn, const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, 
		LS_STAG_CR_PD, _inv_V, LS_STAG_CR_PD, _Mz, LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, 
		const double, _nu, const double, _coefIntegr, const bool, _isFlow);
	LS_STAG_DECLARATION_R_KERNEL_2D(ibFluxHelmZ, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType,
		const coordinate*, _theta, const coordinate*, _alpha, LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _inv_delta_ib_tn, 
        const coordinate*, _delta_ib_nb, const coordinate*, _delta_ib_st, LS_STAG_CR_PC, _Vib, LS_STAG_CR_PC, _VibZY, 
		 LS_STAG_CR_PC, _VibYX, LS_STAG_CR_PD, _inv_V, const int, _N, const int, _NM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(setPlusOutFlowZ, LS_STAG_R_PD, _PlusOutFlow, LS_STAG_CR_PC, _delta,
		const double, _nu, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxNShz, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta,
		const double, V_S, const double, V_N, const int, _N, const int, _NM, const double, _nu);
	LS_STAG_DECLARATION_KERNEL(addFluxTBhz, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, V_B, const double, V_T, const int, _N, const int, _NMK2, const double, _nu);
	LS_STAG_DECLARATION_KERNEL_2D(inFlowFluxhz, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta,
		const double, _V_W, const double, _V_E, const bool, isCavity, const int, _N, const int, _NM, const double, _nu);
/*!
\brief Класс, реализующий дискретный оператор Гельмгольца на z-сетке

Данный класс реализует работу с дискретным аналогом оператора Гельмгольца на z-сетке
*/
	struct HelmholtzDiscreteOperatorZ : HelmholtzDiscreteOperatorBase
	{//m_param_x = My, m_param_y = 1/V
		LS_STAG_KERNEL(addFluxNShz, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double W_S, const double W_N,
			const int _N, const int _NM, const double _nu)
		{
			int numS = i * _NM, numN = numS + _NM - _N;
			coordinate* curD = m_delta + numS;
			double dzz = (curD->z + (curD + _NM)->z);
			double coefS = _nu * W_S * dzz / curD->y;
			double coefN = _nu * W_N * dzz / (m_delta + numN)->y;
			for (int ii = 0; ii < _N; ii++, numS++, numN++, curD++)
			{
				flux[numS] += curD->x * coefS;
				flux[numN] += curD->x * coefN;
			}
		};
		LS_STAG_KERNEL(addFluxTBhz, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta, const double W_B, const double W_T,
			const int _N, const int _NMK2, const double _nu)
		{
			int numB = i * _N, numT = _NMK2 + numB;
			coordinate* curD = m_delta + numB;
			double coefB = _nu * W_B * curD->y / curD->z;
			double coefT = _nu * W_T * curD->y / (m_delta + numT)->z;
			for (int ii = 0; ii < _N; ii++, numB++, numT++, curD++)
			{
				if (_type[numB] != SolidCell) flux[numB] += curD->x * coefB;
				if (_type[numT] != SolidCell) flux[numT] += curD->x * coefT;
			}
		};
		LS_STAG_KERNEL_2D(inFlowFluxhz, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, const double _W_W, const double _W_E,
			const bool isCavity, const int _N, const int _NM, const double _nu)
		{
			int numF = _N * i + _NM * j, _N2 = _N - 1;
			coordinate* curD = _delta + numF;
			double dyz = _nu * curD->y * (curD->z + (curD + _NM)->z);
			_flux[numF] = _W_W * dyz / curD->x;
			if (isCavity)
			{
				numF += _N2;
				_flux[numF] = _W_E * dyz / (curD + _N2)->x;
			}
		};
		LS_STAG_KERNEL_2D(calcHelmZ, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate* _theta,
			const coordinate* _alpha, const coordinate* _delta, const coordinate* _inv_delta_ib_tn,
			const coordinate* _delta_ib_nb, const coordinate* _delta_ib_st, LS_STAG_CR_PD _inv_V, LS_STAG_CR_PD _Mz,
			LS_STAG_CR_PCT _type, const int _N, const int _M, const int _NM, const double _nu, const double _coefIntegr,
			const bool _isFlow)
		{
			for (int num = _N * i + _NM * j, m_end = num + _N; num < m_end; num++)
			{
				int numA2 = _Number[num];
				switch (_type[num])
				{
				case FluidCell:
				{
					const coordinate* curD = _delta + num;
					const coordinate* curT = _theta + num;
					double dy = curD->y, dx = curD->x, dzz = _nu * ((curD + _NM)->z + curD->z);
					double aux = dzz * dy;
					double k_w = aux / (dx + (curT - 1)->z * (curD - 1)->x);
					double k_e = aux / (dx + (curT + 1)->z * (curD + 1)->x);
					aux = dx * dzz;
					double k_s = aux / (dy + (curT - _N)->z * (curD - _N)->y);
					double k_n = aux / (dy + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					dzz = _nu * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - k_b - k_s - k_w - k_e - k_n - k_t;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_n;
					_Cell[numA2 + 6] = k_t * (curT+_NM)->z;
				}
				break;
				case SolidCell:
					_Cell[numA2] = PLUS_TINY;
					break;
				case N_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const  coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (aux2 + wT->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = dx * _nu * (curT->y * dz + tT->y * dzz) * curA->x / (tz * dy);
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - k_s - k_w - (k_b + k_t) * tz - k_e - k_wallN;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_t * tT->z;
				}
				break;
				case S_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const  coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (aux2 + wT->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy);
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_w - k_e - k_n - k_wallS;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_n;
					_Cell[numA2 + 5] = k_t * tT->z;
				}
				break;
				case E_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const  coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (aux2 + sT->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + tT->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_s - k_w - k_e - k_n;
					_Cell[numA2 + 4] = k_n;
					_Cell[numA2 + 5] = k_t * tT->z;
				}
				break;
				case W_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const 	coordinate* tT = curT + _NM; const  coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (aux2 + sT->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + tT->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_s - k_e - k_n;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_n;
					_Cell[numA2 + 5] = k_t * tT->z;
				}
				break;
				case T_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* wT = curT - 1; const  coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx, az = (curA - 1)->y;
					bool existIb = (az < EQ_DOUBLE);
					double k_w = (existIb) ? 0.0 : (aux * (wT->x * dz + (wT + _NM)->x * dzz) * az / (aux2 + wT->z * (curD - 1)->x));
					az = curA->y; bool existIb2 = (az < EQ_DOUBLE);
					double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dz + (curT + _NM)->x * dzz) * az / (aux2 + (curT + 1)->z * (curD + 1)->x));
// tex: s^{ib,z,x} = (dS^{ib,t}_{yz}(i,j,k) + dS^{ib,b}_{yz}(i,j,k+1))/dx^{ib,t} при alpha^y_{i-1,j,k}=0 ИЛИ alpha^y_{i,j,k}=0
					double ib = (existIb || existIb2) ? (_nu * (_delta_ib_st[num].z + _delta_ib_nb[num + _NM].z) * _inv_delta_ib_tn[num].x) : 0.0;
					aux = dx * _nu, aux2 = tz * dy;
					az = (curA - _N)->x; existIb = (az < EQ_DOUBLE);
					double k_s = (existIb) ? 0.0 : (aux * (sT->y * dz + (sT + _NM)->y * dzz) * az / (aux2 + sT->z * (curD - _N)->y));
					az = curA->x; existIb2 = (az < EQ_DOUBLE);
					double k_n = (existIb2) ? 0.0 : (aux * (curT->y * dz + (curT + _NM)->y * dzz) * az / (aux2 + (curT + _N)->z * (curD + _N)->y));
// tex: s^{ib,z,y} = (dS^{ib,t}_{xz}(i,j,k) + dS^{ib,b}_{xz}(i,j,k+1))/dy^{ib,t} при alpha^x_{i,j-1,k}=0 ИЛИ alpha^x_{i,j,k}=0
					if (existIb || existIb2) ib += _nu * (_delta_ib_st[num].y + _delta_ib_nb[num + _NM].y) * _inv_delta_ib_tn[num].y;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_w - k_e - k_n - ib;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_n;
				}
				break;
				case B_Border:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* wT = curT - 1; const  coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx, az = (curA - 1)->y;
					bool existIb = (az < EQ_DOUBLE);
					double k_w = (existIb) ? 0.0 : (aux * (wT->x * dz + (wT + _NM)->x * dzz) * az / (aux2 + wT->z * (curD - 1)->x));
					az = curA->y; bool existIb2 = (az < EQ_DOUBLE);
					double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dz + tT->x * dzz) * az / (aux2 + (curT + 1)->z * (curD + 1)->x));
// tex: s^{ib,z,x} = (dS^{ib,t}_{yz}(i,j,k) + dS^{ib,b}_{yz}(i,j,k+1))/dx^{ib,t} при alpha^y_{i-1,j,k}=0 ИЛИ alpha^y_{i,j,k}=0
					double ib = (existIb || existIb2) ? (_nu * (_delta_ib_st[num].z + _delta_ib_nb[num + _NM].z) * _inv_delta_ib_tn[num].x) : 0.0;
					aux = dx * _nu, aux2 = tz * dy;
					az = (curA - _N)->x; existIb = (az < EQ_DOUBLE);
					double k_s = (existIb) ? 0.0 : (aux * (sT->y * dz + (sT + _NM)->y * dzz) * az / (aux2 + sT->z * (curD - _N)->y));
					az = curA->x; existIb2 = (az < EQ_DOUBLE);
					double k_n = (existIb2) ? 0.0 : (aux * (curT->y * dz + tT->y * dzz) * az / (aux2 + (curT + _N)->z * (curD + _N)->y));
// tex: s^{ib,z,y} = (dS^{ib,t}_{xz}(i,j,k) + dS^{ib,b}_{xz}(i,j,k+1))/dy^{ib,t} при alpha^x_{i,j-1,k}=0 ИЛИ alpha^x_{i,j,k}=0
					if (existIb || existIb2) ib += _nu * (_delta_ib_st[num].y + _delta_ib_nb[num + _NM].y) * _inv_delta_ib_tn[num].y;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - k_s - k_w - k_e - k_n - (k_t + k_wallB) * tz - ib;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_n;
					_Cell[numA2 + 5] = k_t * tT->z;
				}
				break;
				// ---------- Рёбра (4 элемента) ----------
				case TN_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const  coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (aux2 + wT->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_w - k_e - k_wallN;
					_Cell[numA2 + 4] = k_e;
				}
				break;
				case TS_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const  coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (aux2 + wT->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_w - k_e - k_n - k_wallS;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_n;
				}
				break;
				case TE_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* sT = curT - _N; const coordinate* wT = curT - 1;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (aux2 + sT->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + tT->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_w - k_e - k_n;
					_Cell[numA2 + 4] = k_n;
				}
				break;
				case TW_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (aux2 + sT->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + tT->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_e - k_n;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_n;
				}
				break;
				case BN_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* wT = curT - 1; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (aux2 + wT->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - k_s - k_w - k_e - (k_t + k_wallB) * tz - k_wallN;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_t * tT->z;
				}
				break;
				case BS_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx;
					double k_w = aux * ((curT - 1)->x * dz + ((curT - 1) + _NM)->x * dzz) * (curA - 1)->y / (aux2 + (curT - 1)->z * (curD - 1)->x);
					double k_e = aux * (curT->x * dz + tT->x * dzz) * curA->y / (aux2 + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_w;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - k_w - k_e - k_n - (k_t + k_wallB) * tz - k_wallS;
					_Cell[numA2 + 2] = k_e;
					_Cell[numA2 + 3] = k_n;
					_Cell[numA2 + 4] = k_t * tT->z;
				}
				break;
				case BE_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* wT = curT - 1;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * ((curT - _N)->y * dz + ((curT - _N) + _NM)->y * dzz) * (curA - _N)->x / (aux2 + (curT - _N)->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - k_s - k_w - k_e - k_n - (k_t + k_wallB) * tz;
					_Cell[numA2 + 3] = k_n;
					_Cell[numA2 + 4] = k_t * (curT + _NM)->z;
				}
				break;
				case BW_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double aux = dx * _nu, aux2 = tz * dy;
					double k_s = aux * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (aux2 + sT->z * (curD - _N)->y);
					double k_n = aux * (curT->y * dz + tT->y * dzz) * curA->x / (aux2 + (curT + _N)->z * (curD + _N)->y);
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - k_s - k_e - k_n - (k_t + k_wallB) * tz;
					_Cell[numA2 + 2] = k_e;
					_Cell[numA2 + 3] = k_n;
					_Cell[numA2 + 4] = k_t * tT->z;
				}
				break;
				case NE_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* sT = curT - _N; const  coordinate* wT = curT - 1;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_s - k_w - k_e - k_wallN;
					_Cell[numA2 + 4] = k_t * (curT + _NM)->z;
				}
				break;
				case NW_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_s - k_e - k_wallN;
					_Cell[numA2 + 3] = k_e;
					_Cell[numA2 + 4] = k_t * (curT + _NM)->z;
				}
				break;
				case SE_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const  coordinate* wT = curT - 1;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_w - k_e - k_n - k_wallS;
					_Cell[numA2 + 3] = k_n;
					_Cell[numA2 + 4] = k_t * (curT + _NM)->z;
				}
				break;
				case SW_Edge:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - (k_b + k_t) * tz - k_e - k_n - k_wallS;
					_Cell[numA2 + 2] = k_e;
					_Cell[numA2 + 3] = k_n;
					_Cell[numA2 + 4] = k_t * (curT + _NM)->z;
				}
				break;
				// ---------- Углы (3 элемента) ----------
				case TNE_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* sT = curT - _N; const  coordinate* wT = curT - 1;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_w - k_e - k_wallN;
				}
				break;
				case TNW_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const 	coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_s - k_e - k_wallN;
					_Cell[numA2 + 3] = k_e;
				}
				break;
				case TSE_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					const coordinate * wT = curT - 1;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_w - k_e - k_n - k_wallS;
					_Cell[numA2 + 3] = k_n;
				}
				break;
				case TSW_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallT = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num + _NM];
					double k_b = _nu * tz * dx * dx * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - (k_b + k_wallT) * tz - k_e - k_n - k_wallS;
					_Cell[numA2 + 2] = k_e;
					_Cell[numA2 + 3] = k_n;
				}
				break;
				case BNE_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					const coordinate * wT = curT - 1;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					const coordinate * sT = curT - _N;
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = k_w;
					_Cell[numA2 + 2] = -(_coefIntegr * _Mz[num]) - k_s - k_w - k_e - (k_t + k_wallB) * tz - k_wallN;
					_Cell[numA2 + 3] = k_t * (curT + _NM)->z;
				}
				break;
				case BNW_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					const coordinate * sT = curT - _N;
					double k_s = dx * _nu * (sT->y * dz + (sT + _NM)->y * dzz) * (curA - _N)->x / (tz * dy + sT->z * (curD - _N)->y);
					double k_wallN = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_s;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - k_s - k_e - (k_t + k_wallB) * tz - k_wallN;
					_Cell[numA2 + 2] = k_e;
					_Cell[numA2 + 3] = k_t * (curT + _NM)->z;
				}
				break;
				case BSE_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					const 	coordinate * wT = curT - 1;
					double k_w = dy * _nu * (wT->x * dz + (wT + _NM)->x * dzz) * (curA - 1)->y / (tz * dx + wT->z * (curD - 1)->x);
    double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = k_w;
					_Cell[numA2 + 1] = -(_coefIntegr * _Mz[num]) - k_w - k_e - k_n - (k_t + k_wallB) * tz - k_wallS;
					_Cell[numA2 + 2] = k_n;
					_Cell[numA2 + 3] = k_t * (curT + _NM)->z;
				}
				break;
				case BSW_Corner:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double k_e = dy * _nu * (curT->x * dz + (curT + _NM)->x * dzz) * curA->y / (tz * dx + (curT + 1)->z * (curD + 1)->x);
					double k_n = dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy + (curT + _N)->z * (curD + _N)->y);
					double k_wallS = (tz > EQ_DOUBLE) ? (dx * _nu * (curT->y * dz + (curT + _NM)->y * dzz) * curA->x / (tz * dy)) : 0.0;
					dx *= dy;
					double k_wallB = _isFlow ? 0.0 : _nu * tz * dx * dx * _inv_V[num];
					double k_t = _nu * tz * dx * dx * _inv_V[num + _NM];
					_Cell[numA2] = -(_coefIntegr * _Mz[num]) - k_e - k_n - (k_t + k_wallB) * tz - k_wallS;
					_Cell[numA2 + 1] = k_e;
					_Cell[numA2 + 2] = k_n;
					_Cell[numA2 + 3] = k_t * (curT + _NM)->z;
				}
				break;
				default:
				{
					const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curT = _theta + num;
					const coordinate* tT = curT + _NM; const coordinate* wT = curT - 1; const  coordinate* sT = curT - _N;
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, tz = curT->z;
					double aux = dy * _nu, aux2 = tz * dx, az = (curA - 1)->y;
					bool existIb = (az < EQ_DOUBLE);
					double k_w = (existIb) ? 0.0 : (aux * (wT->x * dz + (wT + _NM)->x * dzz) * az / (aux2 + wT->z * (curD - 1)->x));
					az = curA->y; bool existIb2 = (az < EQ_DOUBLE);
					double k_e = (existIb2) ? 0.0 : (aux * (curT->x * dz + tT->x * dzz) * az / (aux2 + (curT + 1)->z * (curD + 1)->x));
// tex: s^{ib,z,x} = (dS^{ib,t}_{yz}(i,j,k) + dS^{ib,b}_{yz}(i,j,k+1))/dx^{ib,t} при alpha^y_{i-1,j,k}=0 ИЛИ alpha^y_{i,j,k}=0
					double ib = (existIb || existIb2) ? (_nu * (_delta_ib_st[num].z + _delta_ib_nb[num + _NM].z) * _inv_delta_ib_tn[num].x) : 0.0;
					aux = dx * _nu, aux2 = tz * dy;
					az = (curA - _N)->x; existIb = (az < EQ_DOUBLE);
					double k_s = (existIb) ? 0.0 : (aux * (sT->y * dz + (sT + _NM)->y * dzz) * az / (aux2 + sT->z * (curD - _N)->y));
					az = curA->x; existIb2 = (az < EQ_DOUBLE);
					double k_n = (existIb2) ? 0.0 : (aux * (curT->y * dz + tT->y * dzz) * az / (aux2 + (curT + _N)->z * (curD + _N)->y));
// tex: s^{ib,z,y} = (dS^{ib,t}_{xz}(i,j,k) + dS^{ib,b}_{xz}(i,j,k+1))/dy^{ib,t} при alpha^x_{i,j-1,k}=0 ИЛИ alpha^x_{i,j,k}=0
					if (existIb || existIb2) ib += _nu * (_delta_ib_st[num].y + _delta_ib_nb[num + _NM].y) * _inv_delta_ib_tn[num].y;
					dx *= dy;
					dzz = _nu * tz * dx * dx;
					double k_t = dzz * _inv_V[num + _NM];
					double k_b = dzz * _inv_V[num];
					_Cell[numA2] = k_b * (curT - _NM)->z;
					_Cell[numA2 + 1] = k_s;
					_Cell[numA2 + 2] = k_w;
					_Cell[numA2 + 3] = -(_coefIntegr * _Mz[num]) - k_s - k_w - (k_b + k_t) * tz - k_e - k_n - ib;
					_Cell[numA2 + 4] = k_e;
					_Cell[numA2 + 5] = k_n;
					_Cell[numA2 + 6] = k_t * tT->z;
				}
				}//switch
			}//for i
		};
		LS_STAG_KERNEL_2D(ibFluxHelmZ, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _SideType, const coordinate* _theta, const coordinate* _alpha,
			LS_STAG_CR_PC _delta, LS_STAG_CR_PC _inv_delta_ib_tn, const coordinate* _delta_ib_nb, const coordinate* _delta_ib_st,
			LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibZY, LS_STAG_CR_PC _VibYX, const double* _inv_V, const int _N, const int _NM,
			const double _nu)
		{
			for (int num = i * _N + j * _NM, bNum = (j) ? _NM : 0, m_end = num + _N; num < m_end; num++)//for(i = 0; i < m_N; i++)
			{
				const CellTypes curType = _type[num];
				bool isCutCell = false;
				if((curType > SolidCell) && (curType < FluidCell))
				isCutCell = true;
				else
				{
					switch(curType)
					{
						case B_Border:
						case T_Border:
					if (_SideType[num].inRange(SolidSide, FluidSide)) isCutCell = true;
					}
				}
				if (isCutCell)
				{
					const coordinate* curT = _theta + num; double tz = curT->z;
					if (tz < EQ_DOUBLE)
						_flux[num] = 0.0;
					else
					{
						const coordinate* curD = _delta + num; const coordinate* curA = _alpha + num; const coordinate* curVib = _Vib + num;
						double dxy = curD->x * curD->y;
						double ib = ((curA->y < EQ_DOUBLE) || ((curA - 1)->y < EQ_DOUBLE)) ?
							((_delta_ib_st[num].z + _delta_ib_nb[num + _NM].z) * _inv_delta_ib_tn[num].x * _VibYX[num].z) : 0.0;
						if ((curA->x < EQ_DOUBLE) || ((curA - _N)->x < EQ_DOUBLE))
							ib += (_delta_ib_st[num].y + _delta_ib_nb[num + _NM].y) * _inv_delta_ib_tn[num].y * _VibZY[num].z;

						_flux[num] = _nu * (ib + tz * dxy * dxy * (((tz - (curT + _NM)->z) * (curVib + _NM)->z * _inv_V[num + _NM]) -
							(((curT - bNum)->z - tz) * curVib->z * _inv_V[num])));
					}
				}// Grid.TypeOfCellX[numX] > 9
				else _flux[num] = 0.0;
			}//for i
		};
		LS_STAG_KERNEL_2D(setPlusOutFlowZ, LS_STAG_R_PD _PlusOutFlow, LS_STAG_CR_PC _delta, const double _nu, const int _N,
			const int _M, const int _NM)
		{
			const int numE = i + j * _M;
			const coordinate* curD = _delta + (j * _NM + (i + 1) * _N - 1);
			_PlusOutFlow[numE] = _nu * (curD->z + (curD + _NM)->z) * curD->y / curD->x;
		};
		using HelmholtzDiscreteOperatorBase::initMeshParams;
		void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			HelmholtzDiscreteOperatorBase::initMeshParams(_cur_mesh, _meshes);
			m_Mz = m_meshes->zMesh.VOF;
			build();
		};
	  private:
		void calcPlusOutFlow()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorZ, setPlusOutFlowZ, m_M, m_K, PlusOutFlow, m_delta, m_nu, m_N, m_M, m_N * m_M)
		}
		void initInFlowFlux()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorZ, inFlowFluxhz, m_M, m_K, flux, m_delta, m_task->W_W,
				m_task->W_E, false, m_N, m_N * m_M, m_nu)
			LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorZ, addFluxNShz, m_K, flux, m_delta, m_task->W_S, m_task->W_N,
				m_N, m_N * m_M, m_nu)
			if(!existOutFlow) LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorZ, addFluxTBhz, m_M, flux, m_type, m_delta, m_task->W_B, m_task->W_T,
				m_N, NMK - m_N * m_M, m_nu)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора Гельмгольца на z-сетке на погруженных границах 
		  */
		void ibFlux()
		{// испр. 31.05.14
			LS_STAG_PAR_FOR_2D_RANGE(HelmholtzDiscreteOperatorZ, ibFluxHelmZ, 1, M2, 0, m_K, flux, m_type, m_baseSideTypes,
				m_theta, m_alpha, m_delta, m_inv_delta_ib_tn, m_delta_ib_nb, m_delta_ib_st, m_Vib, m_VibZY, m_VibYX, 
                m_inv_V, m_N, m_N * m_M, m_nu)
		};
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора Гельмгольца на z-сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(HelmholtzDiscreteOperatorZ, calcHelmZ, m_M, m_K, Number, Cell, m_theta, m_alpha, m_delta, 
				m_inv_delta_ib_tn, m_delta_ib_nb, m_delta_ib_st, m_inv_V, m_Mz, m_type, m_N, m_M, m_N * m_M, m_nu, coefIntegr,
				(m_task->id == kind::TASK_FLOW))
		   addOutFlow();
		};
	  private:
		double *m_Mz;///< массив значений объемов ячеек с z-сетки
	};
}

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_Z_H_
