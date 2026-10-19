/*
* \copyright (C) 2013 Puzikova V.V.
* \file      convective_discrete_operator_y.h
* \author    Puzikova Valeria (vvp) 
* \date      24.10.2013 22:34:44
* \brief     convective_discrete_operator_y.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Y_H_
#define _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Y_H_

// ================================================================== INCLUDES
#include "skew_sym_discrete_operator.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ConvectiveDiscreteOperatorY

Данный файл содержит в себе реализацию класса ConvectiveDiscreteOperatorY
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(baseFluxY, LS_STAG_R_PD, _flux, LS_STAG_R_PC, _delta, const double*, _Ux, 
		LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _U_W, const double, _V_W, const double, _UV_E, const bool, isArtery, 
		const bool, isCavity, const int, _N, const int, _M3, const int, _NM, const int, _vs, const int, _vf, 
		const int, _zs, const int, _zf);
	LS_STAG_DECLARATION_KERNEL_2D(calcConvY, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell, const coordinate*, _theta, 
		coordinate*, _delta, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PCT, _type, 
		const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N2);
	LS_STAG_DECLARATION_KERNEL_2D(ibFluxConvY, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType,
		LS_STAG_CR_PCT, _baseType, coordinate*, _theta, coordinate*, _delta, LS_STAG_CR_PC, _Vib,
		LS_STAG_CR_PC, _VibENT, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PD, _Uib,
		LS_STAG_CR_PD, _Ug, const int, _N, const int, _M3, const int, _NM, const int, baseNM);
	LS_STAG_DECLARATION_KERNEL(addFluxNSy, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		LS_STAG_CR_PD, m_Uy, const double, V_S, const double, V_N, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxTBy, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, VW_B, const double, VW_T, const int, _N, const int, _NMK2);
	LS_STAG_DECLARATION_KERNEL_2D(addFluxTBOutY, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_theta, LS_STAG_CR_PC, m_delta, LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _Uz,
		const int, _N, const int, _NMy, const int, _K, const int, _Nbase, const int, _NMbase, const int, _Kw);
/*!
\brief Класс, реализующий дискретный оператор конвекции на у-сетке

Данный класс реализует работу с дискретным аналогом оператора конвекции на у-сетке
*/
	struct ConvectiveDiscreteOperatorY : SkewSymDiscreteOperator
	{
	LS_STAG_KERNEL(addFluxNSy, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, LS_STAG_CR_PD m_Uy, const double V_S, 
		const double V_N, const int _N, const int _NM)
		{
			int numS = i * _NM;
			coordinate* curD = m_delta + (numS + i * _N);
			double dz = 0.25 * curD->z;
			for(int ii = 0, numN = numS + _NM - _N; ii < _N; ii++, numS++, numN++, curD++)
			{
				double coef = dz * curD->x;
				flux[numS] -= coef * V_S * (m_Uy[numS] + V_S);
				flux[numN] += coef * V_N * (m_Uy[numN] + V_N);
			}
		}
	LS_STAG_KERNEL(addFluxTBy, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta, const double VW_B, const double VW_T,
		const int _N, const int _NMK2)
		{
			int numB = i * _N;
			coordinate* curD = m_delta + numB;
			double dyy = 0.25 * (curD->y + (curD+_N)->y);
			for(int ii = 0, numT = _NMK2 + numB; ii < _N; ii++, numB++, numT++, curD++)
			{
				double coef = curD->x * dyy;
				if (_type[numB] != SolidCell) flux[numB] -= coef * VW_B;
				if (_type[numT] != SolidCell) flux[numT] += coef * VW_T;
			}
		}	
		LS_STAG_KERNEL_2D(addFluxTBOutY, LS_STAG_R_PD flux, LS_STAG_CR_PC m_theta, LS_STAG_CR_PC m_delta, LS_STAG_CR_PD _Uy, LS_STAG_CR_PD _Uz,
			const int _N, const int _NMy, const int _K, const int _Nbase, const int _NMbase, const int _Kw)
		{// конвективный выход на B/T (как элементы матрицы): S_y = -/+ (1/4)*dx*(dy*w_j + dyN*w_{j+1}) с theta.z
			const int numB = i + j * _N;                      // v-ячейка (i,j,0)
			const int numT = numB + _Kw * _NMy;          // v-ячейка (i,j,K-1)
			const coordinate* curD = m_delta + numB;            // v-сетка: N = Nbase --- индексы совпадают
			const coordinate* curT = m_theta + numB;
			// theta.z базовых ячеек B/T: вынос идёт только через жидкую (неусечённую) часть грани
			const double wfB = (_Uz[numB] * curD->y * curT->z + _Uz[numB + _Nbase] * (curD + _Nbase)->y * (curT + _Nbase)->z);
			flux[numB] -= 0.25 * curD->x * _Uy[numB] * wfB;
			const int numToT = _Kw * _NMbase;
			const int wT = numB + numToT - _NMbase;
			const coordinate* curDT = curD + numToT;
			const coordinate* curTT = curT + numToT;
			const double wfT = (_Uz[wT] * curDT->y * curTT->z + _Uz[wT + _Nbase] * (curDT + _Nbase)->y * (curTT + _Nbase)->z);
			flux[numT] += 0.25 * curDT->x * _Uy[numT] * wfT;
		}
		LS_STAG_KERNEL_2D(calcConvY, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate* _theta, coordinate* _delta,
			LS_STAG_CR_PD _Ux, const double* _Uy, const double* _Uz, LS_STAG_CR_PCT _type, const int _N, const int _M, 
			const int _NM, const int baseNM, const int _N2)
		{
			int num = _N * i + j * _NM;// = numY
			int numZ = num + j * _N;// numBase = numZ
			coordinate* curD = _delta + numZ;
			for (int numX = numZ - i - j * _M, m_end = num + _N, numN = numX + _N2, numB = numZ - baseNM; num < m_end; 
				num++, numZ++, numX++, numB++, numN++, curD++)
			{
				const int numA2 = _Number[num];
				switch (_type[num])
				{
				case FluidCell:
				{
					double dy = curD->y, dz  = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz, curV = _Uy[num];
					const coordinate* curT = _theta + numZ; const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
            _Cell[numA2] = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
		    _Cell[numA2+1]   = -0.25 * dxz * (curV + (curT-_N)->y * _Uy[num-_N]); // S
            _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]);                 // W
            _Cell[numA2+3] =  0.25 * dz * (dy * _Ux[numX] + dyy * nT->x * _Ux[numN]);                 // E
            _Cell[numA2+4] =  0.25 * dxz * (curV + nT->y * _Uy[num+_N]);  // N
            _Cell[numA2+5] =  0.25 * dx * (dy * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);  // T
				}
				break;
				case SolidCell:
					_Cell[numA2] = 0.0;
					break;
				case N_Border:
				{
					double dy = curD->y, dz  = curD->z, dx = curD->x, dyy = (curD+_N)->y;
					const coordinate* curT = _theta + numZ; const coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
            _Cell[numA2] = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
		    _Cell[numA2+1]   = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]); // S
            _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]);                 // W
            _Cell[numA2+3] =  0.25 * dz * (dy * curT->x * _Ux[numX] + dyy * nT->x * _Ux[numN]);                 // E
            _Cell[numA2+4] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);  // T
				}
				break;
case S_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
    const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+4] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case E_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
	double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+3] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+4] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case W_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
	double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+4] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case T_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const   coordinate* curT = _theta + numZ; const   coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
	double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+3] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+4] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
}
break;
case B_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+4] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
// ---------- Рёбра (4 элемента) ----------
case TN_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const   coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);      // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+3] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
}
break;
case TS_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const   coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);                // N
}
break;
case TE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
	double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+3] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
}
break;
case TW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
	double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
}
break;
case BN_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N;
    _Cell[numA2]   = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BS_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N;
    _Cell[numA2]   = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+1] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+2] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);            // N
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; double ty = curT->y * _Uy[num];
    _Cell[numA2]   = -0.25 * dxz * (ty + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+2] =  0.25 * dxz * (ty + nT->y * _Uy[num+_N]);                // N
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case NE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;

case NW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case SE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ;  const coordinate* nT = curT + _N;  const coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);            // N
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case SW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+2] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);            // N
    _Cell[numA2+3] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
// ---------- Углы (3 элемента) ----------
case TNE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const   coordinate* curT = _theta + numZ; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);      // S
    _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (curT + _N2)->x * _Ux[numN-1]); // W
}
break;
case TNW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);      // S
    _Cell[numA2+2] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * (curT + _N)->x * _Ux[numN]);     // E
}
break;
case TSE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ; const coordinate* nT = curT + _N; const coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);                // N
}
break;
case TSW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
    _Cell[numA2+1] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+2] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);                // N
}
break;
case BNE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ; const   coordinate* nT = curT + _N;
    _Cell[numA2]   = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+2] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BNW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N;
    _Cell[numA2]   = -0.25 * dx * dz * (curT->y * _Uy[num] + (curT-_N)->y * _Uy[num-_N]);          // S
    _Cell[numA2+1] =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+2] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BSE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N;
    _Cell[numA2]   = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]); // W
    _Cell[numA2+1] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);            // N
    _Cell[numA2+2] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
case BSW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dyy = (curD+_N)->y;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N;
    _Cell[numA2]   =  0.25 * dz * (dy * curT->x * _Ux[numX]   + dyy * nT->x * _Ux[numN]);     // E
    _Cell[numA2+1] =  0.25 * dx * dz * (curT->y * _Uy[num] + nT->y * _Uy[num+_N]);            // N
    _Cell[numA2+2] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);    // T
}
break;
				default:
				{
					double dy = curD->y, dz  = curD->z, dx = curD->x, dyy = (curD+_N)->y, dxz = dx * dz;
					const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + _N; const coordinate* bT = curT - baseNM;
					double curV = curT->y * _Uy[num];
            _Cell[numA2] = -0.25 * dx * (dy * bT->z * _Uz[numB] + dyy * (bT+_N)->z * _Uz[numB+_N]); // B
		    _Cell[numA2+1]   = -0.25 * dxz * (curV + (curT-_N)->y * _Uy[num-_N]); // S
            _Cell[numA2+2] = -0.25 * dz * (dy * (curT-1)->x * _Ux[numX-1] + dyy * (nT-1)->x * _Ux[numN-1]);                 // W
            _Cell[numA2+3] =  0.25 * dz * (dy * curT->x * _Ux[numX] + dyy * nT->x * _Ux[numN]);                 // E
            _Cell[numA2+4] =  0.25 * dxz * (curV + nT->y * _Uy[num+_N]);  // N
            _Cell[numA2+5] =  0.25 * dx * (dy * curT->z * _Uz[numZ] + dyy * nT->z * _Uz[numZ+_N]);  // T
				}
				}//switch
			}// for i 
		}
		LS_STAG_KERNEL_2D(baseFluxY, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, const double* _Ux, LS_STAG_CR_PD _Uy,
			LS_STAG_CR_PD _U_W, const double _V_W, const double _UV_E, const bool isArtery, const bool isCavity, 
			const int _N, const int _M3, const int _NM, const int _vs, const int _vf, const int _zs, const int _zf)
		{
			int numF = _N * i+ j * _NM, N2 = _N - 1;
			coordinate* curD = _delta + (numF + j * _N);
			double dz = curD->z, dy = curD->y, dyy = (curD + _N)->y;
			if((isArtery && ((i <= _vs) || (i >= _vf) || ((j < _zs) || (j > _zf)))))
			    _flux[numF] = 0.0;
			else
			{
				const int numW = i + j * _M3;
				_flux[numF] = -0.25 * (_U_W[numW] * dy + _U_W[numW+1] * dyy) * _V_W * dz;
			}	 
			for (int num = numF + 1, m_end = numF + N2; num < m_end; num++)
				_flux[num] = 0.0;
			numF += N2;
			if (isCavity)
				_flux[numF] = 0.25 * _UV_E * (dy + dyy) * dz;
			else
			{
				const int num2 = numF + j * (_N - _M3) - 1 - i;
				_flux[numF] = 0.25 * _Uy[numF] * (_Ux[num2] * dy + _Ux[num2 + N2] * dyy) * dz;
			}
		}
		LS_STAG_KERNEL_2D(ibFluxConvY, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _baseSideType, LS_STAG_CR_PCT _baseType,
			coordinate* _theta, coordinate* _delta, LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibENT,
			LS_STAG_CR_PD _Ux, const double* _Uy, const double* _Uz, LS_STAG_CR_PD _Uib, LS_STAG_CR_PD _Ug, const int _N,
			 const int _M3, const int _NM, const int baseNM)
		{
			for (int numY = _N * i + j * _NM, m_end = numY + _N; numY < m_end; numY++)
			{
				const CellTypes curType = _type[numY];
				if ((curType > SolidCell) && (curType < FluidCell))
				{
					double cxe = 0.0, cxw = 0.0;
					int num = numY + j * _N;
					switch (_baseType[num])
					{
					case BSE_Heptahedron3:
					case BSW_Heptahedron3:
					case TSE_Heptahedron3:
					case TSW_Heptahedron3:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].y;
						break;
					case BW_Pentahedron:
					case TW_Pentahedron:
					case BE_Pentahedron:
					case TE_Pentahedron:
					case E_Hexahedron4:
					case W_Hexahedron4:
					case T_Hexahedron4:
					case B_Hexahedron4:
					case BSBW_Hexahedron5:
					case BWSW_Hexahedron5:
					case BWNW_Hexahedron5:
					case TWNW_Hexahedron5:
					case TWSW_Hexahedron5:
					case BSBE_Hexahedron5:
					case BNBW_Hexahedron5:
					case BNBE_Hexahedron5:
					case TESE_Hexahedron5:
					case TENE_Hexahedron5:
					case BESE_Hexahedron5:
					case BENE_Hexahedron5:
					case TNTW_Hexahedron5:
					case TNTE_Hexahedron5:
					case TSTW_Hexahedron5:
					case TSTE_Hexahedron5:
					case BW_Heptahedron4:
					case TW_Heptahedron4:
					case BE_Heptahedron4:
					case TE_Heptahedron4:
					case BSW_Heptahedron6:
					case BSE_Heptahedron6:
					case TSW_Heptahedron6:
					case TSE_Heptahedron6:
					case TESE_Heptahedron5:
					case TSTW_Heptahedron5:
					case TSTE_Heptahedron5:
					case BSBW_Heptahedron5:
					case BSBE_Heptahedron5:
					case TWSW_Heptahedron5:
					case BWSW_Heptahedron5:
					case BESE_Heptahedron5:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num].y;
						break;
					case S_Hexahedron4:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->y * curD->x * curD->z * _Uy[numY]) * _Vib[num].y;
					}
						break;
					case NE_Heptahedron4:
					{// 2D NE_Pentagon: 0.25*(Uib-Ug + theta^u*dSy*dSz*Ux)*vib  (тангенциальная Ux, theta^u=theta.x)
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * curD->y * curD->z * _Ux[num-i-j*_M3]) * _VibENT[num].y;
					}
						break;
					case NW_Heptahedron4:
					{// 2D NW_Pentagon: 0.25*(Uib-Ug - theta^u_w*dSy*dSz*Ux_w)*vib
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * curD->y * curD->z * _Ux[num-i-j*_M3 - 1]) * _VibENT[num].y;
					}
						break;
					case SW_Pentahedron:
					case SE_Pentahedron:
					case BS_Pentahedron:
					case TS_Pentahedron:
					case TNTE_Heptahedron5:
					case TNTW_Heptahedron5:
					case BENE_Heptahedron5:
					case TWNW_Heptahedron5:
					case TENE_Heptahedron5:
					case BWNW_Heptahedron5:
					case BNBW_Heptahedron5:
					case BNBE_Heptahedron5:
					case TN_Heptahedron4:
					case BN_Heptahedron4:
					case BSSW_Hexahedron5:
					case BSSE_Hexahedron5:
					case TSSW_Hexahedron5:
					case TSSE_Hexahedron5:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->y * curD->x * curD->z * _Uy[numY]) * _VibENT[num].y;
					}
						break;
					case TNNE_Heptahedron5:
					case TNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dy * curD->z * _Ux[num-i-j*_M3] + 
							curT->z * dy * curD->x * _Uz[num]) * _VibENT[num].y;
					}
						break;
					case BNNE_Heptahedron5:
					case BNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dy * curD->z * _Ux[num-i-j*_M3] -  
							(curT-_NM)->z * dy * curD->x * _Uz[num-baseNM]) * _VibENT[num].y;
					}
						break;
					case TNNW_Heptahedron5:
					case TNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-1)->x * dy * curD->z * _Ux[num-1-i-j*_M3] + 
							curT->z * dy * curD->x * _Uz[num]) * _VibENT[num].y;
					}
						break;
					case BNNW_Heptahedron5:
					case BNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-1)->x * dy * curD->z * _Ux[num-1-i-j*_M3] - 
							(curT-_NM)->z * dy * curD->x * _Uz[num-baseNM]) * _VibENT[num].y;
					}
						break;
					case BSW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y +
							curD->y * (curVib->x * curT->x * dz + curVib->z * curT->z * dx) * _VibENT[num].y);
					}
					break;
					case BSE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y +
							curD->y * (-curVib->x * (curT-1)->x * dz + curVib->z * curT->z * dx) * _VibENT[num].y);
					}
					break;
					case TSW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y +
							curD->y * (curVib->x * curT->x * dz - curVib->z * (curT-baseNM)->z * dx) * _VibENT[num].y);
					}
					break;
					case TSE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y -
							curD->y * (curVib->x * (curT-1)->x * dz + curVib->z * (curT-baseNM)->z * dx) * _VibENT[num].y);
					}
					break;
					}//switch
					num += _N;
					switch (_baseType[num])
					{
					case BNE_Heptahedron3:
					case BNW_Heptahedron3:
					case TNE_Heptahedron3:
					case TNW_Heptahedron3:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].y;
						break;
					case BW_Pentahedron:
					case TW_Pentahedron:
					case BE_Pentahedron:
					case TE_Pentahedron:
					case E_Hexahedron4:
					case W_Hexahedron4:
					case T_Hexahedron4:
					case B_Hexahedron4:
					case BSBW_Hexahedron5:
					case BWSW_Hexahedron5:
					case BWNW_Hexahedron5:
					case TWNW_Hexahedron5:
					case TWSW_Hexahedron5:
					case BSBE_Hexahedron5:
					case BNBW_Hexahedron5:
					case BNBE_Hexahedron5:
					case TESE_Hexahedron5:
					case TENE_Hexahedron5:
					case BESE_Hexahedron5:
					case BENE_Hexahedron5:
					case TNTW_Hexahedron5:
					case TNTE_Hexahedron5:
					case TSTW_Hexahedron5:
					case TSTE_Hexahedron5:
					case BW_Heptahedron4:
					case TW_Heptahedron4:
					case BE_Heptahedron4:
					case TE_Heptahedron4:
					case BNW_Heptahedron6:
					case BNE_Heptahedron6:
					case TNW_Heptahedron6:
					case TNE_Heptahedron6:
					case TENE_Heptahedron5:
					case TNTW_Heptahedron5:
					case TNTE_Heptahedron5:
					case BNBW_Heptahedron5:
					case BNBE_Heptahedron5:
					case TWNW_Heptahedron5:
					case BWNW_Heptahedron5:
					case BENE_Heptahedron5:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num-_N].y;
						break;
					case N_Hexahedron4:
					{
						coordinate* curD = _delta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num - _N)->x * curD->x * curD->z * _Uy[numY]) * _Vib[num].y;
					}
						break;
					case NW_Pentahedron:
					case BN_Pentahedron:
					{// 2D NW_Triangle: 0.25*(((vy+Uy)*dSx_s*thx_s - Ug)*vy - (vx*thy_n*dSy_n + Ug)*vib_n)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y -
							curD->y * (curVib->x * curT->x * dz + curVib->z * curT->z * dx) * _VibENT[num].y);
					}
						break;
					case NE_Pentahedron:
					{// 2D NE_Triangle: 0.25*(((vy+Uy)*dSx_s*thx_s - Ug)*vy + (vx*thy_wn*dSy_n - Ug)*vib_n)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y +
							curD->y * (curVib->x * (curT-1)->x * dz - curVib->z * curT->z * dx) * _VibENT[num].y);
					}
						break;
					case SE_Heptahedron4:
					{// 2D SE_Pentagon: 0.25*(Uib-Ug + theta^u_n*dSy*dSz*Ux_n)*vib  (num уже = северная ячейка; Ux_n по индексу севера)
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * curD->y * curD->z * _Ux[num - (i+1) - j*_M3]) * _VibENT[num].y;
					}
						break;
					case SW_Heptahedron4:
					{// 2D SW_Pentagon: 0.25*(Uib-Ug - theta^u_wn*dSy*dSz*Ux_wn)*vib
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * curD->y * curD->z * _Ux[num - (i+1) - j*_M3 - 1]) * _VibENT[num].y;
					}
						break;
					case TSTE_Heptahedron5:
					case TSTW_Heptahedron5:
					case BESE_Heptahedron5:
					case TWSW_Heptahedron5:
					case TESE_Heptahedron5:
					case BWSW_Heptahedron5:
					case BSBW_Heptahedron5:
					case BSBE_Heptahedron5:
					case TS_Heptahedron4:
					case BS_Heptahedron4:
					case BNNW_Hexahedron5:
					case BNNE_Hexahedron5:
					case TNNW_Hexahedron5:
					case TNNE_Hexahedron5:
					{
						coordinate* curD = _delta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num-_N)->y * curD->x * curD->z * _Uy[numY]) * _VibENT[num-_N].y;
					}
						break;
					case TSSE_Heptahedron5:
					case TSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dy * curD->z * _Ux[num-i-j*_M3] + 
							curT->z * dy * curD->x * _Uz[num]) * _VibENT[num-_N].y;
					}
						break;
					case BSSE_Heptahedron5:
					case BSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dy * curD->z * _Ux[num-i-j*_M3] -  
							(curT-_NM)->z * dy * curD->x * _Uz[num-baseNM]) * _VibENT[num-_N].y;
					}
						break;
					case TSSW_Heptahedron5:
					case TSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-1)->x * dy * curD->z * _Ux[num-1-i-j*_M3] + 
							curT->z * dy * curD->x * _Uz[num]) * _VibENT[num-_N].y;
					}
						break;
					case BSSW_Heptahedron5:
					case BSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-1)->x * dy * curD->z * _Ux[num-1-i-j*_M3] - 
							(curT-_NM)->z * dy * curD->x * _Uz[num-baseNM]) * _VibENT[num-_N].y;
					}
						break;
					case BNW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y -
							curD->y * (curVib->x * curT->x * dz + curVib->z * curT->z * dx) * _VibENT[num-_N].y);
					}
					break;
					case BNE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y +
							curD->y * (curVib->x * (curT-1)->x * dz - curVib->z * curT->z * dx) * _VibENT[num-_N].y);
					}
					break;
					case TNW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y +
							curD->y * (-curVib->x * curT->x * dz + curVib->z * (curT-baseNM)->z * dx) * _VibENT[num-_N].y);
					}
					break;
					case TNE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dz = curD->z;
						cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y +
							curD->y * (curVib->x * (curT-1)->x * dz + curVib->z * (curT-baseNM)->z * dx) * _VibENT[num-_N].y);
					}
					break;
					}//switch
					_flux[numY] = cxe + cxw;
				}// if type > 9
				else if ((curType == B_Border) || (curType == T_Border))
				{
					int num = numY + j * _N;
					const SideTypes szSide = _baseSideType[num].z;
					const SideTypes nzSide = _baseSideType[num + _N].z;
					const bool sCutCell = ((szSide > SolidSide) && (szSide < FluidSide));
					const bool nCutCell = ((nzSide > SolidSide) && (nzSide < FluidSide));
					if(sCutCell || nCutCell)
					{
						double cxe = 0.0, cxw = 0.0;
						if(sCutCell)
						{// S-половинка (грань num): только южные формы среза
							const coordinate* curD = _delta + num;
							switch (szSide)
							{
								case SW_Triangle:
								{// 2D SW_Triangle (общая формула, как SW_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dx = curD->x, dz = curD->z;
									cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y +
										curD->y * (curVib->x * curT->x * dz + curVib->z * curT->z * dx) * _VibENT[num].y);
								}
									break;
								case SE_Triangle:
								{// 2D SE_Triangle (общая формула, как SE_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									// на B-слое (j==0) нижнего соседа нет — продлеваем поле (theta_B := theta)
									const coordinate* curTB = (j > 0) ? (curT - baseNM) : curT;
									double dx = curD->x, dz = curD->z;
									cxe = -0.25 * (((curVib->y + _Uy[numY]) * dx * dz * curT->y + _Ug[num]) * curVib->y -
										curD->y * (curVib->x * (curT-1)->x * dz + curVib->z * curTB->z * dx) * _VibENT[num].y);
								}
									break;
								case NW_Pentagon:
								{// 2D NW_Pentagon (тангенциальная Ux_w), как NW_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * curD->y * curD->z * _Ux[num - i - j*_M3 - 1]) * _VibENT[num].y;
								}
									break;
								case NE_Pentagon:
								{// 2D NE_Pentagon (тангенциальная Ux), как NE_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * curD->y * curD->z * _Ux[num - i - j*_M3]) * _VibENT[num].y;
								}
									break;
								case S_Trapezoid:// как S_Hexahedron4 (простая — VibENT)
									cxe = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num].y;
									break;
								case W_Trapezoid:// как W_Hexahedron4 (с Uy — Vib)
								case E_Trapezoid:// как E_Hexahedron4
									cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->y * curD->x * curD->z * _Uy[numY]) * _Vib[num].y;
									break;
							}
						}
						num += _N;
						if(nCutCell)
						{// N-половинка (грань num+_N): только северные формы среза
							const coordinate* curD = _delta + num;
							switch (nzSide)
							{
								case NW_Triangle:
								{// 2D NW_Triangle (общая формула, как NW_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dx = curD->x, dz = curD->z;
									cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y -
										curD->y * (curVib->x * curT->x * dz + curVib->z * curT->z * dx) * _VibENT[num].y);
								}
									break;
								case NE_Triangle:
								{// 2D NE_Triangle (общая формула, как NE_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dx = curD->x, dz = curD->z;
									cxw = 0.25 * (((curVib->y + _Uy[numY]) * dx * dz * (curT-_N)->y - _Ug[num]) * curVib->y +
										curD->y * (curVib->x * (curT-1)->x * dz - curVib->z * curT->z * dx) * _VibENT[num].y);
								}
									break;
								case SW_Pentagon:
								{// 2D SW_Pentagon (тангенциальная Ux_wn), как SW_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * curD->y * curD->z * _Ux[num - (i+1) - j*_M3 - 1]) * _VibENT[num].y;
								}
									break;
								case SE_Pentagon:
								{// 2D SE_Pentagon (тангенциальная Ux_n), как SE_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * curD->y * curD->z * _Ux[num - (i+1) - j*_M3]) * _VibENT[num].y;
								}
									break;
								case N_Trapezoid:// как N_Hexahedron4 (простая — VibENT)
									cxw = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num - _N].y;
									break;
								case W_Trapezoid:// как E_Hexahedron4 (с Uy — Vib[num])
								case E_Trapezoid:
									cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num - _N)->y * curD->x * curD->z * _Uy[numY]) * _Vib[num].y;
									break;
							}
						}
						_flux[numY] = cxe + cxw;
					}
				}
			}//for i
		}
	  private:
		  /*!
		  Рассчитывает значения источникового члена оператора конвекции на у-сетке на границах расчетной области
		  */
		void baseFluxY()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorY, baseFluxY, m_M, m_K, flux, m_delta, m_Ux, m_Uy, 
				m_task->pU_W, m_task->V_W, m_task->U_E * m_task->V_E, false, 
				false, m_N, M3, m_N * m_M, vs, vf, zs, zf)
			LS_STAG_PAR_FOR(ConvectiveDiscreteOperatorY, addFluxNSy, m_K, flux, m_delta, m_Uy, m_task->V_S, 
				m_task->V_N, m_N, m_N * m_M)
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorY, addFluxTBOutY, m_N, m_M, flux, m_theta, m_delta, m_Uy, m_Uz,
				m_N, m_N * m_M, m_K, m_N, m_N * M3, m_K - 1)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора конвекции на у-сетке на погруженных границах 
		  */
		void ibFlux()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorY, ibFluxConvY, m_M, m_K, flux, m_type, m_baseSideTypes, m_meshes->baseMesh.typeOfCell,
				m_theta, m_delta, m_Vib, m_VibENT, m_Ux, m_Uy, m_Uz, m_Uib, m_Ug, m_N, M3, m_N * m_M, m_N * M3)
		};
	  public:
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора конвекции на у-сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorY, calcConvY, m_M, m_K, Number, Cell, m_theta, m_delta, m_Ux, 
				m_Uy, m_Uz, m_type, m_N, M3, m_N * m_M, m_N * M3, N2)
		};
	};
}

#endif //_LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Y_H_
