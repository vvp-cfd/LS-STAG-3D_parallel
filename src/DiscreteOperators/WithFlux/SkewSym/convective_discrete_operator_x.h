/*
* \copyright (C) 2013 Puzikova V.V.
* \file      convective_discrete_operator_x.h
* \author    Puzikova Valeria (vvp) 
* \date      24.10.2013 22:04:53
* \brief     convective_discrete_operator_x.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_X_H_
#define _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_X_H_

// ================================================================== INCLUDES
#include "skew_sym_discrete_operator.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ConvectiveDiscreteOperatorX

Данный файл содержит в себе реализацию класса ConvectiveDiscreteOperatorX
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(baseFluxX, LS_STAG_R_PD, _flux, LS_STAG_R_PC, _delta, LS_STAG_CR_PD, _Ux, 
		LS_STAG_CR_PD, _U_W, const double, _U_E, const bool, isArtery, const bool, isCavity, const int, _N, const int, _M, 
		const int, _NM, const int, _vs, const int, _vf, const int, _zs, const int, _zf);
	LS_STAG_DECLARATION_KERNEL_2D(calcConvX, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell,
		const coordinate*, _theta, LS_STAG_R_PC, _delta, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz,
		LS_STAG_CR_PCT, _type, const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N3);
	LS_STAG_DECLARATION_KERNEL_2D(ibFlux, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type,
		LS_STAG_CR_P(coord<SideTypes>), _baseSideType, LS_STAG_CR_PCT, _baseType, coordinate*, _theta, LS_STAG_R_PC, _delta, LS_STAG_CR_PC, _Vib,
		LS_STAG_CR_PC, _VibENT, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PD, _Uib,
		LS_STAG_CR_PD, _Ug, const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N3);
	LS_STAG_DECLARATION_KERNEL(addFluxNS, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta, 
		const double, UV_S, const double, UV_N, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxTB, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		const double, UW_B, const double, UW_T, const int, _N, const int, _NMK2);
	LS_STAG_DECLARATION_KERNEL_2D(addFluxTBOut, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_theta, LS_STAG_CR_PC, m_delta, LS_STAG_CR_PD, _Ux, LS_STAG_CR_PD, _Uz,
		const int, _N, const int, _NMx, const int, _K, const int, _Nbase, const int, _NMbase, const int, _Kw, const int, _M2);
/*!
\brief Класс, реализующий дискретный оператор конвекции на х-сетке

Данный класс реализует работу с дискретным аналогом оператора конвекции на х-сетке
*/
	struct ConvectiveDiscreteOperatorX : SkewSymDiscreteOperator
	{
		LS_STAG_KERNEL_2D(calcConvX, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate* _theta, LS_STAG_R_PC _delta,
			LS_STAG_CR_PD _Ux, const double* _Uy, const double* _Uz, LS_STAG_CR_PCT _type, const int _N, const int _M, 
			const int _NM, const int baseNM, const int _N3)
		{
			int num = _N * i + j * _NM;
			int numZ = num + i + j * _M;
			coordinate* curD = _delta + numZ;
			for (int numY = numZ - j * _N3, m_end = num + _N, numS = numY - _N3, numB = numZ - baseNM; num < m_end;
				num++, numZ++, numY++, numS++, numB++, curD++)
			{
				const int numA2 = _Number[num];
				switch (_type[num])
				{
				case FluidCell:
				{
					double dy = curD->y, dz  = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz, curU = _Ux[num];
					const coordinate* curT = _theta + numZ; const coordinate* nT = curT + 1;
					const coordinate* sT = curT - _N3; const coordinate* bT = curT - baseNM;
            _Cell[numA2] = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]); // B
		    _Cell[numA2+1]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS + 1]); // S
            _Cell[numA2+2] = -0.25 * dyz * (curU + (curT - 1)->x * _Ux[num-1]);                 // W
            _Cell[numA2+3] =  0.25 * dyz * (curU + nT->x * _Ux[num+1]);                 // E
            _Cell[numA2+4] =  0.25 * dz * (dx * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);  // N
            _Cell[numA2+5] =  0.25 * dy * (dx * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);  // T
				}
				break;
				case SolidCell:
					_Cell[numA2] = 0.0;
					break;
				case N_Border:
				{
			double dy = curD->y, dz  = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;	
			const coordinate* curT = _theta + numZ; const coordinate* nT = curT + 1; const coordinate* bT = curT - baseNM;
			double tx = curT->x * _Ux[num]; const coordinate* sT = curT - _N3;
        _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB + 1]);   // B
        _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS + 1]);   // S
        _Cell[numA2+2] = -0.25 * dyz * (tx + (curT - 1)->x * _Ux[num - 1]); // W
        _Cell[numA2+3] =  0.25 * dyz * (tx + nT->x * _Ux[num + 1]); // E
        _Cell[numA2+4] =  0.25 * dy * (dx * curT->z * _Uz[numZ]  + dxx * nT->z * _Uz[numZ  + 1]);   // T
				}
				break;
				case S_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const   coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
	double tx = curT->x * _Ux[num];
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dyz * (tx + (curT - 1)->x * _Ux[num-1]);                                       // W
    _Cell[numA2+2] =  0.25 * dyz * (tx +  nT->x * _Ux[num+1]);                                       // E
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+4] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case E_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
	const coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+4] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case W_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const   coordinate* bT = curT - baseNM;
	const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+4] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case T_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const   coordinate* bT = curT - baseNM;
    double tx = curT->x * _Ux[num]; const coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+3] =  0.25 * dyz * (tx + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+4] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
}
break;

case B_Border:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1;
    double tx = curT->x * _Ux[num]; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dyz * (tx + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+4] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
// ---------- Рёбра (4 элемента) ----------
case TN_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const  coordinate* bT = curT - baseNM;
    double tx = curT->x * _Ux[num]; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+3] =  0.25 * dyz * (tx + (curT + 1)->x * _Ux[num+1]);                           // E
}
break;
case TS_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
    double tx = curT->x * _Ux[num];
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dyz * (tx + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
}
break;

case TE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* bT = curT - baseNM; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * (curT + 1)->y * _Uy[numY+1]);      // N
}
break;
case TW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
	const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
}
break;
case BN_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const   coordinate* nT = curT + 1; const  coordinate* sT = curT - _N3;
    double tx = curT->x * _Ux[num];
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dyz * (tx + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case BS_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; double tx = curT->x * _Ux[num];
    _Cell[numA2]   = -0.25 * dyz * (tx + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+1] =  0.25 * dyz * (tx + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case BE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case BW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case NE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const   coordinate* curT = _theta + numZ; const coordinate* bT = curT - baseNM; const coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * (curT + 1)->z * _Uz[numZ+1]);      // T
}
break;
case NW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
	const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case SE_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;

case SW_Edge:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+3] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
// ---------- Углы (3 элемента) ----------
case TNE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const coordinate* bT = curT - baseNM; const  coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
}
break;
case TNW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const coordinate* bT = curT - baseNM; const coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+2] =  0.25 * dy * dz * (curT->x * _Ux[num] + (curT + 1)->x * _Ux[num+1]);                           // E
}
break;
case TSE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * (curT + 1)->y * _Uy[numY+1]);      // N
}
break;
case TSW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const  coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
    _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB+1]);   // B
    _Cell[numA2+1] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+2] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
}
break;
case BNE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const   coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+2] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * (curT + 1)->z * _Uz[numZ+1]);      // T
}
break;
case BNW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const coordinate* curT = _theta + numZ; const   coordinate* nT = curT + 1; const   coordinate* sT = curT - _N3;
    _Cell[numA2]   = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS+1]);   // S
    _Cell[numA2+1] =  0.25 * dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+2] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case BSE_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1;
    _Cell[numA2]   = -0.25 * dy * dz * (curT->x * _Ux[num] + (curT-1)->x * _Ux[num-1]);                     // W
    _Cell[numA2+1] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+2] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
case BSW_Corner:
{
    double dy = curD->y, dz = curD->z, dx = curD->x, dxx = (curD+1)->x;
	const   coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1;
    _Cell[numA2]   =  0.25 *  dy * dz * (curT->x * _Ux[num] + nT->x * _Ux[num+1]);                           // E
    _Cell[numA2+1] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);      // N
    _Cell[numA2+2] =  0.25 * dy * (dx * curT->z * _Uz[numZ] + dxx * nT->z * _Uz[numZ+1]);      // T
}
break;
				default:
				{
			double dy = curD->y, dz  = curD->z, dx = curD->x, dxx = (curD+1)->x, dyz = dy * dz;	
			const 	coordinate* curT = _theta + numZ; const  coordinate* nT = curT + 1; const  coordinate* bT = curT - baseNM;
			double tx = curT->x * _Ux[num]; const  coordinate* sT = curT - _N3;
        _Cell[numA2]   = -0.25 * dy * (dx * bT->z * _Uz[numB] + dxx * (bT+1)->z * _Uz[numB + 1]);   // B
        _Cell[numA2+1] = -0.25 * dz * (dx * sT->y * _Uy[numS] + dxx * (sT+1)->y * _Uy[numS + 1]);   // S
        _Cell[numA2+2] = -0.25 * dyz * (tx + (curT - 1)->x * _Ux[num - 1]); // W
        _Cell[numA2+3] =  0.25 * dyz * (tx + nT->x * _Ux[num + 1]); // E
        _Cell[numA2+4] =  0.25 * dz * (dx * curT->y * _Uy[numY] + dxx * nT->y * _Uy[numY+1]);  // N
        _Cell[numA2+5] =  0.25 * dy * (dx * curT->z * _Uz[numZ]  + dxx * nT->z * _Uz[numZ  + 1]);   // T
				}
				}//switch
			}// for i
		}
		LS_STAG_KERNEL_2D(baseFluxX, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, LS_STAG_CR_PD _Ux, LS_STAG_CR_PD _U_W,
			const double _U_E, const bool isArtery, const bool isCavity, const int _N, const int _M, const int _NM, 
			const int _vs, const int _vf, const int _zs, const int _zf)
		{
			int numF = i * _N + j * _NM, N2 = _N - 1;
			coordinate* curD = _delta + (numF + i + j * _M);
			double dyz = curD->z * curD->y;
			if((isArtery && ((i <= _vs) || (i >= _vf) || ((j < _zs) || (j > _zf)))))
			    _flux[numF] = 0.0;
			else
			{
				double curU_W = _U_W[i+j*_M];
				_flux[numF] = (-0.25 * curU_W * dyz * (curU_W + _Ux[numF]));
			}	 
			for (int num = numF + 1, m_end = numF + N2; num < m_end; num++)
				_flux[num] = 0.0;
			numF += N2;
			if(isCavity) 
			    _flux[numF] = 0.25 * _U_E * dyz * (_U_E + _Ux[numF]);
			else
			{
				double curU = _Ux[numF];
				_flux[numF] = 0.5 * dyz * curU * curU;
			}
		}
	LS_STAG_KERNEL(addFluxNS, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double UV_S, const double UV_N, 
		const int _N, const int _M, const int _NM)
		{
			int numS = i * _NM;
			coordinate* curD = m_delta + (numS + i * _M);
			double dz = 0.25 * curD->z, dx = curD->x; curD++;
			for(int ii = 0, numN = numS + _NM - _N; ii < _N; ii++, numS++, numN++, curD++)
			{
				double dxx = curD->x;
				double coef = dz * (dx + dxx);
				flux[numS] -= coef * UV_S;
				flux[numN] += coef * UV_N;
				dx = dxx;
			}
		}
	LS_STAG_KERNEL(addFluxTB, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta, const double UW_B, const double UW_T,
		const int _N, const int _NMK2)
		{
			int numB = i * _N;
			coordinate* curD = m_delta +  (numB + i);
			double dy = 0.25 * curD->y, dx = curD->x; curD++;
			for(int ii = 0, numT = _NMK2 + numB; ii < _N; ii++, numB++, numT++, curD++)
			{
				double dxx = curD->x;
				double coef = dy * (dx + dxx);
				if(_type[numB] != SolidCell) flux[numB] -= coef * UW_B;
				if (_type[numT] != SolidCell) flux[numT] += coef * UW_T;
				dx = dxx;
			}
		}	
	LS_STAG_KERNEL_2D(addFluxTBOut, LS_STAG_R_PD flux, LS_STAG_CR_PC m_theta, LS_STAG_CR_PC m_delta, LS_STAG_CR_PD _Ux, LS_STAG_CR_PD _Uz,
		const int _N, const int _NMx, const int _K, const int _Nbase, const int _NMbase, const int _Kw, const int _M2)
	{// конвективный выход на B/T (как элементы матрицы): S_x = -/+ (1/4)*dy*(dx*w_i + dxE*w_{i+1}) с theta.z
		const int numB = i + j * _N;                      // u-ячейка (i,j,0): u-сетка Nx = N-1
		const int numT = numB + _Kw * _NMx;          // u-ячейка (i,j,K-1)
		const int base = i + j * _Nbase;                  // базовая ячейка (i,j)
		const coordinate* curD = m_delta + base;
		const coordinate* curT = m_theta + base;
		const int wB = numB + j;                          // w-ячейка (i,j,0): w-сетка Nx = N, шаг по y = N
		// theta.z базовых ячеек B/T: вынос идёт только через жидкую (неусечённую) часть грани
		const double wfB = (_Uz[wB] * curD->x * curT->z + _Uz[wB + 1] * (curD + 1)->x * (curT + 1)->z);
		flux[numB] -= 0.25 * curD->y * _Ux[numB] * wfB;
		const int numToT = _Kw * _NMbase;
		const int wT = wB + numToT - _NMbase;          // w-ячейка (i,j,Kw-1) --- T-грань
		const coordinate* curDT = curD + numToT;
		const coordinate* curTT = curT + numToT;
		const double wfT = (_Uz[wT] * curDT->x * curTT->z + _Uz[wT + 1] * (curDT + 1)->x * (curTT + 1)->z);
		flux[numT] += 0.25 * curDT->y * _Ux[numT] * wfT;
	}
		LS_STAG_KERNEL_2D(ibFlux, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _baseSideType, LS_STAG_CR_PCT _baseType, coordinate* _theta,
			LS_STAG_R_PC _delta, LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibENT, LS_STAG_CR_PD _Ux, const double* _Uy, 
			const double* _Uz, LS_STAG_CR_PD _Uib, LS_STAG_CR_PD _Ug, const int _N, const int _M, const int _NM, 
			const int baseNM, const int _N3)
		{
			for (int numX = _N * i + _NM * j, shiftBase = i + j * _M, shiftY = j * _N3, shiftYN = shiftY + _N3, m_end = numX + _N; numX < m_end; numX++)
			{
				const CellTypes curType = _type[numX];
				if ((curType > SolidCell) && (curType < FluidCell))
				{
					int num = numX + shiftBase;
					double cxe = 0.0, cxw = 0.0;

					switch (_baseType[num])
					{
					case BSE_Heptahedron3:
					case BNE_Heptahedron3:
					case TSE_Heptahedron3:
					case TNE_Heptahedron3:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].x;
						break;
					case BS_Pentahedron:
					case TS_Pentahedron:
					case BN_Pentahedron:
					case TN_Pentahedron:
					case N_Hexahedron4:
					case S_Hexahedron4:
					case T_Hexahedron4:
					case B_Hexahedron4:
					case TSSW_Hexahedron5:
					case BSSW_Hexahedron5:
					case TSSE_Hexahedron5:
					case BSSE_Hexahedron5:
					case BSBW_Hexahedron5:
					case BSBE_Hexahedron5:
					case BNBW_Hexahedron5:
					case BNBE_Hexahedron5:
					case TNTW_Hexahedron5:
					case TNTE_Hexahedron5:
					case TSTW_Hexahedron5:
					case TSTE_Hexahedron5:
					case BNNW_Hexahedron5:
					case TNNW_Hexahedron5:
					case BNNE_Hexahedron5:
					case TNNE_Hexahedron5:
					case BS_Heptahedron4:
					case TS_Heptahedron4:
					case BN_Heptahedron4:
					case TN_Heptahedron4:
					case BSW_Heptahedron6:
					case BNW_Heptahedron6:
					case TSW_Heptahedron6:
					case TNW_Heptahedron6:
					case TNTW_Heptahedron5:
					case TSTW_Heptahedron5:
					case BNNW_Heptahedron5:
					case BSSW_Heptahedron5:
					case BSBW_Heptahedron5:
					case BNBW_Heptahedron5:
					case TSSW_Heptahedron5:
					case TNNW_Heptahedron5:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num].x;
						break;
					case W_Hexahedron4:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->x * curD->y * curD->z * _Ux[numX]) * _Vib[num].x;
					}
						break;
					case NE_Heptahedron4:
					{// 2D NE_Pentagon: 0.25*(Uib-Ug + theta^v*dSx*dSz*Uy)*vib  (тангенциальная Uy, theta^v=theta.y)
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->y * curD->x * curD->z * _Uy[num-shiftY]) * _VibENT[num].x;
					}
						break;
					case SE_Heptahedron4:
					{// 2D SE_Pentagon: 0.25*(Uib-Ug - theta^v_s*dSx*dSz*Uy_s)*vib
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * curD->x * curD->z * _Uy[num-shiftYN]) * _VibENT[num].x;
					}
						break;
					case SW_Pentahedron:
					case BW_Pentahedron:
					{// 2D SW_Triangle (общая формула, как BSW_Tetrahedron)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
							curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case NW_Pentahedron:
					case TW_Pentahedron:
					{// 2D NW_Triangle (как BNW_Tetrahedron)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
							curD->x * (-curVib->y * (curT-_N3)->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case TNTE_Heptahedron5:
					case TSTE_Heptahedron5:
					case BNNE_Heptahedron5:
					case TSSE_Heptahedron5:
					case TNNE_Heptahedron5:
					case BSSE_Heptahedron5:
					case BSBE_Heptahedron5:
					case BNBE_Heptahedron5:
					case TE_Heptahedron4:
					case BE_Heptahedron4:
					case BWSW_Hexahedron5:
					case BWNW_Hexahedron5:
					case TWSW_Hexahedron5:
					case TWNW_Hexahedron5:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->x * curD->y * curD->z * _Ux[numX]) * _VibENT[num].x;
					}
						break;
					case TENE_Heptahedron5:
					case TNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->y * dx * curD->z * _Uy[num-j*_N3] + 
							curT->z * dx * curD->y * _Uz[num]) * _VibENT[num].x;
					}
						break;
					case BENE_Heptahedron5:
					case BNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->y * dx * curD->z * _Uy[num-j*_N3] -  
							(curT-_NM)->z * dx * curD->y * _Uz[num-baseNM]) * _VibENT[num].x;
					}
						break;
					case TESE_Heptahedron5:
					case TSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * dx * curD->z * _Uy[num-(j+1)*_N3] + 
							curT->z * dx * curD->y * _Uz[num]) * _VibENT[num].x;
					}
						break;
					case BESE_Heptahedron5:
					case BSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * dx * curD->z * _Uy[num-(j+1)*_N3] - 
							(curT-_NM)->z * dx * curD->y * _Uz[num-baseNM]) * _VibENT[num].x;
					}
						break;
					case BSW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
							curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case BNW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
							curD->x * (-curVib->y * (curT-_N3)->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case TSW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
							curD->x * (curVib->y * curT->y * dz - curVib->z * (curT-baseNM)->z * dy) * _VibENT[num].x);
					}
					break;
					case TNW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x -
							curD->x * (curVib->y * (curT-_N3)->y * dz + curVib->z * (curT-baseNM)->z * dy) * _VibENT[num].x);
					}
					break;
					}//switch
					num++;
					switch (_baseType[num])
					{
					case BSW_Heptahedron3:
					case BNW_Heptahedron3:
					case TSW_Heptahedron3:
					case TNW_Heptahedron3:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].x;
						break;
					case BS_Pentahedron:
					case TS_Pentahedron:
					case BN_Pentahedron:
					case TN_Pentahedron:
					case N_Hexahedron4:
					case S_Hexahedron4:
					case T_Hexahedron4:
					case B_Hexahedron4:
					case TSSW_Hexahedron5:
					case BSSW_Hexahedron5:
					case TSSE_Hexahedron5:
					case BSSE_Hexahedron5:
					case BSBW_Hexahedron5:
					case BSBE_Hexahedron5:
					case BNBW_Hexahedron5:
					case BNBE_Hexahedron5:
					case TNTW_Hexahedron5:
					case TNTE_Hexahedron5:
					case TSTW_Hexahedron5:
					case TSTE_Hexahedron5:
					case BNNW_Hexahedron5:
					case TNNW_Hexahedron5:
					case BNNE_Hexahedron5:
					case TNNE_Hexahedron5:
					case BS_Heptahedron4:
					case TS_Heptahedron4:
					case BN_Heptahedron4:
					case TN_Heptahedron4:
					case BSE_Heptahedron6:
					case BNE_Heptahedron6:
					case TSE_Heptahedron6:
					case TNE_Heptahedron6:
					case TNTE_Heptahedron5:
					case TSTE_Heptahedron5:
					case BNNE_Heptahedron5:
					case BSSE_Heptahedron5:
					case BSBE_Heptahedron5:
					case BNBE_Heptahedron5:
					case TSSE_Heptahedron5:
					case TNNE_Heptahedron5:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num-1].x;
						break;
					case E_Hexahedron4:
					{// 2D E_Trapezoid: 0.25*(Uib-Ug + theta^u_w*dy*dz*Ux_w)*vib  (num уже = восточная ячейка)
						coordinate* curD = _delta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + (_theta + num - 1)->x * curD->y * curD->z * _Ux[numX]) * _Vib[num].x;
					}
						break;
					case NW_Heptahedron4:
					{// 2D NW_Pentagon: 0.25*(Uib-Ug + theta^v*dx*dz*Uy)*vib  (num уже = восточная ячейка)
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->y * curD->x * curD->z * _Uy[num-shiftY]) * _VibENT[num].x;
					}
						break;
					case SW_Heptahedron4:
					{// 2D SW_Pentagon: 0.25*(Uib-Ug - theta^v_s*dx*dz*Uy_s)*vib
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * curD->x * curD->z * _Uy[num-shiftYN]) * _VibENT[num].x;
					}
						break;
					case SE_Pentahedron:
					case BE_Pentahedron:
					{// 2D SE_Triangle (как BSE_Tetrahedron)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x -
							curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case NE_Pentahedron:
					case TE_Pentahedron:
					{// 2D NE_Triangle (как BNE_Tetrahedron)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x +
							curD->x * (curVib->y * (curT-_N3)->y * dz - curVib->z * curT->z * dy) * _VibENT[num].x);
					}
					break;
					case TNTW_Heptahedron5:
					case TSTW_Heptahedron5:
					case BNNW_Heptahedron5:
					case TSSW_Heptahedron5:
					case TNNW_Heptahedron5:
					case BSSW_Heptahedron5:
					case BSBW_Heptahedron5:
					case BNBW_Heptahedron5:
					case TW_Heptahedron4:
					case BW_Heptahedron4:
					case BESE_Hexahedron5:
					case BENE_Hexahedron5:
					case TESE_Hexahedron5:
					case TENE_Hexahedron5:
					{
						coordinate* curD = _delta + num;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num - 1)->x * curD->y * curD->z * _Ux[numX]) * _VibENT[num-1].x;
					}
						break;
					case TWNW_Heptahedron5:
					case TNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->y * dx * curD->z * _Uy[num-j*_N3] + 
							curT->z * dx * curD->y * _Uz[num]) * _VibENT[num-1].x;
					}
						break;
					case BWNW_Heptahedron5:
					case BNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->y * dx * curD->z * _Uy[num-j*_N3] -  
							(curT-_NM)->z * dx * curD->y * _Uz[num-baseNM]) * _VibENT[num-1].x;
					}
						break;
					case TWSW_Heptahedron5:
					case TSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * dx * curD->z * _Uy[num-(j+1)*_N3] + 
							curT->z * dx * curD->y * _Uz[num]) * _VibENT[num-1].x;
					}
						break;
					case BWSW_Heptahedron5:
					case BSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * dx * curD->z * _Uy[num-(j+1)*_N3] - 
							(curT-_NM)->z * dx * curD->y * _Uz[num-baseNM]) * _VibENT[num-1].x;
					}
						break;
					case BSE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x -
							curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num-1].x);
					}
					break;
					case BNE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x +
							curD->x * (curVib->y * (curT-_N3)->y * dz - curVib->z * curT->z * dy) * _VibENT[num-1].x);
					}
					break;
					case TSE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x -
							curD->x * (curVib->y * curT->y * dz - curVib->z * (curT-baseNM)->z * dy) * _VibENT[num-1].x);
					}
					break;
					case TNE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dy = curD->y, dz = curD->z;
						cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x +
							curD->x * (curVib->y * (curT-_N3)->y * dz + curVib->z * (curT-baseNM)->z * dy) * _VibENT[num-1].x);
					}
					break;
					}// switch
					_flux[numX] = cxe + cxw;
				}//if type >9
				else if ((curType == B_Border) || (curType == T_Border))
				{
					int num = numX + shiftBase;
					const SideTypes wzSide = _baseSideType[num].z;
					const SideTypes ezSide = _baseSideType[num+1].z;
					const bool wCutCell = ((wzSide > SolidSide) && (wzSide < FluidSide));
					const bool eCutCell = ((ezSide > SolidSide) && (ezSide < FluidSide));
					if(wCutCell || eCutCell)
					{
						double cxe = 0.0, cxw = 0.0;
						if(wCutCell)
						{// W-половинка (грань num): только западные формы среза
							const coordinate* curD = _delta + num;
							switch (wzSide)
							{
								case NW_Triangle:
								{// 2D NW_Triangle (общая формула, как NW_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dy = curD->y, dz = curD->z;
									cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
										curD->x * (-curVib->y * (curT-_N3)->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
								}
									break;
								case SW_Triangle:
								{// 2D SW_Triangle (общая формула, как SW_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dy = curD->y, dz = curD->z;
									cxe = -0.25 * (((curVib->x + _Ux[numX]) * dy * dz * curT->x + _Ug[num]) * curVib->x +
										curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
								}
									break;
								case NW_Pentagon:
								{// 2D NW_Pentagon (тангенциальная Uy), как NW_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->y * curD->x * curD->z * _Uy[num-shiftY]) * _VibENT[num].x;
								}
									break;
								case SW_Pentagon:
								{// 2D SW_Pentagon (тангенциальная Uy_s), как SW_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * curD->x * curD->z * _Uy[num-shiftYN]) * _VibENT[num].x;
								}
									break;
								case W_Trapezoid:// как W_Hexahedron4 (с Ux — Vib)
									cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->x * curD->y * curD->z * _Ux[numX]) * _Vib[num].x;
									break;
								case N_Trapezoid:// как N_Hexahedron4 (простая — VibENT)
								case S_Trapezoid:// как S_Hexahedron4
									cxe = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num].x;
									break;
							}
						}
						num++;
						if(eCutCell)
						{// E-половинка (грань num+1): только восточные формы среза
							const coordinate* curD = _delta + num;
							switch (ezSide)
							{
								case NE_Triangle:
								{// 2D NE_Triangle (общая формула, как NE_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dy = curD->y, dz = curD->z;
									cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x +
										curD->x * (curVib->y * (curT-_N3)->y * dz - curVib->z * curT->z * dy) * _VibENT[num].x);
								}
									break;
								case SE_Triangle:
								{// 2D SE_Triangle (общая формула, как SE_Pentahedron)
									const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									double dy = curD->y, dz = curD->z;
									cxw = 0.25 * (((curVib->x + _Ux[numX]) * dy * dz * (curT-1)->x - _Ug[num]) * curVib->x -
										curD->x * (curVib->y * curT->y * dz + curVib->z * curT->z * dy) * _VibENT[num].x);
								}
									break;
								case NE_Pentagon:
								{// 2D NE_Pentagon (тангенциальная Uy), как NE_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->y * curD->x * curD->z * _Uy[num-shiftY]) * _VibENT[num].x;
								}
									break;
								case SE_Pentagon:
								{// 2D SE_Pentagon (тангенциальная Uy_s), как SE_Heptahedron4
									coordinate* curD = _delta + num; coordinate* curT = _theta + num;
									cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT-_N3)->y * curD->x * curD->z * _Uy[num-shiftYN]) * _VibENT[num].x;
								}
									break;
								case E_Trapezoid:// как E_Hexahedron4 (2D: +theta^u_w*dy*dz*Ux_w)
									cxw = 0.25 * (_Uib[num] - _Ug[num] + (_theta + num - 1)->x * curD->y * curD->z * _Ux[numX]) * _Vib[num].x;
									break;
								case N_Trapezoid:// как N_Hexahedron4 (простая — VibENT)
								case S_Trapezoid:// как S_Hexahedron4
									cxw = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num-1].x;
									break;
							}
						}
						_flux[numX] = cxe + cxw;
					}
				}
			}//for i
		}
	  private:
		  /*!
		  Рассчитывает значения источникового члена оператора конвекции на х-сетке на границах расчетной области
		  */
		void baseFluxX()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorX, baseFluxX, m_M, m_K, flux, m_delta, m_Ux, 
				m_task->pU_W, m_task->U_E, false, false, 
				m_N, m_M, m_N * m_M, vs, vf, zs, zf)
			LS_STAG_PAR_FOR(ConvectiveDiscreteOperatorX, addFluxNS, m_K, flux, m_delta, m_task->U_S * m_task->V_S, 
				m_task->U_N * m_task->V_N, m_N, m_M, m_N * m_M)
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorX, addFluxTBOut, m_N, m_M, flux, m_theta, m_delta, m_Ux, m_Uz,
				m_N, m_N * m_M, m_K, N3, N3 * m_M, m_K - 1, M2)
		};
		  /*!
		  Рассчитывает значения источникового члена оператора конвекции на х-сетке на погруженных границах 
		  */
		void ibFlux()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorX, ibFlux, m_M, m_K, flux, m_type, m_baseSideTypes, m_meshes->baseMesh.typeOfCell,
				m_theta, m_delta, m_Vib, m_VibENT, m_Ux, m_Uy, m_Uz, m_Uib, m_Ug, m_N, m_M, m_N * m_M, N3 * m_M, N3)
		};
	  public:
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора конвекции на х-сетке
		  */
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorX, calcConvX, m_M, m_K, Number, Cell, m_theta, m_delta, m_Ux, m_Uy, 
				m_Uz, m_type, m_N, m_M, m_N * m_M, N3 * m_M, N3)
		};
	};
}

#endif //_LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_X_H_
