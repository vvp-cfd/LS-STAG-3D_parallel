/*
* \copyright (C) 2026 Puzikova V.V.
* \file      convective_discrete_operator_z.h
* \author    Puzikova Valeria (vvp) 
* \date      13.06.2026 22:39:29
* \brief     convective_discrete_operator_z.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Z_H_
#define _LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Z_H_

// ================================================================== INCLUDES
#include "skew_sym_discrete_operator.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ConvectiveDiscreteOperatorZ

Данный файл содержит в себе реализацию класса ConvectiveDiscreteOperatorZ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(baseFluxz, LS_STAG_R_PD, _flux, LS_STAG_R_PC, _delta,
		const double*, _Ux, LS_STAG_CR_PD, _Uz, LS_STAG_CR_PD, _U_W, const double, _W_W, const double, _UW_E,
		const bool, isArtery, const bool, isCavity, const int, _N, const int, _M, const int, _NM, const int, _vs,
		const int, _vf, const int, _zs, const int, _zf);
	LS_STAG_DECLARATION_KERNEL_2D(calcConvZ, LS_STAG_CR_PI, _Number, LS_STAG_R_PD, _Cell, const coordinate*, _theta,
		coordinate*, _delta, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PCT, _type,
		const int, _N, const int, _M, const int, _NM, const int, baseNM, const int, _N2);
	LS_STAG_DECLARATION_KERNEL_2D(ibFluxConvZ, LS_STAG_R_PD, _flux, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _baseSideType,
		LS_STAG_CR_PCT, _baseType, coordinate*, _theta, coordinate*, _delta, LS_STAG_CR_PC, _Vib,
		LS_STAG_CR_PC, _VibENT, LS_STAG_CR_PD, _Ux, const double*, _Uy, const double*, _Uz, LS_STAG_CR_PD, _Uib,
		LS_STAG_CR_PD, _Ug, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxNSz, LS_STAG_R_PD, flux, LS_STAG_R_PC, m_delta,
		const double, VW_S, const double, VW_N, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(addFluxTBz, LS_STAG_R_PD, flux, LS_STAG_CR_PCT, _type, LS_STAG_R_PC, m_delta,
		LS_STAG_CR_PD, m_Uz, const double, W_B, const double, W_T, const int, _N, const int, _NMK2);
	LS_STAG_DECLARATION_KERNEL_2D(addFluxTBOutZ, LS_STAG_R_PD, flux, LS_STAG_CR_PC, m_theta, LS_STAG_CR_PC, m_delta, LS_STAG_CR_PD, _Uz,
		const int, _N, const int, NM, const int, _K);
	/*!
	\brief Класс, реализующий дискретный оператор конвекции на z-сетке

	Данный класс реализует работу с дискретным аналогом оператора конвекции на z-сетке
	*/
	struct ConvectiveDiscreteOperatorZ : SkewSymDiscreteOperator
	{
		LS_STAG_KERNEL(addFluxNSz, LS_STAG_R_PD flux, LS_STAG_R_PC m_delta, const double VW_S, const double VW_N,
			const int _N, const int _NM)
		{
			int numS = i * _NM;
			coordinate* curD = m_delta + numS;
			double dzz = 0.25 * (curD->z + (curD + _NM)->z);
			for (int ii = 0, numN = numS + _NM - _N; ii < _N; ii++, numS++, numN++, curD++)
			{
				double coef = curD->x * dzz;
				flux[numS] -= coef * VW_S;
				flux[numN] += coef * VW_N;
			}
		};
		LS_STAG_KERNEL(addFluxTBz, LS_STAG_R_PD flux, LS_STAG_CR_PCT _type, LS_STAG_R_PC m_delta, LS_STAG_CR_PD m_Uz, const double W_B,
			const double W_T, const int _N, const int _NMK2)
		{
			int numB = i * _N;
			coordinate* curD = m_delta + numB;
			double dy = 0.25 * curD->y;
			for (int ii = 0, numT = _NMK2 + numB; ii < _N; ii++, numB++, numT++, curD++)
			{
				double coef = dy * curD->x;
				if (_type[numB] != SolidCell) flux[numB] -= coef * W_B * (m_Uz[numB] + W_B);
				if (_type[numT] != SolidCell) flux[numT] += coef * W_T * (m_Uz[numT] + W_T);
			}
		};
		LS_STAG_KERNEL_2D(addFluxTBOutZ, LS_STAG_R_PD flux, LS_STAG_CR_PC m_theta, LS_STAG_CR_PC m_delta, LS_STAG_CR_PD _Uz,
			const int _N, const int NM, const int _K)
		{// конвективный выход на B/T (как элементы матрицы): S_z = -/+ (1/4)*dx*dy*(theta-взвеш. сумма w по z)
			const int numB = i + j * _N;
			const int numT = numB + (_K - 1) * NM;
			const coordinate* curD = m_delta + numB;
			const coordinate* curT = m_theta + numB;
			const double wB = _Uz[numB];
			flux[numB] -= 0.5 * curD->x * curD->y * curT->z * wB * wB;
			const int numToT = (_K - 1) * NM;
			const coordinate* curDT = curD + numToT;
			const double wT = _Uz[numT];
			flux[numT] += 0.5 * curDT->x * curDT->y * (curT + numToT)->z * wT * wT;
		}
		LS_STAG_KERNEL_2D(calcConvZ, LS_STAG_CR_PI _Number, LS_STAG_R_PD _Cell, const coordinate * _theta, coordinate * _delta,
			LS_STAG_CR_PD _Ux, const double* _Uy, const double* _Uz, LS_STAG_CR_PCT _type, const int _N, const int _M,
			const int _NM, const int _NMy, const int _NMx)
		{
			int num = _N * i + j * _NM;// = numZ
			coordinate* curD = _delta + num;
			for (int numX = num - i - j * _M, numY = num - j * _N, m_end = num + _N, numS = numY - _N, numT = numX + _NMx; num < m_end;
				num++, numX++, numY++, numT++, numS++, curD++)
			{
				const int numA2 = _Number[num];
				switch (_type[num])
				{
				case FluidCell:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy, curW = _Uz[num];
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);                 // W
					_Cell[numA2 + 3] = 0.25 * dy * (dz * _Ux[numX] + dzz * nT->x * _Ux[numT]);                 // E
					_Cell[numA2 + 4] = 0.25 * dx * (dz * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]);  // N
					_Cell[numA2 + 5] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]);  // T
				}
				break;
				case SolidCell:
					_Cell[numA2] = 0.0;
					break;
				case N_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					double curW = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);                 // W
					_Cell[numA2 + 3] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]);                 // E
					_Cell[numA2 + 4] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]);  // T
				}
				break;
				case S_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					double curW = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]); // W
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]);// N
					_Cell[numA2 + 4] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case E_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					double curW = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 4] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case W_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const  coordinate * sT = curT - _N;
					double curW = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 4] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case T_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]); // W
					_Cell[numA2 + 3] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 4] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]);// N
				}
				break;
				case B_Border:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 4] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				// ---------- Рёбра (4 элемента) ----------
				case TN_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const  coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 3] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
				}
				break;
				case TS_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * nT = curT + _NM;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
				}
				break;
				case TE_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
				}
				break;
				case TW_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
				}
				break;
				case BN_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BS_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					_Cell[numA2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 1] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 3] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BE_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const  coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 3] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BW_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 3] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case NE_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					double tz = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (tz + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 3] = 0.25 * dxy * (tz + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case NW_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					double tz = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (tz + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 3] = 0.25 * dxy * (tz + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case SE_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; double tz = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (tz + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 3] = 0.25 * dxy * (tz + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case SW_Edge:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; double tz = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (tz + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 3] = 0.25 * dxy * (tz + nT->z * _Uz[num + _NM]); // T
				}
				break;
				// ---------- Углы (3 элемента) ----------
				case TNE_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * ((curT + _NM) - 1)->x * _Ux[numT - 1]); // W (nT=curT+_NM)
				}
				break;
				case TNW_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * (curT + _NM)->x * _Ux[numT]); // E
				}
				break;
				case TSE_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
				}
				break;
				case TSW_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					_Cell[numA2] = -0.25 * dx * dy * (curT->z * _Uz[num] + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 2] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
				}
				break;
				case BNE_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 2] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BNW_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const  coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					_Cell[numA2] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 1] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 2] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BSE_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					_Cell[numA2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);// W
					_Cell[numA2 + 1] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 2] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				case BSW_Corner:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM;
					_Cell[numA2] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]); // E
					_Cell[numA2 + 1] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]); // N
					_Cell[numA2 + 2] = 0.25 * dx * dy * (curT->z * _Uz[num] + nT->z * _Uz[num + _NM]); // T
				}
				break;
				default:
				{
					double dy = curD->y, dz = curD->z, dx = curD->x, dzz = (curD + _NM)->z, dxy = dx * dy;
					const coordinate * curT = _theta + num; const coordinate * nT = curT + _NM; const coordinate * sT = curT - _N;
					double curW = curT->z * _Uz[num];
					_Cell[numA2] = -0.25 * dxy * (curW + (curT - _NM)->z * _Uz[num - _NM]); // B
					_Cell[numA2 + 1] = -0.25 * dx * (dz * sT->y * _Uy[numS] + dzz * (sT + _NM)->y * _Uy[numS + _NMy]); // S
					_Cell[numA2 + 2] = -0.25 * dy * (dz * (curT - 1)->x * _Ux[numX - 1] + dzz * (nT - 1)->x * _Ux[numT - 1]);                 // W
					_Cell[numA2 + 3] = 0.25 * dy * (dz * curT->x * _Ux[numX] + dzz * nT->x * _Ux[numT]);                 // E
					_Cell[numA2 + 4] = 0.25 * dx * (dz * curT->y * _Uy[numY] + dzz * nT->y * _Uy[numY + _NMy]);  // N
					_Cell[numA2 + 5] = 0.25 * dxy * (curW + nT->z * _Uz[num + _NM]);  // T
				}
				}//switch
			}// for i 
		};
		LS_STAG_KERNEL_2D(baseFluxz, LS_STAG_R_PD _flux, LS_STAG_R_PC _delta, const double* _Ux, LS_STAG_CR_PD _Uz,
			LS_STAG_CR_PD _U_W, const double _W_W, const double _UW_E, const bool isArtery, const bool isCavity,
			const int _N, const int _M, const int _NM, const int _vs, const int _vf, const int _zs, const int _zf)
		{
			int numF = _N * i + _NM * j, _N2 = _N - 1, numW = i + j * _M;
			coordinate* curD = _delta + numF;
			double dz = curD->z, dy = curD->y, dzz = (curD + _NM)->z;
			_flux[numF] = ((isArtery && ((i <= _vs) || (i >= _vf) || ((j < _zs) || (j > _zf))))) ? 0.0 :
				(-0.25 * (_U_W[numW] * dz + _U_W[numW + _M] * dzz) * _W_W * dy);
			for (int num = numF + 1, m_end = numF + _N2; num < m_end; num++)
				_flux[num] = 0.0;
			numF += _N2;
			if (isCavity)
				_flux[numF] = 0.25 * _UW_E * (dz + dzz) * dy;
			else
			{
				const int num2 = numF - numW - 1;
				_flux[numF] = 0.25 * _Uz[numF] * (_Ux[num2] * dz + _Ux[num2 + _N2 * _M] * dzz) * dy;
			}
		};
		LS_STAG_KERNEL_2D(ibFluxConvZ, LS_STAG_R_PD _flux, LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _baseSideType, LS_STAG_CR_PCT _baseType,
			coordinate * _theta, coordinate * _delta, LS_STAG_CR_PC _Vib, LS_STAG_CR_PC _VibENT,
			LS_STAG_CR_PD _Ux, const double* _Uy, const double* _Uz, LS_STAG_CR_PD _Uib, LS_STAG_CR_PD _Ug, const int _N,
			const int _M, const int _NM)
		{
			for (int num = _N * i + _NM * j, m_end = num + _N; num < m_end; num++)
			{
				const CellTypes curType = _type[num];
				if ((curType > SolidCell) && (curType < FluidCell))
				{
					double cxe = 0.0, cxw = 0.0;
					switch (_baseType[num])
					{
					case BSE_Heptahedron3:
					case BSW_Heptahedron3:
					case BNE_Heptahedron3:
					case BNW_Heptahedron3:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].z;
						break;
					case E_Hexahedron4:
					case W_Hexahedron4:
					case N_Hexahedron4:
					case S_Hexahedron4:
					case BSSW_Hexahedron5:
					case BWSW_Hexahedron5:
					case BWNW_Hexahedron5:
					case TWNW_Hexahedron5:
					case TWSW_Hexahedron5:
					case TESE_Hexahedron5:
					case TENE_Hexahedron5:
					case BESE_Hexahedron5:
					case BENE_Hexahedron5:
					case BNNW_Hexahedron5:
					case BNNE_Hexahedron5:
					case TNNW_Hexahedron5:
					case TNNE_Hexahedron5:
					case TSSW_Hexahedron5:
					case TSSE_Hexahedron5:
					case BSSE_Hexahedron5:
					case SW_Heptahedron4:
					case NW_Heptahedron4:
					case SE_Heptahedron4:
					case NE_Heptahedron4:
					case BNW_Heptahedron6:
					case BNE_Heptahedron6:
					case BSW_Heptahedron6:
					case BSE_Heptahedron6:
					case BENE_Heptahedron5:
					case BESE_Heptahedron5:
					case BNNE_Heptahedron5:
					case BNNW_Heptahedron5:
					case BSSE_Heptahedron5:
					case BSSW_Heptahedron5:
					case BWSW_Heptahedron5:
					case BWNW_Heptahedron5:
						cxe = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num].z;
						break;
					case B_Hexahedron4:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->z * curD->x * curD->y * _Uz[num]) * _Vib[num].z;
					}
					break;
					case BW_Pentahedron:
					case BE_Pentahedron:
					case BS_Pentahedron:
					case BN_Pentahedron:
					case TENE_Heptahedron5:
					case TWNW_Heptahedron5:
					case TESE_Heptahedron5:
					case TWSW_Heptahedron5:
					case TNNE_Heptahedron5:
					case TNNW_Heptahedron5:
					case TSSE_Heptahedron5:
					case TSSW_Heptahedron5:
					case TE_Heptahedron4:
					case TW_Heptahedron4:
					case TN_Heptahedron4:
					case TS_Heptahedron4:
					case BSBW_Hexahedron5:
					case BSBE_Hexahedron5:
					case BNBW_Hexahedron5:
					case BNBE_Hexahedron5:
					{
						coordinate* curD = _delta + num;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num)->z * curD->x * curD->y * _Uz[num]) * _VibENT[num].z;
					}
					break;
					case TNTE_Heptahedron5:
					case TNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dz * curD->y * _Ux[num - i - j * _M] +
							curT->y * dz * curD->x * _Uy[num - j * _N]) * _VibENT[num].z;
					}
					break;
					case TSTE_Heptahedron5:
					case TSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxe = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dz * curD->y * _Ux[num - i - j * _M] -
							(curT - _N)->y * dz * curD->x * _Uy[num - (j + 1) * _N]) * _VibENT[num].z;
					}
					break;
					case TNTW_Heptahedron5:
					case TNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * dz * curD->y * _Ux[num - 1 - i - j * _M] +
							curT->y * dz * curD->x * _Uy[num - j * _N]) * _VibENT[num].z;
					}
					break;
					case TSTW_Heptahedron5:
					case TSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxe = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * dz * curD->y * _Ux[num - 1 - i - j * _M] +
							(curT - _N)->y * dz * curD->x * _Uy[num - (j + 1) * _N]) * _VibENT[num].z;
					}
					break;
					case BSW_Tetrahedron:
					case SW_Pentahedron:
					{// 2D SW_Triangle: общая формула (как в x/y)
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y;
						cxe = -0.25 * (((curVib->z + _Uz[num]) * dx * dy * curT->z + _Ug[num]) * curVib->z +
							curD->z * (curVib->x * curT->x * dy + curVib->y * curT->y * dx) * _VibENT[num].z);
					}
					break;
					case BSE_Tetrahedron:
					case SE_Pentahedron:
					{// 2D SE_Triangle: общая формула
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y;
						cxe = -0.25 * (((curVib->z + _Uz[num]) * dx * dy * curT->z + _Ug[num]) * curVib->z +
							curD->z * (-curVib->x * (curT - 1)->x * dy + curVib->y * curT->y * dx) * _VibENT[num].z);
					}
					break;
					case BNW_Tetrahedron:
					case NW_Pentahedron:
					{// 2D NW_Triangle: общая формула
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y;
						cxe = -0.25 * (((curVib->z + _Uz[num]) * dx * dy * curT->z + _Ug[num]) * curVib->z +
							curD->z * (curVib->x * curT->x * dy - curVib->y * (curT - _N)->y * dx) * _VibENT[num].z);
					}
					break;
					case BNE_Tetrahedron:
					case NE_Pentahedron:
					{// 2D NE_Triangle: общая формула
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y;
						cxe = -0.25 * (((curVib->z + _Uz[num]) * dx * dy * curT->z + _Ug[num]) * curVib->z -
							curD->z * (curVib->x * (curT - 1)->x * dy + curVib->y * (curT - _N)->y * dx) * _VibENT[num].z);
					}
					break;
					}//switch
					num += _NM;
					switch (_baseType[num])
					{
					case TSE_Heptahedron3:
					case TSW_Heptahedron3:
					case TNE_Heptahedron3:
					case TNW_Heptahedron3:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _Vib[num].z;
						break;
					case SW_Pentahedron:
					case NW_Pentahedron:
					case NE_Pentahedron:
					case SE_Pentahedron:
					case E_Hexahedron4:
					case W_Hexahedron4:
					case N_Hexahedron4:
					case S_Hexahedron4:
					case BSSW_Hexahedron5:
					case BWSW_Hexahedron5:
					case BWNW_Hexahedron5:
					case TWNW_Hexahedron5:
					case TWSW_Hexahedron5:
					case TESE_Hexahedron5:
					case TENE_Hexahedron5:
					case BESE_Hexahedron5:
					case BENE_Hexahedron5:
					case BNNW_Hexahedron5:
					case BNNE_Hexahedron5:
					case TNNW_Hexahedron5:
					case TNNE_Hexahedron5:
					case TSSW_Hexahedron5:
					case TSSE_Hexahedron5:
					case BSSE_Hexahedron5:
					case SW_Heptahedron4:
					case NW_Heptahedron4:
					case SE_Heptahedron4:
					case NE_Heptahedron4:
					case TNW_Heptahedron6:
					case TNE_Heptahedron6:
					case TSW_Heptahedron6:
					case TSE_Heptahedron6:
					case TENE_Heptahedron5:
					case TESE_Heptahedron5:
					case TNNE_Heptahedron5:
					case TNNW_Heptahedron5:
					case TSSE_Heptahedron5:
					case TSSW_Heptahedron5:
					case TWSW_Heptahedron5:
					case TWNW_Heptahedron5:
						cxw = 0.25 * (_Uib[num] - _Ug[num]) * _VibENT[num - _NM].z;
						break;
					case T_Hexahedron4:
					{
						coordinate* curD = _delta + num; int num2 = num - _NM;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num2)->z * curD->x * curD->y * _Uz[num2]) * _Vib[num].z;
					}
					break;
					case TW_Pentahedron:
					case TE_Pentahedron:
					case TS_Pentahedron:
					case TN_Pentahedron:
					case BENE_Heptahedron5:
					case BWNW_Heptahedron5:
					case BESE_Heptahedron5:
					case BWSW_Heptahedron5:
					case BNNE_Heptahedron5:
					case BNNW_Heptahedron5:
					case BSSE_Heptahedron5:
					case BSSW_Heptahedron5:
					case BE_Heptahedron4:
					case BW_Heptahedron4:
					case BN_Heptahedron4:
					case BS_Heptahedron4:
					case TSTW_Hexahedron5:
					case TSTE_Hexahedron5:
					case TNTW_Hexahedron5:
					case TNTE_Hexahedron5:
					{
						coordinate* curD = _delta + num; int num2 = num - _NM;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (_theta + num2)->z * curD->x * curD->y * _Uz[num2]) * _VibENT[num2].z;
					}
					break;
					case BNBE_Heptahedron5:
					case BNE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dz * curD->y * _Ux[num - i - j * _M] +
							curT->y * dz * curD->x * _Uy[num - j * _N]) * _VibENT[num - _NM].z;
					}
					break;
					case BSBE_Heptahedron5:
					case BSE_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxw = 0.25 * (_Uib[num] - _Ug[num] + curT->x * dz * curD->y * _Ux[num - i - j * _M] -
							(curT - _N)->y * dz * curD->x * _Uy[num - (j + 1) * _N]) * _VibENT[num - _NM].z;
					}
					break;
					case BNBW_Heptahedron5:
					case BNW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * dz * curD->y * _Ux[num - 1 - i - j * _M] +
							curT->y * dz * curD->x * _Uy[num - j * _N]) * _VibENT[num - _NM].z;
					}
					break;
					case BSBW_Heptahedron5:
					case BSW_Heptahedron6:
					{
						coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dz = curD->z;
						cxw = 0.25 * (_Uib[num] - _Ug[num] - (curT - 1)->x * dz * curD->y * _Ux[num - 1 - i - j * _M] +
							(curT - _N)->y * dz * curD->x * _Uy[num - (j + 1) * _N]) * _VibENT[num - _NM].z;
					}
					break;
					case TSW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y; int num2 = num - _NM;
						cxw = 0.25 * (((curVib->z + _Uz[num2]) * dx * dy * (curT - _NM)->z - _Ug[num]) * curVib->z -
							curD->z * (curVib->x * curT->x * dy + curVib->y * curT->y * dx) * _VibENT[num2].z);
					}
					break;
					case TSE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y; int num2 = num - _NM;
						cxw = 0.25 * (((curVib->z + _Uz[num2]) * dx * dy * (curT - _NM)->z - _Ug[num]) * curVib->z -
							curD->z * (-curVib->x * (curT - 1)->x * dy + curVib->y * curT->y * dx) * _VibENT[num2].z);
					}
					break;
					case TNW_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y; int num2 = num - _NM;
						cxw = 0.25 * (((curVib->z + _Uz[num2]) * dx * dy * (curT - _NM)->z - _Ug[num]) * curVib->z -
							curD->z * (curVib->x * curT->x * dy - curVib->y * (curT - _N)->y * dx) * _VibENT[num2].z);
					}
					break;
					case TNE_Tetrahedron:
					{
						const coordinate* curVib = _Vib + num; coordinate* curD = _delta + num; coordinate* curT = _theta + num;
						double dx = curD->x, dy = curD->y; int num2 = num - _NM;
						cxw = 0.25 * (((curVib->z + _Uz[num2]) * dx * dy * (curT - _NM)->z - _Ug[num]) * curVib->z +
							curD->z * (curVib->x * (curT - 1)->x * dy + curVib->y * (curT - _N)->y * dx) * _VibENT[num2].z);
					}
					break;
					}//switch
					num -= _NM;
					_flux[num] = cxe + cxw;
				}// if type > 9
				else if ((curType == B_Border) || (curType == T_Border))
				{
					const SideTypes zSide = _baseSideType[num].z;
					if((zSide > SolidSide) && (zSide < FluidSide))
					{
						const int numT = num + _NM;
						_flux[num] = 0.25 * (_Uib[num] + _Uib[numT] - _Ug[num] - _Ug[numT]) * _VibENT[num].z;
					}
				}
			}//for i
		};
	private:
		/*!
		Рассчитывает значения источникового члена оператора конвекции на у-сетке на границах расчетной области
		*/
		void baseFlux()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorZ, baseFluxz, m_M, m_K, flux, m_delta, m_Ux, m_Uz,
				m_task->pU_W, m_task->W_W, m_task->U_E * m_task->W_E, false,
				false, m_N, m_M, m_N * m_M, vs, vf, zs, zf)
			LS_STAG_PAR_FOR(ConvectiveDiscreteOperatorZ, addFluxNSz, m_K, flux, m_delta, m_task->W_S * m_task->V_S,
					m_task->W_N * m_task->V_N, m_N, m_N * m_M)
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorZ, addFluxTBOutZ, m_N, m_M, flux, m_theta, m_delta, m_Uz,
				m_N, m_N * m_M, m_K)
		};
		/*!
		Рассчитывает значения источникового члена оператора конвекции на у-сетке на погруженных границах
		*/
		void ibFlux()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorZ, ibFluxConvZ, m_M, m_K, flux, m_type, m_baseSideTypes, m_meshes->baseMesh.typeOfCell,
				m_theta, m_delta, m_Vib, m_VibENT, m_Ux, m_Uy, m_Uz, m_Uib, m_Ug, m_N, m_M, m_N * m_M)
		}
	public:
		/*!
		Рассчитывает значения элементов разреженной матрицы дискретного аналога оператора конвекции на у-сетке
		*/
		void calcElem()
		{
			LS_STAG_PAR_FOR_2D(ConvectiveDiscreteOperatorZ, calcConvZ, m_M, m_K, Number, Cell, m_theta, m_delta, m_Ux,
				m_Uy, m_Uz, m_type, m_N, M3, m_N * m_M, m_N * M2, N2 * m_M)
		};
	};
}

#endif //_LS_STAG_CONVECTIVE_DISCRETE_OPERATOR_Y_H_
