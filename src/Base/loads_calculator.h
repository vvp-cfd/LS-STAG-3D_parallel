/*
* \copyright (C) 2013 Puzikova V.V.
* \file      loads_calculator.h
* \author    Puzikova Valeria (vvp) 
* \date      25.10.2013 23:34:20
* \brief     loads_calculator.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LOADS_CALCULATOR_H_
#define _LS_STAG_LOADS_CALCULATOR_H_

// ================================================================== INCLUDES
#include "tstring.h"
#include "LS_STAG_Mesh/cut_cells.h"
#include "MeshFunctions/mesh_function.h"
#include "MeshFunctions/coord_mesh_function.h"
#include "butterworth_low_pass_filter.h"
#include "LS_STAG_Mesh/ls_stag_uniform_mesh.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LoadsCalculator

Данный файл содержит в себе реализацию класса LoadsCalculator
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(calc, LS_STAG_CR_PCT, _type, LS_STAG_CR_PPI, _OnAirfoil, LS_STAG_R_PD, _ForceX,
		LS_STAG_R_PD, _ForceY, LS_STAG_R_PD, _ForceZ, LS_STAG_CR_PC, _Vib, LS_STAG_R_PC, _VibYX, LS_STAG_R_PC, _VibZY,
		LS_STAG_CR_PC, _delta, LS_STAG_R_PC, _inv_delta_ib_ne, LS_STAG_R_PC, _inv_delta_ib_tn, 
		LS_STAG_CR_PC, _delta_ib_max, LS_STAG_CR_PC, _theta, LS_STAG_CR_PD, _Ux, LS_STAG_CR_PD, _Uy, LS_STAG_CR_PD, _Uz,
		LS_STAG_CR_PD, _P, LS_STAG_CR_PD, _invV, LS_STAG_CR_PD, _nu, double, nu, const int, _N, const int, _M, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_R(double, calcMoment, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceX, LS_STAG_CR_PD, _ForceY, LS_STAG_CR_PC, _curC);
	LS_STAG_DECLARATION_KERNEL_R(double, calcMomentX, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceY, LS_STAG_CR_PD, _ForceZ, LS_STAG_CR_PC, _curC);
	LS_STAG_DECLARATION_KERNEL_R(double, calcMomentY, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PPI, _OnAirfoil,
		LS_STAG_CR_PD, _ForceX, LS_STAG_CR_PD, _ForceZ, LS_STAG_CR_PC, _curC);
/*!
\brief Класс, реализующий расчет аэродинамических нагрузок

Данный класс реализует расчет значений компонент аэродинамической силы в усеченных ячейках и стационарных
аэродинамических коэффициентов (лобового сопротивления, подъемной силы, момента)
*/
	struct LoadsCalculator
	{
		LoadsCalculator():
		countCutCells(0),
		C_d(0.0),
		C_y(0.0),
		C_z(0.0),// коэф-т боковой силы
		moment_aer(0.0),
		moment_x_aer(0.0),
		moment_y_aer(0.0),
		F_aer_x(0.0), 
		F_aer_y(0.0),
		F_aer_z(0.0),
		fileExist(false),
		m_nu(NULL),
		m_pNu(NULL)
		{
			;
		}
		~LoadsCalculator()
		{
			if(fileExist)
			{
				cdFile.close();
				cyFile.close();
				czFile.close();
				mzFile.close();
				mxFile.close();
				myFile.close();
			}
		};
		LS_STAG_KERNEL(calc, LS_STAG_CR_PCT _type, LS_STAG_CR_PPI _OnAirfoil, LS_STAG_R_PD _ForceX, LS_STAG_R_PD _ForceY, 
			LS_STAG_R_PD _ForceZ, LS_STAG_CR_PC _Vib, LS_STAG_R_PC _VibYX, LS_STAG_R_PC _VibZY, LS_STAG_CR_PC _delta, 
			LS_STAG_R_PC _inv_delta_ib_ne, LS_STAG_R_PC _inv_delta_ib_tn, LS_STAG_CR_PC _delta_ib_max, 
			LS_STAG_CR_PC _theta, LS_STAG_CR_PD _Ux, LS_STAG_CR_PD _Uy, LS_STAG_CR_PD _Uz, LS_STAG_CR_PD _P,
			LS_STAG_CR_PD _invV, LS_STAG_CR_PD _nu, double nu, const int _N, const int _M, const int _NM)
		{
			const tripletInt* curPoint = _OnAirfoil + i;
			const int num = curPoint->z;
			const int numX = curPoint->x;
			const int numY = curPoint->y; 
			const double VV = _invV[num];
			const double curP = _P[num];
			const coordinate* curD = _delta + num;
			const coordinate* curVib = _Vib + num;
			const coordinate* curDibMax = _delta_ib_max + num;
			if (_nu) nu = _nu[num];
			// Нормальные компоненты силы: (n ΔS) (p - ν du/dx) и т.д.
    // n_x ΔS = (θ^u_{i-1} - θ^u_i) Δy Δz   (аналогично для y, z)
    const coordinate* curT = _theta + num;
    double dx = curD->x, dy = curD->y, dz = curD->z;
	double tx = curT->x, ty = curT->y, tz = curT->z;
	double txx = (curT - 1)->x, tyy = (curT - _N)->y, tzz = (curT - _NM)->z;

    double nx = txx - tx;   // θ^u_{i-1} - θ^u_i
    double ny = tyy - ty;  // θ^v_{j-1} - θ^v_j
    double nz = tzz - tz; // θ^w_{k-1} - θ^w_k

    double area_x = dy * dz;               // Δy Δz
    double area_y = dx * dz;               // Δx Δz
    double area_z = dx * dy;               // Δx Δy

    // Производные скоростей (нормальные) по формулам (normx3d)-(normz3d)
    double dux = (  tx  * _Ux[numX] - txx * _Ux[numX-1] + nx * curVib->x ) * VV * area_x;   // du/dx
    double dvy = (  ty  * _Uy[numY] - tyy * _Uy[numY-_N] + ny * curVib->y ) * VV * area_y;   // dv/dy
    double dwz = (  tz  * _Uz[num] - tzz * _Uz[num-_NM] + nz * curVib->z ) * VV * area_z;   // dw/dz

    // Базовые вклады от давления и нормальных напряжений
    double baseX = nx * area_x * (curP - nu * dux);
    double baseY = ny * area_y * (curP - nu * dvy);
    double baseZ = nz * area_z * (curP - nu * dwz);

    // Касательные Quad члены – инициализируем нулями
    double quad_uy = 0.0, quad_uz = 0.0;   // для силы X
    double quad_vx = 0.0, quad_vz = 0.0;   // для силы Y
    double quad_wx = 0.0, quad_wy = 0.0;   // для силы Z

			switch (_type[num])
			{// ====================== Тетраэдры (8 типов) ======================
case TNW_Tetrahedron: // В-С-З
{
	double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
	double curU = _Ux[numX], curV = _Uy[numY-_N], curW = _Uz[num - _NM]; 
	coordinate* curVibYX = _VibYX + num; 
	coordinate* curVibZY = _VibZY + num;
	coordinate* curDib = _inv_delta_ib_ne + num;
	quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
	quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    curDib -= _N;
	quad_vx = S_yz * ((curVibYX - _N)->y - curV) * curDib->x;
	curDib = _inv_delta_ib_tn + num;
	quad_vz = S_xy * ((curVibZY - _N)->y - curV) * curDib->z;
    curDib -= _NM;
	quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
	quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
    break;
case TNE_Tetrahedron: // В-С-В
{   // Обратные разности для uy, uz, vx, vz, wx, wy
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY-_N], curW = _Uz[num - _NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // ∂u/∂y|^{ib}_{i-1,j,k} 
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * (curDib-1)->y;
    // ∂u/∂z|^{ib}_{i-1,j,k} 
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * (curDib-1)->z;
    // ∂v/∂x|^{ib}_{i,j-1,k} 
    quad_vx = S_yz * ((curVibYX-_N)->y - curV) * (curDib-_N)->x;
    curDib = _inv_delta_ib_tn + num;
    // ∂v/∂z|^{ib}_{i,j-1,k} 
    quad_vz = S_xy * ((curVibZY-_N)->y - curV) * (curDib-_N)->z;
	curDib -= _NM;
    // ∂w/∂x|^{ib}_{i,j,k-1} 
    quad_wx = S_yz * ((curVibYX-_NM)->z - curW) * curDib->x;
    // ∂w/∂y|^{ib}_{i,j,k-1} 
    quad_wy = S_xz * ((curVibZY-_NM)->z - curW) * curDib->y;
}
    break;
	case TSW_Tetrahedron: // В-Ю-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY], curW = _Uz[num - _NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // uy(i,j,k) прямая, uz(i,j,k) прямая
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // vx(i-1,j,k) обратная
    quad_vx = S_yz * (curVibYX->y - curV) * curDib->x;
    curDib = _inv_delta_ib_tn + num;
    // vz(i,j,k) прямая
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // wx(i-1,j,k-1) обратная
	curDib -= _NM;
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // wy(i,j-1,k-1) обратная
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
    break;
}
case TSE_Tetrahedron: // В-Ю-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY], curW = _Uz[num - _NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + (num - 1);
    // uy(i-1,j,k), uz(i-1,j,k) 
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * curDib->y;
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * curDib->z;
    // vx(i,j,k) прямая
    quad_vx = S_yz * (curVibYX->y - curV) * (curDib+1)->x;
    curDib = _inv_delta_ib_tn + num;
    // vz(i,j,k) прямая
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // wx(i,j,k-1) прямая
	curDib -= _NM;
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // wy(i,j-1,k-1) обратная
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
    break;
case BNW_Tetrahedron: // Н-С-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY-_N], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;

    // uy (i,j-1,k) обратная → (u - u_ib)
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // uz (i,j,k) прямая → (u_ib - u)
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // vx (i-1,j,k) обратная → (v - v_ib)
    quad_vx = S_yz * ((curVibYX - _N)->y - curV) * (curDib - _N)->x;
    curDib = _inv_delta_ib_tn + num;
    // vz (i,j,k-1) обратная
    quad_vz = S_xy * ((curVibZY - _N)->y - curV) * (curDib - _N)->z;
    // wx (i-1,j,k) обратная → (w - w_ib)
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // wy (i,j-1,k) обратная → (w - w_ib)
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
    break;
case BNE_Tetrahedron: // Н-С-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY-_N], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // uy (i,j,k) прямая
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * (curDib-1)->y;
    // uz (i,j,k) прямая
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * (curDib-1)->z;
    // vx (i,j-1,k) 
    quad_vx = S_yz * ((curVibYX-_N)->y - curV) * (curDib-_N)->x;
    curDib = _inv_delta_ib_tn + num;
    // vz (i,j,k-1)
    quad_vz = S_xy * ((curVibZY-_N)->y - curV) * (curDib-_N)->z;
    // wx (i,j,k) прямая
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // wy (i,j,k) прямая
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
    break;
case BSW_Tetrahedron: // Н-Ю-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // uy (i,j-1,k) обратная
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // uz (i,j,k) прямая
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // vx (i,j,k) прямая
    quad_vx = S_yz * (curVibYX->y - curV) * curDib->x;
    curDib = _inv_delta_ib_tn + num;
    // vz (i,j,k-1) обратная 
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // wx (i-1,j,k) обратная
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // wy (i,j-1,k) обратная
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
    break;
case BSE_Tetrahedron: // Н-Ю-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + (num - 1);
    // uy (i,j,k) прямая
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * curDib->y;
    // uz (i,j,k) прямая
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * curDib->z;
    // vx (i,j-1,k) прямая со сдвигом
    quad_vx = S_yz * (curVibYX->y - curV) * (curDib+1)->x;
    curDib = _inv_delta_ib_tn + num;
    // vz (i,j,k-1) обратная
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // wx (i,j,k) прямая
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // wy (i,j,k) прямая
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
    break;
	case N_Hexahedron4: // С 4-угольный гексаэдр
case TNNW_Hexahedron5:
case TNNE_Hexahedron5:
case BNNW_Hexahedron5:
case BNNE_Hexahedron5:
{
    double S_xz = curDibMax->y;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // ∂u/∂y|^{ib}_{i-1,j,k} + ∂u/∂y|^{ib}_{i,j,k}
    quad_uy = S_xz * (((curVib - 1)->x - _Ux[numX - 1]) * (curDib - 1)->y + (curVib->x - _Ux[numX]) * curDib->y);
    // ∂w/∂y|^{ib}_{i,j,k-1} +  ∂w/∂y|^{ib}_{i,j,k}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * (((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y + (curVib->z - _Uz[num]) * curDib->y);
}
break;
case S_Hexahedron4: // Ю 4-угольный гексаэдр
case TSSW_Hexahedron5:
case TSSE_Hexahedron5:
case BSSW_Hexahedron5:
case BSSE_Hexahedron5:
{// QuadIbUY: ΔS_xz * (∂u/∂y|^{ib}_{i-1,j-1,k} + ∂u/∂y|^{ib}_{i,j-1,k}) // QuadIbWY: ΔS_xz * (∂w/∂y|^{ib}_{i,j-1,k-1} + ∂w/∂y|^{ib}_{i,j-1,k})
    double S_xz = curDibMax->y;
	coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    quad_uy = S_xz * (((curVib - 1)->x - _Ux[numX - 1]) * (curDib - 1)->y + (curVib->x - _Ux[numX]) * curDib->y);
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * (((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y + (curVib->z - _Uz[num]) * curDib->y);
}
break;
case E_Hexahedron4: // В 4-угольный гексаэдр
case TENE_Hexahedron5:
case TESE_Hexahedron5:
case BENE_Hexahedron5:
case BESE_Hexahedron5:
{// QuadIbVX: ΔS_yz * (∂v/∂x|^{ib}_{i,j-1,k} + ∂v/∂x|^{ib}_{i,j,k})  // QuadIbWX: ΔS_yz * (∂w/∂x|^{ib}_{i,j,k-1} + ∂w/∂x|^{ib}_{i,j,k})
    double S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    quad_vx = S_yz * (((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->x + (curVib->y - _Uy[numY]) * curDib->x);
    curDib = _inv_delta_ib_tn + num;
    quad_wx = S_yz * (((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->x + (curVib->z - _Uz[num]) * curDib->x);
}
break;
case W_Hexahedron4: // З 4-угольный гексаэдр
case TWNW_Hexahedron5:
case TWSW_Hexahedron5:
case BWNW_Hexahedron5:
case BWSW_Hexahedron5:
{// QuadIbVX: ΔS_yz * (∂v/∂x|^{ib}_{i-1,j-1,k} + ∂v/∂x|^{ib}_{i-1,j,k})  // QuadIbWX: ΔS_yz * (∂w/∂x|^{ib}_{i-1,j,k-1} + ∂w/∂x|^{ib}_{i-1,j,k})
    double S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    quad_vx = S_yz * (((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->x + (curVib->y - _Uy[numY] ) * curDib->x);
    curDib = _inv_delta_ib_tn + num;
    quad_wx = S_yz * (((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->x + (curVib->z - _Uz[num]) * curDib->x);
}
break;
case T_Hexahedron4: // Верх 4-угольный гексаэдр
case TNTW_Hexahedron5:
case TNTE_Hexahedron5:
case TSTW_Hexahedron5:
case TSTE_Hexahedron5:
{// QuadIbVZ: ΔS_xy * (∂v/∂z|^{ib}_{i,j-1,k} + ∂v/∂z|^{ib}_{i,j,k}) // (также полагаем наличие QuadIbUZ для ∂u/∂z)
    double S_xy = curDibMax->x;
    coordinate* curVib = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->z + (curVib->y - _Uy[numY]) * curDib->z);
    quad_uz = S_xy * (((curVib - 1)->x - _Ux[numX - 1]) * (curDib - 1)->z + (curVib->x - _Ux[numX]) * curDib->z);
}
break;
case B_Hexahedron4: // Низ 4-угольный гексаэдр
case BNBW_Hexahedron5:
case BNBE_Hexahedron5:
case BSBW_Hexahedron5:
case BSBE_Hexahedron5:
{   // QuadIbVZ: ΔS_xy * (∂v/∂z|^{ib}_{i,j-1,k-1} + ∂v/∂z|^{ib}_{i,j,k-1})
    double S_xy = curDibMax->x;
    coordinate* curVib = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->z + (curVib->y - _Uy[numY]) * curDib->z);
    quad_uz = S_xy * (((curVib - 1)->x - _Ux[numX - 1]) * (curDib - 1)->z + (curVib->x - _Ux[numX]) * curDib->z);
}
break;// ==================== Пентаэдры (12 типов) ====================
case TN_Pentahedron:
{
    double S_xz = curDibMax->y, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // W + E: ∂u/∂y|^{ib}_{i-1,j,k} + ∂u/∂y|^{ib}_{i,j,k}
	double prevU = _Ux[numX - 1], curU = _Ux[numX];
    quad_uy = S_xz * (((curVib - 1)->x - prevU) * (curDib - 1)->y + (curVib->x - curU) * curDib->y);
    // W + E: ∂u/∂z|^{ib}_{i-1,j,k} + ∂u/∂z|^{ib}_{i,j,k}
	curVib = _VibZY + num; 
    quad_uz = S_xy * (((curVib - 1)->x - prevU) * (curDib - 1)->z + (curVib->x - curU) * curDib->z);
    // T: ∂w/∂y|^{ib}_{i,j,k}
    curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y;
    // N: ∂v/∂z|^{ib}_{i,j,k}
    quad_vz = S_xy * ((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->z;
}
break;
case TS_Pentahedron:
{
    double S_xz = curDibMax->y, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // W+E: ∂u/∂y|^{ib}_{i-1,j,k} + ∂u/∂y|^{ib}_{i,j,k} 
	double prevU = _Ux[numX - 1], curU = _Ux[numX];
    quad_uy = S_xz * (((curVib - 1)->x - prevU) * (curDib - 1)->y + (curVib->x - curU) * curDib->y);
    // W+E: ∂u/∂z|^{ib}_{i-1,j,k} + ∂u/∂z|^{ib}_{i,j,k}
    curVib = _VibZY + num;
    quad_uz = S_xy * (((curVib - 1)->x - prevU) * (curDib - 1)->z + (curVib->x - curU) * curDib->z);
    // B: ∂w/∂y|^{ib}_{i,j-1,k-1} 
    curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * (((curVib-_NM)->z - _Uz[num-_NM]) * (curDib-_NM)->y);
    // N: ∂v/∂z|^{ib}_{i,j,k}  (прямая одиночная)
    quad_vz = S_xy * ((curVib->y - _Uy[numY]) * curDib->z);
}
break;
case BN_Pentahedron: // Н-С (низ-север)
{
    double S_xz = curDibMax->y, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num; 
    coordinate* curDib = _inv_delta_ib_ne + num;
    // ∂u/∂y|^{ib}_{i-1,j-1,k} + ∂u/∂y|^{ib}_{i,j-1,k}
	double prevU = _Ux[numX - 1], curU = _Ux[numX];
    quad_uy = S_xz * (((curVib - 1)->x - prevU) * (curDib - 1)->y + (curVib->x - curU) * curDib->y);
    // ∂u/∂z|^{ib}_{i-1,j,k-1} + ∂u/∂z|^{ib}_{i,j,k-1}
    curVib = _VibZY + num;
    quad_uz = S_xy * (((curVib - 1)->x - prevU) * (curDib - 1)->z + (curVib->x - curU) * curDib->z);
    // ∂w/∂y|^{ib}_{i,j,k}  
    curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * ((curVib->z - _Uz[num]) * curDib->y);
    // S: ∂v/∂z|^{ib}_{i,j-1,k-1}  (обратная одиночная)
    quad_vz = S_xy * (((curVib-_N)->y - _Uy[numY-_N]) * (curDib-_N)->z);
}
break;
case BS_Pentahedron: // Н-Ю (низ-юг)
{
    double S_xz = curDibMax->y, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // ∂u/∂y|^{ib}_{i-1,j-1,k} + ∂u/∂y|^{ib}_{i,j-1,k}
	double prevU = _Ux[numX - 1], curU = _Ux[numX];
    quad_uy = S_xz * (((curVib - 1)->x - prevU) * (curDib - 1)->y + (curVib->x - curU) * curDib->y);
    // ∂u/∂z|^{ib}_{i-1,j,k-1} + ∂u/∂z|^{ib}_{i,j,k-1}
    curVib = _VibZY + num; 
    quad_uz = S_xy * (((curVib - 1)->x - prevU) * (curDib - 1)->z + (curVib->x - curU) * curDib->z);
    // S: ∂w/∂y|^{ib}_{i,j-1,k} 
    curDib = _inv_delta_ib_tn + num;
    quad_wy = S_xz * ((curVib->z - _Uz[num]) * curDib->y);
    // N: ∂v/∂z|^{ib}_{i,j,k-1}  (обратная одиночная)
    quad_vz = S_xy * ((curVib->y - _Uy[numY]) * curDib->z);
}
break;
case TE_Pentahedron: // В-В (верх-восток)
{
    double S_yz = curDibMax->z, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDibEN = _inv_delta_ib_ne + num;
    coordinate* curDibTN = _inv_delta_ib_tn + num;
    // E: ∂w/∂x|^{ib}_{i,j,k-1}
    quad_wx = S_yz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDibTN - _NM)->x;
    // E: ∂v/∂x|^{ib}_{i,j-1,k} + ∂v/∂x|^{ib}_{i,j,k}  (прямые суммы)
    double prevV = _Uy[numY - _N], curV  = _Uy[numY];
    quad_vx = S_yz * (((curVib - _N)->y - prevV) * (curDibEN - _N)->x + (curVib->y - curV) * curDibEN->x);
    // T: ∂v/∂z|^{ib}_{i,j-1,k} +  ∂v/∂z|^{ib}_{i,j,k}
    curVib = _VibZY + num; 
    quad_vz = S_xy * (((curVib - _N)->y - prevV) * (curDibTN - _N)->z + (curVib->y - curV) * curDibTN->z);
    // T: ∂u/∂z|^{ib}_{i,j,k}
    quad_uz = S_xy * (((curVib-1)->x - _Ux[numX-1]) * (curDibEN-1)->z);
}
break;
case TW_Pentahedron: // В-З (верх-запад)
{
    double S_yz = curDibMax->z, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDibEN = _inv_delta_ib_ne + num;
    coordinate* curDibTN = _inv_delta_ib_tn + num;
    // E: ∂w/∂x|^{ib}_{i-1,j,k-1}
    quad_wx = S_yz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDibTN - _NM)->x;
    // E: ∂v/∂x|^{ib}_{i-1,j-1,k} + ∂v/∂x|^{ib}_{i-1,j,k}  
    double prevV = _Uy[numY - _N], curV  = _Uy[numY];
    quad_vx = S_yz * (((curVib - _N)->y - prevV) * (curDibEN - _N)->x + (curVib->y - curV) * curDibEN->x);
    // T: ∂v/∂z|^{ib}_{i,j-1,k} +  ∂v/∂z|^{ib}_{i,j,k}
    curVib = _VibZY + num; 
    quad_vz = S_xy * (((curVib - _N)->y - prevV) * (curDibTN - _N)->z + (curVib->y - curV) * curDibTN->z);
    // T: ∂u/∂z|^{ib}_{i-1,j,k}
    quad_uz = S_xy * ((curVib->x - _Ux[numX]) * curDibEN->z);
}
break;
case BE_Pentahedron: // Н-В (низ-восток)
{
    double S_yz = curDibMax->z, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDibEN = _inv_delta_ib_ne + num;
    coordinate* curDibTN = _inv_delta_ib_tn + num;
    // E: ∂w/∂x|^{ib}_{i,j,k}
    quad_wx = S_yz * (curVib->z - _Uz[num]) * curDibTN->x;
    // E: ∂v/∂x|^{ib}_{i,j-1,k} + ∂v/∂x|^{ib}_{i,j,k}  (прямые суммы)
    double prevV = _Uy[numY - _N], curV  = _Uy[numY];
    quad_vx = S_yz * (((curVib - _N)->y - prevV) * (curDibEN - _N)->x + (curVib->y - curV) * curDibEN->x);
    // T: ∂v/∂z|^{ib}_{i,j-1,k-1} +  ∂v/∂z|^{ib}_{i,j,k-1}
    curVib = _VibZY + num; 
    quad_vz = S_xy * (((curVib - _N)->y - prevV) * (curDibTN - _N)->z + (curVib->y - curV) * curDibTN->z);
    // T: ∂u/∂z|^{ib}_{i-1,j,k-1}
    quad_uz = S_xy * (((curVib-1)->x - _Ux[numX-1]) * (curDibEN-1)->z);
}
break;
case BW_Pentahedron: // Н-З (низ-запад)
{
    double S_yz = curDibMax->z, S_xy = curDibMax->x;
    coordinate* curVib = _VibYX + num;
    coordinate* curDibEN = _inv_delta_ib_ne + num;
    coordinate* curDibTN = _inv_delta_ib_tn + num;
    // E: ∂w/∂x|^{ib}_{i-1,j,k}
    quad_wx = S_yz * (curVib->z - _Uz[num]) * curDibTN->x;
    // E: ∂v/∂x|^{ib}_{i-1,j-1,k} + ∂v/∂x|^{ib}_{i-1,j,k}  
    double prevV = _Uy[numY - _N], curV  = _Uy[numY];
    quad_vx = S_yz * (((curVib - _N)->y - prevV) * (curDibEN - _N)->x + (curVib->y - curV) * curDibEN->x);
    // T: ∂v/∂z|^{ib}_{i,j-1,k-1} +  ∂v/∂z|^{ib}_{i,j,k-1}
    curVib = _VibZY + num; 
    quad_vz = S_xy * (((curVib - _N)->y - prevV) * (curDibTN - _N)->z + (curVib->y - curV) * curDibTN->z);
    // T: ∂u/∂z|^{ib}_{i,j,k-1}
    quad_uz = S_xy * ((curVib->x - _Ux[numX]) * curDibEN->z);
}
break;
case NW_Pentahedron: // С-З (север-запад)
{
    double S_xz = curDibMax->y, S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // N: ∂u/∂y|^{ib}_{i,j,k}  
    quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // W: ∂v/∂x|^{ib}_{i-1,j-1,k} 
    quad_vx = S_yz * ((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->x;
    // W: ∂w/∂x|^{ib}_{i-1,j,k-1} + ∂w/∂x|^{ib}_{i-1,j,k}  (обратные суммы)
    curDib = _inv_delta_ib_tn + num;
    double prevW = _Uz[num - _NM], curW = _Uz[num];
    quad_wx = S_yz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->x + (curVib->z - curW) * curDib->x);
    // N: ∂w/∂y|^{ib}_{i,j,k-1} + ∂w/∂y|^{ib}_{i,j,k}  (прямые суммы)
    curVib = _VibZY + num;
    quad_wy = S_xz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->y + (curVib->z - curW) * curDib->y);
}
break;
case NE_Pentahedron: // С-В (север-восток)
{
    double S_xz = curDibMax->y, S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // N: ∂u/∂y|^{ib}_{i-1,j,k}  
    quad_uy = S_xz * ((curVib-1)->x - _Ux[numX-1]) * (curDib-1)->y;
    // W: ∂v/∂x|^{ib}_{i,j-1,k} 
    quad_vx = S_yz * ((curVib - _N)->y - _Uy[numY - _N]) * (curDib - _N)->x;
    // W: ∂w/∂x|^{ib}_{i,j,k-1} + ∂w/∂x|^{ib}_{i,j,k}  (обратные суммы)
    curDib = _inv_delta_ib_tn + num;
    double prevW = _Uz[num - _NM], curW = _Uz[num];
    quad_wx = S_yz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->x + (curVib->z - curW) * curDib->x);
    // N: ∂w/∂y|^{ib}_{i,j,k-1} + ∂w/∂y|^{ib}_{i,j,k}  (прямые суммы)
    curVib = _VibZY + num;
    quad_wy = S_xz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->y + (curVib->z - curW) * curDib->y);
}
break;
case SW_Pentahedron: // Ю-З (юг-запад)
{
    double S_xz = curDibMax->y, S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // N: ∂u/∂y|^{ib}_{i,j-1,k}  
    quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // W: ∂v/∂x|^{ib}_{i-1,j,k} 
    quad_vx = S_yz * (curVib->y - _Uy[numY]) * curDib->x;
    // W: ∂w/∂x|^{ib}_{i-1,j,k-1} + ∂w/∂x|^{ib}_{i-1,j,k}  (обратные суммы)
    curDib = _inv_delta_ib_tn + num;
    double prevW = _Uz[num - _NM], curW = _Uz[num];
    quad_wx = S_yz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->x + (curVib->z - curW) * curDib->x);
    // N: ∂w/∂y|^{ib}_{i,j-1,k-1} + ∂w/∂y|^{ib}_{i,j-1,k}  (прямые суммы)
    curVib = _VibZY + num;
    quad_wy = S_xz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->y + (curVib->z - curW) * curDib->y);
}
break;
case SE_Pentahedron: // Ю-В (юг-восток)
{
    double S_xz = curDibMax->y, S_yz = curDibMax->z;
    coordinate* curVib = _VibYX + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // N: ∂u/∂y|^{ib}_{i-1,j-1,k}  
    quad_uy = S_xz * ((curVib-1)->x - _Ux[numX-1]) * (curDib-1)->y;
    // W: ∂v/∂x|^{ib}_{i,j,k} 
    quad_vx = S_yz * (curVib->y - _Uy[numY]) * curDib->x;
    // W: ∂w/∂x|^{ib}_{i,j,k-1} + ∂w/∂x|^{ib}_{i,j,k}  (обратные суммы)
    curDib = _inv_delta_ib_tn + num;
    double prevW = _Uz[num - _NM], curW = _Uz[num];
    quad_wx = S_yz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->x + (curVib->z - curW) * curDib->x);
    // N: ∂w/∂y|^{ib}_{i,j,k-1} + ∂w/∂y|^{ib}_{i,j,k}  (прямые суммы)
    curVib = _VibZY + num;
    quad_wy = S_xz * (((curVib - _NM)->z - prevW) * (curDib - _NM)->y + (curVib->z - curW) * curDib->y);
}
break;// ==================== 6-угольные гептаэдры ====================
case TNW_Heptahedron6: // В-С-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_ne + (num-1);
    // UY: "В-Ю-В/В-С-В тетраэдр; гептаэдры 6-уг. В-С-З/Н-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * curDib->y;
    // UZ: "В-Ю-В/В-С-В тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k}
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * curDib->z;
    // VX: "В-Ю-З/Н-Ю-З тетраэдр; гептаэдры 6-уг. В-С-З/Н-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
    quad_vx = S_yz * (curVibYX->y - curV) * (curDib+1)->x;
    // VZ: "В-Ю-З/В-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З" → ∂v/∂z|^{ib}_{i,j,k}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // WX: "Н-С-З/Н-Ю-З тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k}
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // WY: "Н-С-В/Н-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З" → ∂w/∂y|^{ib}_{i,j,k}
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
break;
case TNE_Heptahedron6: // В-С-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "В-Ю-З/В-С-З тетраэдр; гептаэдры 6-уг. В-С-В/Н-С-В" → ∂u/∂y|^{ib}_{i,j,k}
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // UZ: "В-Ю-З/В-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k}
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // VX: "В-Ю-В/Н-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/Н-С-В" → ∂v/∂x|^{ib}_{i,j,k}
    quad_vx = S_yz * (curVibYX->y - curV) * curDib->x;
    // VZ: "В-Ю-З/В-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З" → ∂v/∂z|^{ib}_{i,j,k}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    // WX: "Н-С-В/Н-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k}
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // WY: "Н-С-В/Н-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З" → ∂w/∂y|^{ib}_{i,j,k}
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
break;
case TSE_Heptahedron6: // В-Ю-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY-_N], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "Н-Ю-З/Н-С-З тетраэдр; гептаэдры 6-уг. В-Ю-В/Н-Ю-В" → ∂u/∂y|^{ib}_{i,j-1,k}
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // UZ: "В-Ю-З/В-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k}
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // VX: "В-С-В/Н-С-В тетраэдр; гептаэдры 6-уг. В-Ю-В/Н-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
    quad_vx = S_yz * ((curVibYX - _N)->y - curV) * (curDib - _N)->x;
    // VZ: "В-С-В/В-С-З тетраэдр; гептаэдры 6-уг. В-Ю-З/В-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - curV) * (curDib - _N)->z;
    // WX: "Н-С-В/Н-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k}
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // WY: "В-С-В/В-С-З тетраэдр; гептаэдры 6-уг. В-Ю-В/В-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
break;
case BNE_Heptahedron6: // Н-С-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY], curW = _Uz[num-_NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "В-Ю-З/В-С-З тетраэдр; гептаэдры 6-уг. В-С-В/Н-С-В" → ∂u/∂y|^{ib}_{i,j,k}
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // UZ: "Н-Ю-З/Н-С-З тетраэдр; гептаэдры 6-уг. Н-С-В/Н-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // VX: "В-Ю-В/Н-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/Н-С-В" → ∂v/∂x|^{ib}_{i,j,k}
    quad_vx = S_yz * (curVibYX->y - curV) * curDib->x;
    // VZ: "Н-Ю-З/Н-Ю-В тетраэдр; гептаэдры 6-уг. Н-С-В/Н-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    curDib -= _NM;
    // WX: "В-С-В/В-Ю-В тетраэдр; гептаэдры 6-уг. Н-С-В/Н-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // WY: "Н-Ю-В/Н-Ю-З тетраэдр; гептаэдры 6-уг. Н-С-В/Н-С-З" → ∂w/∂y|^{ib}_{i,j-1,k}
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
break;
case TSW_Heptahedron6: // В-Ю-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY-_N], curW = _Uz[num];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "Н-Ю-В/Н-С-В тетраэдр; гептаэдры 6-уг. В-Ю-З/Н-Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * (curDib-1)->y;
    // UZ: "В-Ю-В/В-С-В тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k}
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * (curDib-1)->z;
    // VX: "В-С-З/Н-С-З тетраэдр; гептаэдры 6-уг. В-Ю-З/Н-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
    quad_vx = S_yz * ((curVibYX - _N)->y - curV) * (curDib - _N)->x;
    // VZ: "В-С-В/В-С-З тетраэдр; гептаэдры 6-уг. В-Ю-З/В-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - curV) * (curDib - _N)->z;
    // WX: "Н-С-З/Н-Ю-З тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k}
    quad_wx = S_yz * (curVibYX->z - curW) * curDib->x;
    // WY: "В-С-В/В-С-З тетраэдр; гептаэдры 6-уг. В-Ю-В/В-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
    quad_wy = S_xz * (curVibZY->z - curW) * curDib->y;
}
break;
case BNW_Heptahedron6: // Н-С-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY], curW = _Uz[num-_NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + (num-1);
    // UY: "В-Ю-В/В-С-В тетраэдр; гептаэдры 6-уг. В-С-З/Н-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * curDib->y;
    // UZ: "Н-Ю-В/Н-С-В тетраэдр; гептаэдры 6-уг. Н-С-З/Н-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * curDib->z;
    // VX: "В-Ю-З/Н-Ю-З тетраэдр; гептаэдры 6-уг. В-С-З/Н-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
    quad_vx = S_yz * (curVibYX->y - curV) * (curDib+1)->x;
    // VZ: "Н-Ю-З/Н-Ю-В тетраэдр; гептаэдры 6-уг. Н-С-В/Н-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - curV) * curDib->z;
    curDib -= _NM;
    // WX: "В-С-З/В-Ю-З тетраэдр; гептаэдры 6-уг. Н-С-З/Н-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // WY: "Н-Ю-В/Н-Ю-З тетраэдр; гептаэдры 6-уг. Н-С-В/Н-С-З" → ∂w/∂y|^{ib}_{i,j-1,k}
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
break;
case BSW_Heptahedron6: // Н-Ю-З
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX-1], curV = _Uy[numY-_N], curW = _Uz[num-_NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "Н-Ю-В/Н-С-В тетраэдр; гептаэдры 6-уг. В-Ю-З/Н-Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
    quad_uy = S_xz * ((curVibYX-1)->x - curU) * (curDib-1)->y;
    // UZ: "Н-Ю-В/Н-С-В тетраэдр; гептаэдры 6-уг. Н-С-З/Н-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
    quad_uz = S_xy * ((curVibZY-1)->x - curU) * (curDib-1)->z;
    // VX: "В-С-З/Н-С-З тетраэдр; гептаэдры 6-уг. В-Ю-З/Н-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
    quad_vx = S_yz * ((curVibYX - _N)->y - curV) * (curDib - _N)->x;
    // VZ: "Н-С-В/Н-С-З тетраэдр; гептаэдры 6-уг. Н-Ю-З/Н-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - curV) * (curDib - _N)->z;
	curDib -= _NM;
    // WX: "В-С-З/В-Ю-З тетраэдр; гептаэдры 6-уг. Н-С-З/Н-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // WY: "В-Ю-В/В-Ю-З тетраэдр; гептаэдры 6-уг. Н-Ю-В/Н-Ю-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
break;
case BSE_Heptahedron6: // Н-Ю-В
{
    double S_xy = curDibMax->x, S_xz = curDibMax->y, S_yz = curDibMax->z;
    double curU = _Ux[numX], curV = _Uy[numY-_N], curW = _Uz[num-_NM];
    coordinate* curVibYX = _VibYX + num;
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib   = _inv_delta_ib_ne + num;
    // UY: "Н-Ю-З/Н-С-З тетраэдр; гептаэдры 6-уг. В-Ю-В/Н-Ю-В" → ∂u/∂y|^{ib}_{i,j-1,k}
    quad_uy = S_xz * (curVibYX->x - curU) * curDib->y;
    // UZ: "Н-Ю-З/Н-С-З тетраэдр; гептаэдры 6-уг. Н-С-В/Н-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
    quad_uz = S_xy * (curVibZY->x - curU) * curDib->z;
    // VX: "В-С-В/Н-С-В тетраэдр; гептаэдры 6-уг. В-Ю-В/Н-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
    quad_vx = S_yz * ((curVibYX - _N)->y - curV) * (curDib - _N)->x;
    // VZ: "Н-С-В/Н-С-З тетраэдр; гептаэдры 6-уг. Н-Ю-З/Н-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
    curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - curV) * (curDib - _N)->z;
	curDib -= _NM;
    // WX: "В-С-В/В-Ю-В тетраэдр; гептаэдры 6-уг. Н-С-В/Н-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
    quad_wx = S_yz * ((curVibYX - _NM)->z - curW) * curDib->x;
    // WY: "В-Ю-В/В-Ю-З тетраэдр; гептаэдры 6-уг. Н-Ю-В/Н-Ю-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
    quad_wy = S_xz * ((curVibZY - _NM)->z - curW) * curDib->y;
}
break;// ==================== 5-угольные гептаэдры ====================
case TNTW_Heptahedron5: // В-С-В-З
{
    double S_xy = curDibMax->x, curW = _Uz[num];
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_ne + num - 1;
    // UZ: "В-Ю-В/В-С-В тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З и 5-уг. В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k}
    quad_uz = S_xy * ((curVibZY - 1)->x - _Ux[numX-1]) * curDib->z;
    // VZ: "В-Ю-З/В-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З и 5-уг. В-С-В-В/В-С-В-З/В-С-С-В/В-С-С-З" → ∂v/∂z|^{ib}_{i,j,k}
	curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - _Uy[numY]) * curDib->z;
    // WX: "Н-С-З/Н-Ю-З тетраэдр; гептаэдры 6-уг. В-С-З/В-Ю-З и 5-уг. В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k}
    quad_wx = curDibMax->z * ((_VibYX + num)->z - curW) * curDib->x;
    // WY: "Н-С-В/Н-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-С-З и 5-уг. В-Ю-Ю-В/В-Ю-Ю-З/В-С-В-В/В-С-В-З" → ∂w/∂y|^{ib}_{i,j,k}
    quad_wy = curDibMax->y * (curVibZY->z - curW) * curDib->y;
}
break;
case TNTE_Heptahedron5: // В-С-В-В
{
    double S_xy = curDibMax->x, curW = _Uz[num];
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // UZ: "В-Ю-З/В-С-З тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В и 5-уг. В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k}
    quad_uz = S_xy * (curVibZY->x - _Ux[numX]) * curDib->z;
    // VZ: та же строка, что и TNTW (В-С-В-В) → ∂v/∂z|^{ib}_{i,j,k}
	curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * (curVibZY->y - _Uy[numY]) * curDib->z;
    // WX: "Н-С-В/Н-Ю-В тетраэдр; гептаэдры 6-уг. В-С-В/В-Ю-В и 5-уг. В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k}
    quad_wx = curDibMax->z * ((_VibYX + num)->z - curW) * curDib->x;
    // WY: та же строка, что и TNTW (В-С-В-В) → ∂w/∂y|^{ib}_{i,j,k}
    quad_wy = curDibMax->y * (curVibZY->z - curW) * curDib->y;
}
break;
case TSTW_Heptahedron5: // В-Ю-В-З
{
    double S_xy = curDibMax->x, curW = _Uz[num];
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_ne + num - 1;
    // UZ: в строке с В-Ю-В-З (см. TNTW) → ∂u/∂z|^{ib}_{i-1,j,k}
    quad_uz = S_xy * ((curVibZY - 1)->x - _Ux[numX-1]) * curDib->z;
    // VZ: "В-С-В/В-С-З тетраэдр; гептаэдры 6-уг. В-Ю-З/В-Ю-В и 5-уг. В-Ю-В-В/В-Ю-В-З/В-Ю-Ю-З/В-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k}
	curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib - _N)->z;  // обратная
    // WX: в строке с В-Ю-В-З (см. TNTW) → ∂w/∂x|^{ib}_{i-1,j,k}
    quad_wx = curDibMax->z * ((_VibYX + num)->z - curW) * curDib->x;
    // WY: "Н-Ю-В/Н-Ю-З тетраэдр; гептаэдры 6-уг. Н-С-В/Н-С-З и 5-уг. В-Ю-В-В/В-Ю-В-З/В-С-С-В/В-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k}
    quad_wy = curDibMax->y * (curVibZY->z - curW) * curDib->y;  // обратная
}
break;
case TSTE_Heptahedron5: // В-Ю-В-В
{
    double S_xy = curDibMax->x, curW = _Uz[num];
    coordinate* curVibZY = _VibZY + num;
    coordinate* curDib = _inv_delta_ib_ne + num;
    // UZ: в строке с В-Ю-В-В (см. TNTE) → ∂u/∂z|^{ib}_{i,j,k}
    quad_uz = S_xy * (curVibZY->x - _Ux[numX]) * curDib->z;
    // VZ: в строке с В-Ю-В-В (см. TSTW) → ∂v/∂z|^{ib}_{i,j-1,k}
	curDib = _inv_delta_ib_tn + num;
    quad_vz = S_xy * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib - _N)->z;
    // WX: в строке с В-Ю-В-В (см. TNTE) → ∂w/∂x|^{ib}_{i,j,k}
    quad_wx = curDibMax->z * ((_VibYX + num)->z - curW) * curDib->x;
    // WY: в строке с В-Ю-В-В (см. TSTW) → ∂w/∂y|^{ib}_{i,j-1,k}
    quad_wy = curDibMax->y * (curVibZY->z - curW) * curDib->y;
}
break;
case TNNW_Heptahedron5: // В-С-С-З
    {
        double S_xz = curDibMax->y, curV = _Uy[numY];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-З/В-С-С-В/В-З-Ю-З/Н-З-Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
        quad_uy = S_xz * ((curVib - 1)->x - _Ux[numX-1]) * (curDib - 1)->y;
    // VX: строка "В-С-С-З/Н-С-С-З/В-З-С-З/Н-З-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
        quad_vx = curDibMax->z * (curVib->y - curV) * curDib->x;
		curDib = _inv_delta_ib_tn + num; curVib = _VibZY + num;
    // VZ: строка "В-С-В-В/В-С-В-З/В-С-С-В/В-С-С-З" → ∂v/∂z|^{ib}_{i,j,k}
        quad_vz = curDibMax->x * (curVib->y - curV) * curDib->z;
	// WY:  → ∂w/∂y|^{ib}_{i,j-1,k}
    quad_wy = S_xz * (curVib->z - _Uz[num]) * curDib->y;
    }
    break;
case TWNW_Heptahedron5: // В-З-С-З
    {
        double S_yz = curDibMax->z, curU_prev = _Ux[numX-1];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-З/Н-Ю-Ю-З/В-З-С-З/Н-З-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
        quad_uy = curDibMax->y * ((curVibYX - 1)->x - curU_prev) * (curDib_ne - 1)->y;
    // VX: строка "В-С-С-З/Н-С-С-З/В-З-С-З/Н-З-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
        quad_vx = S_yz * (curVibYX->y - _Uy[numY]) * curDib_ne->x;
    // UZ: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - curU_prev) * (curDib_ne - 1)->z;
    // WX: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k}
        quad_wx = S_yz * (curVibYX->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x;
    }
    break;
case TNNE_Heptahedron5: // В-С-С-В
    {
        double S_xz = curDibMax->y, curV = _Uy[numY];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-В/Н-С-С-В/В-В-Ю-В/Н-В-Ю-В" → ∂u/∂y|^{ib}_{i,j,k}
        quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // VX: строка "В-С-С-В/Н-С-С-В/В-В-С-В/Н-В-С-В" → ∂v/∂x|^{ib}_{i,j,k}
        quad_vx = curDibMax->z * (curVib->y - curV) * curDib->x;
    // VZ: строка "В-С-В-В/В-С-В-З/В-С-С-В/В-С-С-З" → ∂v/∂z|^{ib}_{i,j,k}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * (curVib->y - curV) * curDib->z;
    // WY: строка "В-Ю-В-В/В-Ю-В-З/В-С-С-В/В-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k}
        quad_wy = S_xz * (curVib->z - _Uz[num]) * curDib->y;
    }
    break;
case TENE_Heptahedron5: // В-В-С-В
    {
        double S_yz = curDibMax->z, curU = _Ux[numX];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-В/Н-Ю-Ю-В/В-В-С-В/Н-В-С-В" → ∂u/∂y|^{ib}_{i,j,k}
        quad_uy = curDibMax->y * (curVibYX->x - curU) * curDib_ne->y;
    // UZ: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - curU) * curDib_ne->z;
    // VX: строка "В-С-С-В/Н-С-С-В/В-В-С-В/Н-В-С-В" → ∂v/∂x|^{ib}_{i,j,k}
        quad_vx = S_yz * (curVibYX->y - _Uy[numY]) * curDib_ne->x;
    // WX: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k}
        quad_wx = S_yz * (curVibYX->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x;
    }
    break;
case BNNW_Heptahedron5: // Н-С-С-З
    {
        double S_xz = curDibMax->y, curV = _Uy[numY];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY:  ∂u/∂y|^{ib}_{i-1,j-1,k}
        quad_uy = S_xz * ((curVib-1)->x - _Ux[numX-1]) * (curDib-1)->y;
    // VX: строка "В-С-С-З/Н-С-С-З/В-З-С-З/Н-З-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
        quad_vx = curDibMax->z * (curVib->y - curV) * curDib->x;
    // VZ: строка "Н-С-Н-В/Н-С-Н-З/Н-С-С-В/Н-С-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
	 curVib = _VibZY + num;  curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * (curVib->y - curV) * curDib->z;
    // WY: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-С-С-В/Н-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
        quad_wy = S_xz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y;
    }
    break;
case BNNE_Heptahedron5: // Н-С-С-В
    {
        double S_xz = curDibMax->y, curV = _Uy[numY];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-В/Н-С-С-В/В-В-Ю-В/Н-В-Ю-В" → ∂u/∂y|^{ib}_{i,j,k}
        quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // VX: строка "В-С-С-В/Н-С-С-В/В-В-С-В/Н-В-С-В" → ∂v/∂x|^{ib}_{i,j,k}
        quad_vx = curDibMax->z * (curVib->y - curV) * curDib->x;
    // VZ: строка "Н-С-Н-В/Н-С-Н-З/Н-С-С-В/Н-С-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * (curVib->y - curV) * curDib->z;
    // WY: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-С-С-В/Н-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
        quad_wy = S_xz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y;
    }
    break;
case TESE_Heptahedron5: // В-В-Ю-В
    {
        double S_yz = curDibMax->z, curU = _Ux[numX];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-В/Н-С-С-В/В-В-Ю-В/Н-В-Ю-В" → ∂u/∂y|^{ib}_{i,j-1,k}
        quad_uy = curDibMax->y * (curVibYX->x - curU) * curDib_ne->y;
    // UZ: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - curU) * curDib_ne->z;
    // VX: строка "В-Ю-Ю-В/Н-В-Ю-В/В-В-Ю-В/Н-Ю-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
        quad_vx = S_yz * ((curVibYX-_N)->y - _Uy[numY - _N]) * (curDib_ne-_N)->x;
    // WX: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k}
        quad_wx = S_yz * (curVibYX->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x;
    }
    break;
case TSSE_Heptahedron5: // В-Ю-Ю-В
    {
        double S_xz = curDibMax->y, curV_prev = _Uy[numY - _N];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-В/Н-Ю-Ю-В/В-В-С-В/Н-В-С-В" → ∂u/∂y|^{ib}_{i,j-1,k}
        quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // VX: строка "В-Ю-Ю-В/Н-В-Ю-В/В-В-Ю-В/Н-Ю-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
        quad_vx = curDibMax->z * ((curVib - _N)->y - curV_prev) * (curDib - _N)->x;
    // VZ: строка "В-Ю-В-В/В-Ю-В-З/В-Ю-Ю-З/В-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * ((curVib - _N)->y - curV_prev) * (curDib - _N)->z;
    // WY: строка "В-Ю-Ю-В/В-Ю-Ю-З/В-С-В-В/В-С-В-З" → ∂w/∂y|^{ib}_{i,j,k}
        quad_wy = S_xz * (curVib->z - _Uz[num]) * curDib->y;
    }
    break;
case BENE_Heptahedron5: // Н-В-С-В
    {
        double S_yz = curDibMax->z, curU = _Ux[numX];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-В/Н-Ю-Ю-В/В-В-С-В/Н-В-С-В" → ∂u/∂y|^{ib}_{i,j,k}
        quad_uy = curDibMax->y * (curVibYX->x - curU) * curDib_ne->y;
    // UZ: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - curU) * curDib_ne->z;
    // VX: строка "В-С-С-В/Н-С-С-В/В-В-С-В/Н-В-С-В" → ∂v/∂x|^{ib}_{i,j,k}
        quad_vx = S_yz * (curVibYX->y - _Uy[numY]) * curDib_ne->x;
    // WX: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
        quad_wx = S_yz * ((curVibYX - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x;
    }
    break;
case BESE_Heptahedron5: // Н-В-Ю-В
    {
        double S_yz = curDibMax->z, curU = _Ux[numX];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-В/Н-С-С-В/В-В-Ю-В/Н-В-Ю-В" → ∂u/∂y|^{ib}_{i,j-1,k}
        quad_uy = curDibMax->y * (curVibYX->x - curU) * curDib_ne->y;
    // UZ: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - curU) * curDib_ne->z;
    // VX: строка "В-Ю-Ю-В/Н-В-Ю-В/В-В-Ю-В/Н-Ю-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
        quad_vx = S_yz * ((curVibYX-_N)->y - _Uy[numY - _N]) * (curDib_ne-_N)->x;
    // WX: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
        quad_wx = S_yz * ((curVibYX - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x;
    }
    break;
case BNBE_Heptahedron5: // Н-С-Н-В
    {
        double S_xy = curDibMax->x, curW_km1 = _Uz[num - _NM];
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // UZ: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
        quad_uz = S_xy * (curVibZY->x - _Ux[numX]) * (_inv_delta_ib_ne + num)->z;
    // VZ: строка "Н-С-Н-В/Н-С-Н-З/Н-С-С-В/Н-С-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
        quad_vz = S_xy * (curVibZY->y - _Uy[numY]) * curDib_tn->z;
    // WX: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
	curDib_tn -= _NM;
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - curW_km1) * curDib_tn->x;
    // WY: строка "Н-Ю-Ю-В/Н-С-Н-З/Н-С-Н-В/Н-Ю-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - curW_km1) * curDib_tn->y;
    }
    break;
case TWSW_Heptahedron5: // В-З-Ю-З
    {
        double S_yz = curDibMax->z, curU_prev = _Ux[numX-1];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-З/В-С-С-В/В-З-Ю-З/Н-З-Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
        quad_uy = curDibMax->y * ((curVibYX - 1)->x - curU_prev) * (curDib_ne - 1)->y;
    // UZ: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - curU_prev) * (curDib_ne - 1)->z;
    // VX: строка "В-Ю-Ю-З/Н-З-Ю-З/В-З-Ю-З/Н-Ю-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
        quad_vx = S_yz * ((curVibYX-_N)->y - _Uy[numY-_N]) * (curDib_ne-_N)->x;
    // WX: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k}
        quad_wx = S_yz * (curVibYX->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x;
    }
    break;
case TSSW_Heptahedron5: // В-Ю-Ю-З
    {
        double S_xz = curDibMax->y, curV_prev = _Uy[numY - _N];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-З/Н-Ю-Ю-З/В-З-С-З/Н-З-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
        quad_uy = S_xz * ((curVib - 1)->x - _Ux[numX-1]) * (curDib - 1)->y;
    // VX: строка "В-Ю-Ю-З/Н-З-Ю-З/В-З-Ю-З/Н-Ю-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
        quad_vx = curDibMax->z * ((curVib - _N)->y - curV_prev) * (curDib - _N)->x;
    // VZ: строка "В-Ю-В-В/В-Ю-В-З/В-Ю-Ю-З/В-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * ((curVib - _N)->y - curV_prev) * (curDib - _N)->z;
    // WY: строка "В-Ю-Ю-В/В-Ю-Ю-З/В-С-В-В/В-С-В-З" → ∂w/∂y|^{ib}_{i,j,k}
        quad_wy = S_xz * (curVib->z - _Uz[num]) * curDib->y;
    }
    break;
case BWNW_Heptahedron5: // Н-З-С-З
    {
        double S_yz = curDibMax->z, curU_prev = _Ux[numX-1];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-З/Н-Ю-Ю-З/В-З-С-З/Н-З-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
        quad_uy = curDibMax->y * ((curVibYX - 1)->x - curU_prev) * (curDib_ne - 1)->y;
    // UZ: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - curU_prev) * (curDib_ne - 1)->z;
    // VX: строка "В-С-С-З/Н-С-С-З/В-З-С-З/Н-З-С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
        quad_vx = S_yz * (curVibYX->y - _Uy[numY]) * curDib_ne->x;
    // WX: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
        quad_wx = S_yz * ((curVibYX - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x;
    }
    break;
case BWSW_Heptahedron5: // Н-З-Ю-З
    {
        double S_yz = curDibMax->z, curU_prev = _Ux[numX-1];
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-З/В-С-С-В/В-З-Ю-З/Н-З-Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
        quad_uy = curDibMax->y * ((curVibYX - 1)->x - curU_prev) * (curDib_ne - 1)->y;
    // UZ: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - curU_prev) * (curDib_ne - 1)->z;
    // VX: строка "В-Ю-Ю-З/Н-З-Ю-З/В-З-Ю-З/Н-Ю-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
        quad_vx = S_yz * ((curVibYX-_N)->y - _Uy[numY-_N]) * (curDib_ne-_N)->x;
    // WX: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
        quad_wx = S_yz * ((curVibYX - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x;
    }
    break;
case BNBW_Heptahedron5: // Н-С-Н-З
    {
        double S_xy = curDibMax->x, curW_km1 = _Uz[num - _NM];
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // UZ: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
        quad_uz = S_xy * ((curVibZY - 1)->x - _Ux[numX-1]) * (_inv_delta_ib_ne + num - 1)->z;
    // VZ: строка "Н-С-Н-В/Н-С-Н-З/Н-С-С-В/Н-С-С-З" → ∂v/∂z|^{ib}_{i,j,k-1}
        quad_vz = S_xy * (curVibZY->y - _Uy[numY]) * curDib_tn->z;
    // WX: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - curW_km1) * (curDib_tn - _NM)->x;
    // WY: строка "Н-Ю-Ю-В/Н-С-Н-З/Н-С-Н-В/Н-Ю-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - curW_km1) * (curDib_tn - _NM)->y;
    }
    break;
case BSSE_Heptahedron5: // Н-Ю-Ю-В
    {
        double S_xz = curDibMax->y, curV_prev = _Uy[numY - _N];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-В/Н-Ю-Ю-В/В-В-С-В/Н-В-С-В" → ∂u/∂y|^{ib}_{i,j-1,k}
        quad_uy = S_xz * (curVib->x - _Ux[numX]) * curDib->y;
    // VX: строка "В-Ю-Ю-В/Н-В-Ю-В/В-В-Ю-В/Н-Ю-Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
        quad_vx = curDibMax->z * ((curVib-_N)->y - curV_prev) * (curDib-_N)->x;
    // VZ: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-Ю-Ю-З/Н-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * ((curVib - _N)->y - curV_prev) * (curDib - _N)->z;
    // WY: строка "Н-Ю-Ю-В/Н-С-Н-З/Н-С-Н-В/Н-Ю-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
        quad_wy = S_xz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y;
    }
    break;
case BSSW_Heptahedron5: // Н-Ю-Ю-З
    {
        double S_xz = curDibMax->y, curV_prev = _Uy[numY - _N];
        coordinate* curVib = _VibYX + num;
        coordinate* curDib = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-З/Н-Ю-Ю-З/В-З-С-З/Н-З-С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
        quad_uy = S_xz * ((curVib - 1)->x - _Ux[numX-1]) * (curDib - 1)->y;
    // VX: строка "В-Ю-Ю-З/Н-З-Ю-З/В-З-Ю-З/Н-Ю-Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
        quad_vx = curDibMax->z * ((curVib - _N)->y - curV_prev) * (curDib - _N)->x;
    // VZ: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-Ю-Ю-З/Н-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
	curVib = _VibZY + num; curDib = _inv_delta_ib_tn + num;
        quad_vz = curDibMax->x * ((curVib - _N)->y - curV_prev) * (curDib - _N)->z;
    // WY: строка "Н-Ю-Ю-В/Н-С-Н-З/Н-С-Н-В/Н-Ю-Ю-З" → ∂w/∂y|^{ib}_{i,j,k-1}
        quad_wy = S_xz * ((curVib - _NM)->z - _Uz[num - _NM]) * (curDib - _NM)->y;
    }
    break;
case BSBW_Heptahedron5: // Н-Ю-Н-З
    {
        double S_xy = curDibMax->x, curW_km1 = _Uz[num - _NM];
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // UZ: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
        quad_uz = S_xy * ((curVibZY - 1)->x - _Ux[numX-1]) * (_inv_delta_ib_ne + num - 1)->z;
    // VZ: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-Ю-Ю-З/Н-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
        quad_vz = S_xy * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib_tn - _N)->z;
    // WX: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
	curDib_tn -=_NM;
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - curW_km1) * curDib_tn->x;
    // WY: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-С-С-В/Н-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - curW_km1) * curDib_tn->y;
    }
    break;
case BSBE_Heptahedron5: // Н-Ю-Н-В
    {
        double S_xy = curDibMax->x, curW_km1 = _Uz[num - _NM];
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // UZ: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂u/∂z|^{ib}_{i,j,k-1}
        quad_uz = S_xy * (curVibZY->x - _Ux[numX]) * curDib_tn->z;
    // VZ: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-Ю-Ю-З/Н-Ю-Ю-В" → ∂v/∂z|^{ib}_{i,j-1,k-1}
        quad_vz = S_xy * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib_tn - _N)->z;
    // WX: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В" → ∂w/∂x|^{ib}_{i,j,k-1}
	curDib_tn -=_NM;
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - curW_km1) * curDib_tn->x;
    // WY: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-С-С-В/Н-С-С-З" → ∂w/∂y|^{ib}_{i,j-1,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - curW_km1) * curDib_tn->y;
    }
    break;
// ==================== 4-угольные гептаэдры (Heptahedron4) ====================
case TN_Heptahedron4: // В-С
    {
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // VZ: строка "В-С-В-В/В-С-В-З/В-С-С-В/В-С-С-З, гептаэдры 4-уг. В-С" → ∂v/∂z|^{ib}_{i,j,k}
        quad_vz = curDibMax->x * (curVibZY->y - _Uy[numY]) * curDib_tn->z;
    // WY: строка "В-Ю-Ю-В/В-Ю-Ю-З/В-С-В-В/В-С-В-З, гептаэдры 4-уг. В-С" → ∂w/∂y|^{ib}_{i,j,k}
        quad_wy = curDibMax->y * (curVibZY->z - _Uz[num]) * curDib_tn->y;
    }
    break;
case TE_Heptahedron4: // В-В
    {// UZ: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В, гептаэдры 4-уг. В-В" → ∂u/∂z|^{ib}_{i,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - _Ux[numX]) * (_inv_delta_ib_ne + num)->z;
    // WX: строка "В-С-В-В/В-Ю-В-В/В-В-С-В/В-В-Ю-В, гептаэдры 4-уг. В-В" → ∂w/∂x|^{ib}_{i,j,k}
        quad_wx = curDibMax->z * ((_VibYX + num)->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x;
    }
    break;
case NE_Heptahedron4: // С-В
    {
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-В/Н-Ю-Ю-В/В-В-С-В/Н-В-С-В, гептаэдры 4-уг. С-В" → ∂u/∂y|^{ib}_{i,j,k}
        quad_uy = curDibMax->y * (curVibYX->x - _Ux[numX]) * curDib_ne->y;
    // VX: строка "В-С-С-В/Н-С-С-В/В-В-С-В/Н-В-С-В, гептаэдры 4-уг. С-В" → ∂v/∂x|^{ib}_{i,j,k}
        quad_vx = curDibMax->z * (curVibYX->y - _Uy[numY]) * curDib_ne->x;
    }
    break;
case TW_Heptahedron4: // В-З
    { // UZ: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З, гептаэдры 4-уг. В-З" → ∂u/∂z|^{ib}_{i-1,j,k}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - _Ux[numX - 1]) * (_inv_delta_ib_ne + num - 1)->z;
    // WX: строка "В-С-В-З/В-Ю-В-З/В-З-С-З/В-З-Ю-З, гептаэдры 4-уг. В-З" → ∂w/∂x|^{ib}_{i-1,j,k}
        quad_wx = curDibMax->z * ((_VibYX + num)->z - _Uz[num]) * (_inv_delta_ib_tn + num)->x; // обратная: w - w_ib
    }
    break;
case NW_Heptahedron4: // С-З
    {
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-Ю-Ю-З/Н-Ю-Ю-З/В-З-С-З/Н-З-С-З, гептаэдры 4-уг. С-З" → ∂u/∂y|^{ib}_{i-1,j,k}
        quad_uy = curDibMax->y * ((curVibYX - 1)->x - _Ux[numX-1]) * (curDib_ne - 1)->y;
    // VX: строка "В-С-С-З/Н-С-С-З/В-З-С-З/Н-З-С-З, гептаэдры 4-уг. С-З" → ∂v/∂x|^{ib}_{i-1,j,k}
        quad_vx = curDibMax->z * (curVibYX->y - _Uy[numY]) * curDib_ne->x; // прямая: v_ib - v
    }
    break;
case TS_Heptahedron4: // В-Ю
    {
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // VZ: строка "В-Ю-В-В/В-Ю-В-З/В-Ю-Ю-З/В-Ю-Ю-В, гептаэдры 4-уг. В-Ю" → ∂v/∂z|^{ib}_{i,j-1,k}
        quad_vz = curDibMax->x * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib_tn - _N)->z; // обратная (v - v_ib) с j-1
    // WY: строка "В-Ю-В-В/В-Ю-В-З/В-С-С-В/В-С-С-З, гептаэдры 4-уг. В-Ю" → ∂w/∂y|^{ib}_{i,j-1,k}
        quad_wy = curDibMax->y * (curVibZY->z - _Uz[num]) * curDib_tn->y; 
    }
    break;
case SW_Heptahedron4: // Ю-З
    {
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-З/В-С-С-В/В-З-Ю-З/Н-З-Ю-З, гептаэдры 4-уг. Ю-З" → ∂u/∂y|^{ib}_{i-1,j-1,k}
        quad_uy = -curDibMax->y * ((curVibYX - 1)->x - _Ux[numX-1]) * (curDib_ne - 1)->y; // EXP 26.09.2026: знак по 2D SW_Pentagon (txx<1 => (theta-1)<0)
    // VX: строка "В-Ю-Ю-З/Н-З-Ю-З/В-З-Ю-З/Н-Ю-Ю-З, гептаэдры 4-уг. Ю-З" → ∂v/∂x|^{ib}_{i-1,j-1,k}
        quad_vx = -curDibMax->z * ((curVibYX-_N)->y - _Uy[numY-_N]) * (curDib_ne-_N)->x; // EXP 26.09.2026: знак по 2D SW_Pentagon (tyy<1)
    }
    break;
case SE_Heptahedron4: // Ю-В
    {
        coordinate* curVibYX = _VibYX + num;
        coordinate* curDib_ne = _inv_delta_ib_ne + num;
    // UY: строка "В-С-С-В/Н-С-С-В/В-В-Ю-В/Н-В-Ю-В, гептаэдры 4-уг. Ю-В" → ∂u/∂y|^{ib}_{i,j-1,k}
        quad_uy = -curDibMax->y * (curVibYX->x - _Ux[numX]) * curDib_ne->y; // EXP 26.09.2026: знак по 2D SE_Pentagon
    // VX: строка "В-Ю-Ю-В/Н-В-Ю-В/В-В-Ю-В/Н-Ю-Ю-В, гептаэдры 4-уг. Ю-В" → ∂v/∂x|^{ib}_{i,j-1,k}
        quad_vx = curDibMax->z * ((curVibYX-_N)->y - _Uy[numY - _N]) * (curDib_ne-_N)->x; 
    }
    break;
case BN_Heptahedron4: // Н-С
    {
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // VZ: строка "Н-С-Н-В/Н-С-Н-З/Н-С-С-В/Н-С-С-З, гептаэдры 4-уг. Н-С" → ∂v/∂z|^{ib}_{i,j,k-1}
        quad_vz = curDibMax->x * (curVibZY->y - _Uy[numY]) * curDib_tn->z;
    // WY: строка "Н-Ю-Ю-В/Н-С-Н-З/Н-С-Н-В/Н-Ю-Ю-З, гептаэдры 4-уг. Н-С" → ∂w/∂y|^{ib}_{i,j,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - _Uz[num - _NM]) * (curDib_tn - _NM)->y;
    }
    break;
case BE_Heptahedron4: // Н-В
    { // UZ: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В, гептаэдры 4-уг. Н-В" → ∂u/∂z|^{ib}_{i,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num)->x - _Ux[numX]) * (_inv_delta_ib_ne + num)->z; // обратная (u - u_ib) k-1
    // WX: строка "Н-С-Н-В/Н-Ю-Н-В/Н-В-С-В/Н-В-Ю-В, гептаэдры 4-уг. Н-В" → ∂w/∂x|^{ib}_{i,j,k-1}
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x; // прямая (w_ib - w) k-1
    }
    break;
case BW_Heptahedron4: // Н-З
    {// UZ: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З, гептаэдры 4-уг. Н-З" → ∂u/∂z|^{ib}_{i-1,j,k-1}
        quad_uz = curDibMax->x * ((_VibZY + num - 1)->x - _Ux[numX-1]) * (_inv_delta_ib_ne + num - 1)->z; // обратная (u - u_ib) i-1,k-1
    // WX: строка "Н-С-Н-З/Н-Ю-Н-З/Н-З-С-З/Н-З-Ю-З, гептаэдры 4-уг. Н-З" → ∂w/∂x|^{ib}_{i-1,j,k-1}
        quad_wx = curDibMax->z * ((_VibYX + num - _NM)->z - _Uz[num - _NM]) * (_inv_delta_ib_tn + num - _NM)->x;
    }
    break;
case BS_Heptahedron4: // Н-Ю
    {
        coordinate* curVibZY = _VibZY + num;
        coordinate* curDib_tn = _inv_delta_ib_tn + num;
    // VZ: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-Ю-Ю-З/Н-Ю-Ю-В, гептаэдры 4-уг. Н-Ю" → ∂v/∂z|^{ib}_{i,j-1,k-1}
        quad_vz = curDibMax->x * ((curVibZY - _N)->y - _Uy[numY - _N]) * (curDib_tn - _N)->z; // обратная (v - v_ib) j-1,k-1
    // WY: строка "Н-Ю-Н-В/Н-Ю-Н-З/Н-С-С-В/Н-С-С-З, гептаэдры 4-уг. Н-Ю" → ∂w/∂y|^{ib}_{i,j-1,k-1}
        quad_wy = curDibMax->y * ((curVibZY - _NM)->z - _Uz[num - _NM]) * (curDib_tn - _NM)->y;
    }
    break;
			}//switch// Собираем силы
    _ForceX[i] = baseX - nu * (quad_uy + quad_uz);
    _ForceY[i] = baseY - nu * (quad_vx + quad_vz);
    _ForceZ[i] = baseZ - nu * (quad_wx + quad_wy);
		}
		LS_STAG_KERNEL_R(double, calcMoment, LS_STAG_CR_PC _nodes, LS_STAG_CR_PPI _OnAirfoil, LS_STAG_CR_PD _ForceX,
			LS_STAG_CR_PD _ForceY, LS_STAG_CR_PC _curC)
		{
			const coordinate* curD = _nodes + _OnAirfoil[i].z;
			return ((curD->x - _curC->x) * _ForceY[i] - (curD->y - _curC->y) * _ForceX[i]);
		}
		LS_STAG_KERNEL_R(double, calcMomentX, LS_STAG_CR_PC _nodes, LS_STAG_CR_PPI _OnAirfoil, LS_STAG_CR_PD _ForceY,
			LS_STAG_CR_PD _ForceZ, LS_STAG_CR_PC _curC)
		{// M_x = (y - y_c) F_z - (z - z_c) F_y
			const coordinate* curD = _nodes + _OnAirfoil[i].z;
			return ((curD->y - _curC->y) * _ForceZ[i] - (curD->z - _curC->z) * _ForceY[i]);
		}
		LS_STAG_KERNEL_R(double, calcMomentY, LS_STAG_CR_PC _nodes, LS_STAG_CR_PPI _OnAirfoil, LS_STAG_CR_PD _ForceX,
			LS_STAG_CR_PD _ForceZ, LS_STAG_CR_PC _curC)
		{// M_y = (z - z_c) F_x - (x - x_c) F_z
			const coordinate* curD = _nodes + _OnAirfoil[i].z;
			return ((curD->z - _curC->z) * _ForceX[i] - (curD->x - _curC->x) * _ForceZ[i]);
		}
		void init()
		{// вызываем каждый раз перед заполнением OnAirfoil (по мере его заполнения увеличиваем счётчик countCutCells)
			m_OnAirfoil.set(0);
			countCutCells = 0;
		};
		void resize()
		{// вызываем каждый раз после переопределения  countCutCells (и заполнения OnAirfoil)
			ForcesX.resize(countCutCells);
			ForcesY.resize(countCutCells);
			ForcesZ.resize(countCutCells);
			m_OnAirfoil.resize(countCutCells);
		};
		void set(LSSTAGUniformMesh* _mesh, const double _nu, const double _t, const double _fcoef, const std::tstring &output, const std::tstring &label)
		{// вызываем 1 раз после определения maxCountCutCells
			momentFilter.init(_t, DEFAULT_CUT_OFF_FREQUENCY, DEFAULT_SUPPRESSION_LEVEL);
			momentXFilter.init(_t, DEFAULT_CUT_OFF_FREQUENCY, DEFAULT_SUPPRESSION_LEVEL);
			momentYFilter.init(_t, DEFAULT_CUT_OFF_FREQUENCY, DEFAULT_SUPPRESSION_LEVEL);
			nu = _nu;
			fcoef = _fcoef;
			m_mesh =_mesh;
			m_nodes = m_mesh->nodes;
			m_curC = m_mesh->cur_C;// если нужно считать момент при нескольких профилях --- изменить
			m_delta = m_mesh->delta;
			m_theta = m_mesh->theta;
			m_Vib = m_mesh->v_ib;
			m_VibZY = m_mesh->v_ib_zy;
			m_VibYX = m_mesh->v_ib_yx;
			m_inv_delta_ib_ne = m_mesh->inv_delta_ib_ne;
			m_inv_delta_ib_tn = m_mesh->inv_delta_ib_tn;
			m_delta_ib_max = m_mesh->delta_ib_max;
			m_invV = m_mesh->baseMesh.inv_VOF;
			m_N = m_mesh->baseMesh.N;
			m_M = m_mesh->baseMesh.M;
			m_NM = m_N * m_M;
			m_type = m_mesh->baseMesh.typeOfCell;
			m_OnAirfoil.init(m_mesh->maxCountCutCells, OnAirfoil);
			ForcesX.init(m_mesh->maxCountCutCells, ForceX);
			ForcesY.init(m_mesh->maxCountCutCells, ForceY);
			ForcesZ.init(m_mesh->maxCountCutCells, ForceZ);
			cdFilename = output + _T("_Cd") + label + _T(".dat");
			cyFilename = output + _T("_Cy") + label + _T(".dat");
			czFilename = output + _T("_Cz") + label + _T(".dat");
			mzFilename = output + _T("_Mz") + label + _T(".dat");
			mxFilename = output + _T("_Mx") + label + _T(".dat");
			myFilename = output + _T("_My") + label + _T(".dat");
			cdFile.open(cdFilename.c_str());
			cyFile.open(cyFilename.c_str());
			czFile.open(czFilename.c_str());
			mzFile.open(mzFilename.c_str());
			mxFile.open(mxFilename.c_str());
			myFile.open(myFilename.c_str());
			fileExist = true;
		};
		void initSol()
		{
			m_P = m_mesh->m_P;
			m_Ux = m_mesh->m_Ux;
			m_Uy = m_mesh->m_Uy;
			m_Uz = m_mesh->m_Uz;
		};
		void initGNF(double *_nu){ m_nu = _nu; };
		void initTDF(double *_nu){ m_pNu = _nu; };
		/*!
		Выводит очередные значених аэродинамических коэффициентов в файлы
		*/
		void print()
		{
			cdFile << C_d << "\n";
			cyFile << C_y << "\n";
			czFile << C_z << "\n";
			mzFile << moment_aer << "\n";
			mxFile << moment_x_aer << "\n";
			myFile << moment_y_aer << "\n";
		};
		/*!
		Рассчитывает значения компонент аэродинамической силы для каждой усеченной ячейки и вычисляет значения
		аэродинамических коэффициентов
		*/
		void calculate()//расчёт нагрузок //new-new version --- 19.07.12.
		{
			if(countCutCells < 1) return;// нет усеченных ячеек --- нечего считать (иначе size=0 -> деление на 0 в host-обёртке)
			if(m_pNu) nu = m_pNu[0];
			LS_STAG_PAR_FOR(LoadsCalculator, calc, countCutCells, m_type, OnAirfoil, ForceX, ForceY, ForceZ, m_Vib, m_VibYX, m_VibZY,
				m_delta, m_inv_delta_ib_ne, m_inv_delta_ib_tn, m_delta_ib_max, m_theta, m_Ux, m_Uy, m_Uz, m_P, m_invV, m_nu, nu, m_N, m_M, m_NM)
			F_aer_x = ForcesX.getSum();
			F_aer_y = ForcesY.getSum();
			F_aer_z = ForcesZ.getSum();
			C_d = F_aer_x * fcoef;
			C_y = F_aer_y * fcoef;
			C_z = F_aer_z * fcoef;
			moment();
		};
	  private:
		void printCutCells()
		{// debug -- without GPU support 
			std::ofstream outCC;
			static int dbgDumpCounter = 0;
			char dbgBuf[64];
			sprintf_s(dbgBuf, "onCutCells_%03d.dat", dbgDumpCounter++);
			std::string nameCC = dbgBuf;
			outCC.open(nameCC.c_str());
			tripletInt *curPoint = OnAirfoil;
			for(int i = 0; i < countCutCells; i++, curPoint++)
			{
				int num = curPoint->z;
				int numY = curPoint->y;
				int numX = curPoint->x;//curT = &m_theta[num];
				coordinate* curVib = m_Vib + num;
				outCC << " № ячейки   " << num << " ; " << numX << std::endl;
				outCC << " TypeOfCell[i]   " << m_type[num] << std::endl;
				outCC << " Force.x   " << ForceX[i] << std::endl;
				outCC << " Force.y   " << ForceY[i] << std::endl;
				outCC << " Force.z   " << ForceZ[i] << std::endl;
				outCC << " v_ib.x   " << curVib->x << std::endl;
				outCC << " v_ib.y   " << curVib->y << std::endl;
				outCC << " v_ib.z   " << curVib->z << std::endl;
				outCC << " v_ibZY.x   " << m_VibZY[num].x << std::endl;
				outCC << " v_ibZY.y   " << m_VibZY[num].y << std::endl;
				outCC << " v_ibZY.z   " << m_VibZY[num].z << std::endl;
				outCC << " v_ibYX.x   " << m_VibYX[num].x << std::endl;
				outCC << " v_ibYX.y   " << m_VibYX[num].y << std::endl;
				outCC << " v_ibYX.z   " << m_VibYX[num].z << std::endl;
				outCC << " pressure   " << m_P[num] << std::endl;
				outCC << " Ux   " << m_Ux[numX] << std::endl;
				outCC << " Uy   " << m_Uy[numY] << std::endl;
				outCC << " Uz   " << m_Uz[num] << std::endl;
				outCC << " Ux-N   " << m_Ux[numX-m_N+1] << std::endl;
				outCC << " Ux-1   " << m_Ux[numX-1] << std::endl;
				outCC << " Ux+1   " << m_Ux[numX+1] << std::endl;
				outCC << " Ux+N   " << m_Ux[numX+m_N-1] << std::endl;
				outCC << " Uy-N   " << m_Ux[num-m_N] << std::endl;
				outCC << " Uy-1   " << m_Ux[num-1] << std::endl;
				outCC << " Uy+1   " << m_Ux[num+1] << std::endl;
				outCC << " Uy+N   " << m_Ux[num+m_N] << std::endl;
			}
			outCC.close();
		};
		void moment()
		{
			{
				LS_STAG_PAR_SUM(LoadsCalculator, calcMoment, double, countCutCells, m_nodes, OnAirfoil, ForceX, ForceY, m_curC)
				moment_aer = momentFilter.getValue(LS_STAG_PAR_ACC_VALUE);
			}
			{
				LS_STAG_PAR_SUM(LoadsCalculator, calcMomentX, double, countCutCells, m_nodes, OnAirfoil, ForceY, ForceZ, m_curC)
				moment_x_aer = momentXFilter.getValue(LS_STAG_PAR_ACC_VALUE);
			}
			{
				LS_STAG_PAR_SUM(LoadsCalculator, calcMomentY, double, countCutCells, m_nodes, OnAirfoil, ForceX, ForceZ, m_curC)
				moment_y_aer = momentYFilter.getValue(LS_STAG_PAR_ACC_VALUE);
			}
		};
		ButterworthLowPassFilter momentFilter;///< ФНЧ Баттерворта для сглаживания значений аэродинамического момента
		ButterworthLowPassFilter momentXFilter;///< ФНЧ Баттерворта для сглаживания значений момента M_x
		ButterworthLowPassFilter momentYFilter;///< ФНЧ Баттерворта для сглаживания значений момента M_y
	  public:
		CoordMeshFunction<int> m_OnAirfoil;
		extDouble ForcesX, ForcesY, ForcesZ;///< Сеточные функции со значениями компонент аэродинамической силы на усеченных ячейках
	  private:
		std::tstring cdFilename, cyFilename, czFilename, mzFilename, mxFilename, myFilename;
		std::ofstream cdFile, cyFile, czFile, mzFile, mxFile, myFile;///< Потоки для вывода значений аэродинамических коэффициентов в файлы
		double nu, fcoef;
	  public:
		double moment_aer, moment_x_aer, moment_y_aer, F_aer_x, F_aer_y, F_aer_z, C_d, C_y, C_z;///< Текущие значения аэродинамических коэффициентов
		int countCutCells;///< Количество усеченных ячеек
		tripletInt *OnAirfoil;
		double *m_nu;// for GNF
		double *m_pNu;// for TDF
		double *ForceX, *ForceY, *ForceZ;
	  private:
		LSSTAGUniformMesh* m_mesh;
		coordinate *m_nodes, *m_curC, *m_delta, *m_inv_delta_ib_ne, *m_inv_delta_ib_tn, *m_delta_ib_max, *m_theta, *m_Vib, *m_VibZY, *m_VibYX;
		CellTypes *m_type;
		double *m_invV, *m_P, *m_Ux, *m_Uy, *m_Uz;
		int m_N, m_M, m_NM;
		bool fileExist;
	};
}

#endif //_LS_STAG_LOADS_CALCULATOR_H_
