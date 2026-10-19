/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ls_stag_uniform_mesh.h
* \author    Puzikova Valeria (vvp) 
* \date      23.10.2013 19:39:00
* \brief     ls_stag_uniform_mesh.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_UNIFORM_MESH_H_
#define _LS_STAG_LS_STAG_UNIFORM_MESH_H_

// ================================================================== INCLUDES
#include "math.h"
#include "cut_sides.h"
#include "ls_stag_base_mesh_params.h"
#include "../MeshFunctions/coord_mesh_function.h"
#include "../MeshFunctions/mesh_function.h"
#include "../solution.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LSSTAGUniformMesh

Данный файл содержит в себе реализацию класса LSSTAGUniformMesh
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(corr, LS_STAG_R_PD, _vec, LS_STAG_CR_PCT, _type);
	LS_STAG_DECLARATION_KERNEL_2D(cvu, const CellTypes*, _baseType, LS_STAG_R_PCT, _x_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshM, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(cvv, const CellTypes*, _baseType, LS_STAG_R_PCT, _y_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL(cvw, const CellTypes*, _baseType, LS_STAG_R_PCT, _z_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(cvtxy, const CellTypes*, _y_Type, LS_STAG_R_PCT, _xy_Type, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(cvtxz, const CellTypes*, _z_Type, LS_STAG_R_PCT, _xz_Type, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(cvtyz, const CellTypes*, _z_Type, LS_STAG_R_PCT, _yz_Type, const int, _baseMeshN,
		const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL(initIB, LS_STAG_R_PC, _curV_C, LS_STAG_R_PC, _cur_C);
	LS_STAG_DECLARATION_KERNEL(calcDeltaIB, LS_STAG_CR_PCT, _baseType, LS_STAG_R_PC, _theta, LS_STAG_R_PC, _delta,
		LS_STAG_R_PC, _delta_ib_ws, LS_STAG_R_PC, _delta_ib_en, LS_STAG_R_PC, _delta_ib_nb, LS_STAG_R_PC, _delta_ib_st, LS_STAG_R_PC, _delta_ib_max, LS_STAG_CR_P(coord<SideTypes>), _SideType, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(buildInvDeltaIb, LS_STAG_R_PC, invDib, LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, 
			LS_STAG_R_PC, _delta, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(buildInvBaseDeltaIb, LS_STAG_R_PC, invDibNE, LS_STAG_R_PC, invDibTN, 
		LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(buildDeltaNodes, LS_STAG_R_PC, _nodes, LS_STAG_R_PC, _delta,
		const double, deltaX, const double, deltaY, const double, deltaZ, const double, a, const double, b, 
		const double, c, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL_2D(mx, LS_STAG_R_PD, _x_VOF, LS_STAG_CR_PD, _VOF, const int, _baseMeshN, 
		const int, _baseMeshM, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(my, LS_STAG_R_PD, _y_VOF, LS_STAG_CR_PD, _VOF, const int, _baseMeshN, 
		const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL(mz, LS_STAG_R_PD, _z_VOF, LS_STAG_CR_PD, _VOF, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(v, LS_STAG_R_PD, _VOF, LS_STAG_CR_PCT, _baseType, LS_STAG_R_PC, _theta,
		LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(vXY, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType);
	LS_STAG_DECLARATION_KERNEL(vXZ, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType);
	LS_STAG_DECLARATION_KERNEL(vYZ, LS_STAG_R_PD, _VOF, LS_STAG_R_PD, _weight, LS_STAG_R_PD, _quarterV,
		LS_STAG_CR_PCT, _baseType);
	LS_STAG_DECLARATION_KERNEL_2D(mXY, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(mXZ, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL_2D(mYZ, LS_STAG_R_PD, _VOF, LS_STAG_CR_PD, _quarterV, const int, _baseMeshN,
		const int, _baseMeshNM);
	LS_STAG_DECLARATION_KERNEL(projectionVib, const double, fullV, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PC, _delta,
		LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_yx, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT, const int, _N);
	LS_STAG_DECLARATION_S_KERNEL(calcTypeOfCellBaseBox, LS_STAG_R_PCT, _baseType, const int, xN, const int, yM, const int, zK);
    LS_STAG_DECLARATION_S_KERNEL(calcTypeOfSideBaseBox, LS_STAG_R_P(coord<SideTypes>), _SideType, const int, xN, const int, yM, const int, zK);
	LS_STAG_DECLARATION_S_KERNEL(calcThetaBaseBox, LS_STAG_R_PC, _theta, const int, xN, const int, yM, const int, zK);
	LS_STAG_DECLARATION_S_KERNEL(calcAlphaBaseBox, LS_STAG_R_PC, _alpha, const int, xN, const int, yM, const int, zK);
/*!
\brief Базовый класс, определяющий интерфейс для работы с LS-STAG-сеткой

В базовой стуктуре LSSTAGUniformMesh определен интерфейс для работы с LS-STAG-сеткой. 
Данная структура может использоваться самостоятельно для работы с равномерной  LS-STAG-сеткой. 
*/
	struct LSSTAGUniformMesh
	{
		/// Типы сеток
		enum MeshKinds
		{
			MESH_UNIFORM = 0,///< Равномерная сетка
			MESH_NONUNIFORM, ///< Неравномерная сетка
			MESH_MOVING      ///< Сетка с подвижными погруженными границами
		};
		explicit LSSTAGUniformMesh(MeshKinds _id = MESH_UNIFORM, int _size = sizeof(LSSTAGUniformMesh)):
		id(_id),
		size(_size),
		NMK(0),
		existXYMesh(0),
		maxCountCutCells(0),
		ibCount(0)
		{
			defineIbCount(1);
		}
		LS_STAG_KERNEL(corr, LS_STAG_R_PD _vec, LS_STAG_CR_PCT _type)
		{
			if (_type[i] == SolidCell) _vec[i] = 0.0;
		}
		LS_STAG_KERNEL_2D(cvu, const CellTypes* _baseType, LS_STAG_R_PCT _x_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshN, const int _baseMeshM, const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, numX = num - i - j * _baseMeshM, m_end = numX + _baseMeshN - 1; 
				numX < m_end; numX++, num++)// curTypeX();// 6.06.14
			{
				const CellTypes c = _baseType[num];
				const CellTypes d = _baseType[num + 1];
				_x_Type[numX] = ((_theta + num)->x < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
			}
		}
		LS_STAG_KERNEL_2D(cvv, const CellTypes* _baseType, LS_STAG_R_PCT _y_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshN, const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, numY = num - j * _baseMeshN, m_end = numY + _baseMeshN; 
				numY < m_end; numY++, num++)// curTypeY();
			{
				const CellTypes c = _baseType[num];
				const CellTypes d = _baseType[num + _baseMeshN];
				_y_Type[numY] = ((_theta + num)->y < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
			}
		}
		LS_STAG_KERNEL(cvw, const CellTypes* _baseType, LS_STAG_R_PCT _z_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshNM)
		{
				const CellTypes c = _baseType[i];
				const CellTypes d = _baseType[i + _baseMeshNM];
				_z_Type[i] = ((_theta + i)->z < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
		}
		LS_STAG_KERNEL_2D(cvtxy, const CellTypes* _y_Type, LS_STAG_R_PCT _xy_Type, const int _baseMeshN, const int NM2,
			const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, numXY = num - i - j * NM2, numY = num - j * _baseMeshN,
				m_end = numXY + _baseMeshN - 1; numXY < m_end; numY++, numXY++)// curTypeXY();
			{
				const CellTypes c = _y_Type[numY];
				const CellTypes d = _y_Type[numY + 1];
				_xy_Type[numXY] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL_2D(cvtxz, const CellTypes* _z_Type, LS_STAG_R_PCT _xz_Type, const int _baseMeshN,
			const int _baseMeshM, const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, numX = num - i - j * _baseMeshM, m_end = numX + _baseMeshN - 1; 
				numX < m_end; numX++, num++)// curTypeXY();
			{
				const CellTypes c = _z_Type[num];
				const CellTypes d = _z_Type[num + 1];
				_xz_Type[numX] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL_2D(cvtyz, const CellTypes* _z_Type, LS_STAG_R_PCT _yz_Type, const int _baseMeshN, const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, numY = num - j * _baseMeshN, m_end = numY + _baseMeshN; 
				numY < m_end; numY++, num++)// curTypeXY();
			{
				const CellTypes c = _z_Type[num];
				const CellTypes d = _z_Type[num + _baseMeshN];
				_yz_Type[numY] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL(initIB, LS_STAG_R_PC _curV_C, LS_STAG_R_PC _cur_C)
		{
			(_curV_C + i)->set(0.0);
			(_cur_C + i)->set(0.0);
		}
		LS_STAG_KERNEL_2D(buildDeltaNodes, LS_STAG_R_PC _nodes, LS_STAG_R_PC _delta, const double deltaX, const double deltaY,
			const double deltaZ, const double a, const double b, const double c, const int _N, const int _NM)
		{
			const int startNum = i * _N + j * _NM;
			coordinate* curT = _nodes + startNum;
			coordinate* curD = _delta + startNum;
			for (int k = 0; k < _N; k++, curT++, curD++)
			{
				curD->set(deltaX, deltaY, deltaZ);
				curT->set(a + (k + 0.5) * deltaX, b + (i + 0.5) * deltaY, c + (j + 0.5) * deltaZ);
			};
		}
		LS_STAG_KERNEL(calcDeltaIB, LS_STAG_CR_PCT _baseType, LS_STAG_R_PC _theta, LS_STAG_R_PC _delta,
			LS_STAG_R_PC _delta_ib_ws, LS_STAG_R_PC _delta_ib_en, LS_STAG_R_PC _delta_ib_nb, 
			LS_STAG_R_PC _delta_ib_st, LS_STAG_R_PC _delta_ib_max, LS_STAG_CR_P(coord<SideTypes>) _SideType, const int _N, const int _NM)
		{ 
			coordinate* curD = _delta + i;
			coordinate* curT = _theta + i;
			coordinate* curDws = _delta_ib_ws + i;
			coordinate* curDen = _delta_ib_en + i;
    coordinate* curDnb = _delta_ib_nb + i;
    coordinate* curDst = _delta_ib_st + i;
    coordinate* curDmax = _delta_ib_max + i;
			switch (_baseType[i])
			{
			case SolidCell:// полуячеечные δ (аналог 2D default: 0.5*dx,0.5*dy) — нужны соседним cut-ячейкам в IB-потоках
				curDen->set(0.5 * curD->y * curD->z, 0.5 * curD->x * curD->z, 0.5 * curD->x * curD->y);
				curDws->set(curDen);
				curDnb->set(curDen);
				curDst->set(curDen);
				curDmax->set(curDen);
				break;
			case FluidCell:
				curDen->set(0.0);
				curDws->set(0.0);
				curDnb->set(0.0);
				curDst->set(0.0);
				curDmax->set(0.0);
				break;
		case B_Border:// нижняя z-граница РО: слоя снизу нет — продление поля вниз ((curT - _NM)->z := curT->z)
		case T_Border:// верхняя z-граница РО: элементов сверху нет; формулы используют только доступные (снизу)
			switch (_SideType[i].z)
			{// δ по форме z-стороны (как у соответствующего типа внутри, но без отсутствующего слоя)
				case FluidSide:
case SolidSide:
curDen->set(0.0);// полностью жидкая z-грань => δ=0 (как FluidCell в интерьере)
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
break;
			case N_Trapezoid: // как N_Hexahedron4
			case S_Trapezoid: // как S_Hexahedron4
{
			const double dz = curD->z;
    double Sxz = 0.5 * curD->x * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * curD->y * dz;
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case W_Trapezoid: // как W_Hexahedron4
			case E_Trapezoid: // как E_Hexahedron4
{
			const double dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * curD->x * dz;
    double Syz = 0.5 * curD->y * dz;
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case NW_Triangle: // как NW_Pentahedron
{
			const double dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * curD->x * dz;                 // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * curT->x * curD->y * dz;                         // ϑ^u_{i,j,k}
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case NE_Triangle: // как NE_Pentahedron
{
    double dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * curD->x * dz;                 // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * (curT - 1)->x * curD->y * dz;                 // ϑ^u_{i-1,j,k}
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case SE_Triangle: // как SE_Pentahedron
{
    double dz = curD->z;
    double Sxz = 0.5 * curT->y * curD->x * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * (curT - 1)->x * curD->y * dz;                 // ϑ^u_{i-1,j,k}
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case SW_Triangle: // как SW_Pentahedron
{
    double dz = curD->z;
    double Sxz = 0.5 * curT->y * curD->x * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * curT->x * curD->y * dz;                 // ϑ^u_{i-1,j,k}
    curDen->set(0.0, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
				break;
			case NE_Pentagon: // как NE_Heptahedron4
{
    double dz = curD->z;
    double Sxz_e = (1.0 - curT->y) * curD->x * dz;
    double Sxz_tb = 0.5 * Sxz_e;
    double Syz_n = (1.0 - curT->x) * curD->y * dz;
    double Syz_tb = 0.5 * Syz_n;

    curDen->set(0.0, Sxz_e, Syz_n);
    curDws->set(0.0);
    curDnb->set(0.0, Sxz_tb, Syz_tb);
    curDst->set(curDnb);
    curDmax->set(curDen);
}
				break;
			case NW_Pentagon: // как NW_Heptahedron4
{
    double dz = curD->z;
    double Sxz_w = (1.0 - curT->y) * curD->x * dz;
    double Sxz_tb = 0.5 * Sxz_w;
    double Syz_n = (1.0 - (curT - 1)->x) * curD->y * dz;
    double Syz_tb = 0.5 * Syz_n;

    curDen->set(0.0, 0.0, Syz_n);
    curDws->set(0.0, Sxz_w, 0.0);
    curDnb->set(0.0, Sxz_tb, Syz_tb);
    curDst->set(curDnb);
    curDmax->set(0.0, Sxz_w, Syz_n);
} 
				break;
			case SE_Pentagon: // как SE_Heptahedron4
{
    double dz = curD->z;
    double Sxz_e = (1.0 - (curT - _N)->y) * curD->x * dz;
    double Sxz_tb = 0.5 * Sxz_e;
    double Syz_s = (1.0 - curT->x) * curD->y * dz;
    double Syz_tb = 0.5 * Syz_s;

    curDen->set(0.0, Sxz_e, 0.0);
    curDws->set(0.0, 0.0, Syz_s);
    curDnb->set(0.0, Sxz_tb, Syz_tb);
    curDst->set(curDnb);
    curDmax->set(0.0, Sxz_e, Syz_s);
}
				break;
			case SW_Pentagon: // как SW_Heptahedron4
{
    double dz = curD->z;
    double Sxz_w = (1.0 - (curT - _N)->y) * curD->x * dz;
    double Sxz_tb = 0.5 * Sxz_w;
    double Syz_s = (1.0 - (curT - 1)->x) * curD->y * dz;
    double Syz_tb = 0.5 * Syz_s;

    curDen->set(0.0);
    curDws->set(0.0, Sxz_w, Syz_s);
    curDnb->set(0.0, Sxz_tb, Syz_tb);
    curDst->set(curDnb);
    curDmax->set(curDws);
} 
				break;
			}
				break; 
		case N_Border: 
		case S_Border:
    curDen->set(0.0, 0.5 * curD->x * curD->z, 0.0);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break;
        case E_Border: 
		case W_Border:
    curDen->set(0.0, 0.0, 0.5 * curD->y * curD->z);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break;
        case TN_Edge: 
		case TS_Edge:
        case BN_Edge: 
		case BS_Edge:
    curDen->set(0.5 * curD->x * curD->y, 0.5 * curD->x * curD->z, 0.0);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break; 
		case TE_Edge: 
		case TW_Edge: 
		case BE_Edge: 
		case BW_Edge:
    curDen->set(0.5 * curD->x * curD->y, 0.0, 0.5 * curD->y * curD->z);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break;
        case NE_Edge:
		 case NW_Edge: 
		 case SE_Edge: 
		 case SW_Edge:
    curDen->set(0.0, 0.5 * curD->x * curD->z, 0.5 * curD->y * curD->z);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break;
        case TNE_Corner: 
		case TNW_Corner: 
		case TSE_Corner: 
		case TSW_Corner:
        case BNE_Corner: 
		case BNW_Corner: 
		case BSE_Corner: 
		case BSW_Corner:
    curDen->set(0.5 * curD->x * curD->y, 0.5 * curD->x * curD->z, 0.5 * curD->y * curD->z);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
	break;
			// ========== 4-угольные гексаэдры (6 типов) ==========
case T_Hexahedron4:  // верх
case B_Hexahedron4:  // низ
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case N_Hexahedron4:  // север
case S_Hexahedron4:  // юг
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case E_Hexahedron4:  // восток
case W_Hexahedron4:  // запад
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
			// ========== Тетраэдры ==========
        case TNW_Tetrahedron: // В-С-З
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * (curT - _NM)->z * dx * dy, 0.5 * (curT - _N)->y * dx * dz, 0.5 * curT->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case TNE_Tetrahedron: // В-С-В
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * (curT - _NM)->z * dx * dy, 0.5 * (curT - _N)->y * dx * dz, 0.5 * (curT - 1)->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case TSE_Tetrahedron: // В-Ю-В
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * (curT - _NM)->z * dx * dy, 0.5 * curT->y * dx * dz, 0.5 * (curT - 1)->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case TSW_Tetrahedron: // В-Ю-З
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * (curT - _NM)->z * dx * dy, 0.5 * curT->y * dx * dz, 0.5 * curT->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case BNW_Tetrahedron: // Н-С-З
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * curT->z * dx * dy, 0.5 * (curT - _N)->y * dx * dz, 0.5 * curT->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case BNE_Tetrahedron: // Н-С-В
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * curT->z * dx * dy, 0.5 * (curT - _N)->y * dx * dz, 0.5 * (curT - 1)->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case BSE_Tetrahedron: // Н-Ю-В
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * curT->z * dx * dy, 0.5 * curT->y * dx * dz, 0.5 * (curT - 1)->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
        case BSW_Tetrahedron: // Н-Ю-З
		{
			double dx = curD->x, dy = curD->y, dz = curD->z;
            curDen->set(0.5 * curT->z * dx * dy, 0.5 * curT->y * dx * dz, 0.5 * curT->x * dy * dz);
			curDws->set(curDen);
			curDnb->set(curDen);
			curDst->set(curDen);
				curDmax->set(curDen);
		}
            break;
			// ========== Пентаэдры (12 типов) ==========
case TN_Pentahedron:  // В-С
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;               // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz; // |ϑ^u_{i-1,j,k} - ϑ^u_{i,j,k}|
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;             // ϑ^w_{i,j,k-1}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case TE_Pentahedron:  // В-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz; // |ϑ^v_{i,j-1,k} - ϑ^v_{i,j,k}|
    double Syz = 0.5 * (curT - 1)->x * dy * dz;                 // ϑ^u_{i-1,j,k}
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;               // ϑ^w_{i,j,k-1}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case TW_Pentahedron:  // В-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz; // |ϑ^v_{i,j-1,k} - ϑ^v_{i,j,k}|
    double Syz = 0.5 * curT->x * dy * dz;                         // ϑ^u_{i,j,k}
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;               // ϑ^w_{i,j,k-1}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case TS_Pentahedron:  // В-Ю
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz; // |ϑ^u_{i-1,j,k} - ϑ^u_{i,j,k}|
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;               // ϑ^w_{i,j,k-1}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case NW_Pentahedron:  // С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;                 // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * curT->x * dy * dz;                         // ϑ^u_{i,j,k}
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy; // |ϑ^w_{i,j,k-1} - ϑ^w_{i,j,k}|
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case NE_Pentahedron:  // С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;                 // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * (curT - 1)->x * dy * dz;                 // ϑ^u_{i-1,j,k}
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy; // |ϑ^w_{i,j,k-1} - ϑ^w_{i,j,k}|
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case BN_Pentahedron:  // Н-С
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;                 // ϑ^v_{i,j-1,k}
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz; // |ϑ^u_{i-1,j,k} - ϑ^u_{i,j,k}|
    double Sxy = 0.5 * curT->z * dx * dy;                         // ϑ^w_{i,j,k}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case SE_Pentahedron:  // Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * (curT - 1)->x * dy * dz;                 // ϑ^u_{i-1,j,k}
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy; // |ϑ^w_{i,j,k-1} - ϑ^w_{i,j,k}|
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case BE_Pentahedron:  // Н-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz; // |ϑ^v_{i,j-1,k} - ϑ^v_{i,j,k}|
    double Syz = 0.5 * (curT - 1)->x * dy * dz;                 // ϑ^u_{i-1,j,k}
    double Sxy = 0.5 * curT->z * dx * dy;                         // ϑ^w_{i,j,k}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case SW_Pentahedron:  // Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * curT->x * dy * dz;                         // ϑ^u_{i,j,k}
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy; // |ϑ^w_{i,j,k-1} - ϑ^w_{i,j,k}|
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case BW_Pentahedron:  // Н-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz; // |ϑ^v_{i,j-1,k} - ϑ^v_{i,j,k}|
    double Syz = 0.5 * curT->x * dy * dz;                         // ϑ^u_{i,j,k}
    double Sxy = 0.5 * curT->z * dx * dy;                         // ϑ^w_{i,j,k}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;
case BS_Pentahedron:  // Н-Ю
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;                         // ϑ^v_{i,j,k}
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz; // |ϑ^u_{i-1,j,k} - ϑ^u_{i,j,k}|
    double Sxy = 0.5 * curT->z * dx * dy;                         // ϑ^w_{i,j,k}
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen);
    curDnb->set(curDen);
    curDst->set(curDen);
				curDmax->set(curDen);
}
    break;// ========== 5-угольные гептаэдры (24 типа) ==========
case TNTE_Heptahedron5:  // B-C-B-B (верх, север, верх, восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TNTW_Heptahedron5:  // B-C-B-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    // Oxy: для "В-С-В-З" – (ϑ^w_{i,j,k-1} - ϑ^w_{i,j,k}) (аналогично TNTE)
    double Sxy = 0.5 * (1.0 - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TNNE_Heptahedron5:  // B-C-C-B
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TNNW_Heptahedron5:  // B-C-C-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TSTE_Heptahedron5:  // B-Ю-B-B
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TSTW_Heptahedron5:  // B-Ю-B-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TSSE_Heptahedron5:  // B-Ю-Ю-B
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TSSW_Heptahedron5:  // B-Ю-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BNNE_Heptahedron5:  // H-C-C-B
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BNNW_Heptahedron5:  // H-C-C-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BSSE_Heptahedron5:  // H-Ю-Ю-B
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BSSW_Heptahedron5:  // H-Ю-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (1.0 - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TENE_Heptahedron5:  // B-B-C-B (верх-восток-север-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (1.0 - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TWNW_Heptahedron5:  // B-З-C-З (верх-запад-север-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (1.0 - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TESE_Heptahedron5:  // B-B-Ю-B (верх-восток-юг-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (1.0 - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case TWSW_Heptahedron5:  // B-З-Ю-З (верх-запад-юг-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (1.0 - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BENE_Heptahedron5:  // H-B-C-B (низ-восток-север-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (1.0 - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BWNW_Heptahedron5:  // H-З-C-З (низ-запад-север-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (1.0 - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BESE_Heptahedron5:  // H-B-Ю-B (низ-восток-юг-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (1.0 - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BWSW_Heptahedron5:  // H-З-Ю-З (низ-запад-юг-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (1.0 - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BNBE_Heptahedron5:  // H-C-H-B (низ-север-низ-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BNBW_Heptahedron5:  // H-C-H-З (низ-север-низ-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BSBE_Heptahedron5:  // H-Ю-H-B (низ-юг-низ-восток)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
case BSBW_Heptahedron5:  // H-Ю-H-З (низ-юг-низ-запад)
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
}
    break;
	// ========== 5-угольные гексаэдры (24 типа) ==========
// Группа A (8 типов) – первая строка таблицы 1 (разность ϑ^v_{i,j-1,k} - ϑ^v_{i,j,k})
case BNBE_Hexahedron5: // Н-С-Н-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * curT->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BNBW_Hexahedron5: // Н-С-Н-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * curT->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TNTE_Hexahedron5: // В-С-В-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TNTW_Hexahedron5: // В-С-В-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BENE_Hexahedron5: // Н-В-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (curT - 1)->x * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BWNW_Hexahedron5: // Н-З-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * curT->x * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TENE_Hexahedron5: // В-В-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * (curT - 1)->x * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TWNW_Hexahedron5: // В-З-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * curT->x * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
// Группа B (8 типов) – вторая строка таблицы 1 (разность ϑ^v_{i,j,k} - ϑ^v_{i,j-1,k})
case BSBE_Hexahedron5: // Н-Ю-Н-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * curT->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BSBW_Hexahedron5: // Н-Ю-Н-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * curT->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TSTE_Hexahedron5: // В-Ю-В-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TSTW_Hexahedron5: // В-Ю-В-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * (curT - _NM)->z * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BESE_Hexahedron5: // Н-В-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (curT - 1)->x * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BWSW_Hexahedron5: // Н-З-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * curT->x * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TESE_Hexahedron5: // В-В-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * (curT - 1)->x * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TWSW_Hexahedron5: // В-З-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * curT->x * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
// Группа C (4 типа) – третья строка таблицы 1 (ϑ^v_{i,j,k})
case TSSE_Hexahedron5: // В-Ю-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TSSW_Hexahedron5: // В-Ю-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BSSE_Hexahedron5: // Н-Ю-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} break;
case BSSW_Hexahedron5: // Н-Ю-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * curT->y * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
// Группа D (4 типа) – четвёртая строка таблицы 1 (ϑ^v_{i,j-1,k})
case TNNE_Hexahedron5: // В-С-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case TNNW_Hexahedron5: // В-С-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BNNE_Hexahedron5: // Н-С-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
case BNNW_Hexahedron5: // Н-С-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * (curT - _N)->y * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
				curDmax->set(curDen);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen);
} 
break;
// ========== 4-угольные гептаэдры (12 типов) ==========
// ------------------------------------------------------------
// TE_Heptahedron4  (В-В, Top-East)
// ------------------------------------------------------------
case TE_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz_eq = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz_t = (1.0 - curT->x) * dy * dz;
    double Syz_ns = 0.5 * Syz_t;
    double Sxy_e = (1.0 - curT->z) * dx * dy;
    double Sxy_ns = 0.5 * Sxy_e;

    curDen->set(Sxy_e, Sxz_eq, Syz_ns);
    curDws->set(0.0, Sxz_eq, Syz_ns);
    curDnb->set(Sxy_ns, Sxz_eq, 0.0);
    curDst->set(Sxy_ns, Sxz_eq, Syz_t);
    curDmax->set(Sxy_e, Sxz_eq, Syz_t);
} 
break;
// ------------------------------------------------------------
// TW_Heptahedron4  (В-З, Top-West)
// ------------------------------------------------------------
case TW_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz_eq = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz_t = (1.0 - (curT - 1)->x) * dy * dz;
    double Syz_ns = 0.5 * Syz_t;
    double Sxy_w = (1.0 - curT->z) * dx * dy;
    double Sxy_ns = 0.5 * Sxy_w;

    curDen->set(0.0, Sxz_eq, Syz_ns);
    curDws->set(Sxy_w, Sxz_eq, Syz_ns);
    curDnb->set(Sxy_ns, Sxz_eq, 0.0);
    curDst->set(Sxy_ns, Sxz_eq, Syz_t);
    curDmax->set(Sxy_w, Sxz_eq, Syz_t);
} 
break;
// ------------------------------------------------------------
// BE_Heptahedron4  (Н-В, Bottom-East)
// ------------------------------------------------------------
case BE_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz_eq = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz_b = (1.0 - curT->x) * dy * dz;
    double Syz_ns = 0.5 * Syz_b;
    double Sxy_e = (1.0 - (curT - _NM)->z) * dx * dy;
    double Sxy_ns = 0.5 * Sxy_e;

    curDen->set(Sxy_e, Sxz_eq, Syz_ns);
    curDws->set(0.0, Sxz_eq, Syz_ns);
    curDnb->set(Sxy_ns, Sxz_eq, Syz_b);
    curDst->set(Sxy_ns, Sxz_eq, 0.0);
    curDmax->set(Sxy_e, Sxz_eq, Syz_b);
} 
break;
// ------------------------------------------------------------
// BW_Heptahedron4  (Н-З, Bottom-West)
// ------------------------------------------------------------
case BW_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz_eq = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz_b = (1.0 - (curT - 1)->x) * dy * dz;
    double Sxy_w = (1.0 - (curT - _NM)->z) * dx * dy;
    double Syz_ns = 0.5 * Syz_b;
    double Sxy_ns = 0.5 * Sxy_w;

    curDen->set(0.0, Sxz_eq, Syz_ns);
    curDws->set(Sxy_w, Sxz_eq, Syz_ns);
    curDnb->set(Sxy_ns, Sxz_eq, Syz_b);
    curDst->set(Sxy_ns, Sxz_eq, 0.0);
    curDmax->set(Sxy_w, Sxz_eq, Syz_b);
} 
break;
// ------------------------------------------------------------
// TN_Heptahedron4  (В-С, Top-North)
// ------------------------------------------------------------
case TN_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Syz_eq = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxz_t = (1.0 - curT->y) * dx * dz;
    double Sxz_ew = 0.5 * Sxz_t;
    double Sxy_n = (1.0 - curT->z) * dx * dy;
    double Sxy_ew = 0.5 * Sxy_n;

    curDen->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDws->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDnb->set(Sxy_n, 0.0, Syz_eq);
    curDst->set(0.0, Sxz_t, Syz_eq);
    curDmax->set(Sxy_n, Sxz_t, Syz_eq);
} 
break;
// ------------------------------------------------------------
// TS_Heptahedron4  (В-Ю, Top-South)
// ------------------------------------------------------------
case TS_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Syz_eq = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxz_t = (1.0 - (curT - _N)->y) * dx * dz;
    double Sxz_ew = 0.5 * Sxz_t;
    double Sxy_s = (1.0 - curT->z) * dx * dy;
    double Sxy_ew = 0.5 * Sxy_s;

    curDen->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDws->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDnb->set(0.0, 0.0, Syz_eq);
    curDst->set(Sxy_s, Sxz_t, Syz_eq);
    curDmax->set(curDst);
} 
break;
// ------------------------------------------------------------
// BN_Heptahedron4  (Н-С, Bottom-North)
// ------------------------------------------------------------
case BN_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Syz_eq = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxz_b = (1.0 - curT->y) * dx * dz;
    double Sxz_ew = 0.5 * Sxz_b;
    double Sxy_n = (1.0 - (curT - _NM)->z) * dx * dy;
    double Sxy_ew = 0.5 * Sxy_n;

    curDen->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDws->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDnb->set(Sxy_n, Sxz_b, Syz_eq);
    curDst->set(0.0, 0.0, Syz_eq);
    curDmax->set(curDnb);
} 
break;
// ------------------------------------------------------------
// BS_Heptahedron4  (Н-Ю, Bottom-South)
// ------------------------------------------------------------
case BS_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Syz_eq = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxz_b = (1.0 - (curT - _N)->y) * dx * dz;
    double Sxz_ew = 0.5 * Sxz_b;
    double Sxy_s = (1.0 - (curT - _NM)->z) * dx * dy;
    double Sxy_ew = 0.5 * Sxy_s;

    curDen->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDws->set(Sxy_ew, Sxz_ew, Syz_eq);
    curDnb->set(0.0, Sxz_b, Syz_eq);
    curDst->set(Sxy_s, 0.0, Syz_eq);
    curDmax->set(Sxy_s, Sxz_b, Syz_eq);
} 
break;
// ------------------------------------------------------------
// NE_Heptahedron4  (С-В, North-East)
// ------------------------------------------------------------
case NE_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxy_eq = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    double Sxz_e = (1.0 - curT->y) * dx * dz;
    double Sxz_tb = 0.5 * Sxz_e;
    double Syz_n = (1.0 - curT->x) * dy * dz;
    double Syz_tb = 0.5 * Syz_n;

    curDen->set(Sxy_eq, Sxz_e, Syz_n);
    curDws->set(Sxy_eq, 0.0, 0.0);
    curDnb->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDst->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDmax->set(Sxy_eq, Sxz_e, Syz_n);
} 
break;
// ------------------------------------------------------------
// NW_Heptahedron4  (С-З, North-West)
// ------------------------------------------------------------
case NW_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxy_eq = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    double Sxz_w = (1.0 - curT->y) * dx * dz;
    double Sxz_tb = 0.5 * Sxz_w;
    double Syz_n = (1.0 - (curT - 1)->x) * dy * dz;
    double Syz_tb = 0.5 * Syz_n;

    curDen->set(Sxy_eq, 0.0, Syz_n);
    curDws->set(Sxy_eq, Sxz_w, 0.0);
    curDnb->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDst->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDmax->set(Sxy_eq, Sxz_w, Syz_n);
} 
break;
// ------------------------------------------------------------
// SE_Heptahedron4  (Ю-В, South-East)
// ------------------------------------------------------------
case SE_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxy_eq = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    double Sxz_e = (1.0 - (curT - _N)->y) * dx * dz;
    double Sxz_tb = 0.5 * Sxz_e;
    double Syz_s = (1.0 - curT->x) * dy * dz;
    double Syz_tb = 0.5 * Syz_s;

    curDen->set(Sxy_eq, Sxz_e, 0.0);
    curDws->set(Sxy_eq, 0.0, Syz_s);
    curDnb->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDst->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDmax->set(Sxy_eq, Sxz_e, Syz_s);
} 
break;
// ------------------------------------------------------------
// SW_Heptahedron4  (Ю-З, South-West)
// ------------------------------------------------------------
case SW_Heptahedron4:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxy_eq = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    double Sxz_w = (1.0 - (curT - _N)->y) * dx * dz;
    double Sxz_tb = 0.5 * Sxz_w;
    double Syz_s = (1.0 - (curT - 1)->x) * dy * dz;
    double Syz_tb = 0.5 * Syz_s;

    curDen->set(Sxy_eq, 0.0, 0.0);
    curDws->set(Sxy_eq, Sxz_w, Syz_s);
    curDnb->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDst->set(Sxy_eq, Sxz_tb, Syz_tb);
    curDmax->set(curDws);
} 
break;
// ========== 6-угольные гептаэдры (8 типов) ==========
case TNE_Heptahedron6: // В-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case TNW_Heptahedron6: // В-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case TSE_Heptahedron6: // В-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case TSW_Heptahedron6: // В-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs((curT - _NM)->z - curT->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case BNE_Heptahedron6: // Н-С-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case BNW_Heptahedron6: // Н-С-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs((curT - _N)->y - curT->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case BSE_Heptahedron6: // Н-Ю-В
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs((curT - 1)->x - curT->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
case BSW_Heptahedron6: // Н-Ю-З
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double Sxz = 0.5 * fabs(curT->y - (curT - _N)->y) * dx * dz;
    double Syz = 0.5 * fabs(curT->x - (curT - 1)->x) * dy * dz;
    double Sxy = 0.5 * fabs(curT->z - (curT - _NM)->z) * dx * dy;
    curDen->set(Sxy, Sxz, Syz);
    curDws->set(curDen); curDnb->set(curDen); curDst->set(curDen); curDmax->set(curDen);
} 
break;
// ========== 3-угольные гептаэдры (8 типов) ==========
case TNE_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - curT->y) * dx * dz;
    double val_yz = (1.0 - curT->x) * dy * dz;
    double val_xy = (1.0 - curT->z) * dx * dy;
    curDen->set(val_xy, val_xz, val_yz);
    curDws->set(0.0);
    curDnb->set(val_xy, 0.0, 0.0);
    curDst->set(0.0, val_xz, val_yz); 
	curDmax->set(0.0);
} 
break;
case TNW_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - curT->y) * dx * dz;
    double val_yz = (1.0 - (curT - 1)->x) * dy * dz;
    double val_xy = (1.0 - curT->z) * dx * dy;
    curDen->set(0.0, 0.0, val_yz);
    curDws->set(val_xy, val_xz, 0.0);
    curDnb->set(val_xy, 0.0, 0.0);
    curDst->set(0.0, val_xz, val_yz);
	curDmax->set(0.0);
} 
break;
case TSE_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - (curT - _N)->y) * dx * dz;
    double val_yz = (1.0 - curT->x) * dy * dz;
    double val_xy = (1.0 - curT->z) * dx * dy;
    curDen->set(val_xy, val_xz, 0.0);
    curDws->set(0.0, 0.0, val_yz);
    curDnb->set(0.0);
    curDst->set(val_xy, val_xz, val_yz);
	curDmax->set(0.0);
} 
break;
case TSW_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - (curT - _N)->y) * dx * dz;
    double val_yz = (1.0 - (curT - 1)->x) * dy * dz;
    double val_xy = (1.0 - curT->z) * dx * dy;
    curDen->set(0.0);
    curDws->set(val_xy, val_xz, val_yz);
    curDnb->set(0.0);
    curDst->set(curDws);
	curDmax->set(0.0);
} 
break;
case BNE_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - curT->y) * dx * dz;
    double val_yz = (1.0 - curT->x) * dy * dz;
    double val_xy = (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(val_xy, val_xz, val_yz);
    curDws->set(0.0);
    curDnb->set(curDen);
    curDst->set(0.0);
	curDmax->set(0.0);
} 
break;
case BNW_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - curT->y) * dx * dz;
    double val_yz = (1.0 - (curT - 1)->x) * dy * dz;
    double val_xy = (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(0.0, 0.0, val_yz);
    curDws->set(val_xy, val_xz, 0.0);
    curDnb->set(val_xy, val_xz, val_yz);
    curDst->set(0.0);
	curDmax->set(0.0);
} 
break;
case BSE_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - (curT - _N)->y) * dx * dz;
    double val_yz = (1.0 - curT->x) * dy * dz;
    double val_xy = (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(val_xy, val_xz, 0.0);
    curDws->set(0.0, 0.0, val_yz);
    curDnb->set(0.0, val_xz, val_yz);
    curDst->set(val_xy, 0.0, 0.0);
	curDmax->set(0.0);
} 
break;
case BSW_Heptahedron3:
{
    double dx = curD->x, dy = curD->y, dz = curD->z;
    double val_xz = (1.0 - (curT - _N)->y) * dx * dz;
    double val_yz = (1.0 - (curT - 1)->x) * dy * dz;
    double val_xy = (1.0 - (curT - _NM)->z) * dx * dy;
    curDen->set(0.0);
    curDws->set(val_xy, val_xz, val_yz);
    curDnb->set(0.0, val_xz, val_yz);
    curDst->set(val_xy, 0.0, 0.0);
	curDmax->set(0.0);
} 
break;
			}//switch
		}
		LS_STAG_KERNEL_INLINE double getThetaZ(coordinate*& curA, SideTypes topSideType, const int _N)
		{// внутри, не на границах области
			double tz = 1.0;
			// ========== Вычисление θ^w (верхняя грань) ==========
			switch (topSideType) {
			case FluidSide:
				break;
			case SolidSide:
				tz = 0.0;
				break;
				// Трапеции
			case N_Trapezoid:      // 11 (С)
			case S_Trapezoid:      // 1100 (Ю)
				tz = 0.5 * (curA->y + (curA - 1)->y);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case W_Trapezoid:      // 1001 (З)
			case E_Trapezoid:      // 110 (В)
				tz = 0.5 * (curA->x + (curA - _N)->x);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
				// Пятиугольники
			case NW_Pentagon:      // 1011 (С-З)
				tz = 1.0 - 0.5 * (1.0 - (curA - 1)->y) * (1.0 - curA->x);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case NE_Pentagon:      // 111 (С-В)
				tz = 1.0 - 0.5 * (1.0 - curA->y) * (1.0 - curA->x);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case SW_Pentagon:      // 1101 (Ю-З)
				tz = 1.0 - 0.5 * (1.0 - (curA - 1)->y) * (1.0 - (curA - _N)->x);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case SE_Pentagon:      // 1110 (Ю-В)
				tz = 1.0 - 0.5 * (1.0 - curA->y) * (1.0 - (curA - _N)->x);
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
				// Треугольники
			case NW_Triangle:      // 1 (С-З)
				tz = 0.5 * curA->y * (curA - _N)->x;
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case NE_Triangle:      // 10 (С-В)
				tz = 0.5 * (curA - 1)->y * (curA - _N)->x;
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case SW_Triangle:      // 1000 (Ю-З)
				tz = 0.5 * curA->y * curA->x;
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			case SE_Triangle:      // 100 (Ю-В)
				tz = 0.5 * (curA - 1)->y * curA->x;
				if (tz < THETA_EQ_ZERO) tz = 0.0;
				break;
			}
			return tz;
		}
		LS_STAG_KERNEL_INLINE void tbBuildTheta(coordinate*& curT, LS_STAG_R_PC _theta, LS_STAG_R_PC _alpha,
			LS_STAG_CR_P(coord<SideTypes>) _baseType, const int _N)
		{
			const int num = (int)(curT - _theta);
			coord<SideTypes> curSideTypes = _baseType[num];
			coordinate * curA = _alpha + num; 
			double tz = LSSTAGUniformMesh::getThetaZ(curA, curSideTypes.z, _N); 
			double tx = curA->y, ty = curA->x;
			curT->set((tx < THETA_EQ_ZERO) ? 0.0 : tx, (ty < THETA_EQ_ZERO) ? 0.0 : ty, tz);// alpha
			curT++;
		}
		LS_STAG_KERNEL_INLINE void BuildTheta(coordinate*& curT, LS_STAG_R_PC _theta, LS_STAG_R_PC _alpha, 
			LS_STAG_CR_P(coord<SideTypes>) _baseType, const int _N, const int _NM)
		{// внутри, не на границах области
			const int num = (int)(curT - _theta);
			coord<SideTypes> curSideTypes =  _baseType[num];
			coordinate* curA = _alpha + num;
			double tx = 1.0, ty = 1.0;
			double tz = LSSTAGUniformMesh::getThetaZ(curA, curSideTypes.z, _N);
			// ========== Вычисление θ^u (восточная грань) ==========
switch (curSideTypes.x) {
    case FluidSide:
       break;
    case SolidSide:
	tx = 0.0;
       break;
    // Трапеции
    case N_Trapezoid:      // 11 (С)
    case S_Trapezoid:      // 1100 (Ю)
        tx = 0.5 * (curA->z + (curA - _N)->z);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case W_Trapezoid:      // 1001 (З)
    case E_Trapezoid:      // 110 (В)
        tx = 0.5 * (curA->y + (curA - _NM)->y);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    // Пятиугольники
    case NW_Pentagon:      // 1011 (С-З)
        tx = 1.0 - 0.5 * (1.0 - (curA - _N)->z) * (1.0 - curA->y);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case NE_Pentagon:      // 111 (С-В)
        tx = 1.0 - 0.5 * (1.0 - curA->z) * (1.0 - curA->y);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case SW_Pentagon:      // 1101 (Ю-З)
        tx = 1.0 - 0.5 * (1.0 - (curA - _N)->z) * (1.0 - (curA - _NM)->y);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case SE_Pentagon:      // 1110 (Ю-В)
        tx = 1.0 - 0.5 * (1.0 - curA->z) * (1.0 - (curA - _NM)->y);
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    // Треугольники
    case NW_Triangle:      // 1 (С-З)
        tx = 0.5 * curA->z * (curA - _NM)->y;
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case NE_Triangle:      // 10 (С-В)
        tx = 0.5 * (curA - _N)->z * (curA - _NM)->y;
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case SW_Triangle:      // 1000 (Ю-З)
        tx = 0.5 * curA->z * curA->y;
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
    case SE_Triangle:      // 100 (Ю-В)
        tx = 0.5 * (curA - _N)->z * curA->y;
		if (tx < THETA_EQ_ZERO) tx = 0.0;
        break;
}
// ========== Вычисление θ^v (северная грань) ==========
switch (curSideTypes.y) {
    case FluidSide:
        break;
    case SolidSide:
	ty = 0.0;
       break;
    // Трапеции
    case N_Trapezoid:      // 11 (С)
    case S_Trapezoid:      // 1100 (Ю)
        ty = 0.5 * (curA->z + (curA - 1)->z);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case W_Trapezoid:      // 1001 (З)
    case E_Trapezoid:      // 110 (В)
        ty = 0.5 * (curA->x + (curA - _NM)->x);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    // Пятиугольники
    case NW_Pentagon:      // 1011 (С-З)
        ty = 1.0 - 0.5 * (1.0 - (curA - 1)->z) * (1.0 - curA->x);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case NE_Pentagon:      // 111 (С-В)
        ty = 1.0 - 0.5 * (1.0 - curA->z) * (1.0 - curA->x);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case SW_Pentagon:      // 1101 (Ю-З)
        ty = 1.0 - 0.5 * (1.0 - (curA - 1)->z) * (1.0 - (curA - _NM)->x);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case SE_Pentagon:      // 1110 (Ю-В)
        ty = 1.0 - 0.5 * (1.0 - curA->z) * (1.0 - (curA - _NM)->x);
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    // Треугольники
    case NW_Triangle:      // 1 (С-З)
        ty = 0.5 * curA->z * (curA - _NM)->x;
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case NE_Triangle:      // 10 (С-В)
        ty = 0.5 * (curA - 1)->z * (curA - _NM)->x;
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case SW_Triangle:      // 1000 (Ю-З)
        ty = 0.5 * curA->z * curA->x;
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
    case SE_Triangle:      // 100 (Ю-В)
        ty = 0.5 * (curA - 1)->z * curA->x;
		if (ty < THETA_EQ_ZERO) ty = 0.0;
        break;
}
			curT->set(tx, ty, tz);// alpha
			curT++;
		};
		LS_STAG_KERNEL(v, LS_STAG_R_PD _VOF, LS_STAG_CR_PCT _baseType, LS_STAG_R_PC _theta, LS_STAG_R_PC _alpha, 
			LS_STAG_R_PC _delta, const int _N, const int _NM)
		{
			coordinate curD = _delta[i];
			double fullV = curD.x * curD.y * curD.z;
			double vol = fullV;
			coordinate* curT = _theta + i;
			coordinate* curA = _alpha + i;
			switch (_baseType[i])
			{
			case FluidCell:
				break;
			case B_Border:
			case T_Border:
			vol *= curT->z;
				break;
			case SolidCell:
				vol = 0.0;
				break;
			// ========== Объёмы тетраэдров ==========
case BSW_Tetrahedron: // Н-Ю-З (bottom-south-west)
    vol = curA->x * curT->x * fullV * ONE_THIRD;
    break;
case BSE_Tetrahedron: // Н-Ю-В (bottom-south-east)
    vol = curA->x * (curT - 1)->x * fullV * ONE_THIRD;
    break;
case BNE_Tetrahedron: // Н-С-В (bottom-north-east)
    vol = (curA - _N)->x * (curT - 1)->x * fullV * ONE_THIRD;
    break;
case BNW_Tetrahedron: // Н-С-З (bottom-north-west)
    vol = (curA - _N)->x * curT->x * fullV * ONE_THIRD;
    break;
case TSW_Tetrahedron: // В-Ю-З (top-south-west)
    vol = (curA - _NM)->x * curT->x * fullV * ONE_THIRD;
    break;
case TSE_Tetrahedron: // В-Ю-В (top-south-east)
    vol = (curA - _NM)->x * (curT - 1)->x * fullV * ONE_THIRD;
    break;
case TNE_Tetrahedron: // В-С-В (top-north-east)
    vol = (curA - _N - _NM)->x * (curT - 1)->x * fullV * ONE_THIRD;
    break;
case TNW_Tetrahedron: // В-С-З (top-north-west)
    vol = (curA - _N - _NM)->x * curT->x * fullV * ONE_THIRD;
    break;
// ========== Треугольные гептаэдры (8 типов) ==========
case TNE_Heptahedron3: // В-С-В
    vol = (3.0 - (1.0 - curA->x) * (1.0 - curT->x)) * fullV * ONE_THIRD;
    break;
case TNW_Heptahedron3: // В-С-З
    vol = (3.0 - (1.0 - curA->x) * (1.0 - (curT - 1)->x)) * fullV * ONE_THIRD;
    break;
case TSW_Heptahedron3: // В-Ю-З
    vol = (3.0 - (1.0 - (curA - _N)->x) * (1.0 - (curT - 1)->x)) * fullV * ONE_THIRD;
    break;
case TSE_Heptahedron3: // В-Ю-В
    vol = (3.0 - (1.0 - (curA - _N)->x) * (1.0 - curT->x)) * fullV * ONE_THIRD;
    break;
case BNE_Heptahedron3: // Н-С-В
    vol = (3.0 - (1.0 - (curA - _NM)->x) * (1.0 - curT->x)) * fullV * ONE_THIRD;
    break;
case BNW_Heptahedron3: // Н-С-З
    vol = (3.0 - (1.0 - (curA - _NM)->x) * (1.0 - (curT - 1)->x)) * fullV * ONE_THIRD;
    break;
case BSW_Heptahedron3: // Н-Ю-З
    vol = (3.0 - (1.0 - (curA - _N - _NM)->x) * (1.0 - (curT - 1)->x)) * fullV * ONE_THIRD;
    break;
case BSE_Heptahedron3: // Н-Ю-В
    vol = (3.0 - (1.0 - (curA - _N - _NM)->x) * (1.0 - curT->x)) * fullV * ONE_THIRD;
    break;
// ========== Пентаэдры и четырёхугольные гептаэдры (24 типа) ==========
// ==========  + Четырёхугольные гексаэдры (6 типов) ==========
// Группа 1: через ϑ^w + ϑ^w_{k-1}
case SW_Pentahedron:      // Ю-З (10001000)
case SE_Pentahedron:      // Ю-В (1000100)
case NE_Pentahedron:      // С-В (100010)
case NW_Pentahedron:      // С-З (10001)
case SE_Heptahedron4:     // Ю-В (11101110)
case SW_Heptahedron4:     // Ю-З (11011101)
case NE_Heptahedron4:     // С-В (1110111)
case NW_Heptahedron4:     // С-З (10111011)
    vol = (curT->z + (curT - _NM)->z) * fullV * 0.5;
    break;
// Группа 2: через ϑ^u + ϑ^u_{i-1}
case BS_Pentahedron:      // Н-Ю (11000000)
case BN_Pentahedron:      // Н-С (110000)
case TS_Pentahedron:      // В-Ю (1100)
case TN_Pentahedron:      // В-С (11)
case BS_Heptahedron4:     // Н-Ю (11111100)
case BN_Heptahedron4:     // Н-С (11110011)
case TS_Heptahedron4:     // В-Ю (11001111)
case TN_Heptahedron4:     // В-С (111111)
case T_Hexahedron4:   // В(ерх) (1111)
case B_Hexahedron4:   // Н (11110000)
case N_Hexahedron4:   // С (110011)
case S_Hexahedron4:   // Ю (11001100)
    vol = (curT->x + (curT - 1)->x) * fullV * 0.5;
    break;
// Группа 3: через ϑ^v + ϑ^v_{j-1}
case BW_Pentahedron:      // Н-З (10010000)
case BE_Pentahedron:      // Н-В (1100000)
case TW_Pentahedron:      // В-З (1001)
case TE_Pentahedron:      // В-В (110)
case BW_Heptahedron4:     // Н-З (11111001)
case BE_Heptahedron4:     // Н-В (11110110)
case TW_Heptahedron4:     // В-З (10011111)
case TE_Heptahedron4:     // В-В (1101111)
case E_Hexahedron4:   // В (1100110)
case W_Hexahedron4:   // З (10011001)
    vol = (curT->y + (curT - _N)->y) * fullV * 0.5;
    break;
	// ========== Пятиугольные гептаэдры (24 типа) ==========
// Группа 1: через ϑ^w (θ_{k-1} и θ_k), α^z и α^y
case TENE_Heptahedron5: // В-В-С-В
{
    double theta1 = (curT - _NM)->z;
    double alpha_y = (curA - _NM)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - _N)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = (theta1 + curT->z + term) * fullV * 0.5;
}
    break;
case BESE_Heptahedron5: // Н-В-Ю-В
{
    double theta2 = curT->z;
    double alpha_y = curA->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * curA->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = ((curT - _NM)->z + theta2 + term) * fullV * 0.5;
}
    break;
case BENE_Heptahedron5: // Н-В-С-В
{
    double theta2 = curT->z;
    double alpha_y = curA->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - _N)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = ((curT - _NM)->z + theta2 + term) * fullV * 0.5;
}
    break;
case TESE_Heptahedron5: // В-В-Ю-В
{
    double theta1 = (curT - _NM)->z;
    double alpha_y = (curA - _NM)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * curA->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = (theta1 + curT->z + term) * fullV * 0.5;
}
    break;
case BWSW_Heptahedron5: // Н-З-Ю-З
{
    double theta2 = curT->z;
    double alpha_y = (curA - 1)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - 1)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = ((curT - _NM)->z + theta2 + term) * fullV * 0.5;
}
    break;
case BWNW_Heptahedron5: // Н-З-С-З
{
    double theta2 = curT->z;
    double alpha_y = (curA - 1)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - 1 - _N)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = ((curT - _NM)->z + theta2 + term) * fullV * 0.5;
}
    break;
case TWSW_Heptahedron5: // В-З-Ю-З
{
    double theta1 = (curT - _NM)->z;
    double alpha_y = (curA - 1 - _NM)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - 1)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = (theta1 + curT->z + term) * fullV * 0.5;
}
    break;
case TWNW_Heptahedron5: // В-З-С-З
{
    double theta1 = (curT - _NM)->z;
    double alpha_y = (curA - 1 - _NM)->y;
    double denom = 1.0 - alpha_y;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_y /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - 1 - _N)->z) * alpha_y * alpha_y * ONE_THIRD;
	}
    vol = (theta1 + curT->z + term) * fullV * 0.5;
}
    break;
// Группа 2: через ϑ^u (θ_{i-1} и θ_i), α^x и α^z
case BNNW_Heptahedron5: // Н-С-С-З
{
    double theta2 = curT->x;
    double alpha_z = curA->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * curA->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = ((curT - 1)->x + theta2 + term) * fullV * 0.5;
}
    break;
case BNNE_Heptahedron5: // Н-С-С-В
{
    double theta1 = (curT - 1)->x;
    double alpha_z = (curA - 1)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * curA->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = (theta1 + curT->x + term) * fullV * 0.5;
}
    break;
case TNNW_Heptahedron5: // В-С-С-З
{
    double theta2 = curT->x;
    double alpha_z = curA->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - _NM)->z) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = ((curT - 1)->x + theta2 + term) * fullV * 0.5;
}
    break;
case TNNE_Heptahedron5: // В-С-С-В
{
    double theta1 = (curT - 1)->x;
    double alpha_z = (curA - 1)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - _NM)->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = (theta1 + curT->x + term) * fullV * 0.5;
}
    break;
case BSSE_Heptahedron5: // Н-Ю-Ю-В
{
    double theta1 = (curT - 1)->x;
    double alpha_z = (curA - 1 - _N)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - _N)->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = (theta1 + curT->x + term) * fullV * 0.5;
}
    break;
case BSSW_Heptahedron5: // Н-Ю-Ю-З
{
    double theta2 = curT->x;
    double alpha_z = (curA - _N)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - _N)->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = ((curT - 1)->x + theta2 + term) * fullV * 0.5;
}
    break;
case TSSE_Heptahedron5: // В-Ю-Ю-В
{
    double theta1 = (curT - 1)->x;
    double alpha_z = (curA - 1 - _N)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = (theta1 + curT->x + term) * fullV * 0.5;
}
    break;
case TSSW_Heptahedron5: // В-Ю-Ю-З
{
    double theta2 = curT->x;
    double alpha_z = (curA - _N)->z;
    double denom = 1.0 - alpha_z;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_z /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha_z * alpha_z * ONE_THIRD;
	}
    vol = ((curT - 1)->x + theta2 + term) * fullV * 0.5;
}
    break;
// Группа 3: через ϑ^v (θ_{j-1} и θ_j), α^y и α^x
case TSTW_Heptahedron5: // В-Ю-В-З
{
    double theta2 = curT->y;
    double alpha_x = curA->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * curA->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = ((curT - _N)->y + theta2 + term) * fullV * 0.5;
}
    break;
case TSTE_Heptahedron5: // В-Ю-В-В
{
    double theta2 = curT->y;
    double alpha_x = curA->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - 1)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = ((curT - _N)->y + theta2 + term) * fullV * 0.5;
}
    break;
case TNTW_Heptahedron5: // В-С-В-З
{
    double theta1 = (curT - _N)->y;
    double alpha_x = (curA - _N)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - 1)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = (theta1 + curT->y + term) * fullV * 0.5;
}
    break;
case TNTE_Heptahedron5: // В-С-В-В
{
    double theta1 = (curT - _N)->y;
    double alpha_x = (curA - _N)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * curA->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = (theta1 + curT->y + term) * fullV * 0.5;
}
    break;
case BSBW_Heptahedron5: // Н-Ю-Н-З
{
    double theta2 = curT->y;
    double alpha_x = (curA - _NM)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - _NM)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = ((curT - _N)->y + theta2 + term) * fullV * 0.5;
}
    break;
case BNBW_Heptahedron5: // Н-С-Н-З
{
    double theta1 = (curT - _N)->y;
    double alpha_x = (curA - _N - _NM)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - _NM)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = (theta1 + curT->y + term) * fullV * 0.5;
}
    break;
case BNBE_Heptahedron5: // Н-С-Н-В
{
    double theta1 = (curT - _N)->y;
    double alpha_x = (curA - _N - _NM)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta1) * (3.0 - 2.0 * (curA - 1 - _NM)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = (theta1 + curT->y + term) * fullV * 0.5;
}
    break;
case BSBE_Heptahedron5: // Н-Ю-Н-В
{
    double theta2 = curT->y;
    double alpha_x = (curA - _NM)->x;
    double denom = 1.0 - alpha_x;
	double term = 0.0;
	if(denom > EQ_DOUBLE)
	{
		alpha_x /= denom;
		term = (1.0 - theta2) * (3.0 - 2.0 * (curA - 1 - _NM)->y) * alpha_x * alpha_x * ONE_THIRD;
	}
    vol = ((curT - _N)->y + theta2 + term) * fullV * 0.5;
}
    break;
	// ========== Пятиугольные гексаэдры (24 типа) ==========
// Группа A: через ϑ^w (8 типов)
case BWSW_Hexahedron5:   // Н-З-Ю-З (10011000)
{
    double denom = (curA - _NM)->y;                  // α^y_{i,j,k-1}
    double term = (curT - _NM)->z;
	if(denom > EQ_DOUBLE)
	{
		double coef = (1.0 - denom) / denom;
		denom =  1.0 - (1.0 + 2.0 * (curA - _N)->z) * coef * coef * ONE_THIRD;
		term *= denom;
	}               
    vol = (curT->z + term) * fullV * 0.5;
}
    break;
	case TWNW_Hexahedron5:   // В-З-С-З (11001)
{
    double denom = curA->y;                         // α^y_{i,j,k}
    double term2 = curT->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->z) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _NM)->z + term2) * fullV * 0.5;
} 
break;
case TWSW_Hexahedron5:   // В-З-Ю-З (10001001)
{
    double denom = curA->y;                         // α^y_{i,j,k}
    double term2 = curT->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _N)->z) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _NM)->z + term2) * fullV * 0.5;
} 
break;
case BWNW_Hexahedron5:   // Н-З-С-З (10010001)
{
    double denom = (curA - _NM)->y;                 // α^y_{i,j,k-1}
    double term2 = (curT - _NM)->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->z) * coef * coef * ONE_THIRD);
    }
    vol = (curT->z + term2) * fullV * 0.5;
} 
break;
case TENE_Hexahedron5:   // В-В-С-В (100110)
{
    double denom = (curA - 1)->y;                   // α^y_{i-1,j,k}
    double term2 = curT->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->z) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _NM)->z + term2) * fullV * 0.5;
} 
break;
case TESE_Hexahedron5:   // В-В-Ю-В (1000110)
{
    double denom = (curA - 1)->y;                   // α^y_{i-1,j,k}
    double term2 = curT->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1 - _N)->z) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _NM)->z + term2) * fullV * 0.5;
} 
break;
case BENE_Hexahedron5:   // Н-В-С-В (1100010)
{
    double denom = (curA - 1 - _NM)->y;             // α^y_{i-1,j,k-1}
    double term2 = (curT - _NM)->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->z) * coef * coef * ONE_THIRD);
    }
    vol = (curT->z + term2) * fullV * 0.5;
} 
break;
case BESE_Hexahedron5:   // Н-В-Ю-В (1100100)
{
    double denom = (curA - 1 - _NM)->y;             // α^y_{i-1,j,k-1}
    double term2 = (curT - _NM)->z;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1 - _N)->z) * coef * coef * ONE_THIRD);
    }
    vol = (curT->z + term2) * fullV * 0.5;
} 
break;
// Группа B: через ϑ^u (8 типов)
case TSSW_Hexahedron5:   // В-Ю-Ю-З (10001100)
{
    double denom = (curA - 1)->z;                   // α^z_{i-1,j,k}
    double term2 = (curT - 1)->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->x) * coef * coef * ONE_THIRD);
    }
    vol = (curT->x + term2) * fullV * 0.5;
} 
break;
case TSSE_Hexahedron5:   // В-Ю-Ю-В (1001100)
{
    double denom = curA->z;                         // α^z_{i,j,k}
    double term2 = curT->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->x) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - 1)->x + term2) * fullV * 0.5;
} 
break;
case BSSW_Hexahedron5:   // Н-Ю-Ю-З (11001000)
{
    double denom = (curA - 1)->z;                   // α^z_{i-1,j,k}
    double term2 = (curT - 1)->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _NM)->x) * coef * coef * ONE_THIRD);
    }
    vol = (curT->x + term2) * fullV * 0.5;
} 
break;
case BSSE_Hexahedron5:   // Н-Ю-Ю-В (11000100)
{
    double denom = curA->z;                         // α^z_{i,j,k}
    double term2 = curT->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _NM)->x) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - 1)->x + term2) * fullV * 0.5;
} 
break;
case TNNW_Hexahedron5:   // В-С-С-З (10011)
{
    double denom = (curA - 1 - _N)->z;              // α^z_{i-1,j-1,k}
    double term2 = (curT - 1)->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _N)->x) * coef * coef * ONE_THIRD);
    }
    vol = (curT->x + term2) * fullV * 0.5;
} 
break;
case TNNE_Hexahedron5:   // В-С-С-В (100011)
{
    double denom = (curA - _N)->z;                  // α^z_{i,j-1,k}
    double term2 = curT->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _N)->x) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - 1)->x + term2) * fullV * 0.5;
} 
break;
case BNNW_Hexahedron5:   // Н-С-С-З (110001)
{
    double denom = (curA - 1 - _N)->z;              // α^z_{i-1,j-1,k}
    double term2 = (curT - 1)->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _N - _NM)->x) * coef * coef * ONE_THIRD);
    }
    vol = (curT->x + term2) * fullV * 0.5;
} 
break;
case BNNE_Hexahedron5:   // Н-С-С-В (110010)
{
    double denom = (curA - _N)->z;                  // α^z_{i,j-1,k}
    double term2 = curT->x;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - _N - _NM)->x) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - 1)->x + term2) * fullV * 0.5;
} 
break;
// Группа C: через ϑ^v (8 типов)
case BNBE_Hexahedron5:   // Н-С-Н-В (1110000)
{
    double denom = curA->x;                         // α^x_{i,j,k}
    double term2 = curT->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->y) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _N)->y + term2) * fullV * 0.5;
} 
break;
case BNBW_Hexahedron5:   // Н-С-Н-З (10110000)
{
    double denom = curA->x;                         // α^x_{i,j,k}
    double term2 = curT->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->y) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _N)->y + term2) * fullV * 0.5;
} 
break;
case BSBW_Hexahedron5:   // Н-Ю-Н-З (11010000)
{
    double denom = (curA - _N)->x;                  // α^x_{i,j-1,k}
    double term2 = (curT - _N)->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->y) * coef * coef * ONE_THIRD);
    }
    vol = (curT->y + term2) * fullV * 0.5;
} 
break;
case BSBE_Hexahedron5:   // Н-Ю-Н-В (11100000)
{
    double denom = (curA - _N)->x;                  // α^x_{i,j-1,k}
    double term2 = (curT - _N)->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->y) * coef * coef * ONE_THIRD);
    }
    vol = (curT->y + term2) * fullV * 0.5;
} 
break;
case TNTE_Hexahedron5:   // В-С-В-В (111)
{
    double denom = (curA - _NM)->x;                 // α^x_{i,j,k-1}
    double term2 = curT->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->y) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _N)->y + term2) * fullV * 0.5;
} 
break;
case TNTW_Hexahedron5:   // В-С-В-З (1011)
{
    double denom = (curA - _NM)->x;                 // α^x_{i,j,k-1}
    double term2 = curT->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->y) * coef * coef * ONE_THIRD);
    }
    vol = ((curT - _N)->y + term2) * fullV * 0.5;
} 
break;
case TSTW_Hexahedron5:   // В-Ю-В-З (1101)
{
    double denom = (curA - _N - _NM)->x;            // α^x_{i,j-1,k-1}
    double term2 = (curT - _N)->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * (curA - 1)->y) * coef * coef * ONE_THIRD);
    }
    vol = (curT->y + term2) * fullV * 0.5;
} 
break;
case TSTE_Hexahedron5:   // В-Ю-В-В (1110)
{
    double denom = (curA - _N - _NM)->x;            // α^x_{i,j-1,k-1}
    double term2 = (curT - _N)->y;
    if (denom > EQ_DOUBLE) {
        double coef = (1.0 - denom) / denom;
        term2 *= (1.0 - (1.0 + 2.0 * curA->y) * coef * coef * ONE_THIRD);
    }
    vol = (curT->y + term2) * fullV * 0.5;
} 
break;
// ========== Шестиугольные гептаэдры (8 типов) ==========
case BSW_Heptahedron6:   // Н-Ю-З (11011000)
{
    double theta2 = curT->x;
    double term = 0.0;
    double alpha_z = (curA - _N)->z;               // α^z_{i,j-1,k}
    double denom = 1.0 - alpha_z;
    if (denom > EQ_DOUBLE) {
        alpha_z /= denom;
        term = (3.0 - 2.0 * (curA - _NM)->x) * alpha_z * alpha_z;
    }
	alpha_z = (curA - _NM)->y;              // α^y_{i,j,k-1}
    denom = 1.0 - alpha_z;
    if (denom > EQ_DOUBLE) {
        alpha_z /= denom;
        term += (3.0 - 2.0 * (curA - _N)->x) * alpha_z * alpha_z;
    }

    vol = ((curT - 1)->x + theta2 + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case BSE_Heptahedron6:   // Н-Ю-В (11100100) — 1-θ^u_{i-1,j,k}
{
    double theta2 = (curT - 1)->x;
    double term = 0.0;
    double alpha = (curA - 1 - _N)->z;           // α^z_{i-1,j-1,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * (curA - _N)->x) * alpha * alpha;
    }
    alpha = (curA - 1 - _NM)->y;                // α^y_{i-1,j,k-1}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * (curA - _NM)->x) * alpha * alpha;
    }
    vol = (theta2 + curT->x + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case BNE_Heptahedron6:   // Н-С-В (1110010) — 1-θ^u_{i-1,j,k}
{
    double theta2 = (curT - 1)->x;
    double term = 0.0;
    double alpha = (curA - 1)->z;                // α^z_{i-1,j,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * curA->x) * alpha * alpha;
    }
    alpha = (curA - 1 - _NM)->y;                // α^y_{i-1,j,k-1}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha * alpha;
    }
    vol = (theta2 + curT->x + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case BNW_Heptahedron6:   // Н-С-З (10110001) — 1-θ^u_{i,j,k}
{
    double theta2 = curT->x;
    double term = 0.0;
    double alpha = curA->z;                      // α^z_{i,j,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * curA->x) * alpha * alpha;
    }
    alpha = (curA - _NM)->y;                    // α^y_{i,j,k-1}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha * alpha;
    }
    vol = ((curT - 1)->x + theta2 + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case TSW_Heptahedron6:   // В-Ю-З (10001101) — 1-θ^u_{i,j,k}
{
    double theta2 = curT->x;
    double term = 0.0;
    double alpha = (curA - _N)->z;               // α^z_{i,j-1,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha * alpha;
    }
    alpha = curA->y;                             // α^y_{i,j,k}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * curA->x) * alpha * alpha;
    }
    vol = ((curT - 1)->x + theta2 + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case TSE_Heptahedron6:   // В-Ю-В (1001110) — 1-θ^u_{i-1,j,k}
{
    double theta2 = (curT - 1)->x;
    double term = 0.0;
    double alpha = (curA - 1 - _N)->z;           // α^z_{i-1,j-1,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * (curA - _N - _NM)->x) * alpha * alpha;
    }
    alpha = (curA - 1)->y;                       // α^y_{i-1,j,k}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * curA->x) * alpha * alpha;
    }
    vol = (theta2 + curT->x + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case TNE_Heptahedron6:   // В-С-В (100111) — 1-θ^u_{i-1,j,k}
{
    double theta2 = (curT - 1)->x;
    double term = 0.0;
    double alpha = (curA - 1)->z;                // α^z_{i-1,j,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * (curA - _NM)->x) * alpha * alpha;
    }
    alpha = (curA - 1)->y;                       // α^y_{i-1,j,k}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * (curA - _N)->x) * alpha * alpha;
    }
    vol = (theta2 + curT->x + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
case TNW_Heptahedron6:   // В-С-З (11011) — 1-θ^u_{i,j,k}
{
    double theta2 = curT->x;
    double term = 0.0;
    double alpha = curA->z;                      // α^z_{i,j,k}
    double denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term = (3.0 - 2.0 * (curA - _NM)->x) * alpha * alpha;
    }
    alpha = curA->y;                             // α^y_{i,j,k}
    denom = 1.0 - alpha;
    if (denom > EQ_DOUBLE) {
        alpha /= denom;
        term += (3.0 - 2.0 * (curA - _N)->x) * alpha * alpha;
    }
    vol = ((curT - 1)->x + theta2 + (1.0 - theta2) * term * ONE_THIRD) * fullV * 0.5;
} 
break;
			}//switch
			_VOF[i] = vol + PLUS_TINY;
		}
		LS_STAG_KERNEL(vXY, LS_STAG_R_PD _VOF, LS_STAG_R_PD _weight, LS_STAG_R_PD _quarterV, LS_STAG_CR_PCT _baseType, LS_STAG_CR_P(coord<SideTypes>) _SideType)
		{
			double curWeight = 0.25;
			switch (_baseType[i])
			{
			case SolidCell:
				curWeight = 0.0;
				break;
			case NE_Pentahedron:
			case NW_Pentahedron:
			case SE_Pentahedron:
			case SW_Pentahedron:
			case TNW_Tetrahedron:
			case TNE_Tetrahedron:
			case TSE_Tetrahedron:
			case TSW_Tetrahedron:
			case BNW_Tetrahedron:
			case BNE_Tetrahedron:
			case BSE_Tetrahedron:
			case BSW_Tetrahedron:
				curWeight = ONE_THIRD;
				break;
			case B_Border:
			case T_Border:// только треугольная z-грань даёт вес 1/3
			switch(_SideType[i].z)
			{
				case NW_Triangle:
				case SW_Triangle:
				case NE_Triangle:
				case SE_Triangle:
					curWeight = ONE_THIRD;
			}
				break;
			}//switch
			_weight[i] = curWeight;
			_quarterV[i] = curWeight * _VOF[i];
		}
		LS_STAG_KERNEL(vXZ, LS_STAG_R_PD _VOF, LS_STAG_R_PD _weight, LS_STAG_R_PD _quarterV, LS_STAG_CR_PCT _baseType)
		{
			double curWeight = 0.25;
			switch (_baseType[i])
			{
			case SolidCell:
				curWeight = 0.0;
				break;
			case TE_Pentahedron:
			case TW_Pentahedron:
			case BE_Pentahedron:
			case BW_Pentahedron:
			case TNW_Tetrahedron:
			case TNE_Tetrahedron:
			case TSE_Tetrahedron:
			case TSW_Tetrahedron:
			case BNW_Tetrahedron:
			case BNE_Tetrahedron:
			case BSE_Tetrahedron:
			case BSW_Tetrahedron:
				curWeight = ONE_THIRD;
				break;
			}//switch
			_weight[i] = curWeight;
			_quarterV[i] = curWeight * _VOF[i];
		}
		LS_STAG_KERNEL(vYZ, LS_STAG_R_PD _VOF, LS_STAG_R_PD _weight, LS_STAG_R_PD _quarterV, LS_STAG_CR_PCT _baseType)
		{
			double curWeight = 0.25;
			switch (_baseType[i])
			{
			case SolidCell:
				curWeight = 0.0;
				break;
			case TN_Pentahedron:
			case TS_Pentahedron:
			case BN_Pentahedron:
			case BS_Pentahedron:
			case TNW_Tetrahedron:
			case TNE_Tetrahedron:
			case TSE_Tetrahedron:
			case TSW_Tetrahedron:
			case BNW_Tetrahedron:
			case BNE_Tetrahedron:
			case BSE_Tetrahedron:
			case BSW_Tetrahedron:
				curWeight = ONE_THIRD;
				break;
			}//switch
			_weight[i] = curWeight;
			_quarterV[i] = curWeight * _VOF[i];
		}
		LS_STAG_KERNEL_2D(mx, LS_STAG_R_PD _x_VOF, LS_STAG_CR_PD _VOF, const int _baseMeshN, const int _baseMeshM,
			const int _baseMeshNM)
		{//массы х-ых ячеек
			int num = i * _baseMeshN + j * _baseMeshNM;
			double curV = _VOF[num];
			num++;
			for (int numX = num - 1 - i - j * _baseMeshM, m_end = numX + _baseMeshN - 1; numX < m_end; num++, numX++)
			{
				double eV = _VOF[num];
				_x_VOF[numX] = 0.5 * (curV + eV);
				curV = eV;
			}
		}
		LS_STAG_KERNEL_2D(my, LS_STAG_R_PD _y_VOF, LS_STAG_CR_PD _VOF, const int _baseMeshN, const int _baseMeshNM)
		{// массы y-ых ячеек
			for (int num = i * _baseMeshN + j * _baseMeshNM, numY = num - j * _baseMeshN, m_end = num + _baseMeshN; 
				num < m_end; numY++, num++)
			_y_VOF[numY] = 0.5 * (_VOF[num] + _VOF[num + _baseMeshN]);
		}
		LS_STAG_KERNEL(mz, LS_STAG_R_PD _z_VOF, LS_STAG_CR_PD _VOF, const int _NM)
		{// массы z-ых ячеек
			_z_VOF[i] = 0.5 * (_VOF[i] + _VOF[i + _NM]);
		}
		LS_STAG_KERNEL_2D(mXY, LS_STAG_R_PD _VOF, LS_STAG_CR_PD _quarterV, const int _baseMeshN, const int NM2, const int _baseMeshNM)
		{
			int num = i * _baseMeshN + j * _baseMeshNM;
			int numN = num + _baseMeshN;
			double quarterNW = _quarterV[num];
			double quarterSW = _quarterV[numN];
			num++; numN++;
			for (int numXY = num - 1 - i - j * NM2, m_end = numXY + _baseMeshN - 1; numXY < m_end; numXY++, num++, numN++)
			{
				double quarterNE = _quarterV[num];
				double quarterSE = _quarterV[numN];
				_VOF[numXY] = quarterNE + quarterNW + quarterSE + quarterSW;
				quarterNW = quarterNE;
				quarterSW = quarterSE;
			}
		}
		LS_STAG_KERNEL_2D(mXZ, LS_STAG_R_PD _VOF, LS_STAG_CR_PD _quarterV, const int _baseMeshN, const int _baseMeshM,
			const int _baseMeshNM)
		{
			int num = i * _baseMeshN + j * _baseMeshNM;
			int numN = num + _baseMeshNM;
			double quarterNW = _quarterV[num];
			double quarterSW = _quarterV[numN];
			num++; numN++;
			for (int numX = num - 1 - i - j * _baseMeshM, m_end = numX + _baseMeshN - 1; numX < m_end; numX++, num++, numN++)
			{
				double quarterNE = _quarterV[num];
				double quarterSE = _quarterV[numN];
				_VOF[numX] = quarterNE + quarterNW + quarterSE + quarterSW;
				quarterNW = quarterNE;
				quarterSW = quarterSE;
			}
		}
		LS_STAG_KERNEL_2D(mYZ, LS_STAG_R_PD _VOF, LS_STAG_CR_PD _quarterV, const int _baseMeshN, const int _baseMeshNM)
		{
			for (int num = i * _baseMeshN + j * _baseMeshNM, NM2 = _baseMeshNM + _baseMeshN, m_end = num + _baseMeshN,
				numY = num - j * _baseMeshN; num < m_end; numY++, num++)// curTypeXY();
			{
				_VOF[numY] = _quarterV[num] + _quarterV[num + _baseMeshN] + _quarterV[num + _baseMeshNM] + _quarterV[num + NM2];
			}
		}
		LS_STAG_KERNEL(projectionVib, const double fullV, LS_STAG_CR_PCT _baseType, LS_STAG_CR_P(coord<SideTypes>) _SideType, LS_STAG_CR_PC _delta, LS_STAG_R_PC _alpha,
			LS_STAG_R_PC _v_ib, LS_STAG_R_PC _v_ib_yx, LS_STAG_R_PC _v_ib_zy, LS_STAG_R_PC _v_ib_ENT, const int _N)
		{// вращение вокруг Oz 
			const CellTypes curType = _baseType[i];
			bool isCutCell = false;
			if ((curType > SolidCell) && (curType < FluidCell))
				isCutCell = true;
			else
			{
				switch (curType)
				{
					case B_Border:
					case T_Border:
					{
						const SideTypes zSide = _SideType[i].z;
						if ((zSide > SolidSide) && (zSide < FluidSide)) isCutCell = true;
					}
				}
			}
			if (isCutCell)
			{
				coordinate* curA = _alpha + i;
				const coordinate* curD = _delta + i;
				const double sinb = ((curA - _N)->x - curA->x) * curD->x;
				const double cosb = ((curA - 1)->y - curA->y) * curD->y;
				const double len = sqrt(sinb * sinb + cosb * cosb);
				const double sina = (len > EQ_DOUBLE) ? (fabs(cosb) / len) : 0.0;
				const double cosa = (len > EQ_DOUBLE) ? sqrt(1.0 - sina * sina) : 1.0;// 15.09.12
				const double fs = (cosb > EQ_DOUBLE) ? (-fullV) : ((cosb < MINUS_EQ_DOUBLE) ? fullV : 0.0);// sign_vy * fullV
				const double fw = (sinb > EQ_DOUBLE) ? fullV : ((sinb < MINUS_EQ_DOUBLE) ? (-fullV) : 0.0);// sign_vx * fullV
				coordinate curVel;
				curVel.set(cosa * fw, sina * fs, 0.0);
				_v_ib[i].set(&curVel);
				_v_ib_yx[i].set(&curVel);//18.05.14
				_v_ib_zy[i].set(&curVel);//18.05.14
				_v_ib_ENT[i].set(&curVel);//18.05.14
			}
		};
		LS_STAG_KERNEL_INLINE void checkOscLS(const int i, LS_STAG_R_PD _phi, double& f, const double fwCur, const double fw2)
		{//20.07.12 --- исключаем осциллирующее поведение функции уровня
		//29.09.2026: f<=EQ_DOUBLE --- не осцилляция, а значение на интерфейсе (вырожденный случай:
		//окружность проходит точно через узел). Усреднение соседей сдвигало узел на ~шаг,
		//ломая alpha/theta/VOF в соседних трапециях. См. TA (P_C=0), пары W<->E, S<->N.
			if ((f > EQ_DOUBLE) && (fwCur < 0) && (fw2 < 0))
			{
				f = 0.5 * (fwCur + fw2);
				_phi[i] = f;
			}
		};
		LS_STAG_KERNEL_INLINE void bottomBuildAlpha(coordinate*& curTLoc, LS_STAG_R_PC _alpha, LS_STAG_R_PD _phi, const int _N, const int _NM)
		{
			const int num = (int)(curTLoc - _alpha);
			double f = _phi[num];
			const double fw = _phi[num - 1];
			const double fs = _phi[num - _N];
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fw, _phi[num + 1]);
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fs, _phi[num + _N]);
			const double tx = (f * fw < 0) ? ((f > 0) ? (fw / (fw - f)) : (f / (f - fw))) : ((f > -fw) ? 0.0 : 1.0);
			const double ty = (f * fs < 0) ? ((f > 0) ? (fs / (fs - f)) : (f / (f - fs))) : ((f > -fs) ? 0.0 : 1.0);
			const double tz = (f > 0.0) ? 0.0 : 1.0;// как BuildAlpha при fb := f
			// todo -- можно здесь же и скорректировать, тк оно не параллелится
			curTLoc->set(tx, ty, tz);// alpha
			curTLoc++;
		};
		LS_STAG_KERNEL_INLINE void topBuildAlpha(coordinate*& curTLoc, LS_STAG_R_PC _alpha, LS_STAG_R_PD _phi, const int _N, const int _NM)
		{
			const int num = (int)(curTLoc - _alpha);
			double f = _phi[num];
			const double fw = _phi[num - 1];
			const double fs = _phi[num - _N];
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fw, _phi[num + 1]);
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fs, _phi[num + _N]);
			const double tx = (f * fw < 0) ? ((f > 0) ? (fw / (fw - f)) : (f / (f - fw))) : ((f > -fw) ? 0.0 : 1.0);
			const double ty = (f * fs < 0) ? ((f > 0) ? (fs / (fs - f)) : (f / (f - fs))) : ((f > -fs) ? 0.0 : 1.0);
			const double tz = (f > 0.0) ? 0.0 : 1.0;// как BuildAlpha при fb := f
			// todo -- можно здесь же и скорректировать, тк оно не параллелится
			curTLoc->set(tx, ty, tz);// alpha
			curTLoc++;
		};
		LS_STAG_KERNEL_INLINE void BuildAlpha(coordinate*& curTLoc, LS_STAG_R_PC _alpha, LS_STAG_R_PD _phi, const int _N, const int _NM)
		{
			const int num = (int)(curTLoc - _alpha);
			double f = _phi[num];
			const double fw = _phi[num - 1];
			const double fs = _phi[num - _N];
			const double fb = _phi[num - _NM];
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fw, _phi[num + 1]);
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fs, _phi[num + _N]);
			LSSTAGUniformMesh::checkOscLS(num, _phi, f, fb, _phi[num + _NM]);
			const double tx = (f * fw < 0) ? ((f > 0) ? (fw / (fw - f)) : (f / (f - fw))) : ((f > -fw) ? 0.0 : 1.0);
			const double ty = (f * fs < 0) ? ((f > 0) ? (fs / (fs - f)) : (f / (f - fs))) : ((f > -fs) ? 0.0 : 1.0);
			const double tz = (f * fb < 0) ? ((f > 0) ? (fb / (fb - f)) : (f / (f - fb))) : ((f > -fb) ? 0.0 : 1.0);
            // todo -- можно здесь же и скорректировать, тк оно не параллелится
			curTLoc->set(tx, ty, tz);// alpha
			curTLoc++;
		};
		LS_STAG_KERNEL_INLINE void DefineTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const CellTypes type)
		{
			_baseType[num] = type;
			num++;
		};
		LS_STAG_KERNEL_INLINE SideTypes TypeOfSide(const double phi, const double phiW, 
			const double phiSW, const double phiS)//  определение типов граней
		{
			int number = (phiS  < EQ_DOUBLE) ? 1 : 0;
			if (phiSW < EQ_DOUBLE) number += 10;
			if (phiW < EQ_DOUBLE) number += 100;
			if (phi < EQ_DOUBLE) number += 1000;
			return (SideTypes)number;
		};
		LS_STAG_KERNEL_INLINE void calcInvDeltaIb(const int num, LS_STAG_R_PC invDib, 
			const CellTypes curCellType, const SideTypes zSide, LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM)
		{
			double y = 0.0, z = 0.0, x = 0.0;
			switch(curCellType)
			{
				case B_Border:
				case T_Border:
				{
					if ((zSide < FluidSide) && (zSide > SolidSide))
					{
						switch (zSide)
						{// Triangle -> Pentahedron (только x/y), Trapezoid -> Hexahedron4 (только своя компонента), Pentagon -> нет кейсов (нули)
							case NW_Triangle:
							{// как NW_Pentahedron (x/y)
    coordinate* curA = _alpha + num - _N;
					coordinate* curD = _delta + num;
    double alpha = curA->x;
    double aux = (curA - _NM)->x;    
    double sum = alpha * (alpha + aux) + aux * aux;
    x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
	curA += _N;
			alpha = curA->y;
            aux = (curA - _NM)->y;
			sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
							}
								break;
							case NE_Triangle:
							{// как NE_Pentahedron (x/y)
						coordinate* curD = _delta + num;
								coordinate* pA = _alpha + num - _N;
								double alpha = pA->x, aux = (pA - _NM)->x;
								double sum = alpha * (alpha + aux) + aux * aux;
								x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
								pA--;
								alpha = pA->y; aux = (pA - _NM)->y;
								sum = alpha * (alpha + aux) + aux * aux;
								y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
							}
								break;
							case SW_Triangle:
							{// как SW_Pentahedron (x/y)
						coordinate* curA = _alpha + num;
						coordinate* curD = _delta + num;
								double alpha = curA->x, aux = (curA - _NM)->x;
								double sum = alpha * (alpha + aux) + aux * aux;
								x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
								alpha = curA->y; aux = (curA - _NM)->y;
								sum = alpha * (alpha + aux) + aux * aux;
								y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
							}
								break;
							case SE_Triangle:
							{// как SE_Pentahedron (x/y)
						coordinate* curA = _alpha + num;
						coordinate* curD = _delta + num;
								double alpha = curA->x, aux = (curA - _NM)->x;
								double sum = alpha * (alpha + aux) + aux * aux;
								x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
								curA--;
								alpha = curA->y; aux = (curA - _NM)->y;
								sum = alpha * (alpha + aux) + aux * aux;
								y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
							}
								break;
							case W_Trapezoid:
							{// как W_Hexahedron4 (только x)
						coordinate* curA = _alpha + num;
								double alpha1 = curA->x, alpha2 = (curA - _N - _NM)->x;
								double sum = alpha1 + alpha2;
								double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - _NM)->x * (curA - _N)->x;
								if (denom > EQ_DOUBLE) x = 6.0 * sum / (denom * _delta[num].x);
							}
								break;
							case E_Trapezoid:
							{// как E_Hexahedron4 (только x)
						coordinate* curA = _alpha + num;
								double alpha1 = curA->x, alpha2 = (curA - _N - _NM)->x;
								double sum = alpha1 + alpha2;
								double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - _NM)->x * (curA - _N)->x;
								if (denom > EQ_DOUBLE) x = 6.0 * sum / (denom * _delta[num].x);
							}
								break;
							case N_Trapezoid:
							case S_Trapezoid:
							{// как N/S_Hexahedron4 (только y)
						coordinate* curA = _alpha + num;
								double alpha1 = curA->y, alpha2 = (curA - 1 - _NM)->y;
								double sum = alpha1 + alpha2;
								double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - 1)->y * (curA - _NM)->y;
								if (denom > EQ_DOUBLE) y = 6.0 * sum / (denom * _delta[num].y);
							}
								break;
						}
					}
				}
				break;
				case BSW_Tetrahedron:
				{
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
					double aux = curA->x;
					if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
					aux = curA->y;
					if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
					aux = curA->z;
					if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
				}
				break;
				 case BSE_Tetrahedron: // Н-Ю-В
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
            double aux = curA->x;
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
			curA--;
            aux = curA->y;      // α^y_{i-1,j,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = curA->z;      // α^z_{i-1,j,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case BNW_Tetrahedron: // Н-С-З
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
            double aux = curA->y;                 // α^y_{i,j,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
			curA -= _N;
           aux = curA->x;   // α^x_{i,j-1,k}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
            aux = curA->z;          // α^z_{i,j-1,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case BNE_Tetrahedron: // Н-С-В
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
            double aux = (curA - _N)->x;   // α^x_{i,j-1,k}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
			curA--;
            aux = curA->y;          // α^y_{i-1,j,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = (curA - _N)->z;     // α^z_{i-1,j-1,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case TSW_Tetrahedron: // В-Ю-З
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
            double aux = (curA - _NM)->x;  // α^x_{i,j,k-1}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
            aux = curA->y;                 // α^y_{i,j,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = curA->z;                 // α^z_{i,j,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case TSE_Tetrahedron: // В-Ю-В
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num;
            double aux = (curA - _NM)->x;  // α^x_{i,j,k-1}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
			curA--;
            aux = curA->y;          // α^y_{i-1,j,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = (curA - _NM)->z;    // α^z_{i-1,j,k-1}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case TNW_Tetrahedron: // В-С-З
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num - _N;
            double aux = (curA  - _NM)->x; // α^x_{i,j-1,k-1}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
            aux = curA->y;              // α^y_{i,j-1,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = curA->z;              // α^z_{i,j-1,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        case TNE_Tetrahedron: // В-С-В
        {
					coordinate* curD = _delta + num;
					coordinate* curA = _alpha + num - _N;
            double aux = (curA - _NM)->x; // α^x_{i,j-1,k-1}
            if(aux > EQ_DOUBLE) x = 3.0 / (aux * curD->x);
			curA--;
            aux = curA->y;          // α^y_{i-1,j-1,k}
            if(aux > EQ_DOUBLE) y = 3.0 / (aux * curD->y);
            aux = curA->z;          // α^z_{i-1,j-1,k}
            if(aux > EQ_DOUBLE) z = 3.0 / (aux * curD->z);
        }
            break;
        // ========== 4-угольные гексаэдры ==========
        case E_Hexahedron4: // В (восток)
        case W_Hexahedron4: // З (запад)
        {
					coordinate* curA = _alpha + num;
            double alpha1 = curA->x;                // α^x_{i,j,k}
            double alpha2 = (curA - _N - _NM)->x;   // α^x_{i,j-1,k-1}
            double sum = alpha1 + alpha2;
			double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - _NM)->x * (curA - _N)->x;
            if (denom > EQ_DOUBLE) x = 6.0 * sum / (denom * _delta[num].x);
        }
            break;
        case N_Hexahedron4: // С
        case S_Hexahedron4: // Ю
        {
					coordinate* curA = _alpha + num;
            double alpha1 = curA->y;                // α^y_{i,j,k}
            double alpha2 = (curA - 1 - _NM)->y;   // α^y_{i-1,j,k-1}
            double sum = alpha1 + alpha2;
			double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - 1)->y * (curA - _NM)->y;
            if (denom > EQ_DOUBLE) y = 6.0 * sum / (denom * _delta[num].y);
            } 
            break;
        case T_Hexahedron4: // Верх
        case B_Hexahedron4: // Низ
        {
					coordinate* curA = _alpha + num;
            double alpha1 = curA->z;                // α^z_{i,j,k}
            double alpha2 = (curA - 1 - _N)->z;     // α^z_{i-1,j-1,k}
            double sum = alpha1 + alpha2;
			double denom = 2.0 * sum * sum - alpha1 * alpha2 - (curA - 1)->z * (curA - _N)->z;
            if (denom > EQ_DOUBLE) z = 6.0 * sum / (denom * _delta[num].z);
            }
            break;
        // ========== Пентаэдры ==========
        case SW_Pentahedron: // Ю-З
        {
					coordinate* curA = _alpha + num;
					coordinate* prevA = curA - _NM;
					coordinate* curD = _delta + num;
            double alpha = curA->x;
            double aux = prevA->x;
			double sum = alpha * (alpha + aux) + aux * aux;
            x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
			alpha = curA->y;
            aux = prevA->y;
			sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
        }
        break;
        case SE_Pentahedron: // Ю-В
        {
					coordinate* curA = _alpha + num;
            double alpha = curA->x;
					coordinate* curD = _delta + num;
            double aux = (curA - _NM)->x; // α^x_{i,j,k-1}
			double sum = alpha * (alpha + aux) + aux * aux;
            x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
			curA--;
			alpha = curA->y;
            aux = (curA - _NM)->y;
			sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
        }
        break;
		case NW_Pentahedron: // Северо-Запад
{
    coordinate* curA = _alpha + num - _N;
					coordinate* curD = _delta + num;
    double alpha = curA->x;
    double aux = (curA - _NM)->x;    
    double sum = alpha * (alpha + aux) + aux * aux;
    x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
	curA += _N;
			alpha = curA->y;
            aux = (curA - _NM)->y;
			sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
}
break;
		case NE_Pentahedron: // Северо-Восток
{
    coordinate* curA = _alpha + num - _N;
					coordinate* curD = _delta + num;
    double alpha = curA->x;
    double aux = (curA - _NM)->x;    
    double sum = alpha * (alpha + aux) + aux * aux;
    x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
	curA--;
			alpha = curA->y;
            aux = (curA - _NM)->y;
			sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
}
break;
        case BW_Pentahedron: 
        {
					coordinate* curA = _alpha + num;
					coordinate* prevA = curA - _N;
					coordinate* curD = _delta + num;
            double alpha = curA->x;
            double aux = prevA->x;
			double sum = alpha * (alpha + aux) + aux * aux;
            x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
			alpha = curA->z;
            aux = prevA->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
        }
        break;
        case BE_Pentahedron: 
        {
					coordinate* curA = _alpha + num;
					coordinate* curD = _delta + num;
            double alpha = curA->x;
            double aux = (curA - _N)->x;
			double sum = alpha * (alpha + aux) + aux * aux;
            x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
			curA--;
			alpha = curA->z;
            aux = (curA - _N)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
        }
        break;
		case TW_Pentahedron:
{
    coordinate* curA = _alpha + num - _NM;
					coordinate* curD = _delta + num;
    double alpha = curA->x;
    double aux = (curA - _N)->x;    
    double sum = alpha * (alpha + aux) + aux * aux;
    x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
	curA += _NM;
			alpha = curA->z;
            aux = (curA - _N)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
}
break;
		case TE_Pentahedron: 
{
    coordinate* curA = _alpha + num - _NM;
					coordinate* curD = _delta + num;
    double alpha = curA->x;
    double aux = (curA - _N)->x;    
    double sum = alpha * (alpha + aux) + aux * aux;
    x = 4.0 * sum / ((1.0 + alpha * sum) * curD->x);
	curA += (_NM - 1);
			alpha = curA->z;
            aux = (curA - _N)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
}
break;
        case BS_Pentahedron: 
        {
					coordinate* curA = _alpha + num;
					coordinate* prevA = curA - 1;
					coordinate* curD = _delta + num;
            double alpha = curA->y;
            double aux = prevA->y;
			double sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
			alpha = curA->z;
            aux = prevA->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
        }
        break;
        case BN_Pentahedron: 
        {
					coordinate* curA = _alpha + num;
					coordinate* curD = _delta + num;
            double alpha = curA->y;
            double aux = (curA - 1)->y;
			double sum = alpha * (alpha + aux) + aux * aux;
            y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
			curA -= _N;
			alpha = curA->z;
            aux = (curA - 1)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
        }
        break;
		case TS_Pentahedron:
{
    coordinate* curA = _alpha + num - _NM;
					coordinate* curD = _delta + num;
    double alpha = curA->y;
    double aux = (curA - 1)->y;    
    double sum = alpha * (alpha + aux) + aux * aux;
    y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
	curA += _NM;
			alpha = curA->z;
            aux = (curA - 1)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
}
break;
		case TN_Pentahedron: 
{
    coordinate* curA = _alpha + num - _NM;
					coordinate* curD = _delta + num;
    double alpha = curA->y;
    double aux = (curA - 1)->y;    
    double sum = alpha * (alpha + aux) + aux * aux;
    y = 4.0 * sum / ((1.0 + alpha * sum) * curD->y);
	curA += (_NM - _N);
			alpha = curA->z;
            aux = (curA - 1)->z;
			sum = alpha * (alpha + aux) + aux * aux;
            z = 4.0 * sum / ((1.0 + alpha * sum) * curD->z);
}
break;
			}
			invDib[num].set(x, y, z);
		}
		LS_STAG_KERNEL_INLINE void calcInvBaseDeltaIb(const int num, LS_STAG_R_PC invDibNE,
			LS_STAG_R_PC invDibTN, const CellTypes curCellType, const coord<SideTypes>& curBaseSideType,
			LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM)
		{
			if((curCellType > N_Border) && (curCellType < FluidCell) && (curCellType != SolidCell))
			{
				double yE = 0.0, zE = 0.0, xN = 0.0, zN = 0.0, xT = 0.0, yT = 0.0;
			const int numB = (num < _NM) ? num : (num - _NM);// слой снизу; для k=0 — продление поля вниз
			switch (curBaseSideType.x)
			{
			case FluidSide:
				{
					coordinate* curA = _alpha + num;
					coordinate* curD = _delta + num;
					if ((curA->y < EQ_DOUBLE) || (_alpha[numB].y < EQ_DOUBLE)) zE = 2.0 / curD->z;
					if ((curA->z < EQ_DOUBLE) || ((curA - _N)->z < EQ_DOUBLE)) yE = 2.0 / curD->y;
				}
				break;
			case NE_Triangle:
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[numB].y;
				if (aux > EQ_DOUBLE) yE = 3.0 / (aux * curD->y);
				aux = _alpha[num - _N].z;
				if (aux > EQ_DOUBLE) zE = 3.0 / (aux * curD->z);
			}
			break;
			case NW_Triangle:   // С-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[numB].y;   // α^y_{i,j,k-1}
				if (aux > EQ_DOUBLE) yE = 3.0 / (aux * curD->y);
				aux = _alpha[num].z;         // α^z_{i,j,k}
				if (aux > EQ_DOUBLE) zE = 3.0 / (aux * curD->z);
				break;
			}
			case SE_Triangle:   // Ю-В треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].y;         // α^y_{i,j,k}
				if (aux > EQ_DOUBLE) yE = 3.0 / (aux * curD->y);
				aux = _alpha[num - _N].z;    // α^z_{i,j-1,k}
				if (aux > EQ_DOUBLE) zE = 3.0 / (aux * curD->z);
				break;
			}
			case SW_Triangle:   // Ю-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].y;         // α^y_{i,j,k}
				if (aux > EQ_DOUBLE) yE = 3.0 / (aux * curD->y);
				aux = _alpha[num].z;         // α^z_{i,j,k}
				if (aux > EQ_DOUBLE) zE = 3.0 / (aux * curD->z);
				break;
			}
			case W_Trapezoid:   // В/З трапеция на восточной грани
			case E_Trapezoid:   // В/З трапеция на восточной грани
			{
				double alpha1 = _alpha[numB].y;    // α^y_{i,j,k-1}
				double alpha2 = _alpha[num].y;          // α^y_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) yE = 3.0 * sum_alpha / (denom * _delta[num].y);
			}
			break;
			case S_Trapezoid:   // Ю трапеция на восточной грани
			case N_Trapezoid:   // С трапеция на восточной грани
			{ // Для Δz^{ib,e} – аналогично, с α^z
				double alpha1 = _alpha[num - _N].z;      // α^z_{i,j-1,k}
				double alpha2 = _alpha[num].z;           // α^z_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) zE = 3.0 * sum_alpha / (denom * _delta[num].z);
			}
			break;
			}
			switch (curBaseSideType.y)
			{
			case FluidSide:
				{
					coordinate* curA = _alpha + num;
					coordinate* curD = _delta + num;
					if ((curA->x < EQ_DOUBLE) || (_alpha[numB].x < EQ_DOUBLE)) zN = 2.0 / curD->z;
					if ((curA->z < EQ_DOUBLE) || ((curA - 1)->z < EQ_DOUBLE)) xN = 2.0 / curD->x;
				}
				break;
			case NE_Triangle:  // С-В треугольник (для Δx и Δz)
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[numB].x;          // α^x_{i,j,k-1}
				if (aux > EQ_DOUBLE) xN = 3.0 / (aux * curD->x);
				aux = _alpha[num - 1].z;                   // α^z_{i-1,j,k}
				if (aux > EQ_DOUBLE) zN = 3.0 / (aux * curD->z);
			}
			break;
			case NW_Triangle:  // С-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[numB].x;          // α^x_{i,j,k-1}
				if (aux > EQ_DOUBLE) xN = 3.0 / (aux * curD->x);
				aux = _alpha[num].z;                       // α^z_{i,j,k}
				if (aux > EQ_DOUBLE) zN = 3.0 / (aux * curD->z);
			}
			break;
			case SE_Triangle:  // Ю-В треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].x;                // α^x_{i,j,k}
				if (aux > EQ_DOUBLE) xN = 3.0 / (aux * curD->x);
				aux = _alpha[num - 1].z;                   // α^z_{i-1,j,k}
				if (aux > EQ_DOUBLE) zN = 3.0 / (aux * curD->z);
			}
			break;
			case SW_Triangle:  // Ю-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].x;                // α^x_{i,j,k}
				if (aux > EQ_DOUBLE) xN = 3.0 / (aux * curD->x);
				aux = _alpha[num].z;                       // α^z_{i,j,k}
				if (aux > EQ_DOUBLE) zN = 3.0 / (aux * curD->z);
			}
			break;
			case W_Trapezoid:  // В/З трапеция (для Δx)
			case E_Trapezoid:
			{
				double alpha1 = _alpha[numB].x;       // α^x_{i,j,k-1}
				double alpha2 = _alpha[num].x;             // α^x_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) xN = 3.0 * sum_alpha / (denom * _delta[num].x);
			}
			break;
			case S_Trapezoid:  // С/Ю трапеция (для Δz)
			case N_Trapezoid:
			{
				double alpha1 = _alpha[num - 1].z;         // α^z_{i-1,j,k}
				double alpha2 = _alpha[num].z;             // α^z_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) zN = 3.0 * sum_alpha / (denom * _delta[num].z);
			}
			break;
			}
			switch (curBaseSideType.z)
			{
			case FluidSide:
			{
					coordinate* curA = _alpha + num;
					coordinate* curD = _delta + num;
					if ((curA->x < EQ_DOUBLE) || ((curA - _N)->x < EQ_DOUBLE)) yT = 2.0 / curD->y;
					if ((curA->y < EQ_DOUBLE) || ((curA - 1)->y < EQ_DOUBLE)) xT = 2.0 / curD->x;
				}
				break;
			case NE_Triangle:  // С-В треугольник (для Δx и Δy)
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num - _N].x;           // α^x_{i,j-1,k}
				if (aux > EQ_DOUBLE) xT = 3.0 / (aux * curD->x);
				aux = _alpha[num - 1].y;                   // α^y_{i-1,j,k}
				if (aux > EQ_DOUBLE) yT = 3.0 / (aux * curD->y);
			}
			break;
			case NW_Triangle:  // С-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num - _N].x;           // α^x_{i,j-1,k}
				if (aux > EQ_DOUBLE) xT = 3.0 / (aux * curD->x);
				aux = _alpha[num].y;                       // α^y_{i,j,k}
				if (aux > EQ_DOUBLE) yT = 3.0 / (aux * curD->y);
			}
			break;
			case SE_Triangle:  // Ю-В треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].x;                // α^x_{i,j,k}
				if (aux > EQ_DOUBLE) xT = 3.0 / (aux * curD->x);
				aux = _alpha[num - 1].y;                   // α^y_{i-1,j,k}
				if (aux > EQ_DOUBLE) yT = 3.0 / (aux * curD->y);
			}
			break;
			case SW_Triangle:  // Ю-З треугольник
			{
				coordinate* curD = _delta + num;
				double aux = _alpha[num].x;                // α^x_{i,j,k}
				if (aux > EQ_DOUBLE) xT = 3.0 / (aux * curD->x);
				aux = _alpha[num].y;                       // α^y_{i,j,k}
				if (aux > EQ_DOUBLE) yT = 3.0 / (aux * curD->y);
			}
			break;
			case W_Trapezoid:  // В/З трапеция (для Δx)
			case E_Trapezoid:
			{
				double alpha1 = _alpha[num - _N].x;        // α^x_{i,j-1,k}
				double alpha2 = _alpha[num].x;             // α^x_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) xT = 3.0 * sum_alpha / (denom * _delta[num].x);
			}
			break;
			case S_Trapezoid:  // С/Ю трапеция (для Δy)
			case N_Trapezoid:
			{
				double alpha1 = _alpha[num - 1].y;         // α^y_{i-1,j,k}
				double alpha2 = _alpha[num].y;             // α^y_{i,j,k}
				double sum_alpha = alpha1 + alpha2;
				double denom = (alpha1 * sum_alpha + alpha2 * alpha2);
				if (denom > EQ_DOUBLE) yT = 3.0 * sum_alpha / (denom * _delta[num].y);
			}
			break;
			}
			invDibNE[num].set(xN, yE, zE);
			invDibTN[num].set(xT, yT, zN);
		}
			else 
{
invDibNE[num].set(0.0);
invDibTN[num].set(0.0);
}
		}
		LS_STAG_KERNEL(buildInvDeltaIb, LS_STAG_R_PC invDib, LS_STAG_CR_PCT _CellType, LS_STAG_CR_P(coord<SideTypes>) _SideType, 
			LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM)
		{
			LSSTAGUniformMesh::calcInvDeltaIb(i, invDib, _CellType[i], _SideType[i].z, _alpha, _delta, _N, _NM);
		}
		LS_STAG_KERNEL(buildInvBaseDeltaIb, LS_STAG_R_PC invDibNE, LS_STAG_R_PC invDibTN, LS_STAG_CR_PCT _CellType,
			LS_STAG_CR_P(coord<SideTypes>) _SideType,
			LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM)
		{
			LSSTAGUniformMesh::calcInvBaseDeltaIb(i, invDibNE, invDibTN, _CellType[i], _SideType[i], _alpha, _delta, _N, _NM);
		}
		LS_STAG_KERNEL_INLINE void bottomTypeOfCell(const int numLoc, LS_STAG_R_PD _phi, LS_STAG_R_PCT _baseType,
			LS_STAG_R_P(coord<SideTypes>) _baseSideType, const int _N, const int _NM)//  определение типов ячеек  //
		{
			const double phi = _phi[numLoc];
			const double phiW = _phi[numLoc - 1];
			const int numS = numLoc - _N;
			const double phiS = _phi[numS];
			SideTypes typeT = LSSTAGUniformMesh::TypeOfSide(phi, phiW, _phi[numS - 1], phiS);
			// нижний сосед за границей области отсутствует — продлеваем поле (phi_B := phi, phi_SB := phi_S)
			SideTypes typeE = LSSTAGUniformMesh::TypeOfSide(phi, phiS, phiS, phi);
			SideTypes typeN = LSSTAGUniformMesh::TypeOfSide(phi, phiW, phiW, phi);
			(_baseSideType + numLoc)->set(typeE, typeN, typeT);
			_baseType[numLoc] = (typeT == SolidSide) ? SolidCell : B_Border;
		};
		LS_STAG_KERNEL_INLINE void topTypeOfCell(const int numLoc, LS_STAG_R_PD _phi, LS_STAG_R_PCT _baseType,
			LS_STAG_R_P(coord<SideTypes>) _baseSideType, const int _N, const int _NM)//  определение типов ячеек  //
		{
			const double phi = _phi[numLoc];
			const double phiW = _phi[numLoc - 1];
			const int numS = numLoc - _N;
			const double phiS = _phi[numS];
			SideTypes typeB = LSSTAGUniformMesh::TypeOfSide(phi, phiW, _phi[numS - 1], phiS);
			// верхний сосед за границей области отсутствует — продлеваем поле (phi_T := phi, phi_ST := phi_S)
			SideTypes typeE = LSSTAGUniformMesh::TypeOfSide(phi, phiS, phiS, phi);
			SideTypes typeN = LSSTAGUniformMesh::TypeOfSide(phi, phiW, phiW, phi);
			(_baseSideType + numLoc)->set(typeE, typeN, typeB);
			_baseType[numLoc] = (typeB == SolidSide) ? SolidCell : T_Border;
		};
		LS_STAG_KERNEL_INLINE bool TypeOfCell(const int numLoc, LS_STAG_R_PD _phi, LS_STAG_R_PCT _baseType, 
			LS_STAG_R_P(coord<SideTypes>) _baseSideType, const int _N, const int _NM)//  определение типов ячеек  //
		{
			const double phi = _phi[numLoc];
			const double phiW = _phi[numLoc - 1];
			const int numS = numLoc - _N;
			const int numS2 = numS - 1;
			const int numB = numLoc - _NM;
			const double phiS = _phi[numS];
			const double phiSB = _phi[numS - _NM];
			const double phiB = _phi[numB];
			const double phiBW = _phi[numB - 1];
			SideTypes typeT = LSSTAGUniformMesh::TypeOfSide(phi, phiW, _phi[numS2], phiS);
			SideTypes typeE = LSSTAGUniformMesh::TypeOfSide(phi, phiS, phiSB, phiB);
			SideTypes typeN = LSSTAGUniformMesh::TypeOfSide(phi, phiW, phiBW, phiB);
			(_baseSideType + numLoc)->set(typeE, typeN, typeT);
			int number = ((int)typeT) * 10000;
			typeT = LSSTAGUniformMesh::TypeOfSide(phiB, phiBW, _phi[numS2 - _NM], phiSB);
			number += ((int)typeT);
			const CellTypes c = (CellTypes)number;
			_baseType[numLoc] = c;
			if ((c < FluidCell) && ((c > SolidCell))) return true;
			return false;
		};
		LS_STAG_KERNEL_INLINE void fluidSideRange(coord<SideTypes>*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(FluidSide);
		};
		LS_STAG_KERNEL_INLINE void fluidSides(coord<SideTypes>*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
                LSSTAGUniformMesh::fluidSideRange(curT, _N);
				curT->set(SolidSide, FluidSide, FluidSide);
			}
		};
		LS_STAG_KERNEL_INLINE void nordSidesRange(coord<SideTypes>*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(FluidSide, SolidSide, FluidSide);
		};
		LS_STAG_KERNEL_INLINE void solidSidesRange(coord<SideTypes>*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(SolidSide);
		};
		LS_STAG_KERNEL_INLINE void nordSides(coord<SideTypes>*& curT, const int _N)
		{
            LSSTAGUniformMesh::nordSidesRange(curT, _N);
			curT->set(SolidSide, SolidSide, FluidSide);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void fluidSidesPlanes(coord<SideTypes>*& curT, const int _N, const int _M, 
			const int zf, const int zs = 0)
		{
			for (int k = zs; k < zf; k++)
			{
				LSSTAGUniformMesh::fluidSides(curT, _N, _M);
				LSSTAGUniformMesh::nordSides(curT, _N);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidSidesF(coord<SideTypes>*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs, mN = _N + 1; j < vf; j++)
			{
				LSSTAGUniformMesh::fluidSideRange(curT, mN);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidSidesPlanesF(coord<SideTypes>*& curT, const int _N, const int _M,
			const int zf, const int zs = 0)
		{
			for (int k = zs, M2 = _M + 1; k < zf; k++)
			{
				LSSTAGUniformMesh::fluidSidesF(curT, _N, M2);
			}
		};
		LS_STAG_KERNEL_INLINE void topSidesRange(coord<SideTypes>*& curT, const int _N, const int hs = 0)
		{
				for (int i = hs; i < _N; i++, curT++) curT->set(FluidSide, FluidSide, SolidSide);
		};
		LS_STAG_KERNEL_INLINE void topSides(coord<SideTypes>*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
				LSSTAGUniformMesh::topSidesRange(curT, _N);
				curT->set(SolidSide, FluidSide, SolidSide);
			}
		};
		LS_STAG_KERNEL_INLINE void topNordSidesRange(coord<SideTypes>*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(FluidSide, SolidSide, SolidSide);
		};
		LS_STAG_KERNEL_INLINE void topNordSides(coord<SideTypes>*& curT, const int _N)
		{
            LSSTAGUniformMesh::topNordSidesRange(curT, _N);
			curT->set(SolidSide);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void topSidesPlane(coord<SideTypes>*& curT, const int _N, const int _M)
		{
				LSSTAGUniformMesh::topSides(curT, _N, _M);
				LSSTAGUniformMesh::topNordSides(curT, _N);
		};
		LS_STAG_KERNEL_INLINE void fluidTheta(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
                LSSTAGUniformMesh::fluidAlphaRange(curT, _N);
				curT->set(0.0, 1.0, 1.0);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidThetaF(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs, mN = _N +1; j < vf; j++)
			{
                LSSTAGUniformMesh::fluidAlphaRange(curT, mN);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidThetaPlaneF(coordinate*& curT, const int xN, const int _M, const int vf, const int zs = 0)
		{
			for (int kAdd = zs; kAdd < vf; kAdd++) LSSTAGUniformMesh::fluidThetaF(curT, xN, _M);
		};
		LS_STAG_KERNEL_INLINE void nordTheta(coordinate*& curT, const int _N)
		{
			for (int i = 0; i < _N; i++, curT++) curT->set(1.0, 0.0, 1.0);
			curT->set(0.0, 0.0, 1.0);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void nordThetaRange(coordinate*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(1.0, 0.0, 1.0);
		};
		LS_STAG_KERNEL_INLINE void fluidThetaPlanes(coordinate*& curT, const int _N, const int _M, 
			const int zf, const int zs = 0)
		{
			for (int k = zs; k < zf; k++)
			{
				LSSTAGUniformMesh::fluidTheta(curT, _N, _M);
				LSSTAGUniformMesh::nordTheta(curT, _N);
			}
		};
		LS_STAG_KERNEL_INLINE void solidThetaRange(coordinate*& curT, const int hs, const int hf)
		{
			for (int i = hs; i < hf; i++, curT++) curT->set(0.0);
		};
		LS_STAG_KERNEL_INLINE void topThetaRange(coordinate*& curT, const int _N)
		{
			for (int i = 0; i < _N; i++, curT++) curT->set(1.0, 1.0, 0.0);
		};
		LS_STAG_KERNEL_INLINE void topTheta(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
                LSSTAGUniformMesh::topThetaRange(curT, _N);
                curT->set(0.0, 1.0, 0.0);
			}
		};
		LS_STAG_KERNEL_INLINE void topThetaF(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs, mN = _N + 1; j < vf; j++)
			{
                LSSTAGUniformMesh::topThetaRange(curT, mN);
			}
		};
		LS_STAG_KERNEL_INLINE void topNordTheta(coordinate*& curT, const int _N)
		{
			for (int i = 0; i < _N; i++, curT++) curT->set(1.0, 0.0, 0.0);
			curT->set(0.0);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void topNordThetaRange(coordinate*& curT, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curT++) curT->set(1.0, 0.0, 0.0);
		};
		LS_STAG_KERNEL_INLINE void topThetaPlane(coordinate*& curT, const int _N, const int _M)
		{
				LSSTAGUniformMesh::topTheta(curT, _N, _M);
				LSSTAGUniformMesh::topNordTheta(curT, _N);
		};
		LS_STAG_KERNEL_INLINE void fluidAlphaRange(coordinate*& curA, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curA++) curA->set(1.0);
        }
		LS_STAG_KERNEL_INLINE void fluidAlpha(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
                LSSTAGUniformMesh::fluidAlphaRange(curT, _N);
				curT->set(1.0, 0.0, 0.0);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidAlphaF(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs, mN = _N + 1; j < vf; j++)
			{
                LSSTAGUniformMesh::fluidAlphaRange(curT, mN);
			}
		};
		LS_STAG_KERNEL_INLINE void fluidAlphaPlaneF(coordinate*& curT, const int xN, const int _M, const int zEnd, const int zs = 0)
		{
			for (int kAdd = zs; kAdd < zEnd; kAdd++) LSSTAGUniformMesh::fluidAlphaF(curT, xN, _M);
		}
		LS_STAG_KERNEL_INLINE void nordAlphaRange(coordinate*& curA, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curA++) curA->set(0.0, 1.0, 0.0);
		}
		LS_STAG_KERNEL_INLINE void nordAlpha(coordinate*& curT, const int _N)
		{
			LSSTAGUniformMesh::nordAlphaRange(curT, _N);
			curT->set(0.0);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void fluidAlphaPlanes(coordinate*& curT, const int _N, const int _M, 
			const int zf, const int zs = 0)
		{
			for (int k = zs; k < zf; k++)
			{
				LSSTAGUniformMesh::fluidAlpha(curT, _N, _M);
				LSSTAGUniformMesh::nordAlpha(curT, _N);
			}
		};
        LS_STAG_KERNEL_INLINE void topAlphaRange(coordinate*& curA, const int _N, const int hs = 0)
		{
			for (int i = hs; i < _N; i++, curA++) curA->set(0.0, 0.0, 1.0);
        }
		LS_STAG_KERNEL_INLINE void topAlpha(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs; j < vf; j++, curT++)
			{
                LSSTAGUniformMesh::topAlphaRange(curT, _N);
			curT->set(0.0);
			}
		};
		LS_STAG_KERNEL_INLINE void topAlphaF(coordinate*& curT, const int _N, const int vf, const int vs = 0)
		{
			for (int j = vs, mN = _N + 1; j < vf; j++)
			{
                LSSTAGUniformMesh::topAlphaRange(curT, mN);
			}
		};
		LS_STAG_KERNEL_INLINE void topNordAlpha(coordinate*& curT, const int _N)
		{
			for (int i = 0; i < _N; i++, curT++) curT->set(0.0);
			curT->set(0.0);
			curT++;
		};
		LS_STAG_KERNEL_INLINE void topAlphaPlane(coordinate*& curT, const int _N, const int _M)
		{
				LSSTAGUniformMesh::topAlpha(curT, _N, _M);
				LSSTAGUniformMesh::topNordAlpha(curT, _N);
		};
		LS_STAG_KERNEL_INLINE void RowTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const CellTypes type, const CellTypes type2,
			const CellTypes type3, const int m_vs, const int m_vf, const int _N)
		{
			for (int j = m_vs; j < m_vf; j++)
			{
				LSSTAGUniformMesh::DefineTypeOfCell(num, _baseType, type);
				for (int i = 1; i < _N; i++)
					LSSTAGUniformMesh::DefineTypeOfCell(num, _baseType, type2);
				LSSTAGUniformMesh::DefineTypeOfCell(num, _baseType, type3);
			}
		};
		LS_STAG_KERNEL_INLINE void PlaneTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const CellTypes type, 
			const CellTypes type2, const CellTypes type3, const CellTypes type4, const CellTypes type5,
			const CellTypes type6, const CellTypes type7, const CellTypes type8, const CellTypes type9, 
			const int m_zs, const int m_zf, const int _N, const int _M)
		{
			for (int k = m_zs; k < m_zf; k++)
			{
			LSSTAGUniformMesh::RowTypeOfCell(num, _baseType, type, type2, type3, 0, 1, _N);
			LSSTAGUniformMesh::RowTypeOfCell(num, _baseType, type4, type5, type6, 1, _M, _N);
			LSSTAGUniformMesh::RowTypeOfCell(num, _baseType, type7, type8, type9, 0, 1, _N);
			}
		};
		LS_STAG_KERNEL_INLINE void BottomPlaneTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const int _N, const int _M)
		{
			LSSTAGUniformMesh::PlaneTypeOfCell(num, _baseType, BSW_Corner, BS_Edge, BSE_Corner, BW_Edge, B_Border, BE_Edge, 
				BNW_Corner, BN_Edge, BNE_Corner, 0, 1, _N, _M);
		};
		LS_STAG_KERNEL_INLINE void FluidPlanesTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const int _N, const int _M, const int _K)
		{
			LSSTAGUniformMesh::PlaneTypeOfCell(num, _baseType, SW_Edge, S_Border, SE_Edge, W_Border, FluidCell, E_Border, 
				NW_Edge, N_Border, NE_Edge, 1, _K, _N, _M);
		};
		LS_STAG_KERNEL_INLINE void TopPlaneTypeOfCell(int& num, LS_STAG_R_PCT _baseType, const int _N, const int _M)
		{
			LSSTAGUniformMesh::PlaneTypeOfCell(num, _baseType, TSW_Corner, TS_Edge, TSE_Corner, TW_Edge, T_Border, TE_Edge, 
				TNW_Corner, TN_Edge, TNE_Corner, 0, 1, _N, _M);
		};
		LS_STAG_SINGLE_KERNEL(calcTypeOfSideBaseBox, LS_STAG_R_P(coord<SideTypes>) _SideType, const int xN, const int yM, const int zK)
		{
			coord<SideTypes>* curT = _SideType;
			LSSTAGUniformMesh::fluidSidesPlanes(curT, xN, yM, zK);
			LSSTAGUniformMesh::topSidesPlane(curT, xN, yM);
		}
		LS_STAG_SINGLE_KERNEL(calcTypeOfCellBaseBox, LS_STAG_R_PCT _baseType, const int xN, const int yM, const int zK)
		{
			int i = 0;
			LSSTAGUniformMesh::BottomPlaneTypeOfCell(i, _baseType, xN, yM);
			LSSTAGUniformMesh::FluidPlanesTypeOfCell(i, _baseType, xN, yM, zK);
			LSSTAGUniformMesh::TopPlaneTypeOfCell(i, _baseType, xN, yM);
		}
		LS_STAG_SINGLE_KERNEL(calcThetaBaseBox, LS_STAG_R_PC _theta, const int xN, const int yM, const int zK)
		{
			coordinate* curT = _theta;
			LSSTAGUniformMesh::fluidThetaPlanes(curT, xN, yM, zK);
			LSSTAGUniformMesh::topThetaPlane(curT, xN, yM);
		}
		LS_STAG_SINGLE_KERNEL(calcAlphaBaseBox, LS_STAG_R_PC _alpha, const int xN, const int yM, const int zK)
		{
			coordinate* curT = _alpha;
			LSSTAGUniformMesh::fluidAlphaPlanes(curT, xN, yM, zK);
			LSSTAGUniformMesh::topAlphaPlane(curT, xN, yM);
		}
		virtual void rebuildPrepare(){};//1.06.14
		/*!
		Cтроит сетку перед началом расчета
		*/
		virtual void init()
		{//в baseTaskDescription определили base.N,M, existXYMesh
			baseMesh.calcNodesCount();
			m_VOF = baseMesh.VOF;
			m_baseType = baseMesh.typeOfCell;
			NMK = baseMesh.NMK;
			m_NM2 = baseMesh.N + baseMesh.M - 1;// для ху-сетки 
			m_NM = baseMesh.N * baseMesh.M;
			xMesh.N = baseMesh.N - 1;
			xMesh.M = baseMesh.M;
			xMesh.K = baseMesh.K;
			xMesh.calcNodesCount();
			x_VOF = xMesh.VOF;
			x_Type = xMesh.typeOfCell;
			yMesh.N = baseMesh.N;
			yMesh.M = baseMesh.M - 1;
			yMesh.K = baseMesh.K;
			yMesh.calcNodesCount();
			y_VOF = yMesh.VOF;
			y_Type = yMesh.typeOfCell;
			zMesh.N = baseMesh.N;
			zMesh.M = baseMesh.M;
			zMesh.K = baseMesh.K - 1;
			zMesh.calcNodesCount();
			z_VOF = zMesh.VOF;
			z_Type = zMesh.typeOfCell;
			m_baseSideTypes.init(NMK, baseSideTypes);
			m_theta.init(NMK, theta);
			m_alpha.init(NMK, alpha);
			m_delta.init(NMK, delta);
			m_delta_ib_en.init(NMK, delta_ib_en);
			m_delta_ib_ws.init(NMK, delta_ib_ws);
			m_delta_ib_nb.init(NMK, delta_ib_nb);
			m_delta_ib_st.init(NMK, delta_ib_st);
			m_delta_ib_max.init(NMK, delta_ib_max);
			m_inv_delta_ib_ne.init(NMK, inv_delta_ib_ne); 
			m_inv_delta_ib_tn.init(NMK, inv_delta_ib_tn); 
			m_nodes.init(NMK, nodes);
			m_phi.init(NMK, phi);
			m_v_ib.init(NMK, v_ib);
			m_v_ib_yx.init(NMK, v_ib_yx);
			m_v_ib_zy.init(NMK, v_ib_zy);
			m_v_ib_ENT.init(NMK, v_ib_ENT);
			m_phi.set(0.0);
			m_v_ib.set(0.0);
			m_v_ib_yx.set(0.0);
			m_v_ib_zy.set(0.0);
			m_v_ib_ENT.set(0.0);
			m_U_g.init(NMK, U_g);
			m_U_g.set(0.0);
			if(existXYMesh)
			{// только для вязкоупругих жидкостей + для турбулентных для интерполяций косвенно они нужны
				xyMesh.N = xMesh.N;
				xyMesh.M = yMesh.M;
				xyMesh.K = baseMesh.K;
				xyMesh.calcNodesCount();
				xy_VOF = xyMesh.VOF;
				xy_Type = xyMesh.typeOfCell;
				m_weightXY.init(NMK, weightXY);
				m_quarterVxy.init(NMK, quarterVxy);
				xzMesh.N = xMesh.N;
				xzMesh.M = baseMesh.M;
				xzMesh.K = zMesh.K;
				xzMesh.calcNodesCount();
				xz_VOF = xzMesh.VOF;
				xz_Type = xzMesh.typeOfCell;
				m_weightXZ.init(NMK, weightXZ);
				m_quarterVxz.init(NMK, quarterVxz);
				yzMesh.N = baseMesh.N;
				yzMesh.M = yMesh.M;
				yzMesh.K = zMesh.K;
				yzMesh.calcNodesCount();
				yz_VOF = yzMesh.VOF;
				yz_Type = yzMesh.typeOfCell;
				m_weightYZ.init(NMK, weightYZ);
				m_quarterVyz.init(NMK, quarterVyz);
				m_inv_delta_ib.init(NMK, inv_delta_ib);
			}
		};
		virtual void buildMeshFunctions()
		{
			buildV();
			if(existXYMesh)
			{
				buildV_existXY();
				buildMxy();
				buildV_existXZ();
				buildMxz();
				buildV_existYZ();
				buildMyz();
			}
			buildMx();
			buildMy();
			buildMz();
		};
		void defineIbCount(int _ibCount)
		{
			ibCount = _ibCount;
			m_curV_C.resize(0);// 02.10.2026: init выделяет только при size==0 — сбрасываем, чтобы перевыделить под новый ibCount (ctor звал defineIbCount(1))
			m_cur_C.resize(0);
			m_curV_C.init(ibCount, curV_C);
			m_cur_C.init(ibCount, cur_C);
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, initIB, ibCount, curV_C, cur_C)
			m_curV_C.setCopy(pCurV_C);
			m_cur_C.setCopy(pCur_C);
		};
		void initSol(Solution* cur_sol, Solution* pred_sol, Solution* prev_sol)
		{
			cur_sol->getCopy(m_Ux, m_Uy, m_Uz, m_P);
			pred_sol->getCopy(m_predUx, m_predUy, m_predUz, m_predP);
			prev_sol->getCopy(m_prevUx, m_prevUy, m_prevUz, m_prevP);
		};
		/*!
		Проверяет, имеются ли в расчетной области подвижные погруженные границы
		\return {true в случае наличия движения границ в данный момент расчета, в противном случае false }
		*/
		bool movBody()
		{
			for(int i = 0; i < ibCount; i++)
				if((fabs(pCurV_C[i].x) > EQ_DOUBLE) || (fabs(pCurV_C[i].y) > EQ_DOUBLE) || 
				(fabs(pCurV_C[i].z) > EQ_DOUBLE)) return true;// 07.12.14
			return false;
		};
		/*!
		Перестраивает сетку (вызывается в случае перемещения подвижных границ
		\param deltaX,deltaY,deltaZ Шаги сетки по горизонтали и вертикали, z соответственно
		\param a,b,c Координаты нижнего левого угла расчетной области 
		*/
		void build(const double deltaX, const double deltaY, const double deltaZ, const double a, const double b, const double c)
		{
			minStep = (deltaX < deltaY) ? deltaX : deltaY;
			if (deltaZ < minStep) minStep = deltaZ;
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, buildDeltaNodes, baseMesh.M, baseMesh.K, nodes, delta, deltaX, deltaY, 
				deltaZ, a, b, c, baseMesh.N, m_NM)
		};
		void buildCV_U()//   определение типов КО   // //по x
		{
			LS_STAG_PROFILE("buildCV_U");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, cvu, baseMesh.M, baseMesh.K, m_baseType, x_Type, theta, baseMesh.N, 
				baseMesh.M, m_NM)
		};
		void buildCV_V()//   определение типов КО   // //по y
		{
			LS_STAG_PROFILE("buildCV_V");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, cvv, yMesh.M, baseMesh.K, m_baseType, y_Type, theta, baseMesh.N, m_NM)
		};
		void buildCV_W()//   определение типов КО   // //по z
		{
			LS_STAG_PROFILE("buildCV_W");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, cvw, zMesh.NMK, m_baseType, z_Type, theta, m_NM)
		};
		void buildCV_Txy()//   определение типов КО   // //по xy
		{
			LS_STAG_PROFILE("buildCV_Txy");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, cvtxy, yMesh.M, baseMesh.K, y_Type, xy_Type, baseMesh.N, m_NM2, m_NM)
		};
		void buildCV_Txz()//   определение типов КО   // //по xz
		{
			LS_STAG_PROFILE("buildCV_Txz");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, cvtxz, baseMesh.M, zMesh.K, z_Type, xz_Type, baseMesh.N, 
				baseMesh.M, m_NM)
		};
		void buildCV_Tyz()//   определение типов КО   // //по yz
		{
			LS_STAG_PROFILE("buildCV_Tyz");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, cvtyz, yMesh.M, zMesh.K, z_Type, yz_Type, baseMesh.N, m_NM)
		};
		void buildDelta_ib()
		{
			LS_STAG_PROFILE("buildDeltaIb");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, calcDeltaIB, NMK, m_baseType, theta, delta, delta_ib_ws, delta_ib_en, 
				delta_ib_nb, delta_ib_st, delta_ib_max, baseSideTypes, baseMesh.N, m_NM)
		};
		void buildInvDelta_ib()
		{
			LS_STAG_PROFILE("buildInvDeltaIb");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, buildInvDeltaIb, NMK, inv_delta_ib, m_baseType, baseSideTypes, alpha, 
				delta, baseMesh.N, m_NM)
		};
		void buildInvBaseDelta_ib()
		{
			LS_STAG_PROFILE("buildInvBaseDeltaIb");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, buildInvBaseDeltaIb, NMK, inv_delta_ib_ne, inv_delta_ib_tn,
				m_baseType, baseSideTypes, alpha, delta, baseMesh.N, m_NM)
		};
		void buildTypeOfCellBaseBox()
		{
			LS_STAG_SINGLE(LSSTAGUniformMesh, calcTypeOfCellBaseBox, m_baseType, xMesh.N, yMesh.M, zMesh.K)
			LS_STAG_SINGLE(LSSTAGUniformMesh, calcTypeOfSideBaseBox, baseSideTypes, xMesh.N, yMesh.M, zMesh.K)
		};
		void buildThetaBaseBox()
		{
			LS_STAG_SINGLE(LSSTAGUniformMesh, calcThetaBaseBox, theta, xMesh.N, yMesh.M, zMesh.K)
		};
		void buildAlphaBaseBox()
		{
			LS_STAG_SINGLE(LSSTAGUniformMesh, calcAlphaBaseBox, alpha, xMesh.N, yMesh.M, zMesh.K)
		};
	  private:
		void buildV()
		{
			LS_STAG_PROFILE("buildV");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, v, NMK, m_VOF, m_baseType, theta, alpha, delta, baseMesh.N, m_NM)
			baseMesh.inverseVOF();
		};
		void buildV_existXY()
		{
			LS_STAG_PROFILE("buildVxy");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, vXY, NMK, m_VOF, weightXY, quarterVxy, m_baseType, baseSideTypes)
		};
		void buildV_existXZ()
		{
			LS_STAG_PROFILE("buildVxz");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, vXZ, NMK, m_VOF, weightXZ, quarterVxz, m_baseType)
		};
		void buildV_existYZ()
		{
			LS_STAG_PROFILE("buildVyz");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, vYZ, NMK, m_VOF, weightYZ, quarterVyz, m_baseType)
		};
		void buildMx()
		{
			LS_STAG_PROFILE("buildMx");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, mx, xMesh.M, xMesh.K, x_VOF, m_VOF, baseMesh.N, baseMesh.M, m_NM)
			xMesh.inverseVOF();
		};
		void buildMy()
		{ // массы y-ых ячеек
			LS_STAG_PROFILE("buildMy");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, my, yMesh.M, yMesh.K, y_VOF, m_VOF, baseMesh.N, m_NM)
			yMesh.inverseVOF();
		};
		void buildMz()
		{ // массы z-ых ячеек
			LS_STAG_PROFILE("buildMz");
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, mz, zMesh.NMK, z_VOF, m_VOF, m_NM)
			zMesh.inverseVOF();
		};
		void buildMxy()
		{// массы xy-ячеек
			LS_STAG_PROFILE("buildMxy");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, mXY, xyMesh.M, xyMesh.K, xyMesh.VOF, quarterVxy, baseMesh.N, m_NM2, m_NM)
			xyMesh.inverseVOF();
		}
		void buildMxz()
		{// массы xz-ячеек
			LS_STAG_PROFILE("buildMxz");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, mXZ, xzMesh.M, xzMesh.K, xzMesh.VOF, quarterVxz, baseMesh.N, baseMesh.M, m_NM)
			xzMesh.inverseVOF();
		}
		void buildMyz()
		{// массы yz-ячеек
			LS_STAG_PROFILE("buildMyz");
			LS_STAG_PAR_FOR_2D(LSSTAGUniformMesh, mYZ, yzMesh.M, yzMesh.K, yzMesh.VOF, quarterVyz, baseMesh.N, m_NM)
			yzMesh.inverseVOF();
		}
	  public:
		LSSTAGBaseMeshParams baseMesh, xMesh, yMesh, zMesh, xyMesh, xzMesh, yzMesh;///< Параметры 7 разнесенных сеток
		CoordMeshFunction<double> m_v_ib;
		CoordMeshFunction<double> m_v_ib_yx;///< {u^{ib,E,y}, v^{ib,N,x}, w^{ib,T,x}} -- для расчета Quad^{ib} и вязких членов -- c интерполяцией
		CoordMeshFunction<double> m_v_ib_zy;///< {u^{ib,E,z}, v^{ib,N,z}, w^{ib,T,y}} -- для расчета Quad^{ib} и вязких членов -- c интерполяцией
		CoordMeshFunction<double> m_v_ib_ENT;///< {u^{ib,E}, v^{ib,N}, w^{ib,T}} -- для расчета конвективных членов -- всегда полусумма
		// скорее всего сюда надо еще 6 скоростей в узлах + скорости на гранях
		Array<coordinate> m_theta, m_alpha, m_delta, m_delta_ib_max, m_inv_delta_ib, m_inv_delta_ib_ne, m_inv_delta_ib_tn, m_delta_ib_en, m_delta_ib_ws, m_delta_ib_nb, m_delta_ib_st, m_nodes, m_cur_C, m_curV_C;
		Array<coord<SideTypes>> m_baseSideTypes; ///<  для основной сетки -- типы граней {e,n,t}
		extDouble m_phi, m_U_g;
		coordinate* cur_C;// текущие координаты центра тела
		coordinate* curV_C;// текущая скорость центра тела
		coordinate *pCur_C, *pCurV_C;// на CPU
		coord<SideTypes>* baseSideTypes;
		coordinate* theta;
		coordinate* alpha;
		coordinate* delta;
		coordinate* inv_delta_ib;///< {1/\Delta x^{ib}, 1/\Delta y^{ib}, 1/\Delta z^{ib}} --  для градиента турбулентных величин
		coordinate* inv_delta_ib_ne;///< {1/\Delta x^{ib,n}, 1/\Delta y^{ib,e}, 1/\Delta z^{ib,e}} --  для производных скорости
		coordinate* inv_delta_ib_tn;///< {1/\Delta x^{ib,t}, 1/\Delta y^{ib,t}, 1/\Delta z^{ib,n}} --  для производных скорости
		coordinate* delta_ib_en;///< {\Delta S^{ib,e}_{xy}, \Delta S^{ib,e}_{xz}, \Delta S^{ib,n}_{yz}}
		coordinate* delta_ib_ws;///< {\Delta S^{ib,w}_{xy}, \Delta S^{ib,w}_{xz}, \Delta S^{ib,s}_{yz}}
		coordinate* delta_ib_nb;///< {\Delta S^{ib,n}_{xy}, \Delta S^{ib,b}_{xz}, \Delta S^{ib,b}_{yz}}
		coordinate* delta_ib_st;///< {\Delta S^{ib,s}_{xy}, \Delta S^{ib,t}_{xz}, \Delta S^{ib,t}_{yz}}
		coordinate* delta_ib_max;///< max {delta_ib_en, delta_ib_ws, delta_ib_nb, delta_ib_st}
		coordinate* nodes;
		coordinate *v_ib;//основной шаблон
	    coordinate *v_ib_yx;//основной шаблон //25.09
	    coordinate *v_ib_zy;//основной шаблон //25.09
		coordinate *v_ib_ENT;//основной шаблон
		CellTypes *m_baseType;
		int ibCount;// 17/10/15
	  protected:
		double *m_VOF, *x_VOF, *y_VOF, *z_VOF, *xy_VOF, *xz_VOF, *yz_VOF;
		CellTypes *x_Type, *y_Type, *z_Type, *xy_Type, *xz_Type, *yz_Type;
		coordinate *curT, *curD;
		Array<double> m_weightXY, m_quarterVxy, m_weightXZ, m_quarterVxz, m_weightYZ, m_quarterVyz;
	  public:
		double *weightXY, *quarterVxy, *weightXZ, *quarterVxz, *weightYZ, *quarterVyz;// for xyMesh параметры baseMesh
		double* phi;// level set function
		double *m_Ux, *m_Uy, *m_Uz, *m_P;// указатели на поле скоростей и давление на текущем шаге
		double *m_predUx, *m_predUy, *m_predUz, *m_predP;///< указатели на прогноз поля скоростей и давления на текущем шаге
		double *m_prevUx, *m_prevUy, *m_prevUz, *m_prevP;//</ указатели на поле скоростей и давление на предыдущем шаге
		double* U_g;///< численный поток для подвижных границ // 10.01.2013
		double *U_ib;// указатель на A.S (Flux)
		double base_length;
		double time_coef;// 25.12.13 --- 1/t --- определяем в baseTaskDescription
		double minStep;
		MeshKinds id;
		int NMK, maxCountCutCells, size, existXYMesh;
	  protected:
		int i, j, num, numX, numN;
		int m_NM2, m_NM;
	};
}

#endif //_LS_STAG_LS_STAG_UNIFORM_MESH_H_
