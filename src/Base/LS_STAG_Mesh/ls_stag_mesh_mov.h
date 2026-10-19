/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ls_stag_mesh_mov.h
* \author    Puzikova Valeria (vvp) 
* \date      23.10.2013 20:12:54
* \brief     ls_stag_mesh_mov.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_MESH_MOV_H_
#define _LS_STAG_LS_STAG_MESH_MOV_H_

// ================================================================== INCLUDES
#include "ls_stag_mesh.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LSSTAGMeshMov

Данный файл содержит в себе реализацию класса LSSTAGMeshMov
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(cvuMov, const CellTypes*, _baseType, LS_STAG_R_PCT, _x_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshM, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(cvvMov, const CellTypes*, _baseType, LS_STAG_R_PCT, _y_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(cvw, const CellTypes*, _baseType, LS_STAG_R_PCT, _z_Type, LS_STAG_CR_PC, _theta,
		const int, _baseMeshN, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(cvtxyMov, const CellTypes*, _y_Type, LS_STAG_R_PCT, _xy_Type, const int, _baseMeshN,
		const int, NM2, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(cvtxzMov, const CellTypes*, _z_Type, LS_STAG_R_PCT, _xz_Type, const int, _baseMeshN,
		const int, _baseMeshM, const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(cvtyzMov, const CellTypes*, _z_Type, LS_STAG_R_PCT, _yz_Type, const int, _baseMeshN,
		const int, _baseMeshNM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(buildInvDeltaIb, LS_STAG_R_PC, invDib, LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, 
			LS_STAG_R_PC, _delta, const int, _N, const int, _NM, const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_R_KERNEL_2D(buildInvBaseDeltaIb, LS_STAG_R_PC, invDibNE, LS_STAG_R_PC, invDibTN,
		LS_STAG_CR_PCT, _CellType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_R_PC, _alpha, LS_STAG_R_PC, _delta, const int, _N, const int, _NM,
	    const int, _hsMov, const int, _hfMov);
	LS_STAG_DECLARATION_KERNEL(linVel, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, coordinate*, _v_ib, LS_STAG_CR_PC, _curV_C,
		coordinate*, _v_ib_yx, coordinate*, _v_ib_zy, coordinate*, _v_ib_ENT);
	LS_STAG_DECLARATION_KERNEL(angVel, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, coordinate*, _v_ib, LS_STAG_CR_PC, _nodes, LS_STAG_CR_PC, _delta, LS_STAG_CR_PC, _alpha,
		LS_STAG_CR_PC, P_C, coordinate*, _v_ib_yx, coordinate*, _v_ib_zy, coordinate*, _v_ib_ENT, const double, Omega, const int, _N);
/*!
\brief Класс для работы с LS-STAG-сеткой с подвижными погруженными границами

Данная структура используется для работы с неравномерной LS-STAG-сеткой в случае наличия движения погруженных границ
*/
	struct LSSTAGMeshMov : LSSTAGMesh
	{
		LSSTAGMeshMov():
		LSSTAGMesh(MESH_MOVING, sizeof(LSSTAGMeshMov))
		{
			;
		}
		LS_STAG_KERNEL_2D(cvuMov, const CellTypes* _baseType, LS_STAG_R_PCT _x_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshN, const int _baseMeshM, const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov, 
				numX = num - i - j * _baseMeshM; num < m_end; numX++, num++)// curTypeX();// 6.06.14
			{
				const CellTypes c = _baseType[num];
				const CellTypes d = _baseType[num + 1];
				_x_Type[numX] = ((_theta + num)->x < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
			}
		}
		LS_STAG_KERNEL_2D(cvvMov, const CellTypes* _baseType, LS_STAG_R_PCT _y_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshN, const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov,
				numY = num - j * _baseMeshN; num < m_end; numY++, num++)// curTypeY();
			{
				const CellTypes c = _baseType[num];
				const CellTypes d = _baseType[num + _baseMeshN];
				_y_Type[numY] = ((_theta + num)->y < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
			}
		}
		LS_STAG_KERNEL_2D(cvw, const CellTypes* _baseType, LS_STAG_R_PCT _z_Type, LS_STAG_CR_PC _theta,
			const int _baseMeshN, const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov; num < m_end; num++)// curTypeY();
			{
				const CellTypes c = _baseType[num];
				const CellTypes d = _baseType[num + _baseMeshNM];
				_z_Type[num] = ((_theta + num)->z < EQ_DOUBLE) ? SolidCell : ((d < c) ? d : c);
			}
		}
		LS_STAG_KERNEL_2D(cvtxyMov, const CellTypes* _y_Type, LS_STAG_R_PCT _xy_Type, const int _baseMeshN, const int NM2,
			const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov, 
				numXY = num - i - j * NM2, numY = num - j * _baseMeshN; num < m_end; numY++, numXY++, num++)// curTypeXY();
			{
				const CellTypes c = _y_Type[numY];
				const CellTypes d = _y_Type[numY + 1];
				_xy_Type[numXY] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL_2D(cvtxzMov, const CellTypes* _z_Type, LS_STAG_R_PCT _xz_Type, const int _baseMeshN,
			const int _baseMeshM, const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov, 
				numX = num - i - j * _baseMeshM; num < m_end; numX++, num++)// curTypeXY();
			{
				const CellTypes c = _z_Type[num];
				const CellTypes d = _z_Type[num + 1];
				_xz_Type[numX] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL_2D(cvtyzMov, const CellTypes* _z_Type, LS_STAG_R_PCT _yz_Type, const int _baseMeshN,
			const int _baseMeshNM, const int _hsMov, const int _hfMov)
		{
			for (int numS = i * _baseMeshN + j * _baseMeshNM, num = numS + _hsMov, m_end = numS + _hfMov, 
				numY = num - j * _baseMeshN; num < m_end; numY++, num++)// curTypeXY();
			{
				const CellTypes c = _z_Type[num];
				const CellTypes d = _z_Type[num + _baseMeshN];
				_yz_Type[numY] = (c < d) ? c : d;
			}
		}
		LS_STAG_KERNEL_2D(buildInvDeltaIb, LS_STAG_R_PC invDib, LS_STAG_CR_PCT _CellType, LS_STAG_CR_P(coord<SideTypes>) _SideType,
			LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM, const int _hsMov, const int _hfMov)
		{
			for(int numS = i * _N + j * _NM, num = numS + _hsMov, m_end = numS + _hfMov; num < m_end; num++)
			    LSSTAGUniformMesh::calcInvDeltaIb(num, invDib, _CellType[num], _SideType[num].z, _alpha, _delta, _N, _NM);
		}
		LS_STAG_KERNEL_2D(buildInvBaseDeltaIb, LS_STAG_R_PC invDibNE, LS_STAG_R_PC invDibTN, LS_STAG_CR_PCT _CellType,
			LS_STAG_CR_P(coord<SideTypes>) _SideType,
		    LS_STAG_R_PC _alpha, LS_STAG_R_PC _delta, const int _N, const int _NM, const int _hsMov, const int _hfMov)
		{
			for(int numS = i * _N + j * _NM, num = numS + _hsMov, m_end = numS + _hfMov; num < m_end; num++)
			    LSSTAGUniformMesh::calcInvBaseDeltaIb(num, invDibNE, invDibTN, _CellType[num], _SideType[num], _alpha, _delta, _N, _NM);
		}
		LS_STAG_KERNEL(linVel, LS_STAG_CR_PCT _baseType, LS_STAG_CR_P(coord<SideTypes>) _SideType, coordinate* _v_ib, LS_STAG_CR_PC _curV_C,
			coordinate* _v_ib_yx, coordinate* _v_ib_zy, coordinate* _v_ib_ENT)
		{
			const CellTypes curType = _baseType[i];
			bool isCutCell = false;
			if((curType > SolidCell) && (curType < FluidCell))
			isCutCell = true;
			else
			{
			switch(curType)
			{
				case B_Border:
				case T_Border:
				const SideTypes zSide = _SideType[i].z;
				if((zSide > SolidSide) && (zSide < FluidSide)) isCutCell = true;
			}
			}
			coordinate curVel;
			if (isCutCell)
			{
			curVel.set(_curV_C);
				_v_ib[i].set(&curVel);
			}
			else
			{
				_v_ib[i].set(0.0);
				curVel.set(0.0);
			}
			_v_ib_yx[i].set(&curVel);
			_v_ib_zy[i].set(&curVel);
			_v_ib_ENT[i].set(&curVel);
		}
		LS_STAG_KERNEL(angVel, LS_STAG_CR_PCT _baseType, LS_STAG_CR_P(coord<SideTypes>) _SideType, coordinate* _v_ib, LS_STAG_CR_PC _nodes, LS_STAG_CR_PC _delta, LS_STAG_CR_PC _alpha, LS_STAG_CR_PC P_C,
			coordinate* _v_ib_yx, coordinate* _v_ib_zy, coordinate* _v_ib_ENT, const double Omega, const int _N)
		{// вращение вокруг Oz
			coordinate curVel;	
			const CellTypes curType = _baseType[i];
			bool isCutCell = false;
			if((curType > SolidCell) && (curType < FluidCell))
			isCutCell = true;
			else
			{
			switch(curType)
			{
				case B_Border:
				case T_Border:
				const SideTypes zSide = _SideType[i].z;
				if((zSide > SolidSide) && (zSide < FluidSide)) isCutCell = true;
			}
			}
			if (isCutCell)
			{// cut-cells
				const coordinate* curT = _nodes + i;
				const coordinate* curD = _delta + i;
				const coordinate* curA = _alpha + i;
				const double cosb = curT->x + 0.5 * (curA->x + (curA - 1)->x - 1.0) * curD->x - P_C->x;
				const double sinb = curT->y + 0.5 * (curA->y + (curA - _N)->y - 1.0) * curD->y - P_C->y;
				const double fs = (cosb > EQ_DOUBLE) ? 1.0 : ((cosb < MINUS_EQ_DOUBLE) ? -1.0 : 0.0);
				const double fw = (sinb > EQ_DOUBLE) ? -1.0 : ((sinb < MINUS_EQ_DOUBLE) ? 1.0 : 0.0);
				curVel.set(fw * Omega * fabs(cosb), fs * Omega * fabs(sinb), 0.0);
				_v_ib[i].set(&curVel);
			}
			else 
			{
				curVel.set(0.0);
				_v_ib[i].set(0.0);
			}
			_v_ib_yx[i].set(&curVel);
			_v_ib_zy[i].set(&curVel);
			_v_ib_ENT[i].set(&curVel);
		}
		void init()
		{
			LSSTAGUniformMesh::init();
			m_theta_N.init(NMK, theta_N);
			m_V_N.init(NMK, V_N);
		};
		void rebuildPrepare()
		{
			m_theta_N.copy(theta);
			m_V_N.copy(m_VOF);
			m_v_ib.set(0.0);
			m_v_ib_yx.set(0.0);
			m_v_ib_zy.set(0.0);
			m_v_ib_ENT.set(0.0);
		};
		void buildMeshFunctions()
		{
			LS_STAG_PROFILE("rebuildMeshFunc");
			LSSTAGUniformMesh::buildMeshFunctions();
			m_U_g.copy(m_VOF);
			m_U_g -= V_N;
			m_U_g *= time_coef;
		};
		void rebuildCV_U()//   определение типов КО   // //по x
		{
			LS_STAG_PROFILE("rebuildCV_U");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvuMov, vsMov, vfMov, zsMov, zfMov, m_baseType, x_Type, theta, 
				baseMesh.N, baseMesh.M, m_NM, hsMov, hfMov)
		};
		void rebuildCV_V()//   определение типов КО   // //по y
		{
			LS_STAG_PROFILE("rebuildCV_V");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvvMov, vsMov, vfMov, zsMov, zfMov, m_baseType, y_Type, theta, 
				baseMesh.N, m_NM, hsMov, hfMov)
		};
		void rebuildCV_W()//   определение типов КО   // //по z
		{
			LS_STAG_PROFILE("rebuildCV_W");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvw, vsMov, vfMov, zsMov, zfMov - 1, m_baseType, z_Type, theta, 
				baseMesh.N, m_NM, hsMov, hfMov)
		};
		void rebuildCV_Txy()//   определение типов КО   // //по xy
		{
			LS_STAG_PROFILE("rebuildCV_Txy");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvtxyMov, vsMov, vfMov, zsMov, zfMov, y_Type, xy_Type, baseMesh.N, 
				m_NM2, m_NM, hsMov, hfMov)
		};
		void rebuildCV_Txz()//   определение типов КО   // //по xz
		{
			LS_STAG_PROFILE("rebuildCV_Txz");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvtxzMov, vsMov, vfMov, zsMov, zfMov - 1, z_Type, xz_Type, baseMesh.N, 
				baseMesh.M, m_NM, hsMov, hfMov)
		};
		void rebuildCV_Tyz()//   определение типов КО   // //по yz
		{
			LS_STAG_PROFILE("rebuildCV_Tyz");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, cvtyzMov, vsMov, vfMov, zsMov, zfMov - 1, z_Type, yz_Type, baseMesh.N, 
				m_NM, hsMov, hfMov)
		};
		void rebuildInvDelta_ib()
		{
			LS_STAG_PROFILE("rebuildInvDeltaIb");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, buildInvDeltaIb, vsMov, vfMov, zsMov, zfMov, inv_delta_ib, 
				m_baseType, baseSideTypes, alpha, delta, baseMesh.N, m_NM, hsMov, hfMov)
		};
		void rebuildInvBaseDelta_ib()
		{
			LS_STAG_PROFILE("rebuildInvBaseDeltaIb");
			LS_STAG_PAR_FOR_2D_RANGE(LSSTAGMeshMov, buildInvBaseDeltaIb, vsMov, vfMov, zsMov, zfMov, inv_delta_ib_ne, 
				inv_delta_ib_tn, m_baseType, baseSideTypes, alpha, delta, baseMesh.N, m_NM, hsMov, hfMov)
		};
		void rebuildLinearVelocity()
		{
			LS_STAG_PROFILE("rebuildLinVel");
			LS_STAG_PAR_FOR(LSSTAGMeshMov, linVel, NMK, m_baseType, baseSideTypes, v_ib, curV_C, v_ib_yx, v_ib_zy, v_ib_ENT)
		};
		void rebuildAngleVelocity(const double Omega, const coordinate &P_C)
		{
			LS_STAG_PROFILE("rebuildAngVel");
			LS_STAG_PAR_FOR(LSSTAGMeshMov, angVel, NMK, m_baseType, baseSideTypes, v_ib, nodes, delta, alpha, &P_C, v_ib_yx, v_ib_zy, v_ib_ENT, Omega, baseMesh.N)
		};
	  private:
		CoordMeshFunction<double> m_theta_N;
		extDouble m_V_N;
	  public:
		coordinate* theta_N;// theta с прошлого шага
		double* V_N;// baseMesh.VOF с предыдущего шага
		int hsMov, hfMov, vsMov, vfMov, zsMov, zfMov;
	};
}

#endif //_LS_STAG_LS_STAG_MESH_MOV_H_
