/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 15:07:38
* \brief     flow_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "../../base_task_description.h"
#include "../../../MeshFunctions/coord_mesh_function.h"
#include "../../../LS_STAG_Mesh/ls_stag_mesh.h"
#include "../../../LS_STAG_Mesh/ls_stag_mesh_mov.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowTaskDescription

Данный файл содержит в себе реализацию класса FlowTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcAlphaMovBlock, LS_STAG_R_PC, _alpha, LS_STAG_R_PD, _phi, const int, _hsMov,
		const int, _hfMov, const int, _N, const int, _NM, const int, _zK);
	LS_STAG_DECLARATION_R_KERNEL_2D(calcThetaMovBlock, LS_STAG_R_PC, _theta,  LS_STAG_R_PC, _alpha, 
		LS_STAG_CR_P(coord<SideTypes>), _baseType, const int, _hsMov, const int, _hfMov, const int, _N, const int, _NM, const int, _zK);
	LS_STAG_DECLARATION_R_KERNEL_2D(rebuildMovBlock, const int, hsMov, const int, hfMov, LS_STAG_R_PCT, _baseType, 
		LS_STAG_R_P(coord<SideTypes>), _baseSideType, LS_STAG_R_PD, _phi, LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_R_PI, _cutcellCounts, 
		LS_STAG_R_P(tripletInt*), _airfoilNums, const int, _N, const int, _M, const int, _NM, const int, _zK);
	LS_STAG_DECLARATION_S_KERNEL(calcAlpha, LS_STAG_R_PC, _alpha, const int, hs_loc, const int, hf_loc,
		const int, vs_loc, const int, vf_loc, LS_STAG_R_PD, _phi, const int, xN, const int, yM, const int, zK,
		 const int, N, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(calcTheta, LS_STAG_R_PC, _theta,  LS_STAG_R_PC, _alpha, 
		LS_STAG_CR_P(coord<SideTypes>), _baseType, const int, hs_loc, const int, hf_loc, const int, vs_loc, 
		const int, vf_loc, const int, xN, const int, yM, const int, zK, const int, N, const int, NM);
	LS_STAG_DECLARATION_S_KERNEL(calcTypeOfCells, LS_STAG_R_PD, _phi, LS_STAG_R_PCT, _baseType,
		LS_STAG_R_P(coord<SideTypes>), _baseSideType, const int, hs_loc, const int, hf_loc, const int, vs_loc, const int, vf_loc, const int, xN, const int, yM, const int, zK,
		const int, N, const int, M, const int, NM, LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_R_PI, _cutcellCounts, LS_STAG_R_P(tripletInt*), _airfoilNums);
/*!
\brief Класс, реализующий работу с задачей об обтекании профилей

Данный класс реализует работу с описанием решаемой задачи об обтекании профилей
*/
	struct FlowTaskDescription : BaseTaskDescription
	{
		static const int DEFAULT_NODES_COUNT_IN_UNIFORM_BLOCK = 16;
		/// Тип профиля
		enum AirfoilKind
		{
		   AIRFOIL_UNKNOWN = -1, ///< Тип неизвестен
		   AIRFOIL_SQUARE,       ///< Квадратный профиль
		   AIRFOIL_ELLIPTIC,     ///< Эллиптический профиль
		   AIRFOIL_CIRCLE,       ///< Круговой просиль
		   AIRFOIL_CIRCLES,      ///< Система круговых профилей
		   AIRFOIL_POINTS,       ///< Профиль, заданный поточечно
		   AIRFOIL_ROTOR,        ///< Ротор
		   AIRFOIL_CUBE,         ///< Куб (наследник квадратного цилиндра)
		   AIRFOIL_FINITE_CYLINDER, ///< Цилиндр конечной длины (плоские торцы, наклон в плоскости x-z)
		   AIRFOIL_SPHERE        ///< Сфера (обтекание сферы, D = 2 rad = 1)
		};
	  public:
		explicit FlowTaskDescription(AirfoilKind _id_airfoil = AIRFOIL_UNKNOWN, MotionKinds _id_motion = MOTION_UNKNOWN, 
									LSSTAGUniformMesh* _mesh = (LSSTAGUniformMesh*)0, int _size = sizeof(FlowTaskDescription)):
		  BaseTaskDescription(TASK_FLOW, _id_motion, _mesh, _size),
		  loc_mesh((LSSTAGMesh*)_mesh),
		  id_airfoil(_id_airfoil),
		  zone_y_3(2.125),// --- Cheny & Botella 
		  ANGLE(0.0),//default
		  aa(1.0),// всегда
		  zoneZ(20),//(40),
		  nn(DEFAULT_NODES_COUNT_IN_UNIFORM_BLOCK)//default
		  {
			U_T = 0.0;// x-я скорость набегающего потока на T стенке 
			U_B = 0.0;//  x-я скорость набегающего потока на B стенке
			U_N = V_BASE;// x-я скорость набегающего потока на N стенке 
			U_S = V_BASE;//  x-я скорость набегающего потока на S стенке
			  U_E = V_BASE;//  x-я скорость набегающего потока на E стенке
			  U_W = V_BASE;//  x-я скорость набегающего потока на W стенке
			  V_T = 0.0;// y-я скорость набегающего потока на T стенке 
			  V_B = 0.0;//  y-я скорость набегающего потока на B стенке
			  V_N = 0.0;// y-я скорость набегающего потока на N стенке 
			  V_S = 0.0;//  y-я скорость набегающего потока на S стенке
			  V_E = 0.0;//  y-я скорость набегающего потока на E стенке
			  V_W = 0.0;//  y-я скорость набегающего потока на W стенке
			  W_T = 0.0;// z-я скорость набегающего потока на T стенке 
			  W_B = 0.0;//  z-я скорость набегающего потока на B стенке
			  W_N = 0.0;// z-я скорость набегающего потока на N стенке 
			  W_S = 0.0;//  z-я скорость набегающего потока на S стенке
			  W_E = 0.0;//  z-я скорость набегающего потока на E стенке
			  W_W = 0.0;//  z-я скорость набегающего потока на W стенке
			  P_C.x = 8.0;// default
			  P_C.y = 12.0;// default (при изменении , возможно, нужно изменять блоки сетки )
			  P_C.z = 12.0;// default (при изменении , возможно, нужно изменять блоки сетки )// для цилиндров не сипользуется
			  adaptiveT = 0;// default --- пока схема не отлажена
		  }
		LS_STAG_SINGLE_KERNEL(calcAlpha, LS_STAG_R_PC _alpha, const int hs_loc, const int hf_loc, const int vs_loc,
			const int vf_loc, LS_STAG_R_PD _phi, const int xN, const int yM, const int zK, const int N, const int NM)
		{// _theta --  это alpha
			coordinate* curT = _alpha;
			int rightN = xN - hf_loc;
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
				LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::bottomBuildAlpha(curT, _alpha, _phi, N, NM);
				LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, yM + 1, vf_loc);
			for(int k = 1; k < zK; k++)
			{
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
			    LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::BuildAlpha(curT, _alpha, _phi, N, NM);//окрестность профиля
		     	LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, yM + 1, vf_loc);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
				LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::topBuildAlpha(curT, _alpha, _phi, N, NM);
				LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, yM + 1, vf_loc);
		}
		LS_STAG_SINGLE_KERNEL(calcTheta, LS_STAG_R_PC _theta, LS_STAG_R_PC _alpha, LS_STAG_CR_P(coord<SideTypes>) _baseType, 
		    const int hs_loc, const int hf_loc, const int vs_loc, const int vf_loc, const int xN, 
			const int yM, const int zK, const int N, const int NM)
		{// _theta --  это alpha
			coordinate* curT = _theta;
			int rightN = xN - hf_loc;
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
				LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::tbBuildTheta(curT, _theta, _alpha, _baseType, N);
				LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, yM + 1, vf_loc);
			for(int k = 1; k < zK; k++)
			{
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
                LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::BuildTheta(curT, _theta, _alpha, _baseType, N, NM);//окрестность профиля
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, yM + 1, vf_loc);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs_loc);
			for (int j = vs_loc; j < vf_loc; j++)
			{
				LSSTAGUniformMesh::fluidAlphaRange(curT, hs_loc);
				for (int i = hs_loc; i < hf_loc; i++)
					LSSTAGUniformMesh::tbBuildTheta(curT, _theta, _alpha, _baseType, N);
				LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, yM + 1, vf_loc);
		}
		LS_STAG_SINGLE_KERNEL(calcTypeOfCells, LS_STAG_R_PD _phi, LS_STAG_R_PCT _baseType, LS_STAG_R_P(coord<SideTypes>) _SideType, 
		    const int hs_loc, const int hf_loc, const int vs_loc, const int vf_loc, const int xN, const int yM, const int zK,
			const int N, const int M, const int NM, LS_STAG_CR_PI _cellOnAirfoil,
			LS_STAG_R_PI _cutcellCounts, LS_STAG_R_P(tripletInt*) _airfoilNums)
		{
			int curNum = 0, nordY = yM - vf_loc, rightX = xN - hf_loc;
			coord<SideTypes>* curT = _SideType;
			LSSTAGUniformMesh::fluidSidesF(curT, xN, vs_loc);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BSW_Corner, BS_Edge, BSE_Corner, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, 1, vs_loc, xN);// жидкая область
			for (int jAdd = vs_loc; jAdd < vf_loc; jAdd++)
			{
				LSSTAGUniformMesh::fluidSideRange(curT, hs_loc);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, B_Border, 0, 1, hs_loc - 1);// жидкая область
				for (int i = hs_loc; i < hf_loc; i++, curNum++)
					LSSTAGUniformMesh::bottomTypeOfCell(curNum, _phi, _baseType, _SideType, N, NM);
				curT = _SideType + curNum;
				LSSTAGUniformMesh::fluidSidesF(curT, rightX, 1);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, B_Border, B_Border, BE_Edge, 0, 1, rightX);// жидкая область
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, vf_loc, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BNW_Corner, BN_Edge, BNE_Corner, 0, 1, xN);// N граница РО
			for (int kAdd = 1; kAdd < zK; kAdd++)
			{
				LSSTAGUniformMesh::fluidSidesF(curT, xN, vs_loc);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, SW_Edge, S_Border, SE_Edge, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, 1, vs_loc, xN);// жидкая область
			for (int jAdd = vs_loc; jAdd < vf_loc; jAdd++)
			{
				LSSTAGUniformMesh::fluidSideRange(curT, hs_loc);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
				for (int i = 1; i < hs_loc; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				for (int i = hs_loc; i < hf_loc; i++, curNum++)//окрестность профиля
				{
					if (LSSTAGUniformMesh::TypeOfCell(curNum, _phi, _baseType, _SideType, N, NM))//усечённая ячейка
					{
						const int airNum = LS_STAG_VIRTUAL_INLINE_CALL(CLASS_NAME_AIRFOIL_NUM, getAirfoilNum, curNum, _cellOnAirfoil);
						BaseTaskDescription::addCutCell(airNum, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, N, M);
					}
				}
				curT = _SideType + curNum;
				LSSTAGUniformMesh::fluidSidesF(curT, rightX, 1);
				for (int i = hf_loc; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, vf_loc, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, NW_Edge, N_Border, NE_Edge, 0, 1, xN);// N граница РО			}
		}
			LSSTAGUniformMesh::fluidSidesF(curT, xN, vs_loc);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TSW_Corner, TS_Edge, TSE_Corner, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 1, vs_loc, xN);// жидкая область
			for (int jAdd = vs_loc; jAdd < vf_loc; jAdd++)
			{
				LSSTAGUniformMesh::fluidSideRange(curT, hs_loc);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, T_Border, 0, 1, hs_loc - 1);// жидкая область
				for (int i = hs_loc; i < hf_loc; i++, curNum++)
					LSSTAGUniformMesh::topTypeOfCell(curNum, _phi, _baseType, _SideType, N, NM);
				curT = _SideType + curNum;
				LSSTAGUniformMesh::fluidSidesF(curT, rightX, 1);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, T_Border, T_Border, TE_Edge, 0, 1, rightX);// жидкая область
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, vf_loc, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TNW_Corner, TN_Edge, TNE_Corner, 0, 1, xN);// N граница РО
	}
		LS_STAG_KERNEL_2D(calcAlphaMovBlock, LS_STAG_R_PC _alpha, LS_STAG_R_PD _phi, const int _hsMov, const int _hfMov, const int _N, const int _NM, const int _zK)
		{//_theta --  alpha
			if(j)
			{
			coordinate* curT = _alpha + (i * _N + j * _NM + _hsMov);
			if (j < _zK)
			{
				for (int k = _hsMov; k < _hfMov; k++)//окрестность профиля
				LSSTAGUniformMesh::BuildAlpha(curT, _alpha, _phi, _N, _NM);
			}
			else
			{
				for (int k = _hsMov; k < _hfMov; k++)
				    LSSTAGUniformMesh::topBuildAlpha(curT, _alpha, _phi, _N, _NM);
			}
			}
			else 
			{
			coordinate* curT = _alpha + (i * _N + _hsMov);
				for (int k = _hsMov; k < _hfMov; k++)
					LSSTAGUniformMesh::bottomBuildAlpha(curT, _alpha, _phi, _N, _NM);
			}
		}
		LS_STAG_KERNEL_2D(calcThetaMovBlock, LS_STAG_R_PC _theta, LS_STAG_R_PC _alpha, LS_STAG_CR_P(coord<SideTypes>) _baseType,
		    const int _hsMov, const int _hfMov, const int _N, const int _NM, const int _zK)
		{//_theta --  alpha
			coordinate* curT = _theta + (i * _N + j * _NM + _hsMov);
			if(j && (j < _zK))
			{
				for (int k = _hsMov; k < _hfMov; k++)//окрестность профиля
				LSSTAGUniformMesh::BuildTheta(curT, _theta, _alpha, _baseType, _N, _NM);
			}
			else
			{
				for (int k = _hsMov; k < _hfMov; k++)
				    LSSTAGUniformMesh::tbBuildTheta(curT, _theta, _alpha, _baseType, _N);
			}
		}
		LS_STAG_KERNEL_2D(rebuildMovBlock, const int hsMov, const int hfMov, LS_STAG_R_PCT _baseType,
			LS_STAG_R_P(coord<SideTypes>) _baseSideType, LS_STAG_R_PD _phi, LS_STAG_CR_PI _cellOnAirfoil, LS_STAG_R_PI _cutcellCounts, 
			LS_STAG_R_P(tripletInt*) _airfoilNums, const int _N, const int _M, const int _NM, const int _zK)
		{
			if(j)
			{
			int curNum = i * _N + j * _NM + hsMov;
			if (j < _zK)
			{
				for (int k = hsMov; k < hfMov; k++, curNum++)//окрестность профиля
				{
					if (LSSTAGUniformMesh::TypeOfCell(curNum, _phi, _baseType, _baseSideType, _N, _NM))//усечённая ячейка
					{
						const int airNum = LS_STAG_VIRTUAL_INLINE_CALL(CLASS_NAME_AIRFOIL_NUM, getAirfoilNum, curNum, _cellOnAirfoil);
						BaseTaskDescription::addCutCell(airNum, _cutcellCounts, _airfoilNums, curNum, i, j, _N, _M);
					}
				}
			}
			else
			{
				for (int k = hsMov; k < hfMov; k++, curNum++)
					LSSTAGUniformMesh::topTypeOfCell(curNum, _phi, _baseType, _baseSideType, _N, _NM);
			}
			}
			else 
			{
				for (int curNum = i * _N + hsMov, k = hsMov; k < hfMov; k++, curNum++)
					LSSTAGUniformMesh::bottomTypeOfCell(curNum, _phi, _baseType, _baseSideType, _N, _NM);
			}
		}
		virtual void rebuildMesh()
		{
			LS_STAG_PROFILE("rebuildMesh");
			calcPosition();//5.06.14
			checkExistMov();
			if(existMov)
			{	
				loc_mesh2->rebuildPrepare();
				buildLS(loc_mesh2->zsMov, loc_mesh2->zfMov, loc_mesh2->vsMov, loc_mesh2->vfMov, loc_mesh2->hsMov, loc_mesh2->hfMov);
				{
					LS_STAG_PROFILE("buildAlpha");
					LS_STAG_PAR_FOR_2D_RANGE(FlowTaskDescription, calcAlphaMovBlock, loc_mesh2->vsMov, loc_mesh2->vfMov, 
						loc_mesh2->zsMov, loc_mesh2->zfMov, m_mesh->alpha, m_mesh->phi, loc_mesh2->hsMov, loc_mesh2->hfMov, 
						m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, m_mesh->zMesh.K)// коэффициенты усечения ячеек
				}
				loadsInit();
				{
					LS_STAG_PROFILE("addCutCell");
					LS_STAG_PAR_FOR_2D_RANGE(FlowTaskDescription, rebuildMovBlock, loc_mesh2->vsMov, loc_mesh2->vfMov, loc_mesh2->zsMov, 
						loc_mesh2->zfMov, loc_mesh2->hsMov, loc_mesh2->hfMov, m_mesh->m_baseType, m_mesh->baseSideTypes, m_mesh->phi, 
						getPCellOnAirfoill(), m_cutcellCounts, m_airfoilNums, m_mesh->baseMesh.N, m_mesh->baseMesh.M,
						 m_mesh->baseMesh.N * m_mesh->baseMesh.M, m_mesh->zMesh.K)
				}
				{
					LS_STAG_PROFILE("buildTheta");
					LS_STAG_PAR_FOR_2D_RANGE(FlowTaskDescription, calcThetaMovBlock, loc_mesh2->vsMov, loc_mesh2->vfMov, 
						loc_mesh2->zsMov, loc_mesh2->zfMov, m_mesh->theta, m_mesh->alpha, m_mesh->baseSideTypes, loc_mesh2->hsMov, 
						loc_mesh2->hfMov, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, m_mesh->zMesh.K)//коэффициенты усечения ячеек
				}
				loadsResize();
				loc_mesh2->rebuildCV_U();
				loc_mesh2->rebuildCV_V();
				loc_mesh2->rebuildCV_W();
				if(needXYMesh) 
				{
				loc_mesh2->rebuildCV_Txy();//7.02.14
				loc_mesh2->rebuildCV_Txz();
				loc_mesh2->rebuildCV_Tyz();
				if(useTM) loc_mesh2->rebuildInvDelta_ib();
				}
				loc_mesh2->buildDelta_ib();
				loc_mesh2->rebuildInvBaseDelta_ib();
				buildVib();// пересчитываем скорости на погруженных границах
				loc_mesh2->buildMeshFunctions();
			}
		};
		virtual void load(std::ifstream &in)
		{
			ibCount = 1;
			BaseTaskDescription::load(in);
			loadDiscontin(in);
			in >> tScaleP;
			in >> nn;
			in >> needXYMesh;
			in >> useTM;
			if(useTM || useNNM) needXYMesh = 1; 
		};
	  protected:
		virtual void calcPosition(){};
		virtual void printExtParams(std::ofstream &out)
		{
			 printDiscontin(out);
			 out << tScaleP << "\n";
			 out << nn << "\n";
			 out << needXYMesh << "\n";
			 out << useTM << "\n";
		};
		int getAirfoilKind(){ return (int)id_airfoil; };
		virtual void loadDiscontin(std::ifstream &in){ in >> discontin; };
		virtual void printDiscontin(std::ofstream &out){ out << discontin << "\n"; };
		virtual void calcNodesCount()
		{
			m_mesh->baseMesh.N = (int)floor(7.5 * nn);//номер последнего узла по горизонтали
			m_mesh->baseMesh.M = (int)floor((2.0 * zone_y_3 + 5.0) * nn); //номер последнего узла по вертикали 
			m_mesh->baseMesh.K = zoneZ; //номер последнего узла по вертикали 
			P_C.y = 4.0 * zone_y_3 * aa + 3.5;//Grid.yc = 12.0;
			P_C.z = 4.0* zone_y_3* aa + 3.5;//Grid.yc = 12.0; // для цилиндров не используется
		};
		virtual void buildBaseMesh()
		{
			hs_loc = loc_mesh->hs - 3;
			hf_loc = loc_mesh->hf + 3;
			vs_loc = loc_mesh->vs - 3;
			vf_loc = loc_mesh->vf + 3;
			zs_loc = 0;// loc_mesh->zs - 2;
			zf_loc = loc_mesh->baseMesh.K - zs_loc;// todo loc_mesh->zf + 3;
			loc_mesh->buildNode();
			buildLS(zs_loc, zf_loc, vs_loc, vf_loc, hs_loc, hf_loc);
			BaseTaskDescription::buildBaseMesh();
		};
		virtual void initCenter()
		{
			m_mesh->pCur_C->set(&P_C);
			m_mesh->m_cur_C.sync(m_mesh->pCur_C);
		}; 
		void buildAlpha()
		{
			LS_STAG_SINGLE(FlowTaskDescription, calcAlpha, m_mesh->alpha, hs_loc, hf_loc, vs_loc, vf_loc, m_mesh->phi, 
				m_mesh->xMesh.N, m_mesh->yMesh.M, m_mesh->zMesh.K, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M)
		};
		void buildTheta()
		{
			LS_STAG_SINGLE(FlowTaskDescription, calcTheta, m_mesh->theta, m_mesh->alpha, m_mesh->baseSideTypes, hs_loc, 
				hf_loc, vs_loc, vf_loc, m_mesh->xMesh.N, m_mesh->yMesh.M, m_mesh->zMesh.K, m_mesh->baseMesh.N, 
				m_mesh->baseMesh.N * m_mesh->baseMesh.M)
		};
		void buildTypeOfCells()
		{
			LS_STAG_SINGLE(FlowTaskDescription, calcTypeOfCells, m_mesh->phi, m_mesh->m_baseType, m_mesh->baseSideTypes, hs_loc, 
				hf_loc, vs_loc, vf_loc, m_mesh->xMesh.N, m_mesh->yMesh.M, m_mesh->zMesh.K, m_mesh->baseMesh.N, m_mesh->baseMesh.M,
				m_mesh->baseMesh.N * m_mesh->baseMesh.M, getPCellOnAirfoill(), m_cutcellCounts, m_airfoilNums)
		};
		void initMovMesh()
		{
			loc_mesh->hs = (int)floor(1.125 * nn);
			loc_mesh->hf = m_mesh->baseMesh.N - 3 * loc_mesh->hs; // 25.01.13
			loc_mesh->vs = (int)floor((zone_y_3 + 1) * nn);
			loc_mesh->vf = m_mesh->baseMesh.M - loc_mesh->vs; // 25.01.13
			loc_mesh->zs = 0;// если это тело вытянутое по z (int)floor((zone_y_3 + 1) * nn);
			loc_mesh->zf = m_mesh->baseMesh.K - loc_mesh->zs; // 25.01.13
			loc_mesh2 = (LSSTAGMeshMov*)loc_mesh;
			loc_mesh2->hsMov = loc_mesh->hs - 2;
			loc_mesh2->hfMov = loc_mesh->hf + 3;
			loc_mesh2->vsMov = loc_mesh->vs - 2;
			loc_mesh2->vfMov = loc_mesh->vf + 3;
			loc_mesh2->zsMov = loc_mesh->zs;// - 2;
			loc_mesh2->zfMov = loc_mesh->zf;//  + 3;
		};
		void initFinestMesh()
		{
			m_mesh->maxCountCutCells *= 16; 
			loc_mesh->zone[2].x = (int)floor(0.5 * nn); //число узлов по х , зона2
			loc_mesh->zone[2].y = 12 * nn; //число узлов по у, зона1
			loc_mesh->zone[3].x = 12 * nn; //число узлов по х , зона1
			loc_mesh->zone[4].x = (int)floor(0.5 * nn); //число узлов по х , зона2
			loc_mesh->deltaCurZone[0].x = 0.25 * aa / nn; //шаг сетки по х , зона1
			loc_mesh->deltaCurZone[0].y = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона1	
		};
		virtual void initMeshParams()
		{// это для обтекания цилиндров, для куба/сферы сетка по z д.б. как по у 
			loc_mesh = (LSSTAGMesh*)m_mesh;
			loc_mesh->zone[0].x = (int)floor(0.5 * nn); //число узлов по х, зона4
			loc_mesh->zone[0].y = (int)floor(zone_y_3 * nn); //число узлов по у, зона3
			loc_mesh->zone[0].z = (int)(0.5 * zoneZ); //(int)floor(zone_y_3 * nn); //число узлов по у, зона3
			loc_mesh->zone[1].x = (int)floor(0.375 * nn); //число узлов по х , зона3
			loc_mesh->zone[1].y = nn; //число узлов по  у, зона2
			loc_mesh->zone[1].z = 0;// nn; //число узлов по  у, зона2
			loc_mesh->zone[2].x = (int)floor(0.75 * nn); //число узлов по х , зона2
			loc_mesh->zone[2].y = 3 * nn; //число узлов по у, зона1
			loc_mesh->zone[2].z = 0;// 3 * nn; //число узлов по у, зона1
			loc_mesh->zone[3].x = 2 * nn; //число узлов по х , зона1
			loc_mesh->zone[3].y = nn; //число узлов по  у, зона21
			loc_mesh->zone[3].z = 0;// nn; //число узлов по  у, зона2
			loc_mesh->zone[4].x = (int)floor(0.75 * nn); //число узлов по х , зона2
			loc_mesh->zone[4].y = loc_mesh->zone[0].y; //число узлов по  у, зона3
			loc_mesh->zone[4].z = loc_mesh->zone[0].z; //число узлов по  у, зона3
			loc_mesh->zone[5].x = (int)floor(3.125 * nn); //число узлов по х , зона3
			loc_mesh->zone[5].y = 0; //число узлов по  у, зона4
			loc_mesh->zone[5].z = 0; //число узлов по  у, зона4
			loc_mesh->zone[6].x = 0; //число узлов по х , // для артерии других зон нет
			loc_mesh->zone[6].y = 0; //число узлов по  у, // для артерии других зон нет
			loc_mesh->zone[6].z = 0; //число узлов по  у, // для артерии других зон нет

			loc_mesh->deltaCurZone[0].x = aa / nn; //шаг сетки по х , зона1
			loc_mesh->deltaCurZone[0].y = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона1	
			loc_mesh->deltaCurZone[0].z = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона1	
			loc_mesh->deltaCurZone[1].x = 2.0 * loc_mesh->deltaCurZone[0].x; //шаг сетки по х , зона2
			loc_mesh->deltaCurZone[1].y = loc_mesh->deltaCurZone[1].x; //шаг сетки по у, зона2
			loc_mesh->deltaCurZone[1].z = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона2
			loc_mesh->deltaCurZone[2].x = 4.0 * loc_mesh->deltaCurZone[0].x; //шаг сетки по х , зона3
			loc_mesh->deltaCurZone[2].y = loc_mesh->deltaCurZone[2].x; //шаг сетки по у, зона3
			loc_mesh->deltaCurZone[2].z = loc_mesh->deltaCurZone[1].z;// TEST loc_mesh->deltaCurZone[2].x; //шаг сетки по у, зона3
			loc_mesh->deltaCurZone[3].x = 8.0 * loc_mesh->deltaCurZone[0].x; //шаг сетки по х , зона4
			loc_mesh->deltaCurZone[3].y = 0.0; //шаг сетки по у, зона4
			loc_mesh->deltaCurZone[3].z = 0.0; //шаг сетки по у, зона4

			loc_mesh->hs = (int)floor(2.125 * nn);
			loc_mesh->hf = (int)floor(3.125 * nn);
			loc_mesh->vs = (int)floor((zone_y_3 + 2) * nn);
			loc_mesh->vf = (int)floor((zone_y_3 + 3) * nn);

			initInftyCylMeshZ(loc_mesh);

			m_phi = m_mesh->phi;
		};
		void initInftyCylMeshZ(LSSTAGMesh* _mesh)
		{
			_mesh->maxCountCutCells = nn * nn * zoneZ;//default
			_mesh->zs = 0;// (int)floor((zone_y_3 + 2) * nn);
			_mesh->zf = zoneZ;
			H_BASE = _mesh->deltaCurZone[2].z * (zoneZ - 2);
		};
		LSSTAGMeshMov* loc_mesh2;
		LSSTAGMesh* loc_mesh;
		double* m_phi;
	  public:
		coordinate P_C; ///< координаты центра профиля (в положении равновесия  )
		double zone_y_3;// для сетки --- 2.125 --- Cheny & Botella //task="flow"
		double ANGLE;///< угол атаки
		double aa;//единичный блок
		double dist;
	    int nn;///< число ячеек сетки на характерный размер профиля
		//отношение а/h --- для данного варианта сетки принимает только значения, кратные 4, начиная с 8
		AirfoilKind id_airfoil;///< идентификатор типа профиля
	  private:
		  int hs_loc, hf_loc, vs_loc, vf_loc, zs_loc, zf_loc;
	protected:
		  int zoneZ;
	}; 
}

#endif //_LS_STAG_FLOW_TASK_DESCRIPTION_H_
