/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_square_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 18:15:30
* \brief     flow_square_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_SQUARE_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_SQUARE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowSquareTaskDescription

Данный файл содержит в себе реализацию класса FlowSquareTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(calcAlpha, LS_STAG_R_PC, _alpha, const int, hs, const int, hf,
		const int, vs, const int, vf, const int, xN, const int, yM, const int, zK);
	LS_STAG_DECLARATION_S_KERNEL(calcTheta, LS_STAG_R_PC, _theta, const int, hs, const int, hf,
		const int, vs, const int, vf, const int, xN, const int, yM, const int, zK);
	LS_STAG_DECLARATION_S_KERNEL(calcTypeOfCellsSq, LS_STAG_R_PCT, _baseType, LS_STAG_R_P(coord<SideTypes>), _SideType, 
		const int, hs, const int, hf, const int, vs, const int, vf, const int, zK, const int, xN, const int, yM, const int, N,
		LS_STAG_R_PI, _cutcellCounts, LS_STAG_R_P(tripletInt*), _airfoilNums);
/*!
\brief Класс, реализующий работу с задачей об обтекании квадратного профиля

Данный класс реализует работу с описанием решаемой задачи об обтекании квадратного профиля
*/
	struct FlowSquareTaskDescription : FlowTaskDescription
	{
	  public:
		FlowSquareTaskDescription(AirfoilKind _id_airfoil = AIRFOIL_SQUARE, int _size = sizeof(FlowSquareTaskDescription)):
			FlowTaskDescription(_id_airfoil, MOTION_NO, (LSSTAGUniformMesh*)&mesh, _size)
		{
			discontin = 1;
		}
		LS_STAG_SINGLE_KERNEL(calcAlpha, LS_STAG_R_PC _alpha, const int hs, const int hf, const int vs, const int vf,
			const int xN, const int yM, const int zK)
		{// квадратный цилиндр
			coordinate* curT = _alpha; 
			int sqN = hf - hs, rightN = xN - hf;
			for(int k = -1; k < zK; k++)
			{
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs - 1);
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::nordAlphaRange(curT, sqN);
			LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			for (int j = vs; j < vf; j++)
			{
			LSSTAGUniformMesh::fluidAlpha(curT, hs - 1, 1);
			LSSTAGUniformMesh::solidThetaRange(curT, hs, hf);
			LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, yM+1, vf);
			}
		}
		LS_STAG_SINGLE_KERNEL(calcTheta, LS_STAG_R_PC _theta, const int hs, const int hf, const int vs, const int vf,
			const int xN, const int yM, const int zK)
		{// квадратный цилиндр
			coordinate* curT = _theta; 
			int sqN = hf - hs, rightN = xN - hf;
			for(int k = -1; k < zK; k++)
			{
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs - 1);
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::nordThetaRange(curT, sqN);
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			for (int j = vs; j < vf; j++)
			{
			LSSTAGUniformMesh::fluidTheta(curT, hs - 1, 1);
			LSSTAGUniformMesh::solidThetaRange(curT, hs, hf);
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, yM+1, vf);
			}
		}
		LS_STAG_SINGLE_KERNEL(calcTypeOfCellsSq, LS_STAG_R_PCT _baseType, LS_STAG_R_P(coord<SideTypes>) _SideType, const int hs, const int hf,
			const int vs, const int vf, const int zK, const int xN, const int yM, const int N, LS_STAG_R_PI _cutcellCounts,
			LS_STAG_R_P(tripletInt*) _airfoilNums)
		{
			int curNum = 0, rightN = xN - hf, nordY = yM - vf;
			coord<SideTypes>* curT = _SideType;
			LSSTAGUniformMesh::fluidSidesF(curT, xN, vs - 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BSW_Corner, BS_Edge, BSE_Corner, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, 1, vs - 1, xN);// жидкая область
			LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
		    LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, B_Border, 0, 1, hs - 1);// жидкая область
			for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, BN_Edge);
			LSSTAGUniformMesh::nordSidesRange(curT, hf, hs);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, B_Border, B_Border, BE_Edge, 0, 1, rightN);// жидкая область
			LSSTAGUniformMesh::fluidSidesF(curT,rightN, 1);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, 0, 1, hs - 1);// жидкая область
				LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
				for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, SolidCell);
				LSSTAGUniformMesh::solidSidesRange(curT, hf, hs);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, 0, 1, rightN);// жидкая область
				LSSTAGUniformMesh::fluidSidesF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, B_Border, 0, 1, hs - 1);// жидкая область
			for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, BS_Edge);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, B_Border, B_Border, BE_Edge, 0, 1, rightN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BW_Edge, B_Border, BE_Edge, vf + 1, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, BNW_Corner, BN_Edge, BNE_Corner, 0, 1, xN);// N граница РО
			for (int kAdd = 1; kAdd < zK; kAdd++)
			{
				LSSTAGUniformMesh::fluidSidesF(curT, xN, vs - 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, SW_Edge, S_Border, SE_Edge, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, 1, vs - 1, xN);// жидкая область
			LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
			LSSTAGUniformMesh::nordSidesRange(curT, hf, hs);
			LSSTAGUniformMesh::fluidSidesF(curT,rightN, 1);
			FlowSquareTaskDescription::SquareEdgeTypeOfCell(N_Hexahedron4, vs - 1, kAdd, curNum, hs, hf, xN, yM + 1, _baseType,
				_cutcellCounts, _airfoilNums);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
				LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
				LSSTAGUniformMesh::solidSidesRange(curT, hf, hs);
				LSSTAGUniformMesh::fluidSidesF(curT, rightN, 1);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
				for (int i = 1, hsn = hs - 1; i < hsn; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, N, yM + 1);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Hexahedron4);//E
				for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, SolidCell);
				BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, N, yM + 1);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Hexahedron4);//W
				for (int i = hf + 1; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			FlowSquareTaskDescription::SquareEdgeTypeOfCell(S_Hexahedron4, vf, kAdd, curNum, hs, hf, xN, yM + 1, _baseType, _cutcellCounts, _airfoilNums);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, vf + 1, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, NW_Edge, N_Border, NE_Edge, 0, 1, xN);// N граница РО
		}
		LSSTAGUniformMesh::fluidSidesF(curT, xN, vs - 1);
		LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TSW_Corner, TS_Edge, TSE_Corner, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 1, vs - 1, xN);// жидкая область	
			LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, hs - 1);// жидкая область
			for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, TN_Edge);
			LSSTAGUniformMesh::nordSidesRange(curT, hf, hs);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, rightN);// жидкая область
			LSSTAGUniformMesh::fluidSidesF(curT,rightN, 1);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, hs - 1);// жидкая область
				LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
				for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, SolidCell);
				LSSTAGUniformMesh::solidSidesRange(curT, hf, hs);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, rightN);// жидкая область
				LSSTAGUniformMesh::fluidSidesF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, hs - 1);// жидкая область
			for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, TS_Edge);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, 0, 1, rightN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TW_Edge, T_Border, TE_Edge, vf + 1, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, TNW_Corner, TN_Edge, TNE_Corner, 0, 1, xN);// N граница РО
		}
		LS_STAG_KERNEL_INLINE void SquareEdgeTypeOfCell(const CellTypes type, const int jAdd, const int kAdd, int& curNum, 
			const int hs, const int hf, const int xN, const int M, LS_STAG_R_PCT _baseType, LS_STAG_R_PI _cutcellCounts, LS_STAG_R_P(tripletInt*) _airfoilNums)
		{
			LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
			for (int i = 1; i < hs; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
			for (int i = hs; i < hf; i++)
			{
				BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, xN + 1, M);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, type);
			}
			for (int i = hf; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
			LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
		};
	  protected:
		void initMeshParams()
		{
			FlowTaskDescription::initMeshParams();
			mesh.maxCountCutCells = 4 * (nn + 1) * zoneZ;// check
		};
		void buildBaseMesh()
		{
			mesh.buildNode();
			BaseTaskDescription::buildBaseMesh();
		};
		void buildAlpha()
		{
			LS_STAG_SINGLE(FlowSquareTaskDescription, calcAlpha, mesh.alpha, mesh.hs, mesh.hf, mesh.vs, mesh.vf,
				mesh.xMesh.N, mesh.yMesh.M, mesh.zMesh.K)
		};
		void buildTheta()
		{
			LS_STAG_SINGLE(FlowSquareTaskDescription, calcTheta, mesh.theta, mesh.hs, mesh.hf, mesh.vs, mesh.vf,
				mesh.xMesh.N, mesh.yMesh.M, mesh.zMesh.K)
		};
		void buildTypeOfCells()
		{
			LS_STAG_SINGLE(FlowSquareTaskDescription, calcTypeOfCellsSq, mesh.m_baseType, mesh.baseSideTypes, mesh.hs, 
				mesh.hf, mesh.vs, mesh.vf, mesh.zMesh.K, mesh.xMesh.N, mesh.yMesh.M, mesh.baseMesh.N, m_cutcellCounts, 
				m_airfoilNums)
		};
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
	};
}

#endif //_LS_STAG_FLOW_SQUARE_TASK_DESCRIPTION_H_
