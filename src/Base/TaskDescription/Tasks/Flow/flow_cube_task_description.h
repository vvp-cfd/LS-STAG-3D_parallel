/*
* \copyright (C) 2026 Puzikova V.V.
* \file      flow_cube_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      05.10.2026
* \brief     flow_cube_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CUBE_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CUBE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_square_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCubeTaskDescription

Данный файл содержит в себе реализацию класса FlowCubeTaskDescription ---
обтекание куба: квадратное сечение nn x nn (как у FlowSquareTaskDescription),
высота --- единичная (nn слоёв по z), куб расположен по центру по z.
Тело занимает слои [zs, zf); под кубом (слой zs-1) --- T_Hexahedron4
(срез --- верхняя грань, тело выше), над кубом (слой zf) --- B_Hexahedron4 (срез --- нижняя грань).
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(calcAlphaCube, LS_STAG_R_PC, _alpha, const int, hs, const int, hf, const int, vs,
		const int, vf, const int, xN, const int, yM, const int, zK, const int, zs, const int, zf);
	LS_STAG_DECLARATION_S_KERNEL(calcThetaCube, LS_STAG_R_PC, _theta, const int, hs, const int, hf, const int, vs,
		const int, vf, const int, xN, const int, yM, const int, zK, const int, zs, const int, zf);
	LS_STAG_DECLARATION_S_KERNEL(calcTypeOfCellsCube, LS_STAG_R_PCT, _baseType, LS_STAG_R_P(coord<SideTypes>), _SideType,
		const int, hs, const int, hf, const int, vs, const int, vf, const int, zK, const int, xN, const int, yM,
		const int, N, const int, zs, const int, zf, LS_STAG_R_PI, _cutcellCounts, LS_STAG_R_P(tripletInt*), _airfoilNums);
/*!
\brief Класс, реализующий работу с задачей об обтекании куба

Данный класс реализует работу с описанием решаемой задачи об обтекании куба:
квадратное сечение берётся из FlowSquareTaskDescription, высота --- единичная
(nn слоёв по z), куб расположен по центру по z.
*/
	struct FlowCubeTaskDescription : FlowSquareTaskDescription
	{
	  public:
		FlowCubeTaskDescription():
			FlowSquareTaskDescription(AIRFOIL_CUBE, sizeof(FlowCubeTaskDescription))
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(calcAlphaCube, LS_STAG_R_PC _alpha, const int hs, const int hf, const int vs, const int vf,
			const int xN, const int yM, const int zK, const int zs, const int zf)
		{// куб: жижа (0..zs-2) -> плоскость zs-1 -> куб -> плоскость zf -> жижа
			int sqN = hf - hs, rightN = xN - hf, _M = yM + 1;
			coordinate* curT = _alpha;
			LSSTAGUniformMesh::fluidAlphaPlaneF(curT, xN, _M, zs - 1);
			{// плоскость под кубом (слой zs-1): клетки следа --- (0,0,1) (тело выше)
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::topAlphaRange(curT, sqN);
			LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, _M, vf);
			}
			for (int kAdd = zs; kAdd < zf; kAdd++)
			{// слои тела --- как внутренние слои квадрата
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, vs - 1);
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::nordAlphaRange(curT, sqN);
			LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
			LSSTAGUniformMesh::fluidAlpha(curT, hs - 1, 1);
			LSSTAGUniformMesh::solidThetaRange(curT, hs, hf);
			LSSTAGUniformMesh::fluidAlphaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidAlphaF(curT, xN, _M, vf);
			}
			LSSTAGUniformMesh::fluidAlphaPlaneF(curT, xN, _M, zK + 1, zf);
		}
		LS_STAG_SINGLE_KERNEL(calcThetaCube, LS_STAG_R_PC _theta, const int hs, const int hf, const int vs, const int vf,
			const int xN, const int yM, const int zK, const int zs, const int zf)
		{// куб: жижа (0..zs-2) -> плоскость zs-1 -> куб -> плоскость zf -> жижа
			int sqN = hf - hs, rightN = xN - hf, _M = yM + 1;
			coordinate* curT = _theta;
			LSSTAGUniformMesh::fluidThetaPlaneF(curT, xN, _M, zs - 1);
			{// плоскость под кубом (слой zs-1): клетки следа --- T-апертура 0 (тело выше)
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::topThetaRange(curT, sqN);
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, _M, vf);
			}
			for (int kAdd = zs; kAdd < zf; kAdd++)
			{// слои тела --- как внутренние слои квадрата
			LSSTAGUniformMesh::fluidThetaF(curT, xN, vs - 1);
			LSSTAGUniformMesh::fluidAlphaRange(curT, hs);
			LSSTAGUniformMesh::nordThetaRange(curT, sqN);
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
			LSSTAGUniformMesh::fluidTheta(curT, hs - 1, 1);
			LSSTAGUniformMesh::solidThetaRange(curT, hs, hf);
			LSSTAGUniformMesh::fluidThetaF(curT, rightN, 1);
			}
			LSSTAGUniformMesh::fluidThetaF(curT, xN, _M, vf);
			}
			LSSTAGUniformMesh::fluidThetaPlaneF(curT, xN, _M, zK + 1, zf);
		}
		LS_STAG_SINGLE_KERNEL(calcTypeOfCellsCube, LS_STAG_R_PCT _baseType, LS_STAG_R_P(coord<SideTypes>) _SideType,
			const int hs, const int hf, const int vs, const int vf, const int zK, const int xN, const int yM,
			const int N, const int zs, const int zf, LS_STAG_R_PI _cutcellCounts, LS_STAG_R_P(tripletInt*) _airfoilNums)
		{// куб: жижа (0..zs-2) -> плоскость zs-1 (T_Hexahedron4) -> куб -> плоскость zf (B_Hexahedron4) -> жижа
			const int M = yM + 1;
			int curNum = 0, rightN = xN - hf, nordY = yM - vf, sqN = hf - hs;
			coord<SideTypes>* curT = _SideType;
			// жидкая область: слои 0..zs-2
			LSSTAGUniformMesh::BottomPlaneTypeOfCell(curNum, _baseType, xN, yM);// B граница РО (слой 0)
			LSSTAGUniformMesh::FluidPlanesTypeOfCell(curNum, _baseType, xN, yM, zs - 1);// слои 1..zs-2
			LSSTAGUniformMesh::fluidSidesPlanes(curT, xN, yM, zs - 1);
			// плоскость под кубом (слой zs-1): жидкость + T_Hexahedron4 в следе
			LSSTAGUniformMesh::fluidSidesF(curT, xN, vs);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, SW_Edge, S_Border, SE_Edge, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, 1, vs, xN);// жидкая область
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
				LSSTAGUniformMesh::fluidSideRange(curT, hs);
				LSSTAGUniformMesh::topSidesRange(curT, sqN);
				LSSTAGUniformMesh::fluidSidesF(curT, rightN, 1);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
				for (int i = 1; i < hs; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				for (int i = hs; i < hf; i++)
				{
					BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, zs - 1, N, M);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, T_Hexahedron4);// срез --- верхняя грань (тело выше)
				}
				for (int i = hf; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
			}
			LSSTAGUniformMesh::fluidSidesF(curT, xN, yM - vf + 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, vf, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, NW_Edge, N_Border, NE_Edge, 0, 1, xN);// N граница РО
			// куб: слои [zs, zf) --- как внутренние слои квадрата
			for (int kAdd = zs; kAdd < zf; kAdd++)
			{
				LSSTAGUniformMesh::fluidSidesF(curT, xN, vs - 1);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, SW_Edge, S_Border, SE_Edge, 0, 1, xN);// S граница РО
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, 1, vs - 1, xN);// жидкая область
				LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
				LSSTAGUniformMesh::nordSidesRange(curT, hf, hs);
				LSSTAGUniformMesh::fluidSidesF(curT,rightN, 1);
				FlowSquareTaskDescription::SquareEdgeTypeOfCell(N_Hexahedron4, vs - 1, kAdd, curNum, hs, hf, xN, yM + 1,
					_baseType, _cutcellCounts, _airfoilNums);
				for (int jAdd = vs; jAdd < vf; jAdd++)
				{
					LSSTAGUniformMesh::fluidSides(curT, hs - 1, 1);
					LSSTAGUniformMesh::solidSidesRange(curT, hf, hs);
					LSSTAGUniformMesh::fluidSidesF(curT, rightN, 1);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
					for (int i = 1, hsn = hs - 1; i < hsn; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);
					BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, N, yM + 1);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Hexahedron4);//E
					for (int i = hs; i < hf; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, SolidCell);
					BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, kAdd, N, yM + 1);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Hexahedron4);//W
					for (int i = hf + 1; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
				}
				LSSTAGUniformMesh::fluidSidesPlanesF(curT, xN, nordY, 1);
				FlowSquareTaskDescription::SquareEdgeTypeOfCell(S_Hexahedron4, vf, kAdd, curNum, hs, hf, xN, yM + 1,
					_baseType, _cutcellCounts, _airfoilNums);
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, vf + 1, yM, xN);// жидкая область
				LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, NW_Edge, N_Border, NE_Edge, 0, 1, xN);// N граница РО
			}
			// плоскость над кубом (слой zf): жидкость + B_Hexahedron4 в следе
			LSSTAGUniformMesh::fluidSidesF(curT, xN, vs);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, SW_Edge, S_Border, SE_Edge, 0, 1, xN);// S граница РО
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, 1, vs, xN);// жидкая область
			for (int jAdd = vs; jAdd < vf; jAdd++)
			{
				LSSTAGUniformMesh::fluidSidesF(curT, xN, 1);
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, W_Border);
				for (int i = 1; i < hs; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				for (int i = hs; i < hf; i++)
				{
					BaseTaskDescription::addCutCell(0, _cutcellCounts, _airfoilNums, curNum, jAdd, zf, N, M);
					LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, B_Hexahedron4);// срез --- нижняя грань (тело ниже)
				}
				for (int i = hf; i < xN; i++) LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, FluidCell);//полностью жидкая
				LSSTAGUniformMesh::DefineTypeOfCell(curNum, _baseType, E_Border);//E
			}
			LSSTAGUniformMesh::fluidSidesF(curT, xN, yM - vf + 1);
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, W_Border, FluidCell, E_Border, vf, yM, xN);// жидкая область
			LSSTAGUniformMesh::RowTypeOfCell(curNum, _baseType, NW_Edge, N_Border, NE_Edge, 0, 1, xN);// N граница РО
			// жидкая область сверху: слои zf+1..zK-1 + T граница РО (слой zK)
			LSSTAGUniformMesh::FluidPlanesTypeOfCell(curNum, _baseType, xN, yM, zK - zf);
			LSSTAGUniformMesh::TopPlaneTypeOfCell(curNum, _baseType, xN, yM);
			LSSTAGUniformMesh::fluidSidesPlanes(curT, xN, yM, zK - zf - 1);
			LSSTAGUniformMesh::topSidesPlane(curT, xN, yM);
		}
	  protected:
		void initMeshParams() override
		{// высота --- единичная (nn слоёв), куб по центру по z
			FlowTaskDescription::initMeshParams();
			mesh.maxCountCutCells = 6 * (nn + 1) * nn;// поверхность куба: 4 боковых грани * nn слоёв + 2 граневых плоскости (запас (nn+1))
			mesh.zs = (mesh.baseMesh.K - nn) / 2;
			mesh.zf = mesh.zs + nn;
			H_BASE = L_BASE;// эталонная высота --- ребро куба
		}
		void buildAlpha() override
		{
			LS_STAG_SINGLE(FlowCubeTaskDescription, calcAlphaCube, mesh.alpha, mesh.hs, mesh.hf, mesh.vs, mesh.vf,
				mesh.xMesh.N, mesh.yMesh.M, mesh.zMesh.K, mesh.zs, mesh.zf)
		}
		void buildTheta() override
		{
			LS_STAG_SINGLE(FlowCubeTaskDescription, calcThetaCube, mesh.theta, mesh.hs, mesh.hf, mesh.vs, mesh.vf,
				mesh.xMesh.N, mesh.yMesh.M, mesh.zMesh.K, mesh.zs, mesh.zf)
		}
		void buildTypeOfCells() override
		{
			LS_STAG_SINGLE(FlowCubeTaskDescription, calcTypeOfCellsCube, mesh.m_baseType, mesh.baseSideTypes, mesh.hs,
				mesh.hf, mesh.vs, mesh.vf, mesh.zMesh.K, mesh.xMesh.N, mesh.yMesh.M, mesh.baseMesh.N, mesh.zs, mesh.zf,
				m_cutcellCounts, m_airfoilNums)
		}
	};
}

#endif //_LS_STAG_FLOW_CUBE_TASK_DESCRIPTION_H_
