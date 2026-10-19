/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ls_stag_mesh.h
* \author    Puzikova Valeria (vvp) 
* \date      23.10.2013 15:57:06
* \brief     ls_stag_mesh.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_MESH_H_
#define _LS_STAG_LS_STAG_MESH_H_

// ================================================================== INCLUDES
#include "ls_stag_uniform_mesh.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса LSSTAGMesh

Данный файл содержит в себе реализацию класса LSSTAGMesh
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL_2D(buildNodeX, LS_STAG_R_PC, _nodes, LS_STAG_R_PC, _delta, const double, curDx, 
		const double, curDy, const double, curDz, const double, curShiftX, const double, curShiftY, 
		const double, curShiftZ, const int, curCountX, const int, xStartNum, const int, yStartNum, const int, zStartNum,
		const int, _N, const int, _NM);
/*!
\brief Класс для работы с неравномерной LS-STAG-сеткой

Данная структура используется для работы с неравномерной  LS-STAG-сеткой в случае отсутствия движения погруженных границ. 
*/
	struct LSSTAGMesh : LSSTAGUniformMesh
	{
		static const int ZONES_COUNT = 7;
		static const int DELTA_COUNT = 4;
		static constexpr int ZONE_ORDER_X[ZONES_COUNT] = { 3, 2, 1, 0, 1, 2, 3 };// соответствие x-зон и deltaCurZone
		explicit LSSTAGMesh(MeshKinds _id = MESH_NONUNIFORM, int _size = sizeof(LSSTAGMesh)):
		LSSTAGUniformMesh(_id, _size),
		boundCount(1)
		{
			;
		}
    // Ядро заполнения координат узлов и дельт для одной зоны (теперь трёхмерное)
    LS_STAG_KERNEL_2D(buildNodeX, LS_STAG_R_PC _nodes, LS_STAG_R_PC _delta,
                   const double curDx, const double curDy, const double curDz,
                   const double curShiftX, const double curShiftY, const double curShiftZ, const int curCountX,
                   const int xStartNum, const int yStartNum, const int zStartNum, const int _N, const int _NM)
    {// i – индекс по y внутри зоны (от 0 до curCountY-1), j – по z внутри зоны (от 0 до curCountZ-1)
        const double y = (i + 0.5) * curDy + curShiftY;
        const double z = (j + 0.5) * curDz + curShiftZ;
        for (int k = 0, curNum = (j + zStartNum) * _NM + (i + yStartNum) * _N + xStartNum; k < curCountX; ++k, curNum++) 
		{
            _nodes[curNum].set((k + 0.5) * curDx + curShiftX, y, z);
            _delta[curNum].set(curDx, curDy, curDz);
        }
    }
		void calcZoneWidth()
		{//для артерии уже пересчитали Grid.Zone и delta + boundCount
			zoneWidth[0].x = zone[0].x * deltaCurZone[3].x; //ширина по х, зона4
			zoneWidth[0].y = zone[0].y * deltaCurZone[2].y; //ширина по у, зона3
			zoneWidth[1].x = zone[1].x * deltaCurZone[2].x; //ширина по х , зона3
			zoneWidth[1].y = zone[1].y * deltaCurZone[1].y; //ширина по  у, зона2
			zoneWidth[2].x = zone[2].x * deltaCurZone[1].x; //ширина по х , зона2
			zoneWidth[2].y = zone[2].y * deltaCurZone[0].y; //ширина по у, зона1
			zoneWidth[3].x = zone[3].x * deltaCurZone[0].x; //ширина по х , зона1
			zoneWidth[3].y = zone[3].y * deltaCurZone[1].y; //ширина по у, зона2
			zoneWidth[4].x = zone[4].x * deltaCurZone[1].x; //ширина по х , зона2
			zoneWidth[4].y = zone[4].y * deltaCurZone[2].y; //ширина по у, зона3
			zoneWidth[5].x = zone[5].x * deltaCurZone[2].x; //ширина по х , зона3
			zoneWidth[5].y = zone[5].y * deltaCurZone[1].y; //ширина по у, зона2!
			zoneWidth[6].x = zone[6].x * deltaCurZone[3].x; //ширина по х , зона4
			zoneWidth[6].y = zone[6].y * deltaCurZone[0].y; //ширина по у, зона1!
        // z – по аналогии с y (можно выбрать другое распределение)
        zoneWidth[0].z = zone[0].z * deltaCurZone[2].z;
        zoneWidth[1].z = zone[1].z * deltaCurZone[1].z;
        zoneWidth[2].z = zone[2].z * deltaCurZone[0].z;
        zoneWidth[3].z = zone[3].z * deltaCurZone[1].z;
        zoneWidth[4].z = zone[4].z * deltaCurZone[2].z;
        zoneWidth[5].z = zone[5].z * deltaCurZone[1].z;
        zoneWidth[6].z = zone[6].z * deltaCurZone[0].z;
			minStep = deltaCurZone[0].x;
		};
	  public:
		  void buildNode()// координаты узлов --- середины ячеек (вызываем 1 раз)
		  {
			  calcZoneWidth();
			  const int deltaOrder[2][ZONES_COUNT] = { {2,1,0,1,2,1,0},{0,1,2,1,0,1,2} };
			  const int zoneOrder[2][ZONES_COUNT] = { {0,1,2,3,4,5,6},{6,5,4,3,2,1,0} };
			  // Порядок зон для z – симметричен y
			  const int zoneOrderZ[2][ZONES_COUNT] = { {0,1,2,3,4,5,6}, {6,5,4,3,2,1,0} }; 
			  double curShiftZ = 0.0;
			  int zStartNum = 0;
			  for (int k2 = 0; k2 < boundCount; k2++)
			  {
				  const int num2 = k2 % 2;
				  for (int kz = 0; kz < ZONES_COUNT; curShiftZ += zoneWidth[zoneOrderZ[num2][kz]].z, kz++)
				  {
					  const int curCountZ = zone[zoneOrderZ[num2][kz]].z;
					  if (curCountZ)
					  {
						  const double curDz = deltaCurZone[deltaOrder[num2][kz]].z;
						  double curShiftY = 0.0;
						  int yStartNum = 0;
						  for (int k = 0; k < boundCount; k++)
						  {
							  const int num = k % 2;
							  for (int ky = 0; ky < ZONES_COUNT; curShiftY += zoneWidth[zoneOrder[num][ky]].y, ky++)
							  {
								  const int curCountY = zone[zoneOrder[num][ky]].y;
								  if (curCountY)
								  {
									  const double curDy = deltaCurZone[deltaOrder[num][ky]].y;
									  double curShiftX = 0.0;
									  for (int kx = 0, xStartNum = 0, curCountX = 0; kx < ZONES_COUNT; curShiftX += zoneWidth[kx].x, kx++,
										  xStartNum += curCountX)
									  {// координаты узлов --- середины ячеек --- заполнение одной зоны сетки по у
										  curCountX = zone[kx].x;
										  LS_STAG_PAR_FOR_2D(LSSTAGMesh, buildNodeX, curCountY, curCountZ, nodes, delta,
											  deltaCurZone[ZONE_ORDER_X[kx]].x, curDy, curDz, curShiftX, curShiftY, curShiftZ, curCountX,
											  xStartNum, yStartNum, zStartNum, baseMesh.N, m_NM);
									  }
									  yStartNum += curCountY;
								  }
							  }
						  }
						  zStartNum += curCountZ;
					  }
				  }
			  }
		  };
		/*!
		Индекс базовой ячейки по координате x (раскладка x-зон:
		zoneWidth[k].x = zone[k].x * deltaCurZone[ZONE_ORDER_X[k]].x)
		*/
		int xCellIndex(const double x)
		{
			double cur = 0.0;
			int idx = 0;
			for (int k = 0; k < ZONES_COUNT; k++)
			{
				const int cnt = zone[k].x;
				if (cnt <= 0) continue;
				const double d = deltaCurZone[ZONE_ORDER_X[k]].x;
				const double w = cnt * d;
				if (x < cur + w)
				{
					int i = idx + (int)floor((x - cur) / d);
					if (i < 0) i = 0;
					if (i > baseMesh.N - 1) i = baseMesh.N - 1;
					return i;
				}
				cur += w;
				idx += cnt;
			}
			return baseMesh.N - 1;
		};
		/*!
		Индекс базовой ячейки по координате y (раскладка y-зон:
		zoneWidth[k].y = zone[k].y * deltaCurZone[deltaOrderY[k]].y)
		*/
		int yCellIndex(const double y)
		{
			const int deltaOrderY[ZONES_COUNT] = { 2, 1, 0, 1, 2, 1, 0 };
			double cur = 0.0;
			int idx = 0;
			for (int k = 0; k < ZONES_COUNT; k++)
			{
				const int cnt = zone[k].y;
				if (cnt <= 0) continue;
				const double d = deltaCurZone[deltaOrderY[k]].y;
				const double w = cnt * d;
				if (y < cur + w)
				{
					int i = idx + (int)floor((y - cur) / d);
					if (i < 0) i = 0;
					if (i > baseMesh.M - 1) i = baseMesh.M - 1;
					return i;
				}
				cur += w;
				idx += cnt;
			}
			return baseMesh.M - 1;
		};
		coordinate zoneWidth[ZONES_COUNT];///< Ширина каждой зоны сетки
		coordinate deltaCurZone[DELTA_COUNT];///< Шаги сетки для каждой из зон
		tripletInt zone[ZONES_COUNT];
	  public:
		int boundCount;
		int hs, hf, vs, vf, zs, zf;
	}; 
}

#endif //_LS_STAG_LS_STAG_MESH_H_
