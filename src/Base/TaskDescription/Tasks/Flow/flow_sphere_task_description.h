/*
* \copyright (C) 2026 Puzikova V.V.
* \file      flow_sphere_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      15.10.2026
* \brief     flow_sphere_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_SPHERE_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_SPHERE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowSphereTaskDescription

Данный файл содержит в себе реализацию класса FlowSphereTaskDescription ---
обтекание сферы диаметра D (радиус rad = D/2 = 0.5) с центром в P_C.
Геометрия задаётся точной функцией уровня сферы (phi > 0 внутри):
varphi = rad - sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2) (см. tex, гл.1),
поэтому alpha/theta/типы усечённых ячеек строятся стандартной 3D-машинерией LS-STAG.
Нормировка сил --- по площади миделя S = pi*D^2/4 (H_BASE = pi*rad^2 при L_BASE = D = 1).
Эталон для сверки: Campregher et al. (COBEM 2005) и Cao & Tafti 2018 (Cd(Re), Re = 10..1000).
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetSphere, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM,
		const double, rad);
/*!
\brief Класс, реализующий работу с задачей об обтекании сферы

Данный класс реализует работу с описанием решаемой задачи об обтекании сферы
диаметра D = 2 rad = 1 с центром в P_C (классический 3D-кейс, тело без ориентации).
*/
	struct FlowSphereTaskDescription : FlowTaskDescription
	{
	  public:
		FlowSphereTaskDescription():
			FlowTaskDescription(AIRFOIL_SPHERE, MOTION_NO, (LSSTAGUniformMesh*)&mesh, sizeof(FlowSphereTaskDescription)),
			rad(0.5),
			Lz(8.0),
			useFinestMesh(false)// если сгущаем, то true
		{
			discontin = 0;// 3D-разрывное н.у. неработоспособно (см. FlowCubeTaskDescription)
			// finest-сетка (шаги зон различаются в 8 раз) вызывает стагнацию FGMRES
			// при малом окне рестарта: p-решатель упирается в лимит итераций и расчёт
			// расходится --- расширяем окно ТОЛЬКО для уравнения Пуассона (проверено:
			// 29 итераций против 200+); u/v/w сходятся и с окном 13
			FGMRES_RESTART_NUM_P = useFinestMesh ? DEFAULT_FGMRES_RESTART_NUM_P_FINEST : DEFAULT_FGMRES_RESTART_NUM;// окно Пуассона: finest требует 64, стандартная сетка --- 13
		}
		LS_STAG_KERNEL_2D(calcLevelSetSphere, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double rad)
		{// точная функция уровня сферы: phi > 0 внутри
			int startNum = _NM * j + _N * i;
			int curNum = startNum + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.5 * curD->z);
			for (int end = startNum + hf; curNum < end; curNum++)
			{
				const double qx = curPoint.x - P_C.x;
				const double qy = curPoint.y - P_C.y;
				const double qz = curPoint.z - P_C.z;
				_phi[curNum] = rad - sqrt(qx * qx + qy * qy + qz * qz);
				curD++;
				curPoint.x += curD->x;
			};
		}
		virtual double getDistance(const double x, const double y, const double z) override
		{// расстояние до поверхности сферы (для SA/SST); положительно в жидкости
			const double qx = x - P_C.x;
			const double qy = y - P_C.y;
			const double qz = z - P_C.z;
			return sqrt(qx * qx + qy * qy + qz * qz) - rad;
		};
		void load(std::ifstream &in) override
		{
			FlowTaskDescription::load(in);
			in >> Lz;
		};
	  protected:
		virtual void printExtParams(std::ofstream &out) override
		{
			FlowTaskDescription::printExtParams(out);
			out << Lz << "\n";
		};
		virtual void calcNodesCount() override
		{// z-домен: Lz (в диаметрах); finest-режим --- мелкий шаг 0.25*aa/nn в z (32 ячейки/D)
			const double dz = (useFinestMesh ? 0.25 : 1.0) * aa / nn;
			zoneZ = (int)ceil(Lz * aa / dz);
			if (zoneZ < 4) zoneZ = 4;
			if (zoneZ % 2) zoneZ++;// zone[0].z = zone[4].z = zoneZ/2 --- должно быть целым
			FlowTaskDescription::calcNodesCount();
			if (useFinestMesh)
			{// initFinestMesh меняет суммы зон x/y (как у роторов) --- компенсируем N и M
				mesh.baseMesh.N += (int)floor(9.5 * nn);
				mesh.baseMesh.M += (int)floor(9.0 * nn);
			}
			P_C.z = 0.5 * zoneZ * dz;// центр области по z
			// Сдвиг на пол-ячейки по y и z: при rad = k*dz сфера касалась бы узлов сетки полюсами
			// (phi = 0 в узлах), alpha/theta у приполюсных ячеек скачут и уравнения теряют
			// устойчивость (ср. сдвиг торцов у конечного цилиндра).
			P_C.y += 0.5 * dz;// мелкая y-зона
			P_C.z += 0.5 * dz;
		};
		virtual void initMeshParams() override
		{
			FlowTaskDescription::initMeshParams();
			{// запас на стандартный шаг aa/nn; в finest-режиме initFinestMesh умножит на 16
				mesh.maxCountCutCells = 2 * (int)ceil(4.0 * PI * rad * rad * nn * nn / (aa * aa)) + 10 * nn;
			}
			if (useFinestMesh)
			{// мелкий шаг 0.25*aa/nn в x/y + удлинение мелких зон (как у роторов)
				initFinestMesh();
				// z: initFinestMesh шаг по z не меняет --- измельчаем вручную
				// (zone[0].z и zone[4].z читают deltaCurZone[2].z, см. calcZoneWidth)
				mesh.deltaCurZone[0].z = 0.25 * aa / nn;
				mesh.deltaCurZone[1].z = 0.25 * aa / nn;
				mesh.deltaCurZone[2].z = 0.25 * aa / nn;
			}
			{// x-блок LS --- под тело (радиус rad)
				mesh.hs = mesh.xCellIndex(P_C.x - rad);
				mesh.hf = mesh.xCellIndex(P_C.x + rad) + 1;
				if (mesh.hs < 2) mesh.hs = 2;
				if (mesh.hf > mesh.baseMesh.N - 3) mesh.hf = mesh.baseMesh.N - 3;
			}
			{// y-блок LS
				mesh.vs = mesh.yCellIndex(P_C.y - rad);
				mesh.vf = mesh.yCellIndex(P_C.y + rad) + 1;
				if (mesh.vs < 2) mesh.vs = 2;
				if (mesh.vf > mesh.baseMesh.M - 3) mesh.vf = mesh.baseMesh.M - 3;
			}
			// Нормировка S = pi*rad^2 (мидель сферы): fcoef = 2/(rho*V^2*L_BASE*H_BASE), L_BASE = 1
			H_BASE = PI * rad * rad;
		};
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf) override
		{// phi --- 3D (зависит от z), z-диапазон --- полная область
			LS_STAG_PROFILE("buildLS");
			LS_STAG_PAR_FOR_2D_RANGE(FlowSphereTaskDescription, calcLevelSetSphere, vs, vf, zs, zf, hs, hf, mesh.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad);
		};
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
		double rad;///< радиус сферы (D = 2 rad = 1)
	  public:
		double Lz;///< размер области по z (в диаметрах)
		bool useFinestMesh;///< true --- finest-сетка (шаг тела 0.25*aa/nn = 4*nn ячеек/D, как у роторов); false --- обычная сетка aa/nn
	}; 
}

#endif //_LS_STAG_FLOW_SPHERE_TASK_DESCRIPTION_H_
