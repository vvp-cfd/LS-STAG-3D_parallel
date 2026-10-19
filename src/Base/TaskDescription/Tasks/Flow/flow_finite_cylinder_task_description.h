/*
* \copyright (C) 2026 Puzikova V.V.
* \file      flow_finite_cylinder_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      08.10.2026
* \brief     flow_finite_cylinder_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_FINITE_CYLINDER_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_FINITE_CYLINDER_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowFiniteCylinderTaskDescription

Данный файл содержит в себе реализацию класса FlowFiniteCylinderTaskDescription ---
обтекание цилиндра конечной длины с плоскими торцами. Ось тела лежит в плоскости x-y
и составляет угол ANGLE с осью x: e = (cos ANGLE, sin ANGLE, 0) (как в Cao & Tafti 2018,
Pierson et al. 2019 --- подъёмная сила в y, т.е. наш C_y). ANGLE = pi/2 --- поперечный
обдув (ось вдоль y, свободные торцы), ANGLE = 0 --- осевой. Осевой размах тела уходит
в широкий y-домен (мелкая зона 3D, запас ~10D), а z --- поперечное обтекание: размер Lz
задаётся в диаметрах с запасом ~2-3D с каждой стороны. Геометрия задаётся точной функцией
уровня "capped cylinder" (phi > 0 внутри), поэтому alpha/theta/типы усечённых ячеек
строятся стандартной 3D-машинерией LS-STAG. Длина L = ar * D (D = 2 rad = 1),
нормировка сил --- S = L * D (H_BASE = L).
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetFiniteCylinder, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM,
		const double, rad, const double, halfLength, const double, ex, const double, ey, const double, ez);
/*!
\brief Класс, реализующий работу с задачей об обтекании цилиндра конечной длины

Данный класс реализует работу с описанием решаемой задачи об обтекании цилиндра
конечной длины с плоскими торцами (в т.ч. наклонного --- ось в плоскости x-y)
*/
	struct FlowFiniteCylinderTaskDescription : FlowTaskDescription
	{
	  public:
		FlowFiniteCylinderTaskDescription():
			FlowTaskDescription(AIRFOIL_FINITE_CYLINDER, MOTION_NO, (LSSTAGUniformMesh*)&mesh, sizeof(FlowFiniteCylinderTaskDescription)),
			rad(0.5),
			ex(0.0),
			ey(1.0),
			ez(0.0),
			ar(2.0),
			Lz(5.0)
		{
			ANGLE = PI_2;// по умолчанию --- поперечный обдув (ось вдоль y)
			discontin = 0;// 3D-разрывное н.у. неработоспособно (см. FlowCubeTaskDescription)
			length = ar * 2.0 * rad;
			setAxis();
		}
		LS_STAG_KERNEL_2D(calcLevelSetFiniteCylinder, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double rad, const double halfLength, const double ex,
			const double ey, const double ez)
		{// точная функция уровня конечного цилиндра (плоские торцы): phi > 0 внутри
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
				const double t = qx * ex + qy * ey + qz * ez;// продольная координата (вдоль оси)
				const double rx = qx - t * ex;
				const double ry = qy - t * ey;
				const double rz = qz - t * ez;
				const double r = sqrt(rx * rx + ry * ry + rz * rz);// расстояние до оси
				const double d1 = r - rad;
				const double d2 = fabs(t) - halfLength;
				const double ox = (d1 > 0.0) ? d1 : 0.0;
				const double oy = (d2 > 0.0) ? d2 : 0.0;
				const double inMax = (d1 > d2) ? d1 : d2;
				_phi[curNum] = -(sqrt(ox * ox + oy * oy) + ((inMax < 0.0) ? inMax : 0.0));
				curD++;
				curPoint.x += curD->x;
			};
		}
		virtual double getDistance(const double x, const double y, const double z) override
		{// расстояние до поверхности цилиндра (для SA/SST); отрицательно внутри тела
			const double qx = x - P_C.x;
			const double qy = y - P_C.y;
			const double qz = z - P_C.z;
			const double t = qx * ex + qy * ey + qz * ez;
			const double rx = qx - t * ex;
			const double ry = qy - t * ey;
			const double rz = qz - t * ez;
			const double r = sqrt(rx * rx + ry * ry + rz * rz);
			const double d1 = r - rad;
			const double d2 = fabs(t) - 0.5 * length;
			const double ox = (d1 > 0.0) ? d1 : 0.0;
			const double oy = (d2 > 0.0) ? d2 : 0.0;
			const double inMax = (d1 > d2) ? d1 : d2;
			return -(sqrt(ox * ox + oy * oy) + ((inMax < 0.0) ? inMax : 0.0));
		};
		void load(std::ifstream &in) override
		{
			FlowTaskDescription::load(in);
			in >> ar;
			in >> ANGLE;
			in >> Lz;
			length = ar * 2.0 * rad;
			setAxis();
		};
	  protected:
		void setAxis()
		{// ось единичной длины в плоскости x-y: e = (cos ANGLE, sin ANGLE, 0)
			ex = cos(ANGLE);
			ey = sin(ANGLE);
			ez = 0.0;
		};
		virtual void printExtParams(std::ofstream &out) override
		{
			FlowTaskDescription::printExtParams(out);
			out << ar << "\n";
			out << ANGLE << "\n";
			out << Lz << "\n";
		};
		virtual void calcNodesCount() override
		{// z-домен: Lz (в диаметрах) при равномерном шаге dz = aa/nn (см. initMeshParams базовой задачи)
			zoneZ = (int)ceil(Lz * nn / aa);
			if (zoneZ < 4) zoneZ = 4;
			if (zoneZ % 2) zoneZ++;// zone[0].z = zone[4].z = zoneZ/2 --- должно быть целым
			FlowTaskDescription::calcNodesCount();
			P_C.z = 0.5 * zoneZ * (aa / nn);// центр области по z
			// Торцы при ANGLE=pi/2 (ось вдоль y) сдвигаем на пол-ячейки: при точном совпадении
			// плоскости торца с гранью ячеек phi в углах на торце вырождается (phi = rad - r
			// вместо 0), alpha/theta у приторцевых ячеек скачут и уравнения теряют устойчивость.
			if (fabs(ey) > 1.0 - EQ_DOUBLE) P_C.y += 0.5 * aa / nn;// мелкая y-зона
		};
		virtual void initMeshParams() override
		{
			length = ar * 2.0 * rad;// L = ar * D
			setAxis();
			FlowTaskDescription::initMeshParams();
			{// x-блок LS --- под проекцию тела на x (полудлина (L/2)|ex| + rad*sqrt(1-ex^2))
				const double halfX = 0.5 * length * fabs(ex) + rad * sqrt(1.0 - ex * ex);
				mesh.hs = mesh.xCellIndex(P_C.x - halfX);
				mesh.hf = mesh.xCellIndex(P_C.x + halfX) + 1;
				if (mesh.hs < 2) mesh.hs = 2;
				if (mesh.hf > mesh.baseMesh.N - 3) mesh.hf = mesh.baseMesh.N - 3;
			}
			{// y-блок LS --- под проекцию тела на y (полудлина (L/2)|ey| + rad*sqrt(1-ey^2))
				const double halfY = 0.5 * length * fabs(ey) + rad * sqrt(1.0 - ey * ey);
				mesh.vs = mesh.yCellIndex(P_C.y - halfY);
				mesh.vf = mesh.yCellIndex(P_C.y + halfY) + 1;
				if (mesh.vs < 2) mesh.vs = 2;
				if (mesh.vf > mesh.baseMesh.M - 3) mesh.vf = mesh.baseMesh.M - 3;
			}
			H_BASE = length;// нормировка S = L * D
		};
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf) override
		{// phi --- 3D (зависит от z), z-диапазон --- полная область
			LS_STAG_PROFILE("buildLS");
			LS_STAG_PAR_FOR_2D_RANGE(FlowFiniteCylinderTaskDescription, calcLevelSetFiniteCylinder, vs, vf, zs, zf, hs, hf, mesh.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad,
				0.5 * length, ex, ey, ez)
		};
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
		double rad;///< радиус цилиндра (D = 2 rad = 1)
		double length;///< длина цилиндра L = ar * D
		double ex, ey, ez;///< компоненты единичной оси тела (e = (cos ANGLE, sin ANGLE, 0))
	  public:
		double ar;///< удлинение L/D
		double Lz;///< размер области по z (в диаметрах, шаг dz = aa/nn): поперечное обтекание, запас ~2-3D с каждой стороны
	}; 
}

#endif //_LS_STAG_FLOW_FINITE_CYLINDER_TASK_DESCRIPTION_H_
