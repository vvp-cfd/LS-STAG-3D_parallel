/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circles_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 23:36:40
* \brief     flow_circles_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLES_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLES_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "../flow_task_description.h"
#include "../../../../MeshFunctions/mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleSTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleSTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(cellsOnAir, LS_STAG_R_PI, _cellOnAirfoil, const int, airNum, const int, _xs,
		const int, _xf, const int, _N, const int, _NM);
	LS_STAG_DECLARATION_KERNEL(vib, LS_STAG_CR_PD, _radiusA, LS_STAG_CR_PD, _OmegaF,
		LS_STAG_CR_PI, _cellOnAirfoil, LS_STAG_CR_PCT, _baseType, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PC, _delta, LS_STAG_R_PC, _alpha,
		LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_yx, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT, const int, _N);
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetCircles, const int, hs, const int, hf, LS_STAG_R_PD, _phi, LS_STAG_CR_PC, _nodes, 
		LS_STAG_R_PC, _delta, LS_STAG_CR_PI, _cellOnAirfoil, const coordinate*, _cur_C, const int, _N, const int, _NM, const double*, _radiusA);
/*!
\brief Класс, реализующий работу с задачей об обтекании системы круговых профилей

Данный класс, реализующий работу с описанием решаемой задачи об обтекании системы круговых профилей. Профили могут быть
 неподвижными или вращаться как с постоянной угловой скоростью, так и совершая вращательные колебания
*/
	struct FlowCircleSTaskDescription : FlowTaskDescription
	{
		explicit FlowCircleSTaskDescription(MotionKinds _id_motion = MOTION_NO, int _size = sizeof(FlowCircleSTaskDescription)):
		FlowTaskDescription(AIRFOIL_CIRCLES, _id_motion, (LSSTAGUniformMesh*)&mesh, _size),
		rad(0.5),
		needFineMesh(false)
		{
			radius = 0.5 * L_BASE;
			Omega_0 = 2.0 * V_BASE / L_BASE;
			C_MAX.set(0.0);
			C_MIN.set(10000.0);
			discontin = 1;
			ibCount = 0;// меняем только при загрузке из файла
			#ifdef CLASS_NAME_AIRFOIL_NUM 
				#undef CLASS_NAME_AIRFOIL_NUM
			#endif
			#define CLASS_NAME_AIRFOIL_NUM FlowCircleSTaskDescription
		}
		~FlowCircleSTaskDescription()
		{
			if(ibCount > 0)
			{
				delete[] Ampl;
				delete[] Freq;
			}
		};
		LS_STAG_KERNEL_2D(cellsOnAir, LS_STAG_R_PI _cellOnAirfoil, const int airNum, const int _xs, const int _xf, const int _N, const int _NM)
		{
			for (int k = _xs, num = j * _NM + i * _N + _xs; k < _xf; k++, num++)
			{
				_cellOnAirfoil[num] = airNum;
			}//for k
		}
		LS_STAG_KERNEL(vib, LS_STAG_CR_PD _radiusA, LS_STAG_CR_PD _OmegaF, LS_STAG_CR_PI _cellOnAirfoil,
			LS_STAG_CR_PCT _baseType, LS_STAG_CR_P(coord<SideTypes>) _SideType, LS_STAG_CR_PC _delta, LS_STAG_R_PC _alpha, LS_STAG_R_PC _v_ib, LS_STAG_R_PC _v_ib_yx,
			LS_STAG_R_PC _v_ib_zy, LS_STAG_R_PC _v_ib_ENT, const int _N)
		{
			const int airNum = _cellOnAirfoil[i];
			if (fabs(_OmegaF[airNum]) > EQ_DOUBLE)// нулевая угловая скорость --- скорость IB нулевая
				LSSTAGUniformMesh::_k_projectionVib(i, _radiusA[airNum] * _OmegaF[airNum], _baseType, _SideType,
					_delta, _alpha, _v_ib, _v_ib_yx, _v_ib_zy, _v_ib_ENT, _N);
		}
		LS_STAG_KERNEL_2D(calcLevelSetCircles, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			LS_STAG_CR_PI _cellOnAirfoil, const coordinate* _cur_C, const int _N, const int _NM, const double* _radiusA)
		{
			int startNum = _NM * j + _N * i;
			int curNum = startNum + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = startNum + hf; curNum < end; curNum++)
			{
				coordinate pToC;
				const int airNum = _cellOnAirfoil[curNum];
				pToC.setDiffXY(curPoint, _cur_C[airNum]);
				_phi[curNum] = _radiusA[airNum] - pToC.norm();
				curPoint.x += curD->x;
			};
		}
		// ibCount --- в init считывается +
		virtual void rebuildMesh()// 25.12.13
		{
			for(int i = 0; i < ibCount; i++)
			{
				if(fabs(Freq[i]) > EQ_DOUBLE)//31.05.14
				     pOmegaF[i] = Ampl[i] * Omega_0 * sin(Freq[i] * (*m_curTime));
			}
			m_OmegaF.sync(pOmegaF);
			buildVib();
		};
		virtual void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> ibCount;
			Ampl = new double[ibCount];
			Freq = new double[ibCount];
			m_OmegaF.init(ibCount, OmegaF);
			m_OmegaF.setCopy(pOmegaF);
			m_radiusA.init(ibCount, radiusA);
			m_radiusA.setCopy(pRadiusA);
			m_Centers.init(ibCount, Centers);
			m_Centers.setCopy(pCenters);
			if(Re > 990) needFineMesh = true;// но кооррдинаты в файле уже должны быть с учетом этого
			if(needFineMesh)
			{
				zone_y_3 = 0.75;
				P_C.x = 4.0;
				P_C.y = 6.5;
				P_C.z = 6.5;// для цилиндра не используется
			}
			for(int i = 0; i < ibCount; i++)
			{
				in >> pCenters[i].x;
				in >> pCenters[i].y;
				in >> pCenters[i].z;
				in >> Ampl[i];
				in >> Freq[i];
				in >> pRadiusA[i];
				pOmegaF[i] = Ampl[i] * Omega_0;
				if(pCenters[i].x > C_MAX.x) C_MAX.x = pCenters[i].x;
				if(pCenters[i].x < C_MIN.x) C_MIN.x = pCenters[i].x;
				if(pCenters[i].y > C_MAX.y) C_MAX.y = pCenters[i].y;
				if(pCenters[i].y < C_MIN.y) C_MIN.y = pCenters[i].y;
			}
			m_radiusA.sync(pRadiusA); 
			m_OmegaF.sync(pOmegaF);
			m_Centers.sync(pCenters);
		};
		LS_STAG_KERNEL_INLINE int getAirfoilNum(const int _num, LS_STAG_CR_PI _cellOnAirfoil) { return _cellOnAirfoil ? _cellOnAirfoil[_num] : 0; };
		const int* getPCellOnAirfoill() { return m_cellOnAirfoil; }
	  protected:
		void calcNodesCount()
		{
			if(needFineMesh)
			{
				m_mesh->baseMesh.N = (int)floor((4 + C_MAX.x - C_MIN.x) * 4.0 * nn) + (int)(0.875 * nn) + 
				(int)floor((4.125 + 0.25 * (P_C.x - C_MAX.x)) * nn);//номер последнего узла по горизонтали
				m_mesh->baseMesh.M =  (int)floor((zone_y_3 + 0.25 * (C_MIN.y - P_C.y)) * nn) + (int)(1.5 * nn) + 
				(int)floor((zone_y_3 + 0.25 * (P_C.y - C_MAX.y)) * nn) +
				(int)floor((4 + C_MAX.y - C_MIN.y) * 4.0 * nn); //номер последнего узла по вертикали 
				m_mesh->baseMesh.K = zoneZ;/*(int)floor((zone_y_3 + 0.25 * (C_MIN.z - P_C.z)) * nn) + (int)(1.5 * nn) +
				(int)floor((zone_y_3 + 0.25 * (P_C.z - C_MAX.z)) * nn) +
				(int)floor((4 + C_MAX.z - C_MIN.z) * 4.0 * nn); //номер последнего узла по вертикали */

			}
			else
			{
				m_mesh->baseMesh.N = (int)floor((7.5 + 0.75 * C_MAX.x + 0.125 * P_C.x - 0.875 * C_MIN.x) * nn);//номер последнего узла по горизонтали
				m_mesh->baseMesh.M = (int)floor((2.0 * zone_y_3 + 5.0 + 0.75 * (C_MAX.y - C_MIN.y)) * nn); //номер последнего узла по вертикали 
				m_mesh->baseMesh.K = zoneZ; //(int)floor((2.0 * zone_y_3 + 5.0 + 0.75 * (C_MAX.z - C_MIN.z)) * nn); //номер последнего узла по вертикали 
			}
		};
		virtual void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << ibCount << "\n";
			for(int i = 0; i < ibCount; i++)
			{
				out << pCenters[i] << " " << Ampl[i] << " " << Freq[i] << " " << pRadiusA[i] << "\n";
			}
		};
		virtual void buildVib()
		{
			if(ibCount > 0)
			{
				LS_STAG_PAR_FOR(FlowCircleSTaskDescription, vib, mesh.NMK, radiusA, OmegaF, m_cellOnAirfoil, mesh.m_baseType, mesh.baseSideTypes, mesh.delta,
					mesh.alpha, mesh.v_ib, mesh.v_ib_yx, mesh.v_ib_zy, mesh.v_ib_ENT, mesh.baseMesh.N);
			}
		};
		virtual void buildBaseMesh()
		{
			FlowTaskDescription::buildBaseMesh();
			rebuildMesh();//buildVib();
		};
		/*!
		Вычисляет значения функции уровня в углах ячеек основной сетки
		\param zs номер ряда по z, с которого рассматриваем ячейки
		\param zf номер ряда по z, до которого рассматриваем ячейки
		\param vs номер ряда по вертикали, с которого рассматриваем ячейки
		\param vf номер ряда по вертикали, до которого рассматриваем ячейки
		\param hs номер ряда по горизонтали, с которого рассматриваем ячейки
		\param hf номер ряда по горизонтали, до которого рассматриваем ячейки
		*/
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf)
		{
			LS_STAG_PROFILE("buildLS");
			if (ibCount > 0) LS_STAG_PAR_FOR_2D_RANGE(FlowCircleSTaskDescription, calcLevelSetCircles, vs, vf, zs, zf, hs, hf, mesh.phi,
				m_mesh->nodes, m_mesh->delta, m_cellOnAirfoil, Centers, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, radiusA)
		};
		virtual void initMeshParams()
		{
			FlowTaskDescription::initMeshParams();
			if(needFineMesh)
			{
				loc_mesh->maxCountCutCells *= 16;
				loc_mesh->deltaCurZone[0].x = 0.25 * aa / nn; //шаг сетки по х , зона1
				loc_mesh->deltaCurZone[0].y = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона1
				loc_mesh->deltaCurZone[0].z = loc_mesh->deltaCurZone[0].x; //шаг сетки по у, зона1
				loc_mesh->zone[0].x = 0;//(int)floor((0.5 + 0.125 * (C_MIN.x - P_C.x)) * nn); //число узлов по х, зона4
				loc_mesh->zone[0].y = (int)floor((zone_y_3 + 0.25 * (C_MIN.y - P_C.y)) * nn); //число узлов по у, зона3
				loc_mesh->zone[1].y = (int)(0.75 * nn); //число узлов по  у, зона2
				loc_mesh->zone[0].z = (int)(0.5 * zoneZ); //(int)floor((zone_y_3 + 0.25 * (C_MIN.z - P_C.z)) * nn); //число узлов по у, зона3
				loc_mesh->zone[1].z = 0;// (int)(0.75 * nn); //число узлов по  у, зона2
				loc_mesh->zone[2].x = (int)floor(0.25 * nn); //число узлов по х , зона2
				loc_mesh->zone[2].y = (int)floor((4 + C_MAX.y - C_MIN.y) * 4.0 * nn); //число узлов по у, зона1
				loc_mesh->zone[2].z = 0;// (int)floor((4 + C_MAX.z - C_MIN.z) * 4.0 * nn); //число узлов по у, зона1
				loc_mesh->zone[3].x = (int)floor((4 + C_MAX.x - C_MIN.x) * 4.0 * nn); //число узлов по х , зона1
				loc_mesh->zone[3].y = (int)(0.75 * nn); //число узлов по  у, зона2
				loc_mesh->zone[3].z = 0;// (int)(0.75 * nn); //число узлов по  у, зона2
				loc_mesh->zone[4].x = (int)floor(0.25 * nn); //число узлов по х , зона2
				loc_mesh->zone[4].y = (int)floor((zone_y_3 + 0.25 * (P_C.y - C_MAX.y)) * nn); //число узлов по  у, зона3
				loc_mesh->zone[4].z = loc_mesh->zone[0].z;// (int)floor((zone_y_3 + 0.25 * (P_C.z - C_MAX.z)) * nn); //число узлов по  у, зона3
				loc_mesh->zone[5].x = (int)floor((4.125 + 0.25 * (P_C.x - C_MAX.x)) * nn); //число узлов по х , зона3


				loc_mesh->hs = (int)floor((0.625 + 1 + 0.125 * (C_MIN.x - P_C.x)) *  4 * nn);
				loc_mesh->hf = (int)floor((16.625 + 4.0 * C_MAX.x - 0.125 * P_C.x - 3.875 * C_MIN.x) *nn);
				loc_mesh->vs = (int)floor((zone_y_3 + 0.75 +  0.25 *(C_MIN.y - P_C.y)) * 4 * nn);//
				loc_mesh->vf = (int)floor((zone_y_3 + 16.75 + 4.0 * C_MAX.y - 0.25 * P_C.y - 3.75 * C_MIN.y) *nn);
				loc_mesh->zs = 0;// (int)floor((zone_y_3 + 0.75 + 0.25 * (C_MIN.z - P_C.z)) * 4 * nn);//
				loc_mesh->zf = zoneZ;// (int)floor((zone_y_3 + 16.75 + 4.0 * C_MAX.z - 0.25 * P_C.z - 3.75 * C_MIN.z) * nn);
			}
			else
			{
				loc_mesh->zone[0].x = (int)floor((0.5 + 0.125 * (C_MIN.x - P_C.x)) * nn); //число узлов по х, зона4
				loc_mesh->zone[0].y = (int)floor((zone_y_3 + 0.25 * (C_MIN.y - P_C.y)) * nn); //число узлов по у, зона3
				loc_mesh->zone[2].y = (int)floor((3 + C_MAX.y - C_MIN.y) * nn); //число узлов по у, зона1
				loc_mesh->zone[0].z = (int)(0.5 * zoneZ); // (int)floor((zone_y_3 + 0.25 * (C_MIN.z - P_C.z)) * nn); //число узлов по у, зона3
				loc_mesh->zone[2].z = 0;// (int)floor((3 + C_MAX.z - C_MIN.z) * nn); //число узлов по у, зона1
				loc_mesh->zone[3].x = (int)floor((2 + C_MAX.x - C_MIN.x) * nn); //число узлов по х , зона1
				loc_mesh->zone[4].y = (int)floor((zone_y_3 + 0.25 * (P_C.y - C_MAX.y)) * nn); //число узлов по  у, зона3
				loc_mesh->zone[4].z = loc_mesh->zone[0].z;// (int)floor((zone_y_3 + 0.25 * (P_C.z - C_MAX.z)) * nn); //число узлов по  у, зона3
				loc_mesh->zone[5].x = (int)floor((3.125 + 0.25 * (P_C.x - C_MAX.x)) * nn); //число узлов по х , зона3


				loc_mesh->hs = (int)floor((2.125 + 0.125 * (C_MIN.x - P_C.x)) * nn);
				loc_mesh->hf = (int)floor((3.125 + C_MAX.x - 0.125 * P_C.x - 0.875 * C_MIN.x) * nn);
				loc_mesh->vs = (int)floor((zone_y_3 + 2 + 0.25 * (C_MIN.y - P_C.y)) * nn);
				loc_mesh->vf = (int)floor((zone_y_3 + 3 + C_MAX.y - 0.25 * P_C.y - 0.75 * C_MIN.y) * nn);
				loc_mesh->zs = 0;// (int)floor((zone_y_3 + 2 + 0.25 * (C_MIN.z - P_C.z)) * nn);
				loc_mesh->zf = zoneZ;// (int)floor((zone_y_3 + 3 + C_MAX.z - 0.25 * P_C.z - 0.75 * C_MIN.z) * nn);
			}
			loc_mesh->base_length = L_BASE;
			cellOnAirfoil.init(loc_mesh->NMK, m_cellOnAirfoil);
			cellOnAirfoil.set(0);
			if(needFineMesh)
			{
				int coef = (ibCount > 2) ? 0 : 8;// 15.12.16 
				int coef2 = (ibCount > 2) ? 8 : 24;// 15.12.16 
				// аналитические номера ячеек по у: центры профилей лежат в мелкой зоне (zone[2]),
				// перед ней зоны 0 и 1 с шагами deltaCurZone[2] и deltaCurZone[1]
				const int z0y = loc_mesh->zone[0].y;
				const int z1y = loc_mesh->zone[1].y;
				const double d0y = loc_mesh->deltaCurZone[0].y;
				const double yFineStart = z0y * loc_mesh->deltaCurZone[2].y + z1y * loc_mesh->deltaCurZone[1].y;
				for(int airNumIdx = 0, lastIb = ibCount - 1; airNumIdx < ibCount; airNumIdx++)
				{//выставляем номера профилей
					int loc_hs = loc_mesh->hs - 1 + (int)floor((pCenters[airNumIdx].x - C_MIN.x) * 4 * nn) - coef * nn;
					int loc_hf = loc_hs + coef2 * nn + 3;
					if(loc_hs < loc_mesh->hs) loc_hs = loc_mesh->hs;
					if(loc_hf > loc_mesh->hf) loc_hf = loc_mesh->hf;
					// окно по у --- до середины до соседнего профиля: окна не должны перекрываться,
					// иначе второй профиль перезаписывает первый и цилиндр теряется
					int loc_vs = (airNumIdx > 0)
						? z0y + z1y + (int)floor((0.5 * (pCenters[airNumIdx - 1].y + pCenters[airNumIdx].y) - yFineStart) / d0y)
						: loc_mesh->vs;
					int loc_vf = (airNumIdx < lastIb)
						? z0y + z1y + (int)floor((0.5 * (pCenters[airNumIdx].y + pCenters[airNumIdx + 1].y) - yFineStart) / d0y)
						: loc_mesh->vf;
					if(loc_vs < loc_mesh->vs) loc_vs = loc_mesh->vs;
					if(loc_vf > loc_mesh->vf) loc_vf = loc_mesh->vf;
					int loc_zs = 0;// loc_mesh->zs - 1 + (int)floor((pCenters[airNumIdx].z - C_MIN.z) * 4 * nn) - coef * nn;
					int loc_zf = zoneZ;// loc_zs + coef2 * nn + 3;
					LS_STAG_PAR_FOR_2D_RANGE(FlowCircleSTaskDescription, cellsOnAir, loc_vs, loc_vf, loc_zs, loc_zf, m_cellOnAirfoil, airNumIdx,
						loc_hs, loc_hf, loc_mesh->baseMesh.N, loc_mesh->baseMesh.N * loc_mesh->baseMesh.M)
				}//for i
			}
			else
			{
for(int airNumIdx = 0; airNumIdx < ibCount; airNumIdx++)
			{//выставляем номера профилей
				int loc_hs = loc_mesh->hs - 1 + (int)floor((pCenters[airNumIdx].x - C_MIN.x) * nn);
				int loc_hf = loc_hs + nn + 3;
				if(loc_hs < loc_mesh->hs) loc_hs = loc_mesh->hs;
				if(loc_hf > loc_mesh->hf) loc_hf = loc_mesh->hf;
				int loc_vs = loc_mesh->vs - 1 + (int)floor((pCenters[airNumIdx].y - C_MIN.y) * nn);
				int loc_vf = loc_vs + nn + 3;
				if(loc_vs < loc_mesh->vs) loc_vs = loc_mesh->vs;
				if(loc_vf > loc_mesh->vf) loc_vf = loc_mesh->vf;
				int loc_zs = 0;// mesh.zs - 1 + (int)floor((pCenters[airNumIdx].z - C_MIN.z) * nn);
				int loc_zf = zoneZ;// loc_zs + nn + 3;
				LS_STAG_PAR_FOR_2D_RANGE(FlowCircleSTaskDescription, cellsOnAir, loc_vs, loc_vf, loc_zs, loc_zf, m_cellOnAirfoil, airNumIdx,
					loc_hs, loc_hf, loc_mesh->baseMesh.N, loc_mesh->baseMesh.N * loc_mesh->baseMesh.M)
			}//for i
			}
		};
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
		Array<double> m_radiusA, m_OmegaF;
		Array<coordinate> m_Centers;
		MeshFunction<int> cellOnAirfoil;
		double Omega_0;///< характерная угловая скорость
		double radius;///< характерный радиус
	  protected:
		double rad;
		int* m_cellOnAirfoil;
	  public:
		coordinate C_MAX;///< наибольшие координаты в массиве координат центров профилей
		coordinate C_MIN;///< наименьшие координаты в массиве координат центров профилей
		coordinate *Centers;///< массив координат центров профилей
		coordinate* pCenters;// на CPU
		double *Ampl;///< массив амплитуд угловых скоростей профилей
		double *Freq;///< массив частот вращательных колебаний профилей
		double *OmegaF;///< массив угловых скоростей профилей
		double *radiusA;///< массив радиусов профилей
		double *pRadiusA, *pOmegaF;// на CPU
	  protected:
		bool needFineMesh;
	};
}

#endif //_LS_STAG_FLOW_CIRCLES_TASK_DESCRIPTION_H_
