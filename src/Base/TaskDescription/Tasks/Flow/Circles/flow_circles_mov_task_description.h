/*
* \copyright (C) 2015 Puzikova V.V.
* \file      flow_circles_mov_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2015 17:31:40
* \brief     flow_circles_mov_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLES_MOV_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLES_MOV_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circles_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleSMovTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleSMovTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(build, LS_STAG_R_PC, _v_ib, LS_STAG_R_PC, _v_ib_zy, LS_STAG_R_PC, _v_ib_ENT,
		LS_STAG_R_PC, _v_ib_yx, LS_STAG_CR_PCT, _type, LS_STAG_CR_P(coord<SideTypes>), _SideType, LS_STAG_CR_PI, _OnAirfoil, LS_STAG_CR_PC, _curV_C);
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetCirclesMov, const int, hs, const int, hf, LS_STAG_R_PD, _phi, LS_STAG_CR_PC, _nodes,
		LS_STAG_R_PC, _delta, LS_STAG_CR_PI, _cellOnAirfoil, const coordinate*, _cur_C, const int, _N, const int, _NM, const double, rad);
/*!
\brief Класс, реализующий работу с задачей о вязкоупругих колебаниях системы круговых профилей

Данный класс, реализующий работу с описанием решаемой задачи об обтекании системы круговых профилей. Профили могут быть
 неподвижными или совершать вязкоупругие колебания
*/
	struct FlowCircleSMovTaskDescription : FlowCircleSTaskDescription
	{
		FlowCircleSMovTaskDescription():
		FlowCircleSTaskDescription(MOTION_MOV, sizeof(FlowCircleSMovTaskDescription)),
		MASS(DEFAULT_CIRCLES_MOV_MASS),// default
		COEF_DEMP(DEFAULT_CIRCLES_MOV_COEF_DEMP),//default
		S0(0.24),
		freqRatio(1.0),
		Se(0.24),// default --- резонансная частота
		t_0(1.0),//(0.01),//(1.0),
		fileExist(false)
		{
			m_mesh = (LSSTAGUniformMesh*)&meshMov;// изменяем mesh
			loc_mesh = &meshMov;
			loc_mesh2 = &meshMov;
			Re = 1000.0;
		}
		~FlowCircleSMovTaskDescription()
		{
			if(fileExist)
			{
				for(int i = 0; i < ibCount; i++)
				{
					yFile[i].close();
					vFile[i].close();
					xFile[i].close();
					uFile[i].close();
				}
			}
			if(ibCount)
			{
				delete [] y_n; 
				delete [] y_n_1;
				delete [] x_n; 
				delete [] x_n_1;
				delete [] coef4X;
				delete [] coef4Y;
				delete [] yFilename;
				delete [] vFilename;
				delete [] xFilename;
				delete [] uFilename;
				delete [] yFile;
				delete [] vFile;
				delete [] xFile;
				delete [] uFile;
			}
		};
		LS_STAG_KERNEL(build, LS_STAG_R_PC _v_ib, LS_STAG_R_PC _v_ib_zy, LS_STAG_R_PC _v_ib_ENT, LS_STAG_R_PC _v_ib_yx,
			LS_STAG_CR_PCT _type, LS_STAG_CR_P(coord<SideTypes>) _SideType, LS_STAG_CR_PI _OnAirfoil, LS_STAG_CR_PC _curV_C)
		{
			const CellTypes curType = _type[i];
			coordinate curVel;
			if ((curType > SolidCell) && (curType < FluidCell))
			{
			curVel.set(_curV_C + _OnAirfoil[i]);
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
		LS_STAG_KERNEL_2D(calcLevelSetCirclesMov, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			LS_STAG_CR_PI _cellOnAirfoil, const coordinate* _cur_C, const int _N, const int _NM, const double rad)
		{
			int startNum = _NM * j + _N * i;
			int curNum = startNum + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = startNum + hf; curNum < end; curNum++)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, _cur_C[_cellOnAirfoil[curNum]]);
				_phi[curNum] = rad - pToC.norm();
				curD++;
				curPoint.x += curD->x;
			};
		}
		void checkExistMov()// 1.06.14
		{
			if(!existMov)
			{
				if((*m_curTime) > t_0) existMov = true;
			}
		};
		void initCenter(){};
		void load(std::ifstream &in)
		{
			FlowCircleSTaskDescription::load(in);
			y_n = new double[ibCount]; 
			y_n_1 = new double[ibCount]; 
			x_n = new double[ibCount]; 
			x_n_1 = new double[ibCount];
			coef4X = new double[ibCount];
			coef4Y = new double[ibCount];
			yFilename = new std::tstring[ibCount];
			vFilename = new std::tstring[ibCount];
			xFilename = new std::tstring[ibCount];
			uFilename = new std::tstring[ibCount];
			yFile = new std::ofstream[ibCount];
			vFile = new std::ofstream[ibCount];
			xFile = new std::ofstream[ibCount];
			uFile = new std::ofstream[ibCount];
			meshMov.defineIbCount(ibCount);
			in >> S0;
			in >> freqRatio;
			in >> MASS;
			in >> COEF_DEMP;
			for(int i = 0; i < ibCount; i++)
			{
				pCenters[i].z = 0.0;
				(meshMov.pCur_C + i)->set(pCenters + i);
			}
			meshMov.m_cur_C.sync(meshMov.pCur_C);
		};
		void saveState(std::ofstream &out)
		{
			BaseTaskDescription::saveState(out);
			for(int i = 0; i < ibCount; i++)
			{
				out.write((char*)(y_n + i), sizeof(double));
				out.write((char*)(y_n_1 + i), sizeof(double));
				out.write((char*)(x_n + i), sizeof(double));
				out.write((char*)(x_n_1 + i), sizeof(double));
			}
		};
		void loadState(std::ifstream &in)
		{
			BaseTaskDescription::loadState(in);
			for(int i = 0; i < ibCount; i++)
			{
				in.read((char*)(y_n + i), sizeof(double));
				in.read((char*)(y_n_1 + i), sizeof(double));
				in.read((char*)(x_n + i), sizeof(double));
				in.read((char*)(x_n_1 + i), sizeof(double));
				(meshMov.pCur_C + i)->set(x_n[i], y_n[i], 0.0);
			}
			meshMov.m_cur_C.sync(meshMov.pCur_C);
		};
		void printLoads()
		{
			BaseTaskDescription::printLoads();
			for(int i = 0; i < ibCount; i++)
			{
				yFile[i] << meshMov.pCur_C[i].y << "\n";
				vFile[i] << meshMov.pCurV_C[i].y << "\n";
				xFile[i] << meshMov.pCur_C[i].x << "\n";
				uFile[i] << meshMov.pCurV_C[i].x << "\n";
			}
		};
		/*!
		Вычисляет текущие координаты центров профилей и их скорости (решает уравнения движения)
		*/
		void calcPosition()
		{
			if((*m_curTime) > t_0)
			{
				for(int i = 0; i < ibCount; i++)
				{
					meshMov.pCur_C[i].y = coef4Y[i] + coef3 * m_loads[i].F_aer_y + coef2 * y_n_1[i] + coef * y_n[i];
					double curDev = meshMov.pCur_C[i].y - y_n[i];
					meshMov.pCurV_C[i].y  = curDev / t;
					y_n_1[i] = y_n[i];
					y_n[i] = meshMov.pCur_C[i].y;
					meshMov.pCur_C[i].x = coef4X[i] + coef3 * m_loads[i].F_aer_x + coef2 * x_n_1[i] + coef * x_n[i];
					curDev = meshMov.pCur_C[i].x - x_n[i];
					meshMov.pCurV_C[i].x  = curDev / t;
					x_n_1[i] = x_n[i];
					x_n[i] = meshMov.pCur_C[i].x;
				}
				meshMov.m_cur_C.sync(meshMov.pCur_C);
				meshMov.m_curV_C.sync(meshMov.pCurV_C);
			}
		};
	  private:
		void printExtParams(std::ofstream &out)
		{
			FlowCircleSTaskDescription::printExtParams(out);
			out << S0 << "\n";
			out << freqRatio << "\n";
			out << MASS << "\n";
			out << COEF_DEMP << "\n";
		};
		void buildVib()
		{
			if(ibCount > 0) LS_STAG_PAR_FOR(FlowCircleSMovTaskDescription, build, meshMov.NMK, meshMov.v_ib, meshMov.v_ib_zy, meshMov.v_ib_ENT, 
				meshMov.v_ib_yx, meshMov.m_baseType, meshMov.baseSideTypes, m_cellOnAirfoil, meshMov.curV_C)
		};
		void buildBaseMesh()
		{
			FlowTaskDescription::buildBaseMesh();
		};
		void rebuildMesh()
		{
			FlowTaskDescription::rebuildMesh();
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
			if (ibCount > 0) LS_STAG_PAR_FOR_2D_RANGE(FlowCircleSMovTaskDescription, calcLevelSetCirclesMov, vs, vf, zs, zf, hs, hf, meshMov.phi, 
					m_mesh->nodes, m_mesh->delta, m_cellOnAirfoil, meshMov.cur_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad)
		}; 
		void initMeshParams()
		{
			Se = S0 * freqRatio;
			COEF_DEMP2 = TWO_PI * Se;
			COEF_ELASTICITY = COEF_DEMP2 * COEF_DEMP2 * MASS;
			COEF_DEMP2 *= 2.0;
			COEF_DEMP2 *= MASS;
			COEF_DEMP2 *= COEF_DEMP;// это b  
			double tmp = MASS + MASS;
			coef_aux = 1.0 / (COEF_DEMP2 * t + tmp);
			double sqrT = t * t;
			double tmp2 = tmp - COEF_ELASTICITY * sqrT;
			coef = (tmp2 + tmp2) * coef_aux;
			coef2 = (COEF_DEMP2 * t - tmp) * coef_aux;
			coef3 = (sqrT + sqrT) * coef_aux;
			for(int i = 0; i < ibCount; i++)
			{
				y_n_1[i] = pCenters[i].y;
				y_n[i] = pCenters[i].y;
				coef4Y[i] = coef3 * COEF_ELASTICITY * pCenters[i].y;
				x_n_1[i] = pCenters[i].x;
				x_n[i] = pCenters[i].x;
				coef4X[i] = coef3 * COEF_ELASTICITY * pCenters[i].x;
			}
			FlowCircleSTaskDescription::initMeshParams();
			meshMov.hsMov = meshMov.hs;// - 1 + (int)floor((Centers[0].x - C_MIN.x) * nn) - 2 * nn - 2;
			meshMov.hfMov = meshMov.hf;//s - 1 + (int)floor((Centers[ibCount-1].x - C_MIN.x) * nn) + 5 * nn + 3;
			meshMov.vsMov = meshMov.vs;// - 1 + (int)floor((Centers[0].y - C_MIN.y) * nn) - 2 * nn - 2;
			meshMov.vfMov = meshMov.vf;//s - 1 + (int)floor((Centers[ibCount-1].y - C_MIN.y) * nn) + 5 * nn + 3;
			meshMov.zsMov = meshMov.zs;// - 1 + (int)floor((Centers[0].y - C_MIN.y) * nn) - 2 * nn - 2;
			meshMov.zfMov = meshMov.zf;//s - 1 + (int)floor((Centers[ibCount-1].y - C_MIN.y) * nn) + 5 * nn + 3;
			for(int i = 0; i < ibCount; i++)
			{
				yFilename[i] = outputDirectoryFull + strConv.toTString(i) + _T("_Y_C") + label + _T(".dat");
				vFilename[i] = outputDirectoryFull + strConv.toTString(i) + _T("_V_C") + label + _T(".dat");
				yFile[i].open(yFilename[i].c_str());
				vFile[i].open(vFilename[i].c_str());
				xFilename[i] = outputDirectoryFull + strConv.toTString(i) + _T("_X_C") + label + _T(".dat");
				uFilename[i] = outputDirectoryFull + strConv.toTString(i) + _T("_U_C") + label + _T(".dat");
				xFile[i].open(xFilename[i].c_str());
				uFile[i].open(uFilename[i].c_str());
			}
			fileExist = true;
		};
		LSSTAGMeshMov meshMov;///< неравномерная LS-STAG-сетка с подвижными погруженными границами
		std::tstring *yFilename, *vFilename, *xFilename, *uFilename;
		std::ofstream *yFile, *vFile, *xFile, *uFile;
	  public:
		double t_0;///< время, с которого начинается движение профилей
		double Se;///< безразмерная частота соственных колебаний системы
		double S0, freqRatio;
		double MASS;///< массa профиля
		double COEF_DEMP, COEF_DEMP2;///< коэффициент демпфирования (b)
		double COEF_ELASTICITY;///< коэффициент упругости (с) , ранее брали COEF_FRICTION из роторов
		double *y_n, *y_n_1, *x_n, *x_n_1;///< массив координат центров профилей с предыдущих шагов по времени
	  private:
		double coef, coef2, coef3, coef_aux;
		double *coef4X, *coef4Y;
		bool fileExist;
	};
}

#endif //_LS_STAG_FLOW_CIRCLES_MOV_TASK_DESCRIPTION_H_
