/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_ve_osc_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 23:12:51
* \brief     flow_circle_ve_osc_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_VE_OSC_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_VE_OSC_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circle_mov_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleVEOscTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleVEOscTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с задачей о вязкоупругих колебаниях кругового профиля

Данный класс реализует работу с описанием решаемой задачи о моделировании вязкоупругих колебаний кругового профиля
*/
	struct FlowCircleVEOscTaskDescription : FlowCircleMovTaskDescription
	{
		FlowCircleVEOscTaskDescription():
		FlowCircleMovTaskDescription(OSC_VISCO_ELASTIC, sizeof(FlowCircleVEOscTaskDescription)),
		MASS(DEFAULT_VE_OSC_AIRFOIL_MASS),// default
		COEF_DEMP(DEFAULT_VE_OSC_COEF_DEMP),//default
		Se(0.24),// default --- резонансная частота
		fileExist(false)
		{
			Re = 1000.0;
			t_0 = 1.0;
			V_BASE = 3.0;//
			fcoef = 2.0 / (rho * V_BASE * V_BASE);
			U_N = V_BASE;// x-я скорость набегающего потока на N стенке 
			U_S = V_BASE;//  x-я скорость набегающего потока на S стенке
			U_B = 0.0;// V_BASE;// x-я скорость набегающего потока на T стенке 
			U_T = 0.0;// V_BASE;//  x-я скорость набегающего потока на B стенке
			U_E = V_BASE;//  x-я скорость набегающего потока на E стенке
			U_W = V_BASE;//  x-я скорость набегающего потока на W стенке
			maxDevTM = L_BASE * 0.1;
		}
		~FlowCircleVEOscTaskDescription()
		{
			if(fileExist)
			{
				yFile.close();
				vFile.close();
			}
		};
		void printLoads()
		{
			BaseTaskDescription::printLoads();
			yFile << meshMov.pCur_C->y << "\n";
			vFile << meshMov.pCurV_C->y << "\n";
		};
		void saveState(std::ofstream &out)
		{
			BaseTaskDescription::saveState(out);
			out.write((char*)(&y_n), sizeof(double));
			out.write((char*)(&y_n_1), sizeof(double));
		};
		void loadState(std::ifstream &in)
		{
			BaseTaskDescription::loadState(in);
			in.read((char*)(&y_n), sizeof(double));
			in.read((char*)(&y_n_1), sizeof(double));
		};
		/*!
		Вычисляет текущие координаты центра профиля и его скорость (решает уравнение движения)
		*/
		void calcPosition()
		{
			if((*m_curTime) > t_0)
			{
				meshMov.pCur_C->y = coef4 + coef3 * m_loads->F_aer_y + coef2 * y_n_1 + coef * y_n;
				double curDev = meshMov.pCur_C->y - y_n;
				meshMov.pCurV_C->y  = curDev / t;
				if(useTM && (fabs(curDev) > maxDevTM))
				{
					meshMov.pCurV_C->y = 0.0;
					meshMov.pCur_C->y = y_n;
				}
				y_n_1 = y_n;
				y_n = meshMov.pCur_C->y;
				meshMov.m_cur_C.sync(meshMov.pCur_C);
				meshMov.m_curV_C.sync(meshMov.pCurV_C);
			}
		};
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> Se;
			in >> MASS;
			in >> COEF_DEMP;
		};
	  private:
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << Se << "\n";
			out << MASS << "\n";
			out << COEF_DEMP << "\n";
		};
		void calcNodesCount()
		{
			FlowTaskDescription::calcNodesCount();
			meshMov.baseMesh.N += (int)floor(9.5 * nn);//номер последнего узла по горизонтали
			meshMov.baseMesh.M += (int)floor(9.0 * nn); //номер последнего узла по вертикали 
		};
		void initMeshParams()
		{
			Omega = Se * TWO_PI * V_BASE / L_BASE;
			COEF_ELASTICITY = Omega * Omega * MASS;
			double tmp = MASS + MASS;
			coef_aux = 1.0 / (COEF_DEMP * t + tmp);
			double sqrT = t * t;
			double tmp2 = tmp - COEF_ELASTICITY * sqrT;
			coef = (tmp2 + tmp2) * coef_aux;
			coef2 = (COEF_DEMP * t - tmp) * coef_aux;
			coef3 = (sqrT + sqrT) * coef_aux;
			coef4 = coef3 * COEF_ELASTICITY * P_C.y;
			y_n_1 = P_C.y;
			y_n = P_C.y;
			FlowCircleMovTaskDescription::initMeshParams();
			initFinestMesh();
			yFilename = outputDirectoryFull + _T("Y_C") + label + _T(".dat");
			vFilename = outputDirectoryFull + _T("V_C") + label + _T(".dat");
			yFile.open(yFilename.c_str());
			vFile.open(vFilename.c_str());
			fileExist = true;
		};
		std::tstring yFilename, vFilename;
		std::ofstream yFile, vFile;
	  public:
		double Omega;///< соственная частота системы
		double Se;///< безразмерная частота соственных колебаний системы
		double MASS;///< массa профиля
		double COEF_DEMP;///< коэффициент демпфирования (b)
		double COEF_ELASTICITY;///< коэффициент упругости (с) , ранее брали COEF_FRICTION из роторов
		double y_n, y_n_1;///< ординаты центра профиля с предыдущих шагов по времени
	  private:
		double coef, coef2, coef3, coef4, coef_aux, maxDevTM;
		bool fileExist;
	};
}

#endif //_LS_STAG_FLOW_CIRCLE_VE_OSC_TASK_DESCRIPTION_H_
