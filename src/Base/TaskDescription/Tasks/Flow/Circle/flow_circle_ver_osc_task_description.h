/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_ver_osc_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 23:11:41
* \brief     flow_circle_ver_osc_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_VER_OSC_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_VER_OSC_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circle_mov_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleVerOscTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleVerOscTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с задачей о продольных колебаниях кругового профиля

Данный класс реализует работу с описанием решаемой задачи о моделировании продольных колебаний кругового профиля
*/
	struct FlowCircleVerOscTaskDescription : FlowCircleMovTaskDescription
	{
		FlowCircleVerOscTaskDescription():
		FlowCircleMovTaskDescription(OSC_VERTICAL, sizeof(FlowCircleVerOscTaskDescription)),
		A_MAX(DEFAULT_VER_OSC_A_MAX),
		S_0(DEFAULT_VER_OSC_S_0),
		n_f(DEFAULT_VER_OSC_N_F)
		{
			Re = 185.0;
			t_0 = DEFAULT_VER_OSC_MOV_START;//  0.01;// test 5.06.14 ;
		}
		void initMeshParams()
		{
			Se = n_f * S_0;
			freq = TWO_PI * Se;
			vel_ampl = -freq * A_MAX;
			meshMov.pCur_C->y += A_MAX;
			meshMov.m_cur_C.sync(meshMov.pCur_C);
			FlowCircleMovTaskDescription::initMeshParams();
		};
		/*!
		Вычисляет текущие координаты центра профиля и его скорость
		*/
		void calcPosition()
		{
			if((*m_curTime) > t_0)
			{
				double aux = freq * ((*m_curTime) - t_0);
				meshMov.pCur_C->y = P_C.y + A_MAX * cos(aux);// поперечные колебания  + A_MAX - 200.0 * ((*m_curTime) - t_0); // test 17.11.14
				meshMov.pCurV_C->y = vel_ampl * sin(aux);//-200.0;// test 17.11.14 
				meshMov.m_cur_C.sync(meshMov.pCur_C);
				meshMov.m_curV_C.sync(meshMov.pCurV_C);
			}
		};
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in); 
			in >> n_f;
			in >> S_0;
			in >> A_MAX;
		};
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << n_f << "\n";
			out << S_0 << "\n";
			out << A_MAX << "\n";
		};
		double Se;///< кинематическое число Струхаля
		double S_0, n_f;// n_f = Se / S_0
		double A_MAX;///< максимальная амплитуда колебаний профиля
	  private:
		double freq, vel_ampl;
	};
}

#endif //_LS_STAG_FLOW_CIRCLE_VER_OSC_TASK_DESCRIPTION_H_
