/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_hor_osc_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 23:08:55
* \brief     flow_circle_hor_osc_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_HOR_OSC_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_HOR_OSC_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circle_mov_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleHorOscTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleHorOscTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с задачей о поперечных колебаниях кругового профиля

Данный класс реализует работу с описанием решаемой задачи о моделировании поперечных колебаний кругового профиля
*/
	struct FlowCircleHorOscTaskDescription : FlowCircleMovTaskDescription
	{
		FlowCircleHorOscTaskDescription():
		FlowCircleMovTaskDescription(OSC_HORIZONTAL, sizeof(FlowCircleHorOscTaskDescription)),
		A_MAX(DEFAULT_HOR_OSC_A_MAX),
		Se(DEFAULT_HOR_OSC_SE)
		{
			Re = 100.0;
			t_0 = DEFAULT_HOR_OSC_MOV_START;
		}
		void initMeshParams()
		{
			freq = TWO_PI * Se;
			vel_ampl = -freq * A_MAX;
			meshMov.pCur_C->x += A_MAX;
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
				meshMov.pCur_C->x = P_C.x + A_MAX * cos(aux);// поперечные колебания
				meshMov.pCurV_C->x = vel_ampl * sin(aux);
				meshMov.m_cur_C.sync(meshMov.pCur_C);
				meshMov.m_curV_C.sync(meshMov.pCurV_C);
			}
		};
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in); 
			in >> Se;
			in >> A_MAX;
		};
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << Se << "\n";
			out << A_MAX << "\n";
		};
		double Se;///< кинематическое число Струхаля
		double A_MAX;///< максимальная амплитуда колебаний профиля
	  private:
		double freq, vel_ampl;
	};
}

#endif //_LS_STAG_FLOW_CIRCLE_HOR_OSC_TASK_DESCRIPTION_H_
