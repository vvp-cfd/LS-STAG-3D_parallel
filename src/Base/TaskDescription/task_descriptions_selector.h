/*
* \copyright (C) 2014 Puzikova V.V.
* \file      task_descriptions_selector.h
* \author    Puzikova Valeria (vvp) 
* \date      04.01.2014 12:43:34
* \brief     task_descriptions_selector.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TASK_DESCRIPTIONS_SELECTOR_H_
#define _LS_STAG_TASK_DESCRIPTIONS_SELECTOR_H_

// ================================================================== INCLUDES
#include <cstdlib>
#include <iostream>
#include "Tasks/Flow/Circle/flow_circle_hor_osc_task_description.h"
#include "Tasks/Flow/Circle/flow_circle_rotary_osc_task_description.h"
#include "Tasks/Flow/Circle/flow_circle_ve_osc_task_description.h"
#include "Tasks/Flow/Circle/flow_circle_ver_osc_task_description.h"
#include "Tasks/Flow/Circles/flow_circles_task_description.h"
#include "Tasks/Flow/Circles/flow_circles_mov_task_description.h"
#include "Tasks/Flow/flow_elliptic_task_description.h"
#include "Tasks/Flow/flow_points_task_description.h"
#include "Tasks/Flow/Rotor/flow_rotor_darrieus_task_description.h"
#include "Tasks/Flow/Rotor/flow_rotor_plate_task_description.h"
#include "Tasks/Flow/Rotor/flow_rotor_propeller_task_description.h"
#include "Tasks/Flow/Rotor/flow_rotor_savonius_task_description.h"
#include "Tasks/Flow/flow_square_task_description.h"
#include "Tasks/Flow/flow_cube_task_description.h"
#include "Tasks/Flow/flow_finite_cylinder_task_description.h"
#include "Tasks/Flow/flow_sphere_task_description.h"
#include "Tasks/Flow/Circle/flow_circle_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TaskDescriptionsSelector

Данный файл содержит в себе реализацию класса TaskDescriptionsSelector
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий выбор задачи

Данный класс реализует хранение и загрузку/выгрузку параметров задачи
*/
	struct TaskDescriptionsSelector
	{
		typedef BaseTaskDescription::TaskKinds task;
		typedef BaseTaskDescription::MotionKinds motion;
		typedef BaseTaskDescription::OscKinds osc;
		typedef FlowTaskDescription::AirfoilKind airfoil;
		typedef FlowRotorTaskDescription::RotorKinds rotor;
		TaskDescriptionsSelector():
			m_task(NULL),
		id(task::TASK_UNKNOWN),
		id_motion(motion::MOTION_UNKNOWN),
		id_airfoil(airfoil::AIRFOIL_UNKNOWN),
		id_rotor(rotor::ROTOR_UNKNOWN)
		{
			;
		}
		~TaskDescriptionsSelector()
		{
			if (m_task) delete m_task;
		}
		void load(const std::string &name)
		{
			std::ifstream f;
			// считываем параметры, однозначно определяющие структуру TaskDescription
			f.open(name.c_str());
			f >> aux;
			id = (task)aux;
			f >> aux;
			id_motion = (motion)aux;
			f >> aux;
			id_osc = (osc)aux;
			f >> aux;
			id_airfoil = (airfoil)aux;
			f >> aux;
			id_rotor = (rotor)aux;
			selector();// m_task = указатель на нужную структуру
			if (!m_task)
			{
				std::cerr << "Unsupported task description" << std::endl;
				std::exit(1);
			}
			m_task->load(f);// загружаем параметры расчёта
			f.close();
		};
		void print()
		{
			std::ofstream out;
			m_task->print(out);
			out.close();
		}
	  private:
		void selector()
		{
			switch(id)
			{
			  case task::TASK_FLOW:
				switch(id_airfoil)
				{
				  case airfoil::AIRFOIL_SQUARE:
					m_task = (BaseTaskDescription*)(new FlowSquareTaskDescription);// &square;
					break;
				  case airfoil::AIRFOIL_CUBE:
					m_task = (BaseTaskDescription*)(new FlowCubeTaskDescription);// &cube;
					break;
				  case airfoil::AIRFOIL_FINITE_CYLINDER:
					m_task = (BaseTaskDescription*)(new FlowFiniteCylinderTaskDescription);// &finiteCylinder;
					break;
				  case airfoil::AIRFOIL_SPHERE:
					m_task = (BaseTaskDescription*)(new FlowSphereTaskDescription);// &sphere;
					break;
				  case airfoil::AIRFOIL_ELLIPTIC:
					m_task = (BaseTaskDescription*)(new FlowEllipticTaskDescription);// &elliptic;
					break;
				  case airfoil::AIRFOIL_POINTS:
					m_task = (BaseTaskDescription*)(new FlowPointsTaskDescription);// &points;
					break;
				  case airfoil::AIRFOIL_CIRCLES:
					if(id_motion == motion::MOTION_MOV)
						m_task = (BaseTaskDescription*)(new FlowCircleSMovTaskDescription);// &circleSMov;
					else
						m_task = (BaseTaskDescription*)(new FlowCircleSTaskDescription);// &circleS;
					break;
				  case airfoil::AIRFOIL_ROTOR:
					switch(id_rotor)
					{
					  case rotor::ROTOR_PLATE:
						m_task = (BaseTaskDescription*)(new FlowRotorPlateTaskDescription);// &plate;
						break;
					  case rotor::ROTOR_PROPELLER:
						m_task = (BaseTaskDescription*)(new FlowRotorPropellerTaskDescription);// &propeller;
						break;
					  case rotor::ROTOR_SAVONIUS:
						m_task = (BaseTaskDescription*)(new FlowRotorSavoniusTaskDescription);// &savonius;
						break;
					  case rotor::ROTOR_DARRIEUS:
						m_task = (BaseTaskDescription*)(new FlowRotorDarrieusTaskDescription);// &darrieus;
						break;
					}
					break;
				  case airfoil::AIRFOIL_CIRCLE:
					switch(id_motion)
					{
					  case motion::MOTION_NO:
						m_task = (BaseTaskDescription*)(new FlowCircleTaskDescription);// &circle;
						break;
					  case motion::MOTION_OSC:
						m_task = (BaseTaskDescription*)(new FlowCircleRotaryOscTaskDescription);// &circleRotOsc;
						break;
					  case motion::MOTION_MOV:
						switch(id_osc)
						{
						  case osc::OSC_HORIZONTAL:
							m_task = (BaseTaskDescription*)(new FlowCircleHorOscTaskDescription);// &circleHorOsc;
							break;
						  case osc::OSC_VERTICAL:
							m_task = (BaseTaskDescription*)(new FlowCircleVerOscTaskDescription);// &circleVerOsc;
							break;
						  case osc::OSC_VISCO_ELASTIC:
							m_task = (BaseTaskDescription*)(new FlowCircleVEOscTaskDescription);// &circleVEOsc;
							break;
						}
						break;
					}
					break;
				}
				break;
			}
		};
	  public:
		BaseTaskDescription* m_task;
	  private:
		task id;
		motion id_motion; 
		osc id_osc;
		airfoil id_airfoil;
		rotor id_rotor;
		int aux;
	};
}

#endif //_LS_STAG_TASK_DESCRIPTIONS_SELECTOR_H_
