/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_rotary_osc_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 22:57:01
* \brief     flow_circle_rotary_osc_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_ROTARY_OSC_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_ROTARY_OSC_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circle_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleRotaryOscTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleRotaryOscTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с задачей о вращательных колебаниях кругового профиля

Данный класс реализует работу с описанием решаемой задачи о моделировании вращательных колебаний кругового профиля
*/
	struct FlowCircleRotaryOscTaskDescription : FlowCircleTaskDescription
	{
		FlowCircleRotaryOscTaskDescription():
		FlowCircleTaskDescription(MOTION_OSC, sizeof(FlowCircleRotaryOscTaskDescription)),
		n_f(DEFAULT_TANEDA_N_F),//default Taneda
		S_0(DEFAULT_TANEDA_S_0),//default Taneda
		Taneda(0),
		Omega_0(DEFAULT_ROT_OSC_OMEGA_0),// default low Cd
		Se(DEFAULT_ROT_OSC_SE)// default low Cd
		{
			radius = 0.5 * L_BASE;
		}
		void load(std::ifstream &in)
		{
			FlowTaskDescription::load(in);
			in >> Taneda;
			if(Taneda)
			{
				in >> n_f;
				in >> S_0;
			}
			else
			{
				in >> Se;
				in >> Omega_0;
			}
		};
		void printExtParams(std::ofstream &out)
		{
			FlowTaskDescription::printExtParams(out);
			out << Taneda << "\n";
			if(Taneda)
			{
				out << n_f << "\n";
				out << S_0 << "\n";
			}
			else
			{
				out << Se << "\n";
				out << Omega_0 << "\n";
				Omega = Omega_0;
			}
		};
		void calcNodesCount()
		{
			if(Taneda)
			{
				Se = S_0 * n_f;
				zone_y_3 = 0.75;
				discontin = 0;
				mesh.baseMesh.N = (int)floor(1.1 * nn) + 158;//номер последнего узла по горизонтали
				mesh.baseMesh.M = (int)floor(1.3 * nn) + 58; //номер последнего узла по вертикали
				mesh.baseMesh.K = (int)floor(1.3 * nn) + 58; //номер последнего узла по вертикали
				mesh.pCur_C->set(3.0, 2.5, 2.5);
				mesh.m_cur_C.sync(mesh.pCur_C);
			}
			else FlowTaskDescription::calcNodesCount();
		};
		void initMeshParams()
		{
			mesh.base_length = L_BASE;//22.05.14
			if(Taneda)
			{
				coef = PI_2 * PI * Se;
				coef2 = TWO_PI * Se;
				mesh.zone[0].x = 2; //число узлов по х, зона4
				mesh.zone[0].y = 4; //число узлов по у, зона3
				mesh.zone[0].z = 4; //число узлов по у, зона3
				mesh.zone[1].x = 10; //число узлов по х , зона3
				mesh.zone[1].y = 25; //число узлов по  у, зона2
				mesh.zone[1].z = 25; //число узлов по  у, зона2
				mesh.zone[2].x = (int)floor(0.4 * nn); //число узлов по х , зона2
				mesh.zone[2].y = (int)floor(1.3 * nn); //число узлов по у, зона1
				mesh.zone[2].z = (int)floor(1.3 * nn); //число узлов по у, зона1
				mesh.zone[3].x = 130; //число узлов по х , зона1
				mesh.zone[3].y = 25; //число узлов по  у, зона2
				mesh.zone[3].z = 25; //число узлов по  у, зона2
				mesh.zone[4].x = (int)floor(0.7 * nn); //число узлов по х , зона2
				mesh.zone[4].y = 4; //число узлов по  у, зона3
				mesh.zone[4].z = 4; //число узлов по  у, зона3
				mesh.zone[5].x = 4; //число узлов по х , зона3
				mesh.zone[5].y = 0; //число узлов по  у, зона4
				mesh.zone[5].z = 0; //число узлов по  у, зона4
				mesh.zone[6].x = 12; //число узлов по х , зона4
				mesh.zone[6].y = 0; //число узлов по  у, зона5
				mesh.zone[6].z = 0; //число узлов по  у, зона5

				mesh.deltaCurZone[0].x = 0.005; //шаг сетки по х , зона1
				mesh.deltaCurZone[0].y = aa / nn; //шаг сетки по у, зона1
				mesh.deltaCurZone[0].z = aa / nn; //шаг сетки по у, зона1	
				mesh.deltaCurZone[1].x = mesh.deltaCurZone[0].y; //шаг сетки по х , зона2
				mesh.deltaCurZone[1].y = 0.85 / mesh.zone[1].y; //шаг сетки по у, зона2
				mesh.deltaCurZone[1].z = 0.85 / mesh.zone[1].z; //шаг сетки по у, зона2
				mesh.deltaCurZone[2].x = 0.125; //шаг сетки по х , зона3
				mesh.deltaCurZone[2].y = 0.25; //шаг сетки по у, зона3
				mesh.deltaCurZone[2].z = 0.25; //шаг сетки по у, зона3
				mesh.deltaCurZone[3].x = 0.5; //шаг сетки по х , зона4
				mesh.deltaCurZone[3].y = 0.0; //шаг сетки по у, зона4
				mesh.deltaCurZone[3].z = 0.0; //шаг сетки по у, зона4

				mesh.hs = 12;
				mesh.hf = mesh.hs + ((int)floor(1.1 * nn) + 130);
				mesh.vs = 29;
				mesh.vf = mesh.vs + (int)floor(1.3 * nn);
				FlowTaskDescription::initInftyCylMeshZ(&mesh);
				m_phi = m_mesh->phi;
			}
			else
			{
				coef2 = TWO_PI * Se * V_BASE / L_BASE;
				FlowTaskDescription::initMeshParams();
			}
		};
		void buildBaseMesh()
		{
			FlowTaskDescription::buildBaseMesh();
			rebuildMesh();//buildVib();
		};
		void buildVib()//29.12.13
		{// радиус в каждом слое свой, если это не цилиндр, а шар
			LS_STAG_PAR_FOR(LSSTAGUniformMesh, projectionVib, mesh.NMK, radius * Omega, mesh.m_baseType, mesh.baseSideTypes, mesh.delta, mesh.alpha,
				mesh.v_ib, mesh.v_ib_yx, mesh.v_ib_zy, mesh.v_ib_ENT, mesh.baseMesh.N);
		};
		/*!
		Пересчитывает скорости на погруженной границе (сетка не перестраивается, т.к. погруженная граница не 
		перемещается в пространстве)
		*/
		void rebuildMesh()// 29.12.13
		{
			if(Taneda)
			{
				if(((*m_curTime) < 17.5) || ((*m_curTime) > 33.0))
				{
					double aux = coef2 * (*m_curTime);
					Omega = coef * sin(aux);
				}
				else Omega = 0.0;
			}
			else if(Se > EQ_DOUBLE)
			{
				double aux = coef2 * (*m_curTime);
				Omega = Omega_0 * sin(aux);
			}
			buildVib();
		};
	  private:
		double coef, coef2, radius;
	  public:
		double Se;///< безразмерная частота вращательных колебаний
		double Omega_0;///< амплитуда угловой скорости
		double S_0;///< число Струхаля при обтекании неподвижного профиля, характеризуемом тем же значением числа Re
		double n_f;///< относительная частота вынужденных колебаний 
		// n_f = Se / S_0
		double Omega;///< угловая скорость
		int Taneda;///< 1 --- высокочастотные колебания --- эффект Танеды
		// true => zone_y_3 = 0.75, discontin = 0
	};
}

#endif //_LS_STAG_FLOW_CIRCLE_ROTARY_OSC_TASK_DESCRIPTION_H_
