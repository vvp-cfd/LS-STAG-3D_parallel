/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_rotor_plate_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 19:30:08
* \brief     flow_rotor_plate_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ROTOR_PLATE_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_ROTOR_PLATE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_rotor_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowRotorPlateTaskDescription

Данный файл содержит в себе реализацию класса FlowRotorPlateTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetRotorPlate, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, length,
		const double, width, const double, ANGLE);
/*!
\brief Класс, реализующий работу с задачей об обтекании пластины

Данный класс реализует работу с описанием решаемой задачи об обтекании пластины. Пластина может авторотировать или 
совершать автоколебания
*/
	struct FlowRotorPlateTaskDescription : FlowRotorTaskDescription
	{
	  public:
		FlowRotorPlateTaskDescription():
		FlowRotorTaskDescription(ROTOR_PLATE, sizeof(FlowRotorPlateTaskDescription)),
		length(L_BASE),
		lwRatio(0.1)
		{
			calcWidth();
			numBlade = 1;
		}
		LS_STAG_KERNEL_2D(calcLevelSetRotorPlate, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double length, const double width, const double ANGLE)
		{
			int curNum = _NM * j + _N * i + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = curNum + hf; curNum < end; curNum++, curD++, curPoint.x += curD->x)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				double phi_L, phi_W;
				if ((fabs(ANGLE) < EQ_DOUBLE) || (fabs(ANGLE - PI) < EQ_DOUBLE))
				{
					phi_L = length - fabs(pToC.x);
					phi_W = width - fabs(pToC.y);
				}
				else
				{
					if ((fabs(ANGLE - PI_2) < EQ_DOUBLE) || (fabs(ANGLE - PI_32) < EQ_DOUBLE))
					{
						phi_W = width - fabs(pToC.x);
						phi_L = length - fabs(pToC.y);
					}
					else
					{
						double tgA = tan(ANGLE);
						double ctgA = 1.0 / tgA;
						phi_W = width - fabs(pToC.y - tgA * pToC.x) / sqrt(1.0 + tgA * tgA);
						phi_L = length - fabs(pToC.y + ctgA * pToC.x) / sqrt(1.0 + ctgA * ctgA);
					}
				}
				_phi[curNum] = (((fabs(phi_W) < EQ_DOUBLE) && (phi_L > MINUS_EQ_DOUBLE)) ||
					((fabs(phi_L) < EQ_DOUBLE) && (phi_W > MINUS_EQ_DOUBLE))) ? 0.0 :
					((((fabs(phi_W) < EQ_DOUBLE)) && (phi_L < MINUS_EQ_DOUBLE)) ? phi_L :
						(((fabs(phi_L) < EQ_DOUBLE) && (phi_W < MINUS_EQ_DOUBLE)) ? phi_W :
							(((phi_L * phi_W) > EQ_DOUBLE) ? ((fabs(phi_W) < fabs(phi_L)) ? phi_W : phi_L)
								: ((phi_W < phi_L) ? phi_W : phi_L))));
			}
		}
		void load(std::ifstream &in)
		{
			FlowRotorTaskDescription::load(in);
			in >> lwRatio;
			calcWidth();
		};
		void printExtParams(std::ofstream &out)
		{
			FlowRotorTaskDescription::printExtParams(out);
			out << lwRatio << "\n";
		};
	  private:
		void loadNumBlade(std::ifstream &in){  };
		void printNumBlade(std::ofstream &out){  };
		void calcWidth(){ width = 0.5 * L_BASE * lwRatio; };
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
			LS_STAG_PAR_FOR_2D_RANGE(FlowRotorPlateTaskDescription, calcLevelSetRotorPlate, vs, vf, zs, zf, hs, hf, meshMov.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, length, width, ANGLE)
		}; 
		double length, width;///< длина и ширина пластины
		double lwRatio;
	};
}

#endif //_LS_STAG_FLOW_ROTOR_PLATE_TASK_DESCRIPTION_H_
