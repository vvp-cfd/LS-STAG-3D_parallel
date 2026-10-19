/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_rotor_savonius_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 20:36:13
* \brief     flow_rotor_savonius_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ROTOR_SAVONIUS_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_ROTOR_SAVONIUS_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_rotor_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowRotorSavoniusTaskDescription

Данный файл содержит в себе реализацию класса FlowRotorSavoniusTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetRotorSavonius, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, AxisRadius,
		const int, numBlade, const double, BladeRadius, const double, BladeWidth, const double, step2, const double, step,
		const double, ANGLE, const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor);
/*!
\brief Класс, реализующий работу с задачей об обтекании ротора Савониуса

Данный класс реализует работу с описанием решаемой задачи об обтекании ротора Савониуса. Ротор Савониуса может 
авторотировать или совершать автоколебания
*/
	struct FlowRotorSavoniusTaskDescription : FlowRotorTaskDescription
	{
	  public:
		FlowRotorSavoniusTaskDescription():
		FlowRotorTaskDescription(ROTOR_SAVONIUS, sizeof(FlowRotorSavoniusTaskDescription)),
		rotor_time_aux(DEFAULT_SAVONIUS_ROTOR_TIME_AUX),//default
		RotorRadius(L_BASE),
		rwRatio(0.15), 
		raRatio(0.1)
		{
			BladeRadius = 0.5 * RotorRadius;
			calcWidth();
		};
		LS_STAG_KERNEL_2D(calcLevelSetRotorSavonius, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double AxisRadius, const int numBlade, const double BladeRadius,
			const double BladeWidth, const double step2, const double step, const double ANGLE, const double level,
			const RotorKinds id_rotor)
		{
			int curNum = _NM * j + _N * i + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = curNum + hf; curNum < end; curNum++, curD++, curPoint.x += curD->x)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				double RC = AxisRadius - pToC.norm();//
				if (RC > EQ_DOUBLE)//
					_phi[curNum] = RC;//
				else//
				{
					double phi = 0.0;
					double phiBlade = FlowRotorTaskDescription::CurBlade(pToC, phi, numBlade, step2, step, ANGLE, level, id_rotor);
					double phiBladeNext = phiBlade + step2; // 11.02.13
					if (phiBlade > NEAR_TWO_PI) phiBladeNext -= TWO_PI; // 17.02
					FlowRotorTaskDescription::checkPhiBlade(phiBlade);
					double phiBladeLast = ANGLE + step * (2.0 * numBlade - level);
					double res = 0.0;
					if (((phi < phiBlade) && (phi > (phiBladeNext - TWO_PI))) || ((phi > phiBladeLast) && (ANGLE < PI_2)))
					{
						if ((fabs(phiBlade) < EQ_DOUBLE) || (fabs(phiBlade - PI) < EQ_DOUBLE) || (fabs(phiBlade - TWO_PI) < EQ_DOUBLE))
							res = -fabs(pToC.y);
						else
						{
							if ((fabs(phiBlade - PI_2) < EQ_DOUBLE) || (fabs(phiBlade - PI_32) < EQ_DOUBLE))
								res = -fabs(pToC.x);
							else
							{
								double phi_0 = tan(phiBlade);
								res = -fabs(pToC.y - phi_0 * pToC.x) / sqrt(1.0 + phi_0 * phi_0);
							}
						}
					}
					else
					{
						pToC.diff(cos(phiBlade), sin(phiBlade), 0.0, BladeRadius);
						double aux = BladeRadius - pToC.norm();
						double aux2 = aux - BladeWidth;
						if (aux2 > 0.0)
							res = -aux2;
						else
						{
							if (aux < 0.0)
								res = aux;
							else
							{
								res = fabs(aux);
								if (fabs(aux2) < res) res = fabs(aux2);
							}
						}
					}
					if (res < RC) res = RC;//
					_phi[curNum] = res;
				}
			};
		}
		void load(std::ifstream &in)
		{
			FlowRotorTaskDescription::load(in);
			in >> rotor_time_aux;
			in >> rwRatio;
			in >> raRatio;
			calcWidth();
		};
	  private:  
		void printExtParams(std::ofstream &out)
		{
			FlowRotorTaskDescription::printExtParams(out);
			out << rotor_time_aux << "\n";
			out << rwRatio << "\n";
			out << raRatio << "\n";
		};
		void calcWidth()
		{
			AxisRadius = raRatio * RotorRadius;
			BladeWidth = rwRatio * BladeRadius;
		}
		void initRotorParams()
		{
			FlowRotorTaskDescription::initRotorParams();
			level = 0.2;
			if(numBlade == 1) level = 1.0;
			if(numBlade > 4) level = 0.0;// всё равно при numBlade = 5 и более секторы лопастей перекрываются
		};
		void calcOmega()
		{
			if(((*m_curTime) > rotor_time_aux) && ((*m_curTime) < rotor_time2) && (Omega < TipSpeedRatio)) Omega = TipSpeedRatio;
			FlowRotorTaskDescription::calcOmega();
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
			LS_STAG_PAR_FOR_2D_RANGE(FlowRotorSavoniusTaskDescription, calcLevelSetRotorSavonius, vs, vf, zs, zf, hs, hf, meshMov.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N* m_mesh->baseMesh.M, AxisRadius, numBlade, BladeRadius, BladeWidth, step2,
				step, ANGLE, level, id_rotor)
		};
		double RotorRadius;///< радиус ротора
		double BladeRadius;///< радиус лопасти ротора
		double BladeWidth;///< толщина лопасти
		double AxisRadius;///< радиус оси ротора
		double rwRatio, raRatio;
	  public:
		double rotor_time_aux;// момент времени с которого ротор Савониуса ускоряется (до 0.4)
	};
}

#endif //_LS_STAG_FLOW_ROTOR_SAVONIUS_TASK_DESCRIPTION_H_
