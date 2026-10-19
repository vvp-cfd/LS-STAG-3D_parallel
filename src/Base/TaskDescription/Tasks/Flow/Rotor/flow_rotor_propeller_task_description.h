/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_rotor_propeller_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 19:34:46
* \brief     flow_rotor_propeller_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_ROTOR_PROPELLER_DESCRIPTION_H_
#define _LS_STAG_FLOW_ROTOR_PROPELLER_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_rotor_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowRotorPropellerTaskDescription

Данный файл содержит в себе реализацию класса FlowRotorPropellerTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetRotorPropeller, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, PropellerRadius,
		const double, PropellerWidth, const int, numBlade, const double, L_BASE, const double, step2, const double, step,
		const double, ANGLE, const double, level, const FlowRotorTaskDescription::RotorKinds, id_rotor);
/*!
\brief Класс, реализующий работу с задачей об обтекании пропеллера

Данный класс реализует работу с описанием решаемой задачи об обтекании пропеллера. Пропеллер может авторотировать или 
совершать автоколебания
*/
	struct FlowRotorPropellerTaskDescription : FlowRotorTaskDescription
	{
	  public:
		FlowRotorPropellerTaskDescription():
		FlowRotorTaskDescription(ROTOR_PROPELLER, sizeof(FlowRotorPropellerTaskDescription)),
		h(DEFAULT_BLADE_LENGTH),
		lwRatio(0.1)
		{
			withoutFinestMesh = true;
		}
		LS_STAG_KERNEL_2D(calcLevelSetRotorPropeller, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double PropellerRadius, const double PropellerWidth, const int numBlade,
			const double L_BASE, const double step2, const double step, const double ANGLE, const double level,
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
				const double r = pToC.norm();
				if (r > L_BASE)
					_phi[curNum] = L_BASE - r;
				else
				{
					if (r < PropellerRadius)
						_phi[curNum] = PropellerRadius - r;
					else
					{
						double phi = 0.0;
						double phiBlade = FlowRotorTaskDescription::CurBlade(pToC, phi, numBlade, step2, step, ANGLE, level, id_rotor);
						FlowRotorTaskDescription::checkPhiBlade(phiBlade);
						const double edge = PropellerRadius - r;//aux
						const double edge2 = L_BASE - r;//aux3
						double phi_0 = 0.0;
						if ((fabs(phiBlade) < EQ_DOUBLE) || (fabs(phiBlade - PI) < EQ_DOUBLE) || (fabs(phiBlade - TWO_PI) < EQ_DOUBLE))
							phi_0 = PropellerWidth - fabs(pToC.y);
						else
						{
							if ((fabs(phiBlade - PI_2) < EQ_DOUBLE) || (fabs(phiBlade - PI_32) < EQ_DOUBLE))
								phi_0 = PropellerWidth - fabs(pToC.x);
							else
							{
								const double tgP = tan(phiBlade);
								phi_0 = PropellerWidth - fabs(pToC.y - tgP * pToC.x) / sqrt(1.0 + tgP * tgP);
							}
						}
						if (fabs(edge) < EQ_DOUBLE)
						{
							_phi[curNum] = (phi_0 > EQ_DOUBLE) ? phi_0 : 0.0;
						}
						else
						{
							if (fabs(edge2) < EQ_DOUBLE)
								_phi[curNum] = (phi_0 > MINUS_EQ_DOUBLE) ? 0.0 : phi_0;
							else
							{
								if (fabs(phi_0) < EQ_DOUBLE)
									_phi[curNum] = 0.0;
								else
								{
									double res = 0.0;
									if (phi_0 > EQ_DOUBLE)
									{
										res = -edge;
										if (edge2 < res) res = edge2;
										if (phi_0 < res) res = phi_0;
									}
									else
									{
										res = -edge2;
										if (edge > res) res = edge;
										if (phi_0 > res) res = phi_0;
									}
									_phi[curNum] = res;
								}
							}
						}
					}
				}
			}
		};
		void load(std::ifstream &in)
		{
			FlowRotorTaskDescription::load(in);
			in >> h;
			in >> lwRatio;
		};
	  private:
		void printExtParams(std::ofstream &out)
		{
			FlowRotorTaskDescription::printExtParams(out);
			out << h << "\n";
			out << lwRatio << "\n";
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
			static const double PropellerRadius = L_BASE - h;
			static const double PropellerWidth = 0.5 * lwRatio * h;
			LS_STAG_PAR_FOR_2D_RANGE(FlowRotorPropellerTaskDescription, calcLevelSetRotorPropeller, vs, vf, zs, zf, hs, hf, meshMov.phi,
				m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N* m_mesh->baseMesh.M, PropellerRadius, PropellerWidth, numBlade, L_BASE,
				step2, step, ANGLE, level, id_rotor)
		};
	  public:
		double h;///< длина лопасти
		double lwRatio;
	};
}

#endif //_LS_STAG_FLOW_ROTOR_PROPELLER_DESCRIPTION_H_
