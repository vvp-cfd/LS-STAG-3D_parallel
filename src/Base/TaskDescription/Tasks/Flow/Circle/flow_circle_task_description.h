/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 21:47:27
* \brief     flow_circle_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "../flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_R_KERNEL_2D(calcLevelSetCircle, const int, hs, const int, hf, LS_STAG_R_PD, _phi,
		LS_STAG_CR_PC, _nodes, LS_STAG_R_PC, _delta, const coordinate, P_C, const int, _N, const int, _NM, const double, rad);
/*!
\brief Класс, реализующий работу с задачей об обтекании кругового профиля

Данный класс реализует работу с описанием решаемой задачи об обтекании кругового профиля
*/
	struct FlowCircleTaskDescription : FlowTaskDescription
	{
	  public:
		explicit FlowCircleTaskDescription(MotionKinds _id_motion = MOTION_NO, int _size = sizeof(FlowCircleTaskDescription)):
		  FlowTaskDescription(AIRFOIL_CIRCLE, _id_motion, (LSSTAGUniformMesh*)&mesh, _size),//default mesh
		  rad(0.5)
		  {
			  discontin = (_id_motion != MOTION_MOV) ? 1 : 0;
		  }
		LS_STAG_KERNEL_2D(calcLevelSetCircle, const int hs, const int hf, LS_STAG_R_PD _phi, LS_STAG_CR_PC _nodes, LS_STAG_R_PC _delta,
			const coordinate P_C, const int _N, const int _NM, const double rad)
		{
			int startNum = _NM * j + _N * i;
			int curNum = startNum + hs;
			coordinate curPoint = _nodes[curNum];
			coordinate* curD = _delta + curNum;
			curPoint.add(0.5 * curD->x, 0.5 * curD->y, 0.0);
			for (int end = startNum + hf; curNum < end; curNum++)
			{
				coordinate pToC;
				pToC.setDiffXY(curPoint, P_C);
				_phi[curNum] = rad - pToC.norm();
				curD++;
				curPoint.x += curD->x;
			};
		}
		virtual double getDistance(const double x, const double y, const double z)
		{
			double aux = x - P_C.x;
			double aux2 = y - P_C.y;// для цилиндра не используем double aux3 = z - P_C.z;
			double dist = -rad + sqrt(aux * aux + aux2 * aux2);// +aux3 * aux3);
			return dist;
		};
	  protected:
		  /*!
		  Вычисляет значения функции уровня в углах ячеек основной сетки
		\param zs номер ряда по z, с которого рассматриваем ячейки
		\param zf номер ряда по z, до которого рассматриваем ячейки
		\param vs номер ряда по вертикали, с которого рассматриваем ячейки
		\param vf номер ряда по вертикали, до которого рассматриваем ячейки
		\param hs номер ряда по горизонтали, с которого рассматриваем ячейки
		\param hf номер ряда по горизонтали, до которого рассматриваем ячейки
		*/
		virtual void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf) override
		  {
			  LS_STAG_PROFILE("buildLS");
			  LS_STAG_PAR_FOR_2D_RANGE(FlowCircleTaskDescription, calcLevelSetCircle, vs, vf, zs, zf, hs, hf, mesh.phi,
				  m_mesh->nodes, m_mesh->delta, P_C, m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad)
		  };
		LSSTAGMesh mesh;///< неравномерная LS-STAG-сетка
		double rad;
	}; 
}

#endif //_LS_STAG_FLOW_CIRCLE_TASK_DESCRIPTION_H_
