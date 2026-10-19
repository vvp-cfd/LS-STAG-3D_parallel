/*
* \copyright (C) 2013 Puzikova V.V.
* \file      flow_circle_mov_task_description.h
* \author    Puzikova Valeria (vvp) 
* \date      20.10.2013 23:00:42
* \brief     flow_circle_mov_task_description.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FLOW_CIRCLE_MOV_TASK_DESCRIPTION_H_
#define _LS_STAG_FLOW_CIRCLE_MOV_TASK_DESCRIPTION_H_

// ================================================================== INCLUDES
#include "flow_circle_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FlowCircleMovTaskDescription

Данный файл содержит в себе реализацию класса FlowCircleMovTaskDescription
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с задачей об обтекании движущегося кругового профиля

Данный класс реализует работу с описанием решаемой задачи об обтекании движущегося кругового профиля
*/
	struct FlowCircleMovTaskDescription : FlowCircleTaskDescription
	{
		explicit FlowCircleMovTaskDescription(OscKinds _id_osc = OSC_UNKNOWN, int _size = sizeof(FlowCircleMovTaskDescription)):
		FlowCircleTaskDescription(MOTION_MOV, _size),
		id_osc(_id_osc)
		{
			m_mesh = (LSSTAGUniformMesh*)&meshMov;// изменяем mesh
		}
		void checkExistMov()// 1.06.14
		{
			if(!existMov)
			{
				if((*m_curTime) > t_0) existMov = true;
			}
		};
	  protected:
		void loadDiscontin(std::ifstream &in){  };
		void printDiscontin(std::ofstream &out){  };
		virtual void initMeshParams()
		{
			FlowTaskDescription::initMeshParams();
			initMovMesh();
		};
		int getOscKind(){ return (int)id_osc; };
		void buildVib()
		{
			meshMov.rebuildLinearVelocity();
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
		void buildLS(const int zs, const int zf, const int vs, const int vf, const int hs, const int hf) override
		{
			LS_STAG_PROFILE("buildLS");
			LS_STAG_PAR_FOR_2D_RANGE(FlowCircleTaskDescription, calcLevelSetCircle, vs, vf, zs, zf, hs, hf, m_phi,
				m_mesh->nodes, m_mesh->delta, *(meshMov.pCur_C), m_mesh->baseMesh.N, m_mesh->baseMesh.N * m_mesh->baseMesh.M, rad)
		};
		double getDistance(const double x, const double y, const double z)
		{
			double aux = x - meshMov.pCur_C->x;
			double aux2 = y - meshMov.pCur_C->y;
			double aux3 = z - meshMov.pCur_C->z;
			double dist = -rad + sqrt(aux * aux + aux2 * aux2 + aux3 * aux3);
			return dist;
		};
		LSSTAGMeshMov meshMov;///< неравномерная LS-STAG-сетка c подвижными погруженными границами
	  public:
		double t_0;///< время начала движения профиля 
		OscKinds id_osc;///< идентификатор типа колебаний профиля
	};
}

#endif //_LS_STAG_FLOW_CIRCLE_MOV_TASK_DESCRIPTION_H_
