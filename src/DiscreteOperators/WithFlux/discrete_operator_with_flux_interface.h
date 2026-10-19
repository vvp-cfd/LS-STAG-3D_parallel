/*
* \copyright (C) 2013 Puzikova V.V.
* \file      discrete_operator_with_flux_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      01.10.2013 23:07:53
* \brief     discrete_operator_with_flux_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DISCRETE_OPERATOR_WITH_FLUX_INTERFACE_H_
#define _LS_STAG_DISCRETE_OPERATOR_WITH_FLUX_INTERFACE_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
#include "../../Base/MeshFunctions/mesh_function.h"
#include "../../Base/TaskDescription/Tasks/Flow/flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DiscreteOperatorWithFluxInterface

Данный файл содержит в себе реализацию класса DiscreteOperatorWithFluxInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий интерфейс дискретного оператора с источниковым членом

Данный класс реализует интерфейс дискретного оператора с источниковым членом
*/
	struct DiscreteOperatorWithFluxInterface : DiscreteOperatorInterface
	{	
		typedef BaseTaskDescription::TaskKinds kind;
		typedef BaseTaskDescription::MotionKinds motion;
		typedef FlowTaskDescription::AirfoilKind airfoil;
		DiscreteOperatorWithFluxInterface():
		vs(0),
		vf(0),
		zs(0),
		zf(0),
		baseRecalcFlux(false), 
		ibRecalcFlux(false)
		{
			;
		}
	  protected:
		virtual void setRecalcFluxMarkers(){};
		virtual void ibFlux(){};
	  public:
		virtual void calcFlux(){};
		virtual void recalc(){};
		virtual void recalcTDF(){};
	  protected:
		  using SparseMatrix::set;
		virtual void set(const int _N, const int _M, const int _K)
		{
			SparseMatrix::set(_N, _M, _K);
			N2 = m_N - 1;
			N3 = m_N + 1;// Grid.N для х-й сетки
			M2 = m_M - 1;
			M3 = m_M + 1;
			K2 = m_K - 1;
			K3 = m_K + 1;
			S.init(NMK, flux);
			S.set(0.0);
		};
	  public:
		  using DiscreteOperatorInterface::initMeshParams;
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorInterface::initMeshParams(_cur_mesh, _meshes);
			m_Vib = m_meshes->v_ib;
		};
		void init(BaseTaskDescription* _task, LSSTAGBaseMeshParams* _cur_mesh)
		{
			m_task = _task;
			if(m_task->id == kind::TASK_FLOW)
				flow_task = (FlowTaskDescription*)m_task;
			setRecalcFluxMarkers();
			initMeshParams(_cur_mesh, m_task->m_mesh);//+set
		};
		  /*!
		  Строит разреженную матрицу дискретного оператора и массив значений источникового члена
		  */
		using DiscreteOperatorInterface::build;
		virtual void build()
		{
			DiscreteOperatorInterface::build();
			calcFlux();
		};
		void printFlux(const std::string& name)
		{
			S.print(name + "_flux");
		};
		void print(const std::string &name)
		{
			SparseMatrix::print(name);
			printFlux(name);
		};
		void printNotSparse(const std::string &name)
		{
			SparseMatrix::printNotSparse(name);
			printFlux(name);
		};
		extDouble S; ///< сеточная функция источникового члена
		double* flux;///< указатель на массив значений источникового члена
	  protected:
		BaseTaskDescription* m_task;///< указатель на структуру с параметрами задачии
		FlowTaskDescription* flow_task;
		coordinate *m_Vib, *curVib;///< массивы скоростей на погруженной границе
		int N2, N3, M2, M3, K2, K3, vs, vf, zs, zf;
		bool baseRecalcFlux, ibRecalcFlux;///< флаги пересчета источникового члена
	};
}
#endif //_LS_STAG_DISCRETE_OPERATOR_WITH_FLUX_INTERFACE_H_
