/*
* \copyright (C) 2026 Puzikova V.V.
* \file      discrete_operator_only_flux_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      04.06.2026 12:25:13
* \brief     discrete_operator_only_flux_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DISCRETE_OPERATOR_ONLY_FLUX_INTERFACE_H_
#define _LS_STAG_DISCRETE_OPERATOR_ONLY_FLUX_INTERFACE_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
#include "../../Base/MeshFunctions/mesh_function.h"
#include "../../Base/TaskDescription/Tasks/Flow/flow_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DiscreteOperatorOnlyFluxInterface

Данный файл содержит в себе реализацию класса DiscreteOperatorOnlyFluxInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий интерфейс фиктивного дискретного оператора с источниковым членом

Данный класс реализует интерфейс фиктивного дискретного оператора с источниковым членом
*/
	struct DiscreteOperatorOnlyFluxInterface : DiscreteOperatorInterface
	{	
		DiscreteOperatorOnlyFluxInterface():
		m_operator(NULL)
		{
			;
		}
		virtual void calcFlux(){};// расчет итогового потока 
		using SparseMatrix::set;
		using DiscreteOperatorInterface::initMeshParams;
	  protected:
		void set(const int _N, const int _M, const int _K)
		{
			setSizes(_N, _M, _K);
			S.init(NMK, flux);
			S.set(0.0);
			ibS.init(NMK, ibFlux);
			ibS.set(0.0);
		};
	  public:
		void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			m_cur_mesh = _cur_mesh;
			m_meshes = _meshes;
			m_type = m_cur_mesh->typeOfCell;
			m_delta = m_meshes->delta;
			m_theta = m_meshes->theta;
			set(m_cur_mesh->N, m_cur_mesh->M, m_cur_mesh->K);
		};
		void set(SparseMatrix* _operator) 
		{
			m_operator = _operator;
		} 
		void build()
		{
			calcElem();// пересчет ibFlux при перестроении сетки
		};
		virtual void print(const std::string &name)
		{
			S.print(name+"_flux");
		};
		void printNotSparse(const std::string &name)
		{
			print(name);
		};
		RealValuedMeshFunction S;  ///< сеточная функция итогового источникового члена
	  protected:
		extDouble ibS;
		double* ibFlux;
	    SparseMatrix* m_operator;///< указатель на оператор (без потока, к которму здесь он достраивается)
	public:
		double* flux;///< указатель на массив значений источникового члена
	};
}
#endif //_LS_STAG_DISCRETE_OPERATOR_ONLY_FLUX_INTERFACE_H_
