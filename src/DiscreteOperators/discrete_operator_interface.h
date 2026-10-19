/*
* \copyright (C) 2013 Puzikova V.V.
* \file      discrete_operator_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      01.10.2013 12:00:15
* \brief     discrete_operator_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DISCRETE_OPERATOR_INTERFACE_H_
#define _LS_STAG_DISCRETE_OPERATOR_INTERFACE_H_

// ================================================================== INCLUDES
#include "../Base/SparseMatrix/sparse_matrix.h"
#include "../Base/LS_STAG_Mesh/cut_cells.h"
#include "../Base/TaskDescription/base_task_description.h"
#include "../Base/MeshFunctions/coord_mesh_function.h"
#include "../Base/LS_STAG_Mesh/ls_stag_uniform_mesh.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DiscreteOperatorInterface

Данный файл содержит в себе реализацию класса SmootheDiscreteOperatorInterfacerADLJ
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий интерфейс дискретного оператора

Данный класс реализует интерфейс дискретного оператора
*/
	struct DiscreteOperatorInterface : SparseMatrix
	{
	  protected:
		virtual void buildNum(){};
		  /*!
		  Строит портрет разреженной матрицы дискретного оператора
		  */
		virtual void buildPortait(){};
	  public:
		  /*!
		  Рассчитывает значения элементов разреженной матрицы дискретного оператора
		  */
		virtual void calcElem(){};
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes, CellTypes* _type) { ; }
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			m_cur_mesh = _cur_mesh;
			m_meshes = _meshes;
			m_type = m_cur_mesh->typeOfCell;
			m_baseSideTypes = m_meshes->baseSideTypes;
			m_delta = m_meshes->delta;
			m_theta = m_meshes->theta;
			set(m_cur_mesh->N, m_cur_mesh->M, m_cur_mesh->K);
		};
		  /*!
		  Строит разреженную матрицу дискретного оператора
		  */
		virtual void build()
		{
			buildNum();
			resize();// memory Portrait + Cell
			buildPortait();
			calcElem();
		}
	  protected:
		LSSTAGBaseMeshParams* m_cur_mesh;///< указатель на параметры сетки, на которой строим дискрестный оператор
		LSSTAGUniformMesh* m_meshes;;///< указатель на все разнесенные сетки LS-STAG-сетки
		coordinate* m_delta;///< шаги по пространству --- c основной сетки
		coordinate* m_theta;///< коэффициенты заполнения ячеек --- c основной сетки
		CellTypes* m_type;///< типы ячеек сетки, на которой строим дискрестный оператор
		coord<SideTypes>* m_baseSideTypes;///< типы граней основной сетки
	};
}
#endif //_LS_STAG_DISCRETE_OPERATOR_INTERFACE_H_
