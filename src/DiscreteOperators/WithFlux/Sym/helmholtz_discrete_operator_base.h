/*
* \copyright (C) 2013 Puzikova V.V.
* \file      helmholtz_discrete_operator_base.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 22:10:54
* \brief     helmholtz_discrete_operator_base.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_BASE_H_
#define _LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_BASE_H_

// ================================================================== INCLUDES
#include "sym_discrete_operator.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса HelmholtzDiscreteOperatorBase

Данный файл содержит в себе реализацию класса HelmholtzDiscreteOperatorBase
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(addPlusOutFlow, LS_STAG_R_PD, _Cell, LS_STAG_CR_PI, _numOutFlow, LS_STAG_CR_PD, _PlusOutFlow);
	LS_STAG_DECLARATION_S_KERNEL(numPlusOutFlow, LS_STAG_R_PI, numOutFlow, LS_STAG_CR_PI, Number, const int, m_N, const int, m_M, const int, m_K);
/*!
\brief Класс, реализующий интерфейс дискретного оператора Гельмгольца

Данный класс реализует интерфейс дискретного аналога оператора Гельмгольца
*/
	struct HelmholtzDiscreteOperatorBase : SymDiscreteOperator
	{
		HelmholtzDiscreteOperatorBase():
		existOutFlow(false)
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(numPlusOutFlow, LS_STAG_R_PI numOutFlow, LS_STAG_CR_PI Number, const int m_N, const int m_M, const int m_K)
		{//PlusOutFlow --- отдельно по x и y
			int num = m_N - 1;
			numOutFlow[0] = Number[num] + 1; num += m_N;
			for (int j = 1; j < m_M; j++, num += m_N)
				numOutFlow[j] = Number[num] + 2;
			for(int k = 1, numE = m_M; k < m_K; k++)
			{
				numOutFlow[numE] = Number[num] + 2;
				num += m_N; numE++;
			for (int j = 1; j < m_M; j++, num += m_N, numE++)
				numOutFlow[numE] = Number[num] + 3;
			}
		}
		LS_STAG_KERNEL(addPlusOutFlow, LS_STAG_R_PD _Cell, LS_STAG_CR_PI _numOutFlow, LS_STAG_CR_PD _PlusOutFlow)
		{
			_Cell[_numOutFlow[i]] += _PlusOutFlow[i];
		}
		void recalc()
		{
			calcFlux();// здесь ibFlux, + OutFlow --- поэтому строго перед calcElem
			calcElem();
		};
		void recalcTDF()
		{
			initNu();
			recalc();
			if(existOutFlow) calcPlusOutFlow();
			initInFlowFlux();
		}
		void build()
		{
			buildNum();
			resize();// memory Portrait + Cell
			buildPortait();
			if(existOutFlow) initNumPlusOutFlow();
			recalc();// + OutFlow --- поэтому строго перед calcElem
		}
		using DiscreteOperatorWithFluxInterface::initMeshParams;
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorWithFluxInterface::initMeshParams(_cur_mesh, _meshes);
			initNu();
			coefIntegr = m_task->COEF_INTEGR_TIME;
			m_inv_V = m_meshes->baseMesh.inv_VOF;
			m_alpha = m_meshes->alpha;
			m_delta_ib_en = m_meshes->delta_ib_en;
			m_delta_ib_ws = m_meshes->delta_ib_ws;
			m_delta_ib_nb = m_meshes->delta_ib_nb;
			m_delta_ib_st = m_meshes->delta_ib_st;
			m_inv_delta_ib_ne = m_meshes->inv_delta_ib_ne;
			m_inv_delta_ib_tn = m_meshes->inv_delta_ib_tn;
				m_VibYX = m_meshes->v_ib_yx;
				m_VibZY = m_meshes->v_ib_zy;
				m_Vib = m_meshes->v_ib;
			initFlux();// вклад от стенок РО, на которых г.у. I рода (соот-но строим один раз)
			if(m_task->id == kind::TASK_FLOW) 
			{
				existOutFlow = true;
				initPlusOutFlow();//7.11.14
			}
		};
	  private:
		void initNu()
		{
			m_nu = m_task->nu;
		}
	  protected:
		virtual void initInFlowFlux(){};
		virtual void calcPlusOutFlow(){};
		virtual void initPlusOutFlow()
		{
			m_PlusOutFlow.init(m_M * m_K, PlusOutFlow);
			m_numOutFlow.init(m_M * m_K, numOutFlow);
			calcPlusOutFlow();
		};
		void initNumPlusOutFlow()
		{
			LS_STAG_SINGLE(HelmholtzDiscreteOperatorBase, numPlusOutFlow, numOutFlow, Number, m_N, m_M, m_K)
		};
		void addOutFlow()
		{
			if(existOutFlow)
			{
				LS_STAG_PAR_FOR(HelmholtzDiscreteOperatorBase, addPlusOutFlow, m_M * m_K, Cell, numOutFlow, PlusOutFlow)
			}
		};
		void initFlux()
		{
			initInFlowFlux();
		};
		Array<double> m_PlusOutFlow;
		Array<int> m_numOutFlow;
		double m_nu, coefIntegr;
		double *m_inv_V, *PlusOutFlow;
		int *numOutFlow;
		coordinate *m_inv_delta_ib_ne, *m_inv_delta_ib_tn, *m_alpha;
		coordinate *m_delta_ib_en, *m_delta_ib_ws, *m_delta_ib_nb, *m_delta_ib_st;
		coordinate *m_VibYX, *m_VibZY, *m_Vib;
		bool existOutFlow;
	};
};

#endif //_LS_STAG_HELMHOLTZ_DISCRETE_OPERATOR_BASE_H_
