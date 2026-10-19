/*
* \copyright (C) 2026 Puzikova V.V.
* \file      gradient_discrete_operator_flux_base.h
* \author    Puzikova Valeria (vvp) 
* \date      04.06.2036 21:50:25
* \brief     gradient_discrete_operator_flux_base.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_BASE_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_BASE_H_

// ================================================================== INCLUDES
#include "../discrete_operator_only_flux_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorFluxBase

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorFluxBase
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numBuildG, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NM);
/*!
\brief Класс, реализующий интерфейс дискретного оператора градиента

Данный класс реализует интерфейс дискретного аналога оператора градиента
*/
	struct GradientDiscreteOperatorFluxBase : DiscreteOperatorOnlyFluxInterface
	{
		GradientDiscreteOperatorFluxBase():
		DiscreteOperatorOnlyFluxInterface(),
		NM(0),
		needSecondFlux(false)
		{
			;
		}
		void calcFlux()// расчет итогового потока 
		{// (M)^{-1}*G * func + g^{ib} * (func - ibValue)
			S.setProductScaleP(ibFlux);
			if(needSecondFlux) S2.setProductScalePM(ibFlux, *ibValue);
		};
		void setSquare() { S.setSquare(); }
		using DiscreteOperatorOnlyFluxInterface::initMeshParams;
		using DiscreteOperatorOnlyFluxInterface::print;
		using DiscreteOperatorOnlyFluxInterface::set;
		void init(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes, double*& copyFlux, double *_func)
		{
			m_func = _func;
			initMeshParams(_cur_mesh, _meshes);
			copyFlux = flux;
		}
		void init(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes, double*& copyFlux, double*& copyFlux2, 
			double *_func, double *_func2, double* _ibValue)
			{
			ibValue = _ibValue;
			m_func2 = _func2;
			needSecondFlux = true;
			init(_cur_mesh, _meshes, copyFlux, _func);
			copyFlux2 = flux2;
			}
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorOnlyFluxInterface::initMeshParams(_cur_mesh, _meshes);// здесь же set
			m_inv_delta_ib = m_meshes->inv_delta_ib;
			build();
		};
		void set(SparseMatrix* _operator) 
		{
			DiscreteOperatorOnlyFluxInterface::set(_operator);
			S.setMultipliers(m_operator, m_func, m_cur_mesh->inv_VOF);
			if(needSecondFlux) S2.setMultipliers(m_operator, m_func2, m_cur_mesh->inv_VOF);
		} 
		void print(const std::string &name)
		{
			DiscreteOperatorOnlyFluxInterface::print(name);
			if(needSecondFlux) S2.print(name+"_flux2");
		};
	  private:
	    void set(const int _N, const int _M, const int _K)
		{
			DiscreteOperatorOnlyFluxInterface::set(_N, _M, _K);
			NM = m_N * m_M;
			if(needSecondFlux)
			{
			S2.init(NMK, flux2);
			S2.set(0.0);
			}
		};
		RealValuedMeshFunction S2;  ///< сеточная функция второго итогового источникового члена
		double* ibValue;
		bool needSecondFlux;
	protected:
		coordinate* m_inv_delta_ib;
		double *m_func, *m_func2;///< указатели на сеточные функции, от которых считаем градиент
	public:
		double* flux2;///< указатель на массив значений второго источникового члена
		protected:
		int NM;
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_FLUX_BASE_H_
