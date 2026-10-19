/*
* \copyright (C) 2013 Puzikova V.V.
* \file      gradient_discrete_operator_base.h
* \author    Puzikova Valeria (vvp) 
* \date      27.10.2013 13:00:49
* \brief     gradient_discrete_operator_base.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_GRADIENT_DISCRETE_OPERATOR_BASE_H_
#define _LS_STAG_GRADIENT_DISCRETE_OPERATOR_BASE_H_

// ================================================================== INCLUDES
#include "../discrete_operator_interface.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса GradientDiscreteOperatorBase

Данный файл содержит в себе реализацию класса GradientDiscreteOperatorBase
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(numBuildG, LS_STAG_R_PI, Number, LS_STAG_CR_PCT, m_type, const int, NMK);
/*!
\brief Класс, реализующий интерфейс дискретного оператора градиента

Данный класс реализует интерфейс дискретного аналога оператора градиента
*/
	struct GradientDiscreteOperatorBase : DiscreteOperatorInterface
	{
		LS_STAG_SINGLE_KERNEL(numBuildG, LS_STAG_R_PI Number, LS_STAG_CR_PCT m_type, const int NMK)
		{
			for (int i = 0, num = 0; i < NMK; i++, Number[i] = num)
			{
				if (m_type[i] == SolidCell)
					num++;
				else
					num += 2;
			}
		}
		using DiscreteOperatorInterface::initMeshParams;
		virtual void initMeshParams(LSSTAGBaseMeshParams* _cur_mesh, LSSTAGUniformMesh* _meshes)
		{
			DiscreteOperatorInterface::initMeshParams(_cur_mesh, _meshes);
			setBaseN();
			build();
		};
	  private:
		  /*!
		  Оценивает максимальное количество ненулевых элементов матрицы дискретного аналога оператора градиента
		  */
		void estimateNeedMemory()
		{
			maxElemCount = 2 * NMK;// max кол-во эл-в м-цы (если нет твёрдых ячеек)
		};
		void buildNum()
		{
			LS_STAG_SINGLE(GradientDiscreteOperatorBase, numBuildG, Number, m_type, NMK)
		};
	  protected:
		virtual void setBaseN(){ baseN = m_N; };
		int baseN;
	};
}

#endif //_LS_STAG_GRADIENT_DISCRETE_OPERATOR_BASE_H_
