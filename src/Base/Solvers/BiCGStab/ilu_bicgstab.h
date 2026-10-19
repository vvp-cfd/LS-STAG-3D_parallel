/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ilu_bicgstab.h
* \author    Puzikova Valeria (vvp) 
* \date      30.10.2013 18:36:54
* \brief     ilu_bicgstab.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ILU_BICGSTAB_H_
#define _LS_STAG_ILU_BICGSTAB_H_

// ================================================================== INCLUDES
#include "bicgstab_interface.h"
#include "../../SparseMatrix/Decomposition/ilu_decomposition_7p_stencil_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ILUBiCGStab

Данный файл содержит в себе реализацию класса ILUBiCGStab
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод BiCGStab c ILU-предобуславливанием 

Данный класс реализует ILU-предобуславливание итерационного метода BiCGStab (построение предобуславливателя, 
его перестроение и применение предобуславливателя к СЛАУ)
*/
	struct ILUBiCGStab : BiCGStabInterface
	{
		ILUBiCGStab()
		{
#ifdef GENERAL_SOLVER
			ILU = new ILUDecompositionSparseMatrix;
#else 
			ILU = new ILUDecomposition7pStencilMatrix;// 3d ILUDecomposition5pStencilMatrix; // 2d
#endif
		}
		~ILUBiCGStab()
		{
			delete ILU;
		}
		void rebuildPreconditioner()
		{
			LS_STAG_PROFILE("rebuildILU");
			ILU->rebuild();// 24.11.14 decompose();
		};
		void printPreconditioner(const std::string &name){ ILU->print(name); };
	  private:
		void setPreconditioner(SparseMatrix* _A, BaseTaskDescription* _task)
		{
			ILU->set(_A);
		};
		void preconditioner(const double* _b, double* _x)
		{
			LS_STAG_PROFILE("ILU");
			ILU->solve(_b, _x);
		};
		ILUDecompositionSparseMatrix* ILU;///< ILU-декомпозиция матрицы СЛАУ для предобуславливания
	};
}

#endif //_LS_STAG_ILU_BICGSTAB_H_
