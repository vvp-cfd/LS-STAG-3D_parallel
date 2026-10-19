/*
* \copyright (C) 2015 Puzikova V.V.
* \file      ilu_fgmres.h
* \author    Puzikova Valeria (vvp) 
* \date      08.08.2015 3:05:16
* \brief     ilu_fgmres.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_ILU_FGMRES_H_
#define _LS_STAG_ILU_FGMRES_H_

// ================================================================== INCLUDES
#include "fgmres_interface.h"
#include "../../SparseMatrix/Decomposition/ilu_decomposition_7p_stencil_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ILUFGMRES

Данный файл содержит в себе реализацию класса ILUFGMRES
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод FGMRES c ILU-предобуславливанием 

Данный класс реализует ILU-предобуславливание итерационного метода FGMRES (построение предобуславливателя, 
его перестроение и применение предобуславливателя к СЛАУ)
*/
	struct ILUFGMRES : FGMRESInterface
	{
		ILUFGMRES()
		{
#ifdef GENERAL_SOLVER
			ILU = new ILUDecompositionSparseMatrix;
#else 
			ILU = new ILUDecomposition7pStencilMatrix;// 3d ILUDecomposition5pStencilMatrix; // 2d
#endif
		}
		~ILUFGMRES()
		{
			delete ILU;
		}
		void rebuildPreconditioner()
		{
			LS_STAG_PROFILE("rebuildILU");
			ILU->rebuild();// 24.11.14 decompose();
		};
		void printPreconditioner(const std::string& name) { ILU->print(name); };
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

#endif //_LS_STAG_ILU_FGMRES_H_
