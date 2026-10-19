/*
* \copyright (C) 2015 Puzikova V.V.
* \file      mg_fgmres.h
* \author    Puzikova Valeria (vvp) 
* \date      08.08.2015 3:11:32
* \brief     mg_fgmres.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_MG_FGMRES_H_
#define _LS_STAG_MG_FGMRES_H_

// ================================================================== INCLUDES
#include "fgmres_interface.h"
#include "../../CoarseGrid/multigrid_solver.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса MGFGMRES

Данный файл содержит в себе реализацию класса MGFGMRES
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод FGMRES c многосеточным предобуславливанием 

Данный класс реализует многосеточное предобуславливание итерационного метода FGMRES (построение предобуславливателя, 
его перестроение и применение предобуславливателя к СЛАУ)
*/
	struct MGFGMRES : FGMRESInterface
	{
		void rebuildPreconditioner()
		{
			LS_STAG_PROFILE("rebuildMG");
			MG.rebuild();
		};
		void printPreconditioner(const std::string &name){ MG.print(name); };
	  private:
		void setPreconditioner(SparseMatrix* _A, BaseTaskDescription* _task)
		{
			MG.set(_A, _task);
		};
		void preconditioner(const double* _b, double* _x)
		{
			LS_STAG_PROFILE("MG");
			MG.solve(_b, _x);
		};
		MultiGridSolver MG;///< Многосеточный решатель для предобуславливания
	};
}

#endif //_LS_STAG_MG_FGMRES_H_
