/*
* \copyright (C) 2013 Puzikova V.V.
* \file      mg_bicgstab.h
* \author    Puzikova Valeria (vvp) 
* \date      30.10.2013 20:29:13
* \brief     mg_bicgstab.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_MG_BICGSTAB_H_
#define _LS_STAG_MG_BICGSTAB_H_

// ================================================================== INCLUDES
#include "bicgstab_interface.h"
#include "../../CoarseGrid/multigrid_solver.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса MGBiCGStab

Данный файл содержит в себе реализацию класса MGBiCGStab
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод BiCGStab c многосеточным предобуславливанием 

Данный класс реализует многосеточное предобуславливание итерационного метода BiCGStab (построение предобуславливателя, 
его перестроение и применение предобуславливателя к СЛАУ)
*/
	struct MGBiCGStab : BiCGStabInterface
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

#endif //_LS_STAG_MG_BICGSTAB_H_
