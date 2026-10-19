/*
* \copyright (C) 2019 Puzikova V.V.
* \file      ddm_ilu_fgmres.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2019 22:39:05
* \brief     ddm_ilu_fgmres.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DDM_ILU_FGMRES_H_
#define _LS_STAG_DDM_ILU_FGMRES_H_

// ================================================================== INCLUDES
#include "ddm_interface.h"
#include "../FGMRES/ilu_fgmres.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DDMILUFGMRES

Данный файл содержит в себе реализацию класса DDMILUFGMRES
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод DDM + FGMRES + ILU

Данный класс реализует метод декомпозиции области с использованием метода FGMRES с ILU-предобуславливанием 
для решения задач в подобластях
*/
	struct DDMILUFGMRES : DDMInterface
	{
	  private:
		virtual IterativeSolverInterface* getNewSolver(){ return ((IterativeSolverInterface*)(new ILUFGMRES)); } ;
	};
};

#endif //_LS_STAG_DDM_ILU_FGMRES_H_
