/*
* \copyright (C) 2019 Puzikova V.V.
* \file      ddm_mg_fgmres.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2019 22:57:40
* \brief     ddm_mg_fgmres.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DDM_MG_FGMRES_H_
#define _LS_STAG_DDM_MG_FGMRES_H_

// ================================================================== INCLUDES
#include "ddm_interface.h"
#include "../FGMRES/mg_fgmres.h"
#include "../FGMRES/ilu_fgmres.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DDMMGFGMRES

Данный файл содержит в себе реализацию класса DDMMGFGMRES
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод DDM + FGMRES + MG

Данный класс реализует метод декомпозиции области с использованием метода FGMRES с многосеточным предобуславливанием 
для решения задач в подобластях
*/
	struct DDMMGFGMRES : DDMInterface
	{
	  private:
		virtual IterativeSolverInterface* getNewSolver(){ return ((IterativeSolverInterface*)(new MGFGMRES)); } ;
	};
};

#endif //_LS_STAG_DDM_MG_FGMRES_H_
