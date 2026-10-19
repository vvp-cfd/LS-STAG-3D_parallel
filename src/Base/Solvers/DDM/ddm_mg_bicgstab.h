/*
* \copyright (C) 2019 Puzikova V.V.
* \file      ddm_mg_bicgstab.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2019 23:01:45
* \brief     ddm_mg_bicgstab.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DDM_MG_BICGSTAB_H_
#define _LS_STAG_DDM_MG_BICGSTAB_H_

// ================================================================== INCLUDES
#include "ddm_interface.h"
#include "../BiCGStab/mg_bicgstab.h"
#include "../BiCGStab/ilu_bicgstab.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DDMMGBiCGStab

Данный файл содержит в себе реализацию класса DDMMGBiCGStab
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод DDM + BiCGStab + MG

Данный класс реализует метод декомпозиции области с использованием метода BiCGStab с многосеточным предобуславливанием 
для решения задач в подобластях
*/
	struct DDMMGBiCGStab : DDMInterface
	{
	  private:
		virtual IterativeSolverInterface* getNewSolver(){ return ((IterativeSolverInterface*)(new MGBiCGStab)); } ;
	};
};

#endif //_LS_STAG_DDM_MG_BICGSTAB_H_
