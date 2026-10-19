/*
* \copyright (C) 2019 Puzikova V.V.
* \file      ddm_ilu_bicgstab.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2019 22:59:15
* \brief     ddm_ilu_bicgstab.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DDM_ILU_BICGSTAB_H_
#define _LS_STAG_DDM_ILU_BICGSTAB_H_

// ================================================================== INCLUDES
#include "ddm_interface.h"
#include "../BiCGStab/ilu_bicgstab.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DDMILUBiCGStab

Данный файл содержит в себе реализацию класса DDMILUBiCGStab
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий метод DDM + BiCGStab + ILU

Данный класс реализует метод декомпозиции области с использованием метода BiCGStab с ILU-предобуславливанием 
для решения задач в подобластях
*/
	struct DDMILUBiCGStab : DDMInterface
	{
	  private:
		virtual IterativeSolverInterface* getNewSolver(){ return ((IterativeSolverInterface*)(new ILUBiCGStab)); } ;
	};
};

#endif //_LS_STAG_DDM_ILU_BICGSTAB_H_
