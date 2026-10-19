/*
* \copyright (C) 2019 Puzikova V.V.
* \file      interface_block_solver_mid.h
* \author    Puzikova Valeria (vvp) 
* \date      09.08.2019 23:21:25
* \brief     interface_block_solver_mid.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_INTERFACE_BLOCK_SOLVER_MID_H_
#define _LS_STAG_INTERFACE_BLOCK_SOLVER_MID_H_

// ================================================================== INCLUDES
#include "interface_block_solver_last.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса InterfaceBlockSolverMid

Данный файл содержит в себе реализацию класса InterfaceBlockSolverMid
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий решение задачи на промежуточных интерфейсных границах

Данный класс реализует прямой и обратный ход метода Гаусса для блоков интерфейсной матрицы, соответствующих 
промежуточным интерфейсным границам (т.е. не первой и не последней)
*/
	struct InterfaceBlockSolverMid : InterfaceBlockSolverLast
	{
		void directGauss()
		{// для средних интерфейсов
			InterfaceBlockSolverLast::directGauss();
			LU.solve();
		}
		void inverseGauss()
		{// для средних интерфейсов
			prepareIG();
			InterfaceBlockSolverLast::inverseGauss();
		}
	};
}
#endif //_LS_STAG_INTERFACE_BLOCK_SOLVER_MID_H_