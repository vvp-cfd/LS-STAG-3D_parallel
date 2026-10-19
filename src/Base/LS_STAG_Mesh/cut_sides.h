/*
* \copyright (C) 2026 Puzikova V.V.
* \file      cut_sides.h
* \author    Puzikova Valeria (vvp) 
* \date      30.05.2026 23:08:23
* \brief     cut_sides.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_CUT_SIDES_H_
#define _LS_STAG_CUT_SIDES_H_

// ================================================================== INCLUDES
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием перечисления SideTypes

Данный файл содержит в себе перечисление типов усеченных двумерных граней 
трехмерных ячеек SideTypes
*/
// ===========================================================================

namespace LS_STAG
{	
	/// Типы усеченных граней
	enum SideTypes
	{
		UnknownSide = -1,///< Неизвестен
		SolidSide,     ///< 0 --- полностью твёрдая грань 
		NW_Triangle = 1,
		NE_Triangle = 10,
		SW_Triangle = 1000, 
		SE_Triangle = 100, 
		N_Trapezoid = 11,
		S_Trapezoid = 1100,
		E_Trapezoid = 110, 
		W_Trapezoid = 1001, 
		NW_Pentagon = 1011,
		NE_Pentagon = 111,
		SW_Pentagon = 1101,
		SE_Pentagon = 1110,
		FluidSide = 1111,     ///< 1111 --- полностью жидкая грань
	};
}
#endif //_LS_STAG_CUT__SIDES_H_
