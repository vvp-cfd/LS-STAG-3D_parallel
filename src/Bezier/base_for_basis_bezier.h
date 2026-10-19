/*
* \copyright (C) 2013 Puzikova V.V.
* \file      base_for_basis_bezier.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 20:36:04
* \brief     base_for_basis_bezier.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BASE_FOR_BASIS_BEZIER_H_
#define _LS_STAG_BASE_FOR_BASIS_BEZIER_H_

// ================================================================== INCLUDES
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BaseForBasisBezier

Данный файл содержит в себе реализацию класса BaseForBasisBezier
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с базисной матрицей Безье

Данный класс реализует хранение и построение базисной матрицы Безье для последующего расчета кривых Безье
*/
	struct BaseForBasisBezier
	{
		static const int LIMIT_BEZIER = 30;///< максимальное число точек в каждой из частей профиля
		static const int LB2 = LIMIT_BEZIER + 1;
		BaseForBasisBezier()
		{
			for(int j = 0; j < LIMIT_BEZIER; j++)
			{
				BaseBezier[0][j] = 1;
				j++;
				BaseBezier[0][j] = -1;
			}
			BaseBezier[0][LIMIT_BEZIER] = 1;
			for(int j = 0; j < LB2; j++)
				BaseBezier[1][j] = (j - LIMIT_BEZIER) * BaseBezier[0][j];
			for(int j = 0, i3 = LIMIT_BEZIER - 1; j < i3; j++)
				BaseBezier[LIMIT_BEZIER-j][j] = 1;
			for(int j = 0, i2 = LIMIT_BEZIER, i3 = i2 - 2; j < i3; j++, i2--)
				BaseBezier[i2-1][j] = -i2;
			for(int j = 2, i2 = LIMIT_BEZIER - j; j < LB2; j++, i2--)
			{	
				for(int i = i2 - 2; i > -1; i--)
					BaseBezier[j][i] = -BaseBezier[j][i+1] * (LIMIT_BEZIER - i) / (i2 - i);
				for(int i = LB2 - j; i < LB2; i++)
					BaseBezier[j][i] = 0;
			};
		}
		int BaseBezier[LB2][LB2];///< базисная матрица Безье // TODO  to constant memory 
	};
}

#endif //_LS_STAG_BASE_FOR_BASIS_BEZIER_H_
