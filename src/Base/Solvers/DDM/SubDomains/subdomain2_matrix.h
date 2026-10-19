/*
* \copyright (C) 2019 Puzikova V.V.
* \file      subdomain2_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      13.04.2019 00:37:39
* \brief     subdomain2_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SUBDOMAIN2_MATRIX_H_
#define _LS_STAG_SUBDOMAIN2_MATRIX_H_

// ================================================================== INCLUDES
#include "subdomain_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса Subdomain2Matrix

Данный файл содержит в себе реализацию класса Subdomain2Matrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий матрицу СЛАУ в промежуточной подобласти

Данный класс реализует построение по разреженной матрице СЛАУ в полной расчетной области построение матрицы СЛАУ в 
промежуточной подобласти (имеющей две соседние) и ее хранение
*/
	struct Subdomain2Matrix : SubdomainMatrix
	{
		/*!
		Строит портрет матрицы в подобласти, у которой две соседние, и заполняет ее по исходной матрице 
		СЛАУ для полной расчетной области
		*/
		void rebuild()
		{
			LS_STAG_PROFILE("rebuildSubdom");
			elemCount = 0;
			for(int i = 0, j = startNum + 1, num = m_A_Num[startNum], num2 = m_A_Num[j]; i < NM; i++, j++, num = num2, num2 = m_A_Num[j], Number[i] = elemCount)
			{// 1-й ряд без B
				if(num2 > (num + 1))
				{
					for(int k = num + 1; k < num2; k++, elemCount++)
					{
						Portrait[elemCount] = m_A_Portrait[k] - startNum;
						Cell[elemCount] = m_A_Cell[k];
					}
				}
				else
				{// на границу попала твердая ячейка
					Portrait[elemCount] = i;
					Cell[elemCount] = PLUS_TINY;//0.0;
					elemCount++;
				}
			}
			for(int i = NM, j = bNum + 1, num = m_A_Num[bNum], num2 = m_A_Num[j]; i < nNum; i++, j++, Number[i] = elemCount, num = num2, num2 = m_A_Num[j])
			{// берем матрицу как есть
				for(int k = num; k < num2; k++, elemCount++)
				{
					Portrait[elemCount] = m_A_Portrait[k] - startNum;
					Cell[elemCount] = m_A_Cell[k];
				}
			}
			for(int i = nNum, j = bNum2 + 1, num = m_A_Num[bNum2], num2 = m_A_Num[j]; i < NMK; i++, j++, num = num2, num2 = m_A_Num[j], Number[i] = elemCount)
			{// последний ряд без T
				if(num2 > (num + 1))
				{
					for(int k = num, end  = num2 - 1; k < end; k++, elemCount++)
					{
						Portrait[elemCount] = m_A_Portrait[k] - startNum;
						Cell[elemCount] = m_A_Cell[k];
					}
				}
				else
				{// на границу попала твердая ячейка
					Portrait[elemCount] = i;
					Cell[elemCount] = PLUS_TINY;//0.0;
					elemCount++;
				}
			}
		}
	  private:
		void initBNum()
		{ 
			bNum = startNum + NM; 
			bNum2 = startNum + NMK - NM;
		}
		int bNum2;
	};
}
#endif //_LS_STAG_SUBDOMAIN2_MATRIX_H_
