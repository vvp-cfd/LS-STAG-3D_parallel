/*
* \copyright (C) 2019 Puzikova V.V.
* \file      subdomain_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      12.04.2019 21:00:39
* \brief     subdomain_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SUBDOMAIN_MATRIX_H_
#define _LS_STAG_SUBDOMAIN_MATRIX_H_

// ================================================================== INCLUDES
#include "../../../SparseMatrix/sparse_matrix.h"
#include "../../../profiler.h"
#include "../../../ls_stag_real_constants.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SubdomainMatrix

Данный файл содержит в себе реализацию класса SubdomainMatrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий матрицу СЛАУ в подобласти

Данный класс реализует построение по разреженной матрице СЛАУ в полной расчетной области построение матрицы СЛАУ в 
подобласти и ее хранение
*/
	struct SubdomainMatrix : SparseMatrix
	{
		void build(SparseMatrix* _A, const int _startNum, const int _N, const int _m)
		{
			NM = _N;
			set(_A->m_N, _A->m_M, _m);
			m_A_Num = _A->Number;
			m_A_Portrait = _A->Portrait;
			m_A_Cell = _A->Cell;
			startNum = _startNum;
			nNum = NMK - _N;
			initBNum();
			rebuild();
		}
		/*!
		Строит портрет матрицы в подобласти, у которой только одна соседняя, и заполняет ее по исходной матрице 
		СЛАУ для полной расчетной области
		*/
		virtual void rebuild()
		{
			LS_STAG_PROFILE("rebuildSubdom");
			elemCount = 0;
			if(startNum)
			{// изменяем только первую плоскость
				for(int i = 0, j = startNum + 1, num = m_A_Num[startNum], num2 = m_A_Num[j]; i < NM; i++, j++, num = num2, num2 = m_A_Num[j], Number[i] = elemCount)
				{// 1-я плоскость без B
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
				for(int i = NM, j = bNum + 1, num = m_A_Num[bNum], num2 = m_A_Num[j]; i < NMK; i++, j++, Number[i] = elemCount, num = num2, num2 = m_A_Num[j])
				{// берем матрицу как есть
					for(int k = num; k < num2; k++, elemCount++)
					{
						Portrait[elemCount] = m_A_Portrait[k] - startNum;
						Cell[elemCount] = m_A_Cell[k];
					}
				}
			}
			else
			{// изменяем только последний ряд
				for(int i = 0, j = startNum + 1, num = m_A_Num[startNum], num2 = m_A_Num[j]; i < nNum; i++, j++, Number[i] = elemCount, num = num2, num2 = m_A_Num[j])
				{// берем матрицу как есть
					for(int k = num; k < num2; k++, elemCount++)
					{
						Portrait[elemCount] = m_A_Portrait[k] - startNum;
						Cell[elemCount] = m_A_Cell[k];
					}
				}
				for(int i = nNum, j = bNum + 1, num = m_A_Num[bNum], num2 = m_A_Num[j]; i < NMK; i++, j++, num = num2, num2 = m_A_Num[j], Number[i] = elemCount)
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
		}
	  private:
		  /*!
		  Оценивает максимальное количество ненулевых элементов матрицы
		  */
		void estimateNeedMemory()
		{
			maxElemCount = 7 * NMK - 2 * NM;// max кол-во эл-в м-цы (если нет твёрдых ячеек)
		};
	  protected:
		virtual void initBNum(){ bNum = (startNum) ? (startNum + NM) : (startNum + NMK - NM); }
		int* m_A_Num;
		int* m_A_Portrait;
		double* m_A_Cell;
		int startNum;///< номер ячейки в глобальной области, являющейся 0-й для подобласти
		int bNum, nNum, NM;
	};
}
#endif //_LS_STAG_SUBDOMAIN_MATRIX_H_
