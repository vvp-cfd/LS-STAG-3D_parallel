/*
* \copyright (C) 2019 Puzikova V.V.
* \file      interface_block_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      04.08.2019 22:24:45
* \brief     interface_block_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_INTERFACE_BLOCK_MATRIX_H_
#define _LS_STAG_INTERFACE_BLOCK_MATRIX_H_

// ================================================================== INCLUDES
#include "dense_matrix.h"
#include "../../../SparseMatrix/sparse_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса InterfaceBlockMatrix

Данный файл содержит в себе реализацию класса InterfaceBlockMatrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий интерфейсную матрицу (как блочную)

Данный класс реализует построение и хранение интерфейсной матрицы из метода декомпозиции области в формате блочной 
трехдиагональной матрицы
*/
	struct InterfaceBlockMatrix
	{
	  InterfaceBlockMatrix():
		iCount(0),
		iCount2(0)
		{
			;
		}
		~InterfaceBlockMatrix()
		{
			if(iCount > 0) delete[] D;
			if(iCount2 > 0)
			{
				delete[] L;
				delete[] U;
			}
		}
		void build(SparseMatrix* _A, const int _domainNM, const int _N, const int _m)
		{
			iCount = _m;
			iCount2 = _m - 1;
			m_N = _N;
			if(iCount) D = new DenseMatrix[iCount];
			if(iCount2) 
			{
				L = new DenseMatrix[iCount2];
				U = new DenseMatrix[iCount2];
			}
			for(int i = 0; i < iCount; i++)
				D[i].set(m_N);
			for(int i = 0; i < iCount2; i++)
			{
				L[i].set(m_N);
				U[i].set(m_N);
			}
			m_A_Num = _A->Number;
			m_A_Portrait = _A->Portrait;
			m_A_Cell = _A->Cell;
			NM = m_N * m_N -1;
			domainNM = _domainNM;
			rebuild();
		}
		/*!
		Заполняет три диагонали (нижнюю, главную и верхнюю) из блоков главной диагонали интерфейсной матрицы 
		по исходной матрице СЛАУ для полной расчетной области (остальные элементы заполняются в решателях) и инициализирует 
		нулями блоки из нижней и верхней поддиагоналей 
		*/
		void rebuild()
		{// заполняем только 5 диагоналей (s,w,p,e,n), остальное --- из решателей идет
			for(int i = 0; i < iCount2; i++)
			{
				L[i].init();
				U[i].init();
			}
			for (int i = 0, g = domainNM, startNum = g; i < iCount; i++, g += domainNM, startNum = g)
        {
            D[i].init();
            double* curD = D[i].Cell;
            // обходим все строки плоскости
            for (int j = 0, num = m_A_Num[g], start = num + 1, num2 = m_A_Num[g+1], end = num2 - 1; j < NM; j += m_N, g++, 
				num = num2, num2 = m_A_Num[g+1], start = num + 1, end = num2 - 1)
            {
				if(start == end) start--;
                for (int k = start; k < end; k++) {
                    int locCol = m_A_Portrait[k] - startNum;
                    curD[j + locCol] = m_A_Cell[k];
                }
			}
            }
        }
	  void print(const std::string &name)
	  {
			std::string fileNameL = name + "_L_";
			std::string fileNameD = name + "_D_";
			std::string fileNameU = name + "_U_";
			for(int i = 0; i < iCount; i++)
				D[i].print(fileNameD + std::to_string(i));
			for(int i = 0; i < iCount2; i++)
			{
				L[i].print(fileNameL + std::to_string(i));
				U[i].print(fileNameU + std::to_string(i));
			}
	  };
		int* m_A_Num;
		int* m_A_Portrait;
		double* m_A_Cell;
	  DenseMatrix *L, *D, *U;///< Массивы блоков, из которых состоит интерфейсная матрица
	  int iCount;///< Число блоков на главной диагонали
	  int iCount2;///< Число блоков на верхней и нижней поддиагоналях
	  int m_N, domainNM, NM;
	};
}
#endif //_LS_STAG_INTERFACE_BLOCK_MATRIX_H_