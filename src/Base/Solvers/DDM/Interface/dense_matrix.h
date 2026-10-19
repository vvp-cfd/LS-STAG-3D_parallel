/*
* \copyright (C) 2019 Puzikova V.V.
* \file      dense_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      03.08.2019 06:15:45
* \brief     dense_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_DENSE_MATRIX_H_
#define _LS_STAG_DENSE_MATRIX_H_

// ================================================================== INCLUDES
#include <fstream>
#include <string>
#include "../../../profiler.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса DenseMatrix

Данный файл содержит в себе реализацию класса DenseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с заполненной матрицей

Данный класс реализует хранение заполненной матрицы и работу с ней, в т.ч. матрично-векторное произведение
*/
	struct DenseMatrix
	{
	  DenseMatrix():
		m_N(0),
		elemCount(0)
		{
			;
		}
		~DenseMatrix()
		{
			if(elemCount > 0) delete[] Cell;
		}
	  void init()
	  {// + инициализируем нулями
		  memset(Cell, 0, capacity);
	  }
	  void set(const int _N)
	  {
		  m_N = _N;
		  elemCount = m_N * m_N;
		  Cell = new double[elemCount];
		  capacity = elemCount * sizeof(double);
	  };
	  void setMultProd(double *_mult, double *_prod)
	  {
		  m_mult = _mult;
		  m_prod = _prod;
	  }
	  /*!
	  Вычисляет матрично-векторное произведение
	  */
	  void product()
	  {
		  LS_STAG_PROFILE("GEMM");
		  for(int i = 0, num = 0; i < m_N; i++)
		  {
			  double aux = 0.0;
			  for(int j = 0; j < m_N; j++, num++)
				  aux += Cell[num] * m_mult[j];
			  m_prod[i] = aux;
		  }
	  }
	  void print(const std::string &name)
	  {
		    std::ofstream f;
			std::string fileName = name + ".dat";
			f.open(fileName.c_str());
			for(int i = 0, num = 0; i < m_N; i++)
			{
				for(int j = 0; j < m_N; j++, num++)
					f << Cell[num] << " ";
				f << std::endl;
			}
			f.close();
	  };
	  double *m_mult, *m_prod;
	  double* Cell; ///< элементы матрицы
	  int m_N, elemCount, capacity;
	};
}
#endif //_LS_STAG_DENSE_MATRIX_H_
