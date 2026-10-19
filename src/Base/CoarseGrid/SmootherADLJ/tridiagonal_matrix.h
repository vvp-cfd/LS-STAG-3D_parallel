/*
* \copyright (C) 2013 Puzikova V.V.
* \file      tridiagonal_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      04.10.2013 14:44:42
* \brief     tridiagonal_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TRIDIAGONAL_MATRIX_H_
#define _LS_STAG_TRIDIAGONAL_MATRIX_H_

// ================================================================== INCLUDES
#include "../../MeshFunctions/mesh_function.h"
#include "../../ls_stag_real_constants.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TriDiagonalMatrix

Данный файл содержит в себе реализацию класса TriDiagonalMatrix
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_S_KERNEL(triLU, LS_STAG_R_PD, d, LS_STAG_R_PD, c, LS_STAG_R_PD, e, const int, NM);
/*!
\brief Класс, реализующий хранение и разложение трехдиагональной матрицы

Данный класс реализует хранение трехдиагональной матрицы, а также содержит метод для построение ее LU-разложения
*/
	struct TriDiagonalMatrix
	{
		TriDiagonalMatrix():
		m_N(0),
		m_M(0),
		NM(0),
		m_K(0),
		NMK(0),
		NMK2(0)
		{
			;
		}
		LS_STAG_SINGLE_KERNEL(triLU, LS_STAG_R_PD d, LS_STAG_R_PD c, LS_STAG_R_PD e, const int NM)
		{
			if (NM)
			{// 25/07/15
				double maxAbsA = 0.0;
				for (int i = 0; i < NM; i++)
				{
					double md = fabs(d[i]);
					if (md > maxAbsA) maxAbsA = md;
				}
				double tiny = EQ_TINY2_DOUBLE* maxAbsA;
				if(tiny < EQ_TINY_DOUBLE) tiny = EQ_TINY_DOUBLE;
				const double absTiny = 1.0e-12;// 01.10.2026: линия целиком твёрдая (диагонали ~PLUS_TINY) — сингулярна, d=0
				double aux = d[0];
				d[0] =  (fabs(aux) < tiny) ? 0.0 : (1.0 / aux);// Дирихле-строка (PLUS_TINY): развязываем линию, p=0 в твёрдой ячейке
				for (int i = 1, j = 0; i < NM; i++, j++)
				{
					e[j] *= d[j];// 25/07/15 /= d[j];
					double aux = d[i] - (c[j] * e[j]);// 25/07/15 d[i] -= c[j] * e[j];
					d[i] =  (fabs(aux) < tiny) ? 0.0 : (1.0 / aux);// Дирихле-строка (PLUS_TINY): развязываем линию, p=0 в твёрдой ячейке
				}
			}
		}
		virtual void init(const int _N, const int _M, const int _K)
		{
			m_N = _N;
			m_M = _M;
			m_K = _K;
			NM = m_N * m_M;
			NMK = NM * m_K;
			NMK2 = NMK - 1;
			LDiagonal.init(NMK2, c);
			Diagonal.init(NMK, d);
			UDiagonal.init(NMK2, e);
		};
		/*!
		Строит LU-разложение трехдиагональной матрицы (хранится в самой матрице)
		*/
		void buildLU()
		{
			LS_STAG_PROFILE("buildTriLU");
			LS_STAG_SINGLE(TriDiagonalMatrix, triLU, d, c, e, NMK)
		};
		void setNull()
		{
			LDiagonal.set(0.0);
			Diagonal.set(0.0);
			UDiagonal.set(0.0);
		};
		void print(const std::string &name)
		{
			LDiagonal.print(name+"_L");
			Diagonal.print(name+"_D");
			UDiagonal.print(name+"_U");
		};
		extDouble LDiagonal, Diagonal, UDiagonal;///< Сеточные функции, соответствующие нижней, главной и верхней диагоналям матрицы
		double* c;// L diagonal
		double* d;// diagonal
		double* e;// U diagonal
		int m_N;///< Число узлов сетки, на которой построена матрица, по горизонтали
		int m_M;///< Число узлов сетки, на которой построена матрица, по вертикали
		int m_K;///< Число узлов сетки, на которой построена матрица, по z
		int NMK;///< Длина главной диагонали
		int NMK2, NM;
	};
}

#endif //_LS_STAG_TRIDIAGONAL_MATRIX_H_
