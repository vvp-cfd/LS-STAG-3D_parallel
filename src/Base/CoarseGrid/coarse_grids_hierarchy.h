/*
* \copyright (C) 2013 Puzikova V.V.
* \file      coarse_grids_hierarchy.h
* \author    Puzikova Valeria (vvp) 
* \date      29.10.2013 12:04:31
* \brief     coarse_grids_hierarchy.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_COARSE_GRIDS_HIERARCHY_H_
#define _LS_STAG_COARSE_GRIDS_HIERARCHY_H_

// ================================================================== INCLUDES
#include "coarse_grid.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса CoarseGridsHierarchy

Данный файл содержит в себе реализацию класса CoarseGridsHierarchy
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с иерархией сеточных уровней

Данный класс реализует построение, хранение и перестроение иерархии сеточных уровней для данной сетки
*/
	struct CoarseGridsHierarchy
	{
		CoarseGridsHierarchy():
		numMG(0),
		numLevMG(0)
		{
			;
		}
		~CoarseGridsHierarchy()
		{
			if (numMG > 0)
			{
				delete[] grids;
				numMG = 0;
			}
		};
		void set(SparseMatrix* _A)
		{
			m_A = _A;
			GridLevels();
			grids = new CoarseGrid[numMG];
			numLevMG = numMG;// default
			lastLevel = numLevMG - 1;
			grids[0].set(m_A);
			for(int i = 1; i < numMG; i++)
				grids[i].set((SparseMatrix*)&grids[i-1].AH);
		};
		virtual void rebuild()// 02.12.14
		{
			LS_STAG_PROFILE("rebuildCoarseGrids");
			for(int i = 0; i < numMG; i++)
				grids[i].rebuild();
		};
		void print(const std::string& name)
		{
			for (int i = 0; i < numLevMG; i++)
				grids[i].AH.print(name + std::to_string(i));
		}
	  private:
		  /*!
		  Определяет максимальное число сеточных уровней для данной сетки
		  */
		void GridLevels()
		{
			double aux = min(min(m_A->m_N, m_A->m_M), m_A->m_K) * 1.0;
			while(aux > 2.0)
			{
				aux *= 0.5;
				numMG++;
			}
		};
	  protected:
		SparseMatrix* m_A;///< Указатель на матрицу системы на самой подробной сетке
	  public:
		CoarseGrid* grids;///< Массив сеточных уровней
		int numMG;///< Число сеточных уровней в иерархии
		int numLevMG;///< Сколько сеточных уровней из иерархии используем
		int lastLevel;///< Номер последнего уровня
	};
};

#endif //_LS_STAG_COARSE_GRIDS_HIERARCHY_H_
