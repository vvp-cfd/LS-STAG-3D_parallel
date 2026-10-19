/*
* \copyright (C) 2019 Puzikova V.V.
* \file      shur_interface_task.h
* \author    Puzikova Valeria (vvp) 
* \date      07.04.2019 23:55:03
* \brief     shur_interface_task.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SHUR_INTERFACE_TASK_H_
#define _LS_STAG_SHUR_INTERFACE_TASK_H_

// ================================================================== INCLUDES
#include "../Interface/interface_block_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ShurInterfaceTask

Данный файл содержит в себе реализацию класса ShurInterfaceTask
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий интерфейсную задачу

Данный класс реализует хранение интерфейсной задачи и работу по пересылкам правой части, решения, начальным приближениям 
между интерфейсами и глобальной областью
*/
	struct ShurInterfaceTask
	{
		ShurInterfaceTask():
		N(0),
		domainNM(0),
		p(0)
		{
			;
		}
		~ShurInterfaceTask()
		{
			if(N) 
			{
				delete[] sol;
				delete[] right;
				N = 0;
			}
		}
		void set(SparseMatrix* _A, double* x, double* b, const int _m, const int _p)
		{// _m --- подобласть + граница, _p --- число подобластей
			N = _A->m_N * _A->m_M; // L
			domainNM = N * (_m - 1);
			p = _p - 1;// границ на 1 меньше чем подобластей
			NM = N * p;
			m_A = _A;
			m_solG = x;
			m_rightG = b;
			sol = new double[NM];
			right = new double[NM];
			initInterface();
			build();
		}
		void build()
		{// построение матриц
			Ass.build(m_A, domainNM, N, p);
		}
		void rebuild()
		{// перестроение матриц
			Ass.rebuild();// ! еще заполнить блоки с помощью решателей
		}
		/*!
		Копирует из глобальной области в интерфейсы векторы правой части и начального приближения
		*/
		void initInterface()
		{// копирование решения и правой части из глобального вектора
			LS_STAG_PROFILE("initInterface");
			for(int i = 0, num = 0, numG = domainNM; i < p; i++, numG += domainNM)
			{
				for(int j = 0; j < N; j++, num++, numG++)
				{
					sol[num] = m_solG[numG];
					right[num] = m_rightG[numG];
				}
			}
		}
		/*!
		Копирует из интерфейсов в глобальную область вектор решения
		*/
		void sendInterface()
		{// запись решения в глобальный вектор
			LS_STAG_PROFILE("sendInterface");
			for(int i = 0, num = 0, numG = domainNM; i < p; i++, numG += domainNM)
			{
				for(int j = 0; j < N; j++, num++, numG++)
				{
					m_solG[numG] = sol[num];// m_rightG[numG] = right[num];
				}
			}
		}
		SparseMatrix* m_A;///< Указатель на матрицу СЛАУ в полной расчетной области
		                 // нужно для построения и перестроения задачи
		double *m_solG, *m_rightG;///< указатели на глобальные векторы решения и правой части
	  public:
		InterfaceBlockMatrix Ass;///< матрица интерфейсной задачи
		double *sol, *right;///< векторы решения и правой части для интерфейсной подзадачи
		int N, p;///<  число ячеек по горизонтали и по вертикали (интерфейсные ряды только) 
		int domainNM;///<  число ячеек в подобласти 
		int NM;
	};
};

#endif //_LS_STAG_SHUR_INTERFACE_TASK_H_
