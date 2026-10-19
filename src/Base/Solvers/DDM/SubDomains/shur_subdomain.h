/*
* \copyright (C) 2019 Puzikova V.V.
* \file      shur_subdomain.h
* \author    Puzikova Valeria (vvp) 
* \date      07.04.2019 21:19:27
* \brief     shur_subdomain.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SHUR_SUBDOMAIN_H_
#define _LS_STAG_SHUR_SUBDOMAIN_H_

// ================================================================== INCLUDES
#include "subdomain2_matrix.h"
#include "../Interface/interface_block_matrix.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ShurSubdomain

Данный файл содержит в себе реализацию класса ShurSubdomain
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий задачу в подобласти

Данный класс реализует хранение задачи в подобласти (первой или последней)
и работу по пересылкам правой части, решения, начальным приближениям 
между подобластью, интерфейсами и глобальной областью
*/
	struct ShurSubdomain
	{
		ShurSubdomain():
		N(0),
		m(0),
		Aps(NULL),
		Asp(NULL),
		App(NULL)
		{
			;
		}
		virtual ~ShurSubdomain()
		{
			if(Aps) 
			{
				delete[] Aps;
				delete[] Asp;
				delete[] sol;
				delete[] right;
				Aps = NULL;
			}
			if(App) delete App;
		}
		void set(SparseMatrix* _A, InterfaceBlockMatrix* _Ass,// 10/08/19 SparseMatrix* _Ass, 
			  double* x, double* b, double* bs, double* xs, const int _num, const int _m)
		{
			m_num = _num;// номер подобласти
			N = _A->m_N * _A->m_M;
			m_A_Num = _A->Number;
			m_A_Cell = _A->Cell;
			startNum = _m * N * m_num;
			m = _m - 1;// 1 на интерфейс
			m_A = _A;
			initBsp(bs, xs, _Ass);// 10/08/19 
			sol = new double[NM];
			right = new double[NM];
			m_sol = x + startNum;// поскольку на подобласти режем по горизонтали
			m_right = b + startNum;
			build();
		}
		/*!
		Строит матрицы (СЛАУ и вспомогательные) в подобласти
		*/
		void build()
		{// построение матриц
			createSP();
			App = getSubdomainMatrix();
			App->build(m_A, startNum, N, m);
			rebuildSP(); // Asp, Aps
		}
		/*!
		Перестраивает матрицы (СЛАУ и вспомогательные) в подобласти
		*/
		void rebuild()
		{// перестроение матриц
			LS_STAG_PROFILE("rebuildSubdomain");
			App->rebuild();
			rebuildSP(); // Asp, Aps
		}
		/*!
		Копирует из глобальной области в данную подобласти векторы правой части и начального приближения
		*/
		void initSubDomain()
		{
			LS_STAG_PROFILE("initSubdomain");
			for(int i = 0; i < NM; i++)
			{
				right[i] = m_right[i];
				sol[i] = m_sol[i];
			}
		}
		/*!
		Копирует из данной подобласти в глобальную область вектор решения
		*/
		void sendSubDomain()
		{
			LS_STAG_PROFILE("sendSubdomain");
			for(int i = 0; i < NM; i++)
				m_sol[i] = sol[i];
		}
		/*!
		Корректирует правую часть в подобласти с учетом рассчитанных интерфейсных неизвестных
		*/
		virtual void recalcRHS()
		{// пересчет правой части в подобласти с учетом рассчитанных интерфейсных неизвестных
			LS_STAG_PROFILE("recalcRHS");
			for(int i = 0, j = startAsp; i < N; i++, j++)
				right[j] -= Aps[i] * m_xsp[i];
		}
		/*!
		Корректирует правую часть на интерфейсе по решению в подобласти с исходной правой частью
		*/
		virtual void recalcInterfaceBHS()
		{// пересчет правой части на интерфейсе по решению в подобласти с исходной правой частью
			LS_STAG_PROFILE("recalcInterfaceBHS");
			for(int i = 0, j = startAsp; i < N; i++, j++)
				m_bsp[i] -= Asp[i] * sol[j];
		}
		/*!
		Корректирует правую часть для расчета интерфейсной матрицы на южной границе
		*/
		virtual void setSouthRHS(const int colNum)
		{// правая часть для расчета интерфейсной матрицы на южной границе
			LS_STAG_PROFILE("setSouthRHS");
			setBoundaryRHS(colNum);
			right[colNum] = Aps[colNum];
		}
		/*!
		Корректирует правую часть для расчета интерфейсной матрицы на северной границе
		*/
		virtual void setNorthRHS(const int colNum)
		{// правая часть для расчета интерфейсной матрицы на северной границе
			LS_STAG_PROFILE("setNorthRHS");
			setBoundaryRHS(colNum);
			right[startAsp + colNum] = Aps[colNum];
		}
		void getInterface()//10/08/19
		{
			LS_STAG_PROFILE("getInterface");
			for(int i = 0, j = startAsp, cur_css = m_colNum; i < N; i++, j++, cur_css += N)
				m_AsS_Cell[cur_css] -= (Asp[i] * sol[j]);
		}
		virtual void getSouthInterface()//10/08/19
		{
			getInterface();
		}
		virtual void getNorthInterface()//10/08/19
		{
			getInterface();
		}
	  protected:
		virtual SubdomainMatrix* getSubdomainMatrix() { return (new SubdomainMatrix); }
		void setBoundaryRHS(const int colNum)
		{
			LS_STAG_PROFILE("setBoundaryRHS");
			for(int i = 0; i < NM; i++)
			{
				right[i] = 0.0;
				sol[i] = 0.0;
			}
			m_colNum = colNum;//10.08.19 m_colNum = num_css + colNum;
		}
		/*!
		Строит и перестраивает матрицы \f$A_{p,s}\f$ и \f$A_{s,p}\f$ 
		*/
		virtual void rebuildSP()
		{
			LS_STAG_PROFILE("rebuildSP");
			if(m_num)
			{// не первая подобласть
				for(int i = 0, j = startAspG + 1, num = m_A_Num[startAspG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
				{// Asp --- T c предыдущего граничного ряда
					if(num2 > (num + 1))
						Asp[i] = m_A_Cell[num2-1];// T всегда последняя
					else 
						Asp[i] = 0.0;// граница может попасть на твердые ячейки
				}
				for(int i = 0, j = startApsG + 1, num = m_A_Num[startApsG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
				{// Aps --- B c первого ряда в подобласти
					if(num2 > (num + 1))
						Aps[i] = m_A_Cell[num];// B всегда первая
					else 
						Aps[i] = 0.0;// граница может попасть на твердые ячейки
				}
			}
			else
		    {// первая подобласть
				for(int i = 0, j = startAspG + 1, num = m_A_Num[startAspG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
				{// Asp --- B cо следующего ряда
					if(num2 > (num + 1))
						Asp[i] = m_A_Cell[num];// B всегда первая
					else 
						Asp[i] = 0.0;// граница может попасть на твердые ячейки
				}
				for(int i = 0, j = startApsG + 1, num = m_A_Num[startApsG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
				{// Aps --- T c последнего ряда в подобласти
					if(num2 > (num + 1))
						Aps[i] = m_A_Cell[num2-1];// T всегда последняя
					else 
						Aps[i] = 0.0;// граница может попасть на твердые ячейки
				}
			}
		}
		virtual void createSP()
		{
			Aps = new double[N];
			Asp = new double[N];
		}
		virtual void initBsp(double* bs, double* xs, InterfaceBlockMatrix* _Ass)// 10/08/19
		{
			if(m_num)
			{
				num_css = (m_num - 1) * N;
				NM = m_A->NMK - startNum; // L * K_p
				m_AsS_Cell = _Ass->D[m_num - 1].Cell;// 10/08/19
				m_bsp = bs + num_css;// последняя подобласть
				m_xsp = xs + num_css;
				m = m_A->m_K - (m + 1) * m_num;// в последней подобласти число строк м.б. чуть больше
				startAsp = 0;// первый блок граничный
				startAspG = startNum - N;
				startApsG = startNum;
			}
			else
			{
				num_css = 0;
				if(_Ass) m_AsS_Cell = _Ass->D[0].Cell;// 10/08/19
				NM = N * m;
				m_bsp = bs;// первая подобласть
				m_xsp = xs;
				startAsp = NM - N;// последний блок граничный
				startApsG = startNum + startAsp;
				startAspG = startApsG + N;
			}
		}
		SparseMatrix* m_A;///< Указатель на матрицу СЛАУ в полной расчетной области
		                 // нужно для построения и перестроения задачи
		int* m_A_Num;
		double* m_A_Cell;
		double *m_AsS_Cell;// 10/08/19 int* m_Ass_Num;
		double *m_sol, *m_right, *m_bsp, *m_xsp;///< указатель на  векторы решения и правой части в подобласти 
		                        // (пока просто внутри глобального вектора)
		int startNum;///< номер ячейки в глобальной области, являющейся 0-й для подобласти
		int startAsp;///< номер внутри решения в подобласти, соответствующий \f$A_{s,p}\f$ --- актуально для N границы
		int startAspG;///< номер  ячейки в глобальной области, открывающей граничный ряд
		int startApsG;///< номер  ячейки в глобальной области, открывающей граничный ряд в подобласти
		int m_colNum;///< номер столбца, по которому считаем в данный момент столбец интерфейсной матрицы
		int num_css;
	  public:
		SubdomainMatrix* App;///< матрицы задачи в подобласти
		double *Aps, *Asp;///< матрицы \f$A_{p,s}\f$ и \f$A_{s,p}\f$ в данном случае просто вектор
		double *sol, *right;//собственные векторы решения и правой части
		int N, m;///< число ячеек по горизонтали и по вертикали N = L, m = K_p
		int NM;
		int m_num;///< номер подобласти 
	};
};

#endif //_LS_STAG_SHUR_SUBDOMAIN_H_
