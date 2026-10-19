/*
* \copyright (C) 2019 Puzikova V.V.
* \file      shur_subdomain_2.h
* \author    Puzikova Valeria (vvp) 
* \date      07.04.2019 23:43:15
* \brief     shur_subdomain_2.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SHUR_SUBDOMAIN_2_H_
#define _LS_STAG_SHUR_SUBDOMAIN_2_H_

// ================================================================== INCLUDES
#include "shur_subdomain.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса ShurSubdomain2

Данный файл содержит в себе реализацию класса ShurSubdomain2
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий задачу в промежуточной подобласти

Данный класс реализует хранение задачи в промежуточной подобласти (имеющей две соседние)
и работу по пересылкам правой части, решения, начальным приближениям 
между подобластью, интерфейсами и глобальной областью
*/
	struct ShurSubdomain2 : ShurSubdomain
	{
		ShurSubdomain2():
		ShurSubdomain(),
		Aps2(NULL),
		Asp2(NULL)
		{
			;
		}
		~ShurSubdomain2()
		{
			if(Aps2) 
			{
				delete[] Aps2;
				delete[] Asp2;
				Aps2 = NULL;
			}
		}
		/*!
		Корректирует правую часть в подобласти с учетом рассчитанных интерфейсных неизвестных
		*/
		void recalcRHS()
		{// пересчет правой части в подобласти с учетом рассчитанных интерфейсных неизвестных
			LS_STAG_PROFILE("recalcRHS");
			for(int i = 0, j = startAsp; i < N; i++, j++)
			{
				right[i] -= (Aps[i] * m_xsp[i]);
				right[j] -= (Aps2[i] * m_xsp2[i]);
			}
		}
		/*!
		Корректирует правую часть на интерфейсе по решению в подобласти с исходной правой частью
		*/
		void recalcInterfaceBHS()
		{// пересчет правой части на интерфейсе по решению в подобласти с исходной правой частью
			LS_STAG_PROFILE("recalcInterfaceBHS");
			for(int i = 0, j = startAsp; i < N; i++, j++)
			{
				m_bsp[i] -= Asp[i] * sol[i];
				m_bsp2[i] -= Asp2[i] * sol[j];
			}
		}
		/*!
		Корректирует правую часть для расчета интерфейсной матрицы на северной границе
		*/
		void setNorthRHS(const int colNum)
		{// правая часть для расчета интерфейсной матрицы на южной границе
			LS_STAG_PROFILE("setNorthRHS");
			setBoundaryRHS(colNum);
			right[startAsp + colNum] = Aps2[colNum];//2
		}
		void getSouthInterface()//10/08/19
		{
			LS_STAG_PROFILE("getSouthInterface");
			for(int i = 0, j = startAsp, cur_css = m_colNum; i < N; i++, j++, cur_css += N)
			{
				m_AsS_Cell[cur_css] -= (Asp[i] * sol[i]);//N
				m_AsS2_Cell[cur_css] -= (Asp2[i] * sol[j]);//S
			}
		}
		void getNorthInterface()//10/08/19
		{
			LS_STAG_PROFILE("getNorthInterface");
			for(int i = 0, j = startAsp, cur_css = m_colNum; i < N; i++, j++, cur_css += N)
			{
				m_AsN_Cell[cur_css] -= (Asp[i] * sol[i]);//N
				m_AsN2_Cell[cur_css] -= (Asp2[i] * sol[j]);//S
			}
		}
	  private:
		SubdomainMatrix* getSubdomainMatrix() { return ((SubdomainMatrix*)(new Subdomain2Matrix)); }
		/*!
		Строит и перестраивает матрицы \f$A_{p,s}\f$ и \f$A_{s,p}\f$ 
		*/
		void rebuildSP()
		{
			LS_STAG_PROFILE("rebuildSP");
			for(int i = 0, j = startAspG2 + 1, num = m_A_Num[startAspG2], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
			{// Asp --- Т c предыдущего граничного ряда
				if(num2 > (num + 1))
					Asp[i] = m_A_Cell[num2-1];// T всегда последняя
				else 
					Asp[i] = 0.0;// граница может попасть на твердые ячейки
			}
			for(int i = 0, j = startNum + 1, num = m_A_Num[startNum], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
			{// Aps --- B c первого ряда в подобласти
				if(num2 > (num + 1))
					Aps[i] = m_A_Cell[num];// B всегда первая
				else 
					Aps[i] = 0.0;// граница может попасть на твердые ячейки
			}
			for(int i = 0, j = startAspG + 1, num = m_A_Num[startAspG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
			{// Asp --- B cо следующего ряда
				if(num2 > (num + 1))
					Asp2[i] = m_A_Cell[num];// B всегда первая
				else 
					Asp2[i] = 0.0;// граница может попасть на твердые ячейки
			}
			for(int i = 0, j = startApsG + 1, num = m_A_Num[startApsG], num2 = m_A_Num[j]; i < N; i++, j++, num = num2, num2 = m_A_Num[j])
			{// Aps --- T c последнего ряда в подобласти
				if(num2 > (num + 1))
					Aps2[i] = m_A_Cell[num2-1];// T всегда последняя
				else 
					Aps2[i] = 0.0;// граница может попасть на твердые ячейки
			}
		}
		void createSP()
		{
			ShurSubdomain::createSP();
			Aps2 = new double[N];
			Asp2 = new double[N];
		}
		void initBsp(double* bs, double* xs, InterfaceBlockMatrix* _Ass)// 10/08/19
		{
			int m_num2 = m_num - 1;
			num_cssS = m_num2 * N;
			m_AsS_Cell = _Ass->D[m_num2].Cell;// 10/08/19
			m_AsN_Cell = _Ass->U[m_num2].Cell;// 10/08/19
			m_AsS2_Cell = _Ass->L[m_num2].Cell;// 10/08/19
			m_AsN2_Cell = _Ass->D[m_num].Cell;// 10/08/19
			NM = N * m;
			m_bsp = bs + num_cssS;// S граница 
			m_xsp = xs + num_cssS;// S граница 
			m_bsp2 = m_bsp + N;// N граница 
			m_xsp2 = m_xsp + N;// N граница 
			startAsp = NM - N;// последний блок
			startApsG = startNum + startAsp;
			startAspG = startApsG + N;
			startAspG2 = startNum - N;
		}
		double *m_AsN_Cell, *m_AsN2_Cell, *m_AsS2_Cell;// 10/08/19 
		int num_cssS; //10/08/19 , num_cssN;
		int startAspG2;// наоборот, для S границы
	  public:
		double *Aps2, *Asp2;///< матрицы Aps и Asp в данном случае просто вектор --- N граница
		double *m_bsp2, *m_xsp2;
	};
};

#endif //_LS_STAG_SHUR_SUBDOMAIN_2_H_
