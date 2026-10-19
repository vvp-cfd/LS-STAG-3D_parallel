/*
* \copyright (C) 2015 Puzikova V.V.
* \file      fgmres_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      08.08.2015 3:03:16
* \brief     fgmres_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_FGMRES_INTERFACE_H_
#define _LS_STAG_FGMRES_INTERFACE_H_

// ================================================================== INCLUDES
#include "../iterative_solver_interface.h"
#include "../../MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса FGMRESInterface

Данный файл содержит в себе реализацию класса FGMRESInterface
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(upd, LS_STAG_CR_PD, _curZ, LS_STAG_R_PD, _sol, const double, _curY);
/*!
\brief Класс, реализующий метод FGMRES

Данный класс реализует итерационный метод FGMRES (метод обобщенных минимальных невязок) для решения СЛАУ 
(в т.ч. с использованием предобуславливания)
*/
	struct FGMRESInterface : IterativeSolverInterface
	{
		FGMRESInterface():
		IterativeSolverInterface(),
		restartNum(0),
		max_size(0),
			m_size(0)
		{
			;
		}
		virtual ~FGMRESInterface()
		{
			if(max_size)
			{
				delete [] c;
				delete [] s;
				delete [] y;
				delete [] H;
				max_size = 0;
			}
		};
		LS_STAG_KERNEL(upd, LS_STAG_CR_PD _curZ, LS_STAG_R_PD _sol, const double _curY)
		{
			_sol[i] += (_curZ[i] * _curY);
		}
		void set(SparseMatrix* _A, double* _rightM, double *_sol, BaseTaskDescription* _task)
		{
			IterativeSolverInterface::set(_A, _rightM, _sol, _task);
			restartNum = _task->FGMRES_RESTART_NUM;
			max_size = restartNum * NM;
			m_v.init(max_size, v);
			m_z.init(max_size, z);
			c = new double[restartNum];
			s = new double[restartNum];
			y = new double[restartNum];
			m_size = restartNum * (restartNum + 1) / 2;
			H = new double[m_size];
			w.init(NM);
			w.setOperator(m_A);
		};
		/*!
		Решает СЛАУ методом FGMRES
		\return it число итераций
		*/
		int solve()
		{
			it = 0;
				norm_r = w.setProductWithResidual(m_sol, m_rightM);
				if(norm_r > eps)
				{
					it = 1;
					double normW = norm_r;
					double curH = norm_r;
					for( ; ; )//  рестарты
					{
						double* curV = v;
						double* curZ = z;
						for(int j = 0, hColStart = 0, hNum = 0; j < restartNum; j++, hColStart = hNum, it++, curV += NM, curZ += NM)
						{
							w.copyTo(curV, normW);
							preconditioner(curV, curZ);
							w.setProduct(curZ);
							double* curV2 = v;
							for(int i = 0, m_end = j + 1; i < m_end; i++, curV2 += NM, hNum++)
							{
								double curH = w * curV2;
								w.subtract(curV2, curH);
								H[hNum] = curH;
							}
							curH = H[hColStart];
							for(int i = 0, curHNum = hColStart, nextHNum = curHNum + 1; i < j; i++, curHNum = nextHNum, nextHNum++)
							{// пересчитываем предыдущие строки в данном столбце с учетом уже найденных косинусов и синусов
								double curC = c[i];
								double curS = s[i];
								double nextH = H[nextHNum];
								H[curHNum] = curC * curH - curS * nextH;
								double newNextH = curS * curH + curC * nextH;// H[nextHNum] = newNextH;
								curH = newNextH;
							} 
							normW = w.getNorm();// if(normW < PLUS_TINY) то norm_r точно окажется < eps
							double curC = curH * curH + normW * normW;
							double curA = sqrt(curC);
							curC = 1.0 / curA;
							double curS = -normW * curC;
							curC *= curH;
							H[hNum-1] = curA;
							y[j] = curC * norm_r;
							norm_r *= curS;// ! через 1 итерацию < 0
							curH = fabs(norm_r);
							if((curH > eps) && (it < max_it))
							{
								c[j] = curC;
								s[j] = curS;
							}
							else
							{// построение системы векторов закончено
								norm_r = curH;
								solUpdate(j, hNum - 1);
								return it;
							}
						}
						// естественный конец цикла без сходимости: фиксируем решение и перезапускаемся с истинной невязки
						solUpdate(restartNum - 1, m_size - 1); // по тождеству Saad ||b-Ax|| = |gamma| --- норму не пересчитываем
						w.setProductDiff(m_sol, m_rightM);// w = r --- направление 1-го базисного вектора след. цикла
						normW = curH;
					}
				}
				else return 0;
		}
		void print(const std::string &name)
		{
			double normAux = w.setProductWithResidual(m_sol, m_rightM);
			w.print(name+"_res");
			w.copy(m_sol);
			w.print(name+"_x");
			w.copy(m_rightM);
			w.print(name+"_b");
		};
	  private:
		void solUpdate(int lastBasisVectorNum, int lastHNum)
		{
			LS_STAG_PROFILE("solUpd");
			double* curZ = z + lastBasisVectorNum * NM;
			for(int k = lastBasisVectorNum, hNum = lastHNum; k > -1; k--, hNum--, curZ -= NM)
			{
				double curY = y[k] / H[hNum];
				for(int i = k - 1; i > -1; i--)
				{
					hNum--;
					double aux = curY * H[hNum];
					y[i] -= aux;
				}
				LS_STAG_PAR_FOR(FGMRESInterface, upd, NM, curZ, m_sol, curY)
			}
		};
		Array<double> m_v, m_z;
		RealValuedMeshFunction w;
		double *v, *z;// матрицы из векторов, хранятся по столбцам
		double *c, *s;///< векторы косинусов и синусов (синусы с обратным знаком)
		double *y;    ///< правая часть и решение задачи оптимизации
		double* H;///< верхний треугольник матрицы Хессенберга, хранится по столбцам
		int restartNum;///< через сколько итераций делать рестарт
		int max_size, m_size;
	};
};

#endif //_LS_STAG_FGMRES_INTERFACE_H_
