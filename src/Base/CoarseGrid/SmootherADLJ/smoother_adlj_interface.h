/*
* \copyright (C) 2013 Puzikova V.V.
* \file      smoother_adlj_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      17.10.2013 14:16:18
* \brief     smoother_adlj_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SMOOTHER_ADLJ_INTERFACE_H_
#define _LS_STAG_SMOOTHER_ADLJ_INTERFACE_H_

// ================================================================== INCLUDES
#include "../../ls_stag_real_constants.h"
#include "tridiagonal_matrix.h"
#include "../../MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SmootherADLJInterface

Данный файл содержит в себе реализацию класса SmootherADLJInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий работу c линией Якоби из метода ADLJ

Данный базовый класс реализует выполнение сглаживающей итерации с релаксацией по одной линии Якоби из метода ADLJ
*/
	struct SmootherADLJInterface : TriDiagonalMatrix
	{
		SmootherADLJInterface():
		omega(DEFAULT_OMEGA_3D)
		{
			;
		}
		/*!
		Решает СЛАУ с трехдиагональной матрицей методом прогонки 
		(т.е. при помощи алгоритма Томаса, TDMA --- Tridiagonal matrix algorithm ), при этом учитывает направление 
		(горизонтальное или вертикальное)
		*/
		virtual void TDMA(){};
		virtual void setSolver(SparseMatrix* _AH, extDouble *_r, extDouble *_c, RealValuedMeshFunction *_q, extDouble *_ch, const int _n_min)
		{// матрица со следующего грубого уровня !
			m_AH = _AH;
			m_r = _r;
			m_c = _c;
			m_q = _q;
			m_ch = _ch;
			N_MIN = _n_min;
			m_c->getCopy(solution);
			m_r->getCopy(right_member);
			m_q->getCopy(m_qV);
		};
		void preSolve()
		{// предсглаживающая итерация
			omegaTDMA(); 
			calcRes();
		}
		/*!
		Выполняет сглаживающую итерацию с релаксацией по одной линии Якоби из метода ADLJ 
		*/
		void solve()
		{
			calcRes();
			omegaTDMA();
			m_ch->operator+=(solution);
		};
	private:
		void calcRes()
		{
			m_q->setProduct(m_AH);// setMultiplier в MG сделан
			m_r->operator-=(m_qV);
		}
		void omegaTDMA()
		{
			TDMA();
			m_c->operator*=(omega);
		}
	public:
		double omega;      ///< Параметр релаксации
		SparseMatrix* m_AH;///< Указатель на разреженную матрицу, по которой строился сглаживатель
		RealValuedMeshFunction *m_q;
		extDouble *m_r, *m_c, *m_ch;
		double *solution, *right_member;
	  protected:
		double* m_qV;
		int N_MIN;
	};
}

#endif //_LS_STAG_SMOOTHER_ADLJ_INTERFACE_H_
