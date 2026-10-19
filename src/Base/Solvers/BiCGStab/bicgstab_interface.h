/*
* \copyright (C) 2013 Puzikova V.V.
* \file      bicgstab_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      30.10.2013 7:56:21
* \brief     bicgstab_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BICGSTAB_INTERFACE_H_
#define _LS_STAG_BICGSTAB_INTERFACE_H_

// ================================================================== INCLUDES
#include "../iterative_solver_interface.h"
#include "../../MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BiCGStabInterface

Данный файл содержит в себе реализацию класса BiCGStabInterface
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(solUpd, LS_STAG_CR_PD, _s, LS_STAG_R_PD, _y, LS_STAG_R_PD, _z,
		LS_STAG_R_PD, _sol, LS_STAG_R_PD, _t, LS_STAG_R_PD, _r, const double, alpha, const double, w);
	LS_STAG_DECLARATION_KERNEL_R(double, pUpd, LS_STAG_R_PD, _p, LS_STAG_CR_PD, _nu_v, LS_STAG_CR_PD, _r,
		const double, _coef, const double, _beta);
/*!
\brief Класс, реализующий метод BiCGStab

Данный класс реализует итерационный метод BiCGStab (метод бисопряженных градиентов со стабилизацией) для решения СЛАУ 
(в т.ч. с использованием предобуславливания)
*/
	struct BiCGStabInterface : IterativeSolverInterface
	{
		LS_STAG_KERNEL(solUpd, LS_STAG_CR_PD _s, LS_STAG_R_PD _y, LS_STAG_R_PD _z, LS_STAG_R_PD _sol,
			LS_STAG_R_PD _t, LS_STAG_R_PD _r, const double alpha, const double w)
		{
			_y[i] *= alpha;
			_z[i] *= w;
			double curSum = _y[i] + _z[i];
			_sol[i] += curSum;
			_t[i] *= w;
			_r[i] = _s[i] - _t[i];
		}
		LS_STAG_KERNEL_R(double, pUpd, LS_STAG_R_PD _p, LS_STAG_CR_PD _nu_v, LS_STAG_CR_PD _r,
			const double _coef, const double _beta)
		{
			double curVal = (_p[i] - _nu_v[i] * _coef) * _beta;
			double curElem = _r[i];
			_p[i] = curElem + curVal;
			return curElem * curElem;
		}
		void set(SparseMatrix* _A, double* _rightM, double *_sol, BaseTaskDescription* _task)
		{
			IterativeSolverInterface::set(_A, _rightM, _sol, _task);
				nu_v.init(NM, m_nu_v);
				t.init(NM, m_t);
				y.init(NM, m_y);
				z.init(NM, m_z);
				r.init(NM, m_r);
				p.init(NM, m_p);
				s.init(NM, m_s);
				rr.init(NM);
				z.setMultipliers(m_A, m_sol);
				nu_v.setMultipliers(m_A, m_y);
				t.setMultipliers(m_A, m_z);
		};
		/*!
		Решает СЛАУ методом BiCGStab
		\return it число итераций
		*/
		int solve()
		{
			it = 1;
				initResidual();
				norm_r = r.getNorm();
				if(norm_r > eps) 
				{
					rr.copy(m_r);
					p.copy(m_r);
					do{
						preconditioner(m_p, m_y);
						nu_v.setProduct();
						double vs = rr * m_r;
						double coef = rr * m_nu_v;
						alpha = vs / coef;
						s.setDiffWithScale(m_r, m_nu_v, alpha);//nu_v *= alpha; s.setDiff(m_r, m_nu_v);
						preconditioner(m_s, m_z);
						t.setProduct();
						w = t * m_s;
						coef = t * m_t;
						w /= coef;
						solUpdate();
						coef = w / alpha;
						vs *= coef;
						double beta = rr * m_r;
						beta /= vs;
						{
							LS_STAG_PROFILE("updP");
							LS_STAG_PAR_SUM(BiCGStabInterface, pUpd, double, NM, m_p, m_nu_v, m_r, coef, beta)
							norm_r = sqrt(LS_STAG_PAR_ACC_VALUE);
						}
						it++;
					}while((norm_r > eps) && (it < max_it));
				}
			return it;
		};
		void print(const std::string &name)
		{
		    RealValuedMeshFunction w;
			w.init(NM);
			w.setOperator(m_A);
			double normAux = w.setProductWithResidual(m_sol, m_rightM);
			w.print(name+"_res");
			w.copy(m_sol);
			w.print(name+"_x");
			w.copy(m_rightM);
			w.print(name+"_b");
		};
	  protected:
		virtual void initResidual()// будет перегружена для AS_BiCGStab
		{
			z.setProduct();
			r.setDiff(m_rightM, m_z);
		};
	  private:
		void solUpdate()
		{
			LS_STAG_PROFILE("solUpd");
			LS_STAG_PAR_FOR(BiCGStabInterface, solUpd, NM, m_s, m_y, m_z, m_sol, m_t, m_r, alpha, w)
		};
		RealValuedMeshFunction z, nu_v, t;
		extDouble r, rr, s, p, y;
		double alpha, w;
		double *m_p, *m_y, *m_s, *m_z, *m_r, *m_nu_v, *m_t; 
	};
}

#endif //_LS_STAG_BICGSTAB_INTERFACE_H_
