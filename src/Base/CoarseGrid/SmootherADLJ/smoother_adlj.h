/*
* \copyright (C) 2013 Puzikova V.V.
* \file      smoother_adlj.h
* \author    Puzikova Valeria (vvp) 
* \date      04.10.2013 16:16:48
* \brief     smoother_adlj.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SMOOTHER_ADLJ_H_
#define _LS_STAG_SMOOTHER_ADLJ_H_

// ================================================================== INCLUDES
#include "smoother_adlj_x.h"
#include "smoother_adlj_y.h"
#include "smoother_adlj_z.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SmootherADLJ

Данный файл содержит в себе реализацию класса SmootherADLJ
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(build, const int*, _Number, LS_STAG_CR_PI, _Portrait, LS_STAG_CR_PD, _Cell,
		LS_STAG_CR_PI, _ADLJ_yvo, LS_STAG_CR_PI, _ADLJ_zvo, LS_STAG_R_PD, _y_c, LS_STAG_R_PD, _x_c, LS_STAG_R_PD, _z_c, 
		LS_STAG_R_PD, _y_d, LS_STAG_R_PD, _x_d, LS_STAG_R_PD, _z_d,
		LS_STAG_R_PD, _y_e, LS_STAG_R_PD, _x_e, LS_STAG_R_PD, _z_e, const int, _N, const int, _NM);
/*!
\brief Класс, реализующий сглаживатель ADLJ

Данный класс реализует построение по разреженной матрице СЛАУ сглаживателя ADLJ --- 
Alternating Damped Line Jacobi (переменная релаксированная линия Якоби) и выполнение сглаживающей итерации с релаксацией
*/
struct SmootherADLJ
{
LS_STAG_KERNEL(build, const int* _Number, LS_STAG_CR_PI _Portrait, LS_STAG_CR_PD _Cell, LS_STAG_CR_PI _ADLJ_yvo, 
	LS_STAG_CR_PI _ADLJ_zvo, LS_STAG_R_PD _y_c, LS_STAG_R_PD _x_c, LS_STAG_R_PD _z_c, LS_STAG_R_PD _y_d, 
	LS_STAG_R_PD _x_d, LS_STAG_R_PD _z_d, LS_STAG_R_PD _y_e, LS_STAG_R_PD _x_e, LS_STAG_R_PD _z_e,
	const int _N, const int _NM)
		{
			for (int j = _Number[i], numA = _Number[i + 1], num = _ADLJ_yvo[i], numZ = _ADLJ_zvo[i],
				iw = i - 1, is = i - _N, ib = i - _NM, it = i + _NM, in = i + _N, ie = i + 1; j < numA; j++)
			{
				const int AP = _Portrait[j];
				if (AP == ib)
					_z_c[numZ - 1] = _Cell[j];
				else
				{
					if (AP == is)
						_y_c[num - 1] = _Cell[j];
					else
					{
						if (AP == iw)
							_x_c[iw] = _Cell[j];
						else
						{
							if (AP == i)
							{
								const double cell = _Cell[j];
								_z_d[numZ] = cell;
								_y_d[num] = cell;
								_x_d[i] = cell;
							}
							else
							{
								if (AP == ie)
									_x_e[i] = _Cell[j];
								else
								{
									if(AP == in)
										_y_e[num] = _Cell[j];
									else
									{
										if(AP == it) _z_e[numZ] = _Cell[j];
									}
								}
							}
						}
					}
				}
			}// for j
		}
    void set(SparseMatrix* _A) {
        m_A = _A;
        m_N = m_A->m_N;
        m_M = m_A->m_M;
        m_K = m_A->m_K;
		m_NM = m_N * m_M;
        m_NMK = m_A->NMK;
        m_Cell = m_A->Cell;
        m_Number = m_A->Number;
        m_Portrait = m_A->Portrait;
        ADLJ_x.init(m_N, m_M, m_K);          // M_x: w, p, e
        ADLJ_y.init(m_N, m_M, m_K);          // M_y: s, p, n
        ADLJ_z.init(m_N, m_M, m_K);          // M_z: b, p, t (в z-упорядочении)
        // инициализируем указатели на диагонали
        x_c = ADLJ_x.c; x_d = ADLJ_x.d; x_e = ADLJ_x.e;
        y_c = ADLJ_y.c; y_d = ADLJ_y.d; y_e = ADLJ_y.e;
        z_c = ADLJ_z.c; z_d = ADLJ_z.d; z_e = ADLJ_z.e;
        rebuild();
    }

    void setSolver(extDouble* _r, extDouble* _c, RealValuedMeshFunction* _q,
                   extDouble* _ch, const int _n_min) {
        ADLJ_x.setSolver(m_A, _r, _c, _q, _ch, _n_min);
        ADLJ_y.setSolver(m_A, _r, _c, _q, _ch, _n_min);
        ADLJ_z.setSolver(m_A, _r, _c, _q, _ch, _n_min);
    }

	void preSolve() {
		LS_STAG_PROFILE("preADLJ");
		ADLJ_z.preSolve();
		ADLJ_y.preSolve();
		ADLJ_x.preSolve();
	}

    void solve() {
        LS_STAG_PROFILE("ADLJ");
        ADLJ_z.solve();
        ADLJ_y.solve();
        ADLJ_x.solve();
    }

    void rebuild() {
        LS_STAG_PROFILE("rebuildADLJ");
        ADLJ_x.setNull(); ADLJ_y.setNull(); ADLJ_z.setNull();
        {
            LS_STAG_PROFILE("buildADLJ");
            LS_STAG_PAR_FOR(SmootherADLJ, build, m_NMK, m_Number, m_Portrait, m_Cell,
                ADLJ_y.vo, ADLJ_z.vo, y_c, x_c, z_c, y_d, x_d, z_d, y_e, x_e, z_e, m_N, m_NM)
        }
        ADLJ_x.buildLU();
        ADLJ_y.buildLU();
        ADLJ_z.buildLU();
    }

    void print(const std::string& name) {
        ADLJ_x.print(name + "_ADLJ_x");
        ADLJ_y.print(name + "_ADLJ_y");
        ADLJ_z.print(name + "_ADLJ_z");
    }

    SmootherADLJx ADLJ_x;
    SmootherADLJy ADLJ_y;
    SmootherADLJz ADLJ_z;
    SparseMatrix* m_A;
private:
    double *m_Cell, *x_c, *x_d, *x_e, *y_c, *y_d, *y_e, *z_c, *z_d, *z_e;
    int *m_Number, *m_Portrait;
    int m_N, m_M, m_K, m_NMK, m_NM;
};
}

#endif //_LS_STAG_SMOOTHER_ADLJ_H_
