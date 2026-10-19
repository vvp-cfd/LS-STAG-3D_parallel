/*
* \copyright (C) 2013 Puzikova V.V.
* \file      coarse_grid.h
* \author    Puzikova Valeria (vvp) 
* \date      04.10.2013 14:35:51
* \brief     coarse_grid.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_COARSE_GRID_H_
#define _LS_STAG_COARSE_GRID_H_

// ================================================================== INCLUDES
#include "../SparseMatrix/transposed_sparse_matrix.h"
#include "../SparseMatrix/product_sparse_matrix.h"
#include "SmootherADLJ/smoother_adlj.h"
#include "TransferOperators/transfer_operators_odd_odd_odd.h"
#include "TransferOperators/transfer_operators_even_even_odd.h"
#include "TransferOperators/transfer_operators_even_odd_odd.h"
#include "TransferOperators/transfer_operators_odd_even_odd.h"
#include "TransferOperators/transfer_operators_odd_odd_even.h"
#include "TransferOperators/transfer_operators_even_even_even.h"
#include "TransferOperators/transfer_operators_even_odd_even.h"
#include "TransferOperators/transfer_operators_odd_even_even.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса CoarseGrid

Данный файл содержит в себе реализацию класса CoarseGrid
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий сеточный уровень многосеточного метода

Данный класс реализует построение и хранение сеточного уровня многосеточного метода по матрице системы с предыдущего 
подробного сеточного слоя (в частности, строятся матрица системы на грубой сетке, дискретные аналоги операторов проекции и 
 продолжения между предыдущим и текущим сеточными слоями, сглаживатель ADLJ)
*/
struct CoarseGrid
{
    CoarseGrid() : 
		baseTransfer(NULL),
		oddN(false), 
		oddM(false), 
		oddK(false) 
	{
		;
	}
	~CoarseGrid()
	{
		if (baseTransfer) 
		{ 
			delete baseTransfer; 
			baseTransfer = NULL; 
		}
	}
    // Вспомогательные функции
    void nodesCount(bool &odd, int &nCoarse, const int nFine) {
        if (fmod(nFine, 2.0) > EQ_DOUBLE) {
            odd = true;
            nCoarse = (int)floor((nFine + 1) * 0.5);
        } else {
            nCoarse = (int)floor(nFine * 0.5);
        }
    }

    void initTransfer() {
        // Выбор одного из 8 вариантов по трём флагам
        if (oddN) {
            if (oddM) {
                 baseTransfer = (oddK) ? ((TransferOperatorsInterface*)(new TransferOperatorsOddOddOdd)) :
				  ((TransferOperatorsInterface*)(new TransferOperatorsOddOddEven));
            } else {
                baseTransfer = (oddK) ? ((TransferOperatorsInterface*)(new TransferOperatorsOddEvenOdd)) :
				((TransferOperatorsInterface*)(new TransferOperatorsOddEvenEven));
            }
        } else {
            if (oddM) {
                baseTransfer = (oddK) ? ((TransferOperatorsInterface*)(new TransferOperatorsEvenOddOdd)) :
				((TransferOperatorsInterface*)(new TransferOperatorsEvenOddEven));
            } else {
                baseTransfer = (oddK) ? ((TransferOperatorsInterface*)(new TransferOperatorsEvenEvenOdd)) :
				((TransferOperatorsInterface*)(new TransferOperatorsEvenEvenEven));
            }
        }
        baseTransfer->build(&P, &R, NN, MM, KK, N, M, K);  // + сразу строятся
    }

    void set(SparseMatrix* _A) {
        m_A = _A;
        NN = m_A->m_N;
        MM = m_A->m_M;
        KK = m_A->m_K;                  // добавлено
        nodesCount(oddN, N, NN);
        nodesCount(oddM, M, MM);
        nodesCount(oddK, K, KK);        // добавлено
        NMK = N * M * K;                // общее число ячеек на грубой сетке
        initTransfer();
        V.set((SparseMatrix*)&R, m_A);
        AH.set((SparseMatrix*)&V, &P);
        ADLJ.set(m_A);                  // теперь SmootherADLJ3D
    }

    void rebuild() {
        ADLJ.rebuild();
        V.resize();
        AH.resize();
    }
    TransposedSparseMatrix R;
    SparseMatrix P;
    ProductSparseMatrix AH;
    ProductSparseMatrix V;
    SmootherADLJ ADLJ;                // <-- трёхмерный сглаживатель
    SparseMatrix* m_A;
    TransferOperatorsInterface* baseTransfer;
    int N, M, K;          // размеры грубой сетки
    int NMK;               // N * M * K
    int NN, MM, KK;        // размеры подробной сетки (предыдущей)
    bool oddN, oddM, oddK; // чётности на подробной сетке
};
}

#endif //_LS_STAG_COARSE_GRID_H_
