/*
* \copyright (C) 2013 Puzikova V.V.
* \file      sparse_matrix.h
* \author    Puzikova Valeria (vvp) 
* \date      01.10.2013 11:34:59
* \brief     sparse_matrix.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SPARSE_MATRIX_H_
#define _LS_STAG_SPARSE_MATRIX_H_

// ================================================================== INCLUDES
#include <fstream>
#include <string>
#include "../MeshFunctions/array.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса SparseMatrix

Данный файл содержит в себе реализацию класса SparseMatrix
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс для работы с разреженными матрицами

Данный базовый класс реализует работу с разреженными матрицами. 
Матрицы хранятся в формате CSR (Compressed Sparse Row).
*/struct SparseMatrix
{
    SparseMatrix()
        : m_N(0), m_M(0), m_K(0), NMK(0),
          elemCount(0), maxElemCount(0),
          Cell(nullptr), Portrait(nullptr), Number(nullptr)
	{
		;
	}

    // ---------- 2D-совместимость ----------
    virtual void set(const int _N, const int _M)
    {
        set(_N, _M, 1);          // вырожденный 3D: K = 1
    }

    // ---------- 3D-версия ----------
    virtual void set(const int _N, const int _M, const int _K)
    {
        setSizes(_N, _M, _K);
        init();
        estimateNeedMemory();
        reserve();
    }

    // Ручное задание размеров и NMK (например, для нерегулярных сеток)
    void set(const int _N, const int _M, const int _K, const int _NMK)
    {
        m_N = _N;
        m_M = _M;
        m_K = _K;
        NMK = _NMK;
        init();
    }

    void setCopy(SparseMatrix* _matrix)
    {
        m_N = _matrix->m_N;
        m_M = _matrix->m_M;
        m_K = _matrix->m_K;
        NMK = _matrix->NMK;
        init();
    }

    virtual void set(SparseMatrix* _Matrix) {}   // оставлено для перегрузки

    // ---------- Печать ----------
    virtual void print(const std::string &name)
    {
        std::ofstream f, g;
        int *cpuNumber, *cpuPortrait;
        double* cpuCell;
        m_Cell.setCopy(cpuCell);
        m_Number.setCopy(cpuNumber);
        m_Portrait.setCopy(cpuPortrait);
        std::string fileName = name + ".dat";
        f.open(fileName.c_str());
        fileName = name + "_P.dat";
        g.open(fileName.c_str());
        for (int i = 0; i < NMK; ++i) {
            int num3 = cpuNumber[i+1];
            int num2 = cpuNumber[i];
            int num = num3 - num2;
            g << num << " ";
            for (int j = num2; j < num3; ++j) {
                f << cpuCell[j] << " ";
                g << cpuPortrait[j] << " ";
            }
            f << std::endl;
            g << std::endl;
        }
        f.close();
        g.close();
    }

    virtual void printNotSparse(const std::string &name)
    {
        std::ofstream f;
        int* cpuNumber, *cpuPortrait;
        double* cpuCell;
        m_Cell.setCopy(cpuCell);
        m_Number.setCopy(cpuNumber);
        m_Portrait.setCopy(cpuPortrait);
        std::string fileName = name + ".dat";
        f.open(fileName.c_str());
        for (int i = 0; i < NMK; ++i) {
            for (int k = 0, num3 = cpuNumber[i+1], j = cpuNumber[i]; k < NMK; ++k) {
                if (k == cpuPortrait[j]) {
                    f << cpuCell[j] << " ";
                    if (j < num3) j++;
                } else
                    f << 0 << "  ";
            }
            f << std::endl;
        }
        f.close();
    }

    // Ядра индексации (без изменений, но зависят от N, M, K)
    LS_STAG_KERNEL_INLINE void setPC(int& num, LS_STAG_R_PD _Cell, LS_STAG_R_PI _Portrait,
                                     const int number, const double cell)
    {
        _Portrait[num] = number;
        _Cell[num] = cell;
        num++;
    }

    LS_STAG_KERNEL_INLINE void NumWidth2(const int m_a, const int m_b, int& num, int* _Number, const int _N2)
    {
        for (int i = 0; i < _N2; i++) {
            _Number[num + 1] = _Number[num] + m_a;
            num++;
            _Number[num + 1] = _Number[num] + m_b;
            num++;
        }
        _Number[num + 1] = _Number[num] + m_a;
        num++;
    }

    LS_STAG_KERNEL_INLINE void NumWidth4(const int m_a, const int m_b, int& num, int* _Number, const int _N2)
    {
        _Number[num + 1] = _Number[num] + m_a;
        num++;
        for (int i = 0; i < _N2; i++) {
            _Number[num + 1] = _Number[num] + m_b;
            num++;
            _Number[num + 1] = _Number[num] + m_b;
            num++;
        }
        _Number[num + 1] = _Number[num] + m_a;
        num++;
    }

    LS_STAG_KERNEL_INLINE void NumWidth_1(const int nn, const int m_a, int& num, int* _Number)
    {
        for (int i = 0; i < nn; i++) {
            _Number[num + 1] = _Number[num] + m_a;
            num++;
        }
    }
	
    LS_STAG_KERNEL_INLINE void NumWidth_2(const int _N2, const int m_a, const int m_b, int& num, int* _Number)
    {
        for (int i = 0; i < _N2; i++) {
            _Number[num + 1] = _Number[num] + m_a;
            num++;
            _Number[num + 1] = _Number[num] + m_b;
            num++;
        }
    }

protected:
    virtual void estimateNeedMemory() {}
    void init()
    {
        m_Number.init(NMK + 1, Number);
        m_Number.setElem(0, 0);
        elemCount = 0;
    }
    void setSizes(const int _N, const int _M, const int _K)
    {
        m_N = _N;
        m_M = _M;
        m_K = _K;
        NMK = m_N * m_M * m_K;
    }
public:
    void reserve()
    {
        if (maxElemCount > 0) {
            m_Cell.init(maxElemCount, Cell);
            m_Portrait.init(maxElemCount, Portrait);
        }
    }

    virtual void resize()
    {
        elemCount = m_Number.getElem(NMK);
    }

    void finResize()
    {
        resize();
        maxElemCount = elemCount;
        reserve();
    }

    Array<double> m_Cell;
    Array<int> m_Portrait, m_Number;
    double* Cell;
    int* Portrait;
    int* Number;
    int m_N, m_M, m_K;   // размеры по трём осям
    int NMK;              // общее число ячеек = N * M * K
    int elemCount;
    int maxElemCount;
};
}
#endif //_LS_STAG_SPARSE_MATRIX_H_
