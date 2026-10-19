/*
* \copyright (C) 2014 Puzikova V.V.
* \file      solution.h
* \author    Puzikova Valeria (vvp) 
* \date      03.01.2014 21:31:58
* \brief     solution.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_SOLUTION_H_
#define _LS_STAG_SOLUTION_H_

// ================================================================== INCLUDES
#include "MeshFunctions/real_valued_mesh_function.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса Solution

Данный файл содержит в себе реализацию класса Solution
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий работу с полями скоростей и давления

Данный класс реализует хранение полей скоростей и давления и некоторые базовые операции с ними
*/
struct Solution {
    void init(const int Nx, const int Ny, const int Nz, const int N) {
        u.init(Nx);
        v.init(Ny);
        w.init(Nz);   // <-- инициализация z-компоненты
        p.init(N);
    }

    // ---- Получение указателей на сырые данные ----
    // Добавлен выходной параметр для w.
    void getCopy(double*& m_u, double*& m_v, double*& m_w, double*& m_p) {
        u.getCopy(m_u);
        v.getCopy(m_v);
        w.getCopy(m_w);   // <-- получение данных w
        p.getCopy(m_p);
    }

    // ---- Копирование полного решения ----
    // Копируются u, v, w и p.
    void copy(Solution* source) {
        u.copyR(&(source->u));
        v.copyR(&(source->v));
        w.copyR(&(source->w));   // <-- копирование w
        p.copyR(&(source->p));
    }

    // ---- Присвоение одного и того же значения всем полям ----
    void set(const double value) {
        u.set(value);
        v.set(value);
        w.set(value);   // <-- присвоение значения w
        p.set(value);
    }

    // ---- Печать полей ----
    // В имя файла добавляется "_w".
    void print(const std::string &name) {
        u.print(name + "_u");
        v.print(name + "_v");
        w.print(name + "_w");   // <-- печать w
        p.print(name + "_p");
    }

    // ---- Сохранение состояния ----
    void loadState(std::ifstream& in) {
        u.loadState(in);
        v.loadState(in);
        w.loadState(in);   // <-- загрузка w
        p.loadState(in);
    }

    void saveState(std::ofstream& out) {
        u.saveState(out);
        v.saveState(out);
        w.saveState(out);   // <-- сохранение w
        p.saveState(out);
    }
    // ---- Data members (расширены для 3D) ----
    RealValuedMeshFunction u;  ///< x-компонента скорости
    RealValuedMeshFunction v;  ///< y-компонента скорости
    RealValuedMeshFunction w;  ///< z-компонента скорости (новая)
    RealValuedMeshFunction p;  ///< давление

};

} // namespace LS_STAG

#endif //_LS_STAG_SOLUTION_H_