/*
* \copyright (C) 2013 Puzikova V.V.
* \file      ls_stag_base_mesh_params.h
* \author    Puzikova Valeria (vvp) 
* \date      23.10.2013 18:53:14
* \brief     ls_stag_base_mesh_params.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_BASE_MESH_PARAMS_H_
#define _LS_STAG_LS_STAG_BASE_MESH_PARAMS_H_

// ================================================================== INCLUDES
#include "../MeshFunctions/array.h"
#include "cut_cells.h"
#include "../ls_stag_real_constants.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием классов BaseMeshParams и LSSTAGBaseMeshParams

Данный файл содержит в себе реализацию классов BaseMeshParams и LSSTAGBaseMeshParams
*/
// ===========================================================================

namespace LS_STAG
{
	LS_STAG_DECLARATION_KERNEL(inv, LS_STAG_CR_PD, _VOF, LS_STAG_R_PD, _inv_VOF);
/*!
\brief Класс, реализующий работу с основными параметрами сетки

Вспомогательный класс BaseMeshParams реализует работу с основными параметрами прямоугольной структурированной сетки
*/
struct BaseMeshParams
{
    BaseMeshParams()
        : N(0), M(0), K(0), NMK(0)
    { }

    /*!
    Вычисляет число ячеек сетки как произведение числа ячеек по трём направлениям.
    */
    virtual void calcNodesCount()
    {
        NMK = N * M * K;
    }

    int N;    ///< число ячеек сетки по горизонтали (Ox)
    int M;    ///< число ячеек сетки по вертикали (Oy)
    int K;    ///< число ячеек сетки по глубине (Oz) – добавлено для 3D
    int NMK;  ///< общее число ячеек сетки (N * M * K)
};

/*!
\brief Класс, реализующий работу с основными параметрами LS-STAG-сетки в 3D

Класс LSSTAGBaseMeshParams реализует работу с основными параметрами прямоугольной
структурированной сетки, входящей в число разнесенных сеток LS-STAG-сетки.
*/
struct LSSTAGBaseMeshParams : BaseMeshParams
{
    LSSTAGBaseMeshParams()
        : BaseMeshParams()
        , typeOfCell(NULL)
        , VOF(NULL)
        , inv_VOF(NULL)
    { }

    void calcNodesCount() override
    {
        BaseMeshParams::calcNodesCount();      // NMK = N * M * K
        m_typeOfCell.init(NMK, typeOfCell);
        m_VOF.init(NMK, VOF);
        m_inv_VOF.init(NMK, inv_VOF);
    }

    LS_STAG_KERNEL(inv, LS_STAG_CR_PD _VOF, LS_STAG_R_PD _inv_VOF)
    {
        _inv_VOF[i] = (_VOF[i] < EQ_TINY_DOUBLE) ? 0.0 : (1.0 / _VOF[i]);
    }

    void inverseVOF()
    {
        LS_STAG_PAR_FOR(LSSTAGBaseMeshParams, inv, NMK, VOF, inv_VOF)
    }

    Array<CellTypes> m_typeOfCell;
    Array<double> m_VOF;
    Array<double> m_inv_VOF;
    CellTypes* typeOfCell;   ///< Типы ячеек сетки
    double* VOF;             ///< Объемы ячеек сетки
    double* inv_VOF;         ///< 1 / VOF
};
}

#endif //_LS_STAG_LS_STAG_BASE_MESH_PARAMS_H_
