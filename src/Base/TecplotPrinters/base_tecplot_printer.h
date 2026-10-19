/*
* \copyright (C) 2016 Puzikova V.V.
* \file      tecplot_printer.h
* \author    Puzikova Valeria (vvp) 
* \date      05.04.2016 15:22:42
* \brief     base_tecplot_printer.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_BASE_TECPLOT_PRINTER_H_
#define _LS_STAG_BASE_TECPLOT_PRINTER_H_
// ================================================================== INCLUDES
#include "../TaskDescription/Tasks/Flow/Rotor/flow_rotor_task_description.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса BaseTecplotPrinter

Данный файл содержит в себе реализацию класса BaseTecplotPrinter
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Базовый класс, реализующий печать в Tecplot Data Format

Данный базовый класс реализует работу по выводу рассчитанных полей неизвестных в Tecplot Data Format.
Можно задать шаг по времени, с которым необходимо выводить данные. Номер шага по времени, на котором осуществляется печать,
отражается в имени файла
*/
	struct BaseTecplotPrinter
	{
        BaseTecplotPrinter()
        {
            ;
        }
		    void init(BaseTaskDescription* _task, int* _it,
              ArrayInterface<double>* _u, ArrayInterface<double>* _v,
              ArrayInterface<double>* _w,          // <-- добавлено
              ArrayInterface<double>* _p, double periodPrint, bool _needAverage)
    {
        u = _u; v = _v; w = _w; p = _p;
        m_it  = _it;
        m_task = _task;
        needAverage = _needAverage;

        m_N  = m_task->m_mesh->baseMesh.N;
        m_M  = m_task->m_mesh->baseMesh.M;
        m_K  = m_task->m_mesh->baseMesh.K;   // <-- новая размерность
        m_NM = m_N * m_M;
        m_NMK = m_NM * m_K;

        m_xN = m_task->m_mesh->xMesh.N;
        m_yM = m_task->m_mesh->yMesh.M;
        m_zK = m_task->m_mesh->zMesh.K;      // <-- размер по z

        type  = &(m_task->m_mesh->baseMesh.m_typeOfCell);
        theta = &(m_task->m_mesh->m_theta);
        delta = &(m_task->m_mesh->m_delta);

        ItPrint   = (int)floor(periodPrint / m_task->t);
        ItCounter = 1;
        defineMeshStr();

        // заголовок теперь включает K
        title = "I=" + auxStr + ", J=" + auxStr2 + ", K=" + auxStr3 + ", F=POINT";
        if (needAverage) fullFileName = filename + _T("Average.dat");
    }

    void print()
    {
        if (ItPrint)
        {
            if ((ItCounter == ItPrint) || ((*m_it) == 1))
            {
                if (!needAverage) fullFileName = filename + NumberIter(*m_it) + _T(".dat");
                printField();
                ItCounter = 1;
            }
            else ItCounter++;
        }
    }

protected:
    virtual void defineMeshStr() {}
    virtual void printField() {}

    std::tstring NumberIter(const int it)
	{
		std::tstring st = strConv.toTString(it);
		if (it < 100000) st = _T("0") + st;
		if (it < 10000) st = _T("0") + st;
		if (it < 1000) st = _T("0") + st;
		if (it < 100) st = _T("0") + st;
		if (it < 10) st = _T("0") + st;
		return st;
	};
    StringConvertor strConv;
    std::tstring filename, fullFileName;
    std::string title, auxStr, auxStr2, auxStr3;  // <-- auxStr3 для K
    double t_y, t_x, t_z, t_Ux, t_Uy, t_Uz;      // <-- добавлены t_z, t_Uz
    int ItPrint, ItCounter;
    int m_N, m_M, m_K, m_xN, m_yM, m_zK, m_NM, m_NMK;         // <-- m_K, m_zK
    BaseTaskDescription* m_task;
    LSSTAGBaseMeshParams* m_baseMesh;
    ArrayInterface<double> *u, *v, *w, *p;        // <-- w
    double *m_u, *m_v, *m_w, *m_p;                // <-- m_w
    Array<coordinate> *theta, *delta;
    coordinate *curT, *curD, *m_theta, *m_delta;
    int* m_it;
    Array<CellTypes>* type;
    CellTypes* m_type;
    bool needAverage;
	};
};

#endif //_LS_STAG_BASE_TECPLOT_PRINTER_H_
