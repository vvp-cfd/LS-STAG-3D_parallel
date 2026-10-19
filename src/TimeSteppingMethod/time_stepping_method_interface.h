/*
* \copyright (C) 2013 Puzikova V.V.
* \file      time_stepping_method_interface.h
* \author    Puzikova Valeria (vvp) 
* \date      29.12.2013 23:14:06
* \brief     time_stepping_method_interface.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_TIME_STEPPING_METHOD_INTERFACE_H_
#define _LS_STAG_TIME_STEPPING_METHOD_INTERFACE_H_

// ================================================================== INCLUDES
#include <iostream>
#include "../Base/timer.h"
#include "../Base/profiler.h"
#include "../Base/TaskDescription/task_descriptions_selector.h"
#include "../Base/set_solutions.h"
#include "../DiscreteEquations/poisson_discrete_equation.h"
#include "../DiscreteEquations/helmholtz_discrete_equation_x.h"
#include "../DiscreteEquations/helmholtz_discrete_equation_y.h"
#include "../DiscreteEquations/helmholtz_discrete_equation_z.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса TimeSteppingMethodInterface

Данный файл содержит в себе реализацию класса TimeSteppingMethodInterface
*/
// ===========================================================================

namespace LS_STAG
{
/*!
\brief Класс, реализующий управление вычислительным процессом

Данный класс реализует управление вычислительным процессом при численном моделировании нестационарных процессов
*/
/*!
\brief Класс, реализующий управление вычислительным процессом в 3D

Данный класс реализует управление вычислительным процессом при численном моделировании нестационарных
трёхмерных течений. Добавлено третье уравнение Гельмгольца для w-компоненты скорости (Oz).
*/
struct TimeSteppingMethodInterface
{
    void computation(const int _nodesCount, const int _domainCount = 1)
    {
        PROFILER_RESET();
        LS_STAG_PROFILE("computation");
#ifdef LS_STAG_USE_CUDA
        domainCount = 0;
#else
        domainCount = _domainCount;
        initParallel(_nodesCount);
#endif
        init();
#ifdef LS_STAG_USE_TBB
        logfile << "Число потоков TBB: " << _nodesCount << std::endl;
        std::cout << "TBB threads number: " << _nodesCount << std::endl;
#endif
        while (curTime < m_task->T_end)
        {
            step();
            it++;
            curTime += m_task->t;
        }
        close();
        closeParallel();
    }
	
private:
    void initParallel(int _nodesCount)
    {
#ifndef LS_STAG_USE_CUDA
#ifdef LS_STAG_USE_TBB
        tbb::task_scheduler_init initTBB(_nodesCount);
#elif LS_STAG_USE_CILK
        StringConvertor strConv;
        std::string nodesNumber = strConv.toString(_nodesCount);
        __cilkrts_set_param("nworkers", nodesNumber.c_str());
#else
        omp_set_num_threads(_nodesCount);
#endif
#ifdef LS_STAG_USE_MKL
        mkl_set_num_threads(_nodesCount);
#endif
#endif
    }

    void closeParallel()
    {
#ifdef LS_STAG_USE_CILK
        __cilkrts_end_cilk();
#endif
    }

    void init()
    {
        LS_STAG_PROFILE("init");
        startSol.setTime();
        tasks.load(std::string(CONFIG_PATH) + "/TaskDescription.dat");
        m_task = tasks.m_task;

#ifdef _WIN32
        std::tstring sep = _T("\\");
#else
        std::tstring sep = "/";
#endif
        outputDirectoryFull = m_task->outputDirectory + sep + startSol.timeInfoString.labelT;
#ifdef _WIN32
        CreateDirectory(outputDirectoryFull.c_str(), NULL);
#else
        mkdir(outputDirectoryFull.c_str(), 0777);
#endif
        outputDirectoryFull += sep;
        m_task->outputDirectoryFull = outputDirectoryFull;
        logfileName = outputDirectoryFull + startSol.timeInfoString.labelT + _T(".log");
        logfile.open(logfileName.c_str());
        logfile << "Время начала расчёта: " << startSol.timeInfoString.date << std::endl;
        std::cout << "Start: " << startSol.timeInfoString.date << std::endl;
        profileName = outputDirectoryFull + _T("profile.log");
        profile.open(profileName.c_str());

#ifndef LS_STAG_USE_CUDA
#ifndef LS_STAG_USE_TBB
#ifdef LS_STAG_USE_CILK
        int theadsNum = __cilkrts_get_nworkers();
        logfile << "Число потоков Cilk Plus: " << theadsNum << std::endl;
        std::cout << "Cilk Plus threads number: " << theadsNum << std::endl;
#else
        int theadsNum;
#pragma omp parallel
        {
#pragma omp master
            {
                theadsNum = omp_get_num_threads();
            }
        }
        logfile << "Число потоков OpenMP: " << theadsNum << std::endl;
        std::cout << "OpenMP threads number: " << theadsNum << std::endl;
#endif
#endif
#endif

        it = 1;
        curTime = 0.0;
        m_task->NUM_SUBDOMAINS = domainCount;
        logfile << "Число подобластей: " << domainCount << std::endl;
        std::cout << "subdomains number: " << domainCount << std::endl;

        m_task->init(&curTime);
        sol.init(m_task, &it);
        m_task->initSol(&(sol.current), &(sol.prediction), &(sol.previous));

        if (m_task->continueSolveTask)
        {
            sol.load(std::string(CONFIG_PATH) + "/continueSolveTask.dat");
            m_task->rebuildMesh();
        }
		
        pEq.set(m_task, &(sol.prediction.p));
        m_task->m_mesh->U_ib = pEq.m_A->flux;

        if (m_task->continueSolveTask)
        {
            it++;
            curTime += m_task->t;
        }
		
        uEq.set(m_task, &(sol.prediction.u));
        vEq.set(m_task, &(sol.prediction.v));
        wEq.set(m_task, &(sol.prediction.w));      // новое уравнение для w

        sol.set((SparseMatrix*)(uEq.m_grad), (SparseMatrix*)(vEq.m_grad),
                (SparseMatrix*)(wEq.m_grad), uEq.coef);   // передаём три градиента

        tasks.print();
    }
    void step()
    {
        LS_STAG_PROFILE("step");
        {
            LS_STAG_PROFILE("preSolve");
            sol.init();
            m_task->rebuildMesh();
            pEq.m_A->calcFlux();
        }
        {
            LS_STAG_PROFILE("AllSolve");
            uEq.solve();
            vEq.solve();
            wEq.solve();      // решение уравнения для w
            pEq.solve();
        }
        {
            LS_STAG_PROFILE("postSolve");
            sol.correction();   // в correction уже добавлена обработка w
            m_task->printLoads();   // нагрузки (в 3D должны учитывать z-компоненту)

            sol.print();
            {
                LS_STAG_PROFILE("end");
                logfile << it << " " << curTime << " ";
                logfile << uEq.solver.it << " " << vEq.solver.it << " " << wEq.solver.it << " " << pEq.solver.it << " ";
                logfile << uEq.solver.norm_r << " " << vEq.solver.norm_r << " " << wEq.solver.norm_r << " " << pEq.solver.norm_r << "\n";
                {
                    std::cout << it << " " << curTime << " ";
                    std::cout << uEq.solver.it << " " << vEq.solver.it << " " << wEq.solver.it << " " << pEq.solver.it << std::endl;
                    if (it % FRAME_STEP_PROFILE_REPORT == 0) PROFILER_REPORT(profile);
                }
            }
        }
    }

    void close()
    {
        LS_STAG_PROFILE("close");
        endSol.setTime();
        logfile << "Время завершения расчёта: " << endSol.timeInfoString.date << std::endl;
        std::cout << "End: " << endSol.timeInfoString.date << std::endl;
        endSol -= startSol.rawtime;
        logfile << "Продолжительность расчёта: " << endSol.timeInfoString.label << " секунд ( ";
        logfile << endSol.timeInfoString.date << " )\n";
        std::cout << "Total: " << endSol.timeInfoString.label << " sec ( " << endSol.timeInfoString.date << " )\n";
        logfile.close();
        profile.close();
        PROFILER_STAT_REPORT(outputDirectoryFull + _T("profileStat"), m_task->m_mesh->baseMesh.NMK);
    }

    TaskDescriptionsSelector tasks;
    PoissonDiscreteEquation pEq;       ///< разностный аналог уравнения Пуассона
    HelmholtzDiscreteEquationX uEq;    ///< уравнение Гельмгольца на x-сетке
    HelmholtzDiscreteEquationY vEq;    ///< уравнение Гельмгольца на y-сетке
    HelmholtzDiscreteEquationZ wEq;    ///< новое уравнение Гельмгольца на z-сетке
    SetSolutions sol;
    Timer startSol, endSol;
    std::tstring outputDirectoryFull, logfileName, profileName;
    std::ofstream logfile, profile;
    BaseTaskDescription* m_task;
    double curTime;
    int it;
    int domainCount;

};
}

#endif //_LS_STAG_TIME_STEPPING_METHOD_INTERFACE_H_
