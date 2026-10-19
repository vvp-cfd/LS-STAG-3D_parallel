/*
* \copyright (C) 2024 Puzikova V.V.
* \file      profiler.h
* \author    Puzikova Valeria (vvp)
* \date      21.02.2024 10:27:05
* \brief     profiler.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_PROFILER_H_
#define _LS_STAG_PROFILER_H_

// ================================================================== INCLUDES
#include <assert.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include "tstring.h"
#ifdef LS_STAG_USE_CUDA
#include <cuda_runtime_api.h>
#include <helper_timer.h>
#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif (CUDART_VERSION >= 10000)
#include <nvtx3/nvToolsExt.h>
#elif (CUDART_VERSION < 10000)
#include <nvToolsExt.h>
#error Unknown CUDART_VERSION!
#endif
#else
#include <omp.h>
#endif
#include "ls_stag_parallel.h"
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с описанием класса Profiler

Данный файл содержит в себе реализацию класса Profiler
*/
// ===========================================================================


namespace LS_STAG
{
	/*!
	\brief Класс для профилирования

	Данный класс реализует профилировку по размеченным секциям кода.
	*/
    class Profiler 
    {
    private:
        struct _Timer 
        {
            bool running;
            double time;
            int invocations;
#ifdef LS_STAG_USE_CUDA
            StopWatchInterface* timer;
            _Timer() : running(false), time(0), invocations(0)
            {
                sdkCreateTimer(&timer);
                start();
            }
            double duration() { return sdkGetTimerValue(&timer) / 1000.0; }
            void stop()
            {
                assert(running);
                sdkStopTimer(&timer);
                time += duration();
                running = false;
                sdkResetTimer(&timer);
            }
            void start()
            {
                if (running) {
                    std::cerr << "ERROR: timer already running. Reported timings will be inaccurate." << std::endl;
                    stop();
                }
                assert(!running);
                running = true;
                ++invocations;
                sdkStartTimer(&timer);
            }
#else
            std::chrono::time_point<std::chrono::steady_clock> startTime;
            _Timer() : running(false), time(0), invocations(0)
            {
                start();
            }
            double duration() { return (std::chrono::steady_clock::now() - startTime).count() / 1000000000.0; }
            void stop()
            {
                assert(running);
                time += duration();
                running = false;
            }
            void start()
            {
                if (running) {
                    std::cerr << "ERROR: timer already running. Reported timings will be inaccurate." << std::endl;
                    stop();
                }
                assert(!running);
                running = true;
                ++invocations;
                startTime = std::chrono::steady_clock::now();
            }
#endif
            // gets elapsed time (even if currently running)
            double elapsed() { return running ? time + duration() : time; }
        };
        typedef std::map<std::string, _Timer> TimerMap;
        struct _Section : public _Timer 
        {
            TimerMap timers;
            _Section() : _Timer() { ; }
            void startTimer(std::string name)
            {
                auto it = timers.find(name);
                if (it != timers.end()) {
                    if (it->second.running) std::cerr << "ERROR : timer " << name << " already started. " << std::endl;
                    it->second.start();
                }
                else timers[name] = _Timer();
            }
            using _Timer::duration;
            using _Timer::start;
            void stop()
            {
                _Timer::stop();
                // Also stop all our sub-timers...
                for (auto& entry : timers) {
                    if (entry.second.running) {
                        std::cerr << "WARNING: stopping timer " << entry.first
                            << " implicitly in enclosing section's stop()" << std::endl;
                        entry.second.stop();
                    }
                }
            }
            void start(const std::string& name)
            {
                auto lb = timers.lower_bound(name);
                if ((lb == timers.end()) || (lb->first != name))
                    timers.emplace_hint(lb, name, _Timer());
                else
                    lb->second.start(); // The full section timer must be started too...
            }
            void stop(const std::string& name) { timers.at(name).stop(); }
            void report(std::ostream& os)
            {
                for (auto& entry : timers)
                    os << displayName(entry.first) << '\t' << entry.second.elapsed() << '\t'
                    << entry.second.invocations << std::endl;
            }
        };
        typedef std::map<std::string, _Section> SectionMap;
        typedef SectionMap::iterator SectionIterator;
        typedef SectionMap::const_iterator SectionConstIterator;
        SectionMap m_sections;
        std::list<std::string> m_sectionStack;
        static std::string displayName(const std::string name)
        {
            size_t levels = 0;
            for (char c : name)
                if (c == ':')
                    ++levels;
            if (levels == 0) return "-> " + name;
            std::string result(4 * levels, ' ');
            result += "|-> ";
            result.append(name, name.rfind(':') + 1, std::string::npos);
            return result;
        }
        static std::string lastName(const std::string name)
        {
            const size_t startPos = name.rfind(':');
            if (startPos == std::string::npos) return name;
            return name.substr(startPos + 1);
        }
    public:
        Profiler() { reset(); }
        void startSection(const std::string& name)
        {
            std::string label = (!m_sectionStack.empty()) ? (m_sectionStack.back() + ':' + name) : name;
            m_sectionStack.push_back(label);
            auto lb = m_sections.lower_bound(label);
            if ((lb == m_sections.end()) || (lb->first != label))
                m_sections.emplace_hint(lb, label, _Section());
            else
                lb->second.start(); // The full section timer must be started too..
#ifdef LS_STAG_USE_CUDA
            {
                nvtxRangePushA(name.c_str());
            }
#endif
        }
        void stopSection(const std::string& name)
        {
            assert(!m_sectionStack.empty());
#ifdef LS_STAG_USE_CUDA
            nvtxRangePop();
#endif
            std::string currentLabel = m_sectionStack.back();
            m_sectionStack.pop_back();
            std::string label;
            label = (!m_sectionStack.empty()) ? (m_sectionStack.back() + ':' + name) : name;
            if (label != currentLabel) {
                std::cerr << "ERROR: sections must be stopped in the reverse of the order they were started." << std::endl;
                std::cerr << "(Expected " << currentLabel << ", but got " << label << ")" << std::endl;
            }
            m_sections.at(label).stop();
        }
        void start(std::string timer)
        {
            std::string sectionName;
            if (!m_sectionStack.empty()) {
                sectionName = m_sectionStack.back();
                timer = sectionName + ':' + timer;
            }
            m_sections.at(sectionName).start(timer);
        }
        void stop(std::string timer)
        {
            std::string sectionName;
            if (!m_sectionStack.empty()) {
                sectionName = m_sectionStack.back();
                timer = sectionName + ':' + timer;
            }
            m_sections.at(sectionName).stop(timer);
        }
        // Remove all timers and restart the global section
        void reset()
        {
            m_sections.clear();
            m_sectionStack.clear();
            // (re)start global section
            m_sections.insert(std::make_pair("", _Section()));
        }
        void report(std::ostream& os)
        {
            for (SectionIterator it = m_sections.begin(); it != m_sections.end(); ++it) {
                if (it->first != "") { // Skip global section... i
                    os << displayName(it->first) << "\t" << it->second.elapsed() << "\t"
                        << it->second.invocations << std::endl;
                    it->second.report(os);
                }
            }
            // TODO add percentage
            m_sections.at("").report(os);
            os << "Full time\t" << m_sections.at("").elapsed() << std::endl;
            os << "========" << std::endl;
        }
        void statReport(const std::tstring& name, const int NM)
        {
            const SectionIterator updIt = std::next(m_sections.begin());
            const std::tstring filename = name + _T(".csv");
            std::ofstream file; 
            file.open(filename.c_str());
            file << "NM," << lastName(updIt->first) << " (ms)";
            for (SectionIterator it = std::next(updIt); it != m_sections.end(); ++it)
                file << "," << lastName(it->first) << " (%)";
            file << std::endl;
            const double updFullTimeSec = updIt->second.elapsed();
            file << NM << ","
                << (1000.0 * updFullTimeSec / updIt->second.invocations); // sec to msec
            for (SectionIterator it = std::next(updIt); it != m_sections.end(); ++it)
                file << "," << (100.0 * it->second.elapsed() / updFullTimeSec);
            file << std::endl;
            file.close();
        }
    };
#ifdef LS_STAG_PROFILING
    void PROFILER_RESET();
    void PROFILER_REPORT();
    void PROFILER_REPORT(std::ostream& str);
    void PROFILER_STAT_REPORT(const std::tstring& name, const int NM);
    struct ScopedProfiler
    {
        ScopedProfiler(const std::string& name) : m_name(name)
        {
            if (!inParallel()) start();
        }
        ~ScopedProfiler()
        {
            if (!inParallel()) end();
        }
    private:
        void start();
        void end();
        bool inParallel()
        {
#if !defined(LS_STAG_USE_TBB) && !defined(LS_STAG_USE_CILK) && !defined(LS_STAG_USE_CUDA)
            return (omp_get_thread_num() > 0);//omp_in_parallel(); // actual for DDM in MT
#endif
            return false;
    }
        std::string m_name;
    };
    #define LS_STAG_PROFILE(label) ScopedProfiler __profile(label)
#else
    ///stubs for compilation
    inline void PROFILER_RESET() {}
    inline void PROFILER_REPORT() {}
    inline void PROFILER_REPORT(std::ostream& out) {}
    inline void PROFILER_STAT_REPORT(const std::tstring, const int) { ; }
    /**
    * Stub function if we do not need benchmarking
    */
    struct ScopedProfiler
    {
        ScopedProfiler(std::string) { ; }
    };
    #define LS_STAG_PROFILE(label)
#endif	
}

#endif //_LS_STAG_PROFILER_H_
