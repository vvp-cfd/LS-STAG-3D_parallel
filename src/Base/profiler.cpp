#include "profiler.h"

namespace LS_STAG {
#ifdef LS_STAG_PROFILING
    std::unique_ptr<Profiler> g_profiler = std::make_unique<Profiler>();
    void PROFILER_RESET() { g_profiler->reset(); }
    void PROFILER_REPORT() { g_profiler->report(std::cout); }
    void PROFILER_REPORT(std::ostream& str) { g_profiler->report(str); }
    void PROFILER_STAT_REPORT(const std::tstring& name, const int NM)
    {
        g_profiler->statReport(name, NM);
    }
    void ScopedProfiler::start() { g_profiler->startSection(m_name); }
    void ScopedProfiler::end() { g_profiler->stopSection(m_name); }
#endif
}