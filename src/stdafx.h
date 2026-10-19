// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma warning(disable: 4482)

#pragma once

#define _CRT_SECURE_NO_WARNINGS

#if !defined(__linux__)

#include "targetver.h"

#include <stdio.h>
#endif

#ifndef LS_STAG_USE_CUDA
#ifdef LS_STAG_USE_TBB
#include <tbb\task_scheduler_init.h>
#include <tbb\parallel_for.h>
#include <tbb\parallel_reduce.h>
#else
#ifdef LS_STAG_USE_CILK// 30/07/15
#include <cilk\cilk.h>
#include <cilk\reducer_opadd.h>
#include <cilk\cilk_api.h>
#else
#include <omp.h>
#endif
#endif
#ifdef LS_STAG_USE_MKL
#include <mkl.h>
#endif
#endif

#if defined(_MSC_VER)
#if _MSC_VER >= 1400
#define __restrict__ \
        __restrict
#else /* _MSC_VER >= 1400 */
#if !defined(__restrict__)
#define __restrict__
#endif
#endif
#endif

// TODO: reference additional headers your program requires here
