/*
* \copyright (C) 2024 Puzikova V.V.
* \file      ls_stag_parallel.h
* \author    Puzikova Valeria (vvp)
* \date      02.01.2024 13:21:15
* \brief     ls_stag_parallel.h
* \fmtstyle  3S
*/

#ifndef _LS_STAG_LS_STAG_PARALLEL_H_
#define _LS_STAG_LS_STAG_PARALLEL_H_

// ================================================================== INCLUDES
#ifdef LS_STAG_USE_CUDA
#include <helper_math.h>
#include <helper_cuda.h>
#include "device_launch_parameters.h"
#endif
// ================================================================== SYNOPSIS
/*!
\file
\brief Заголовочный файл с используемыми для распараллеливания макросами, функциями и т.д.

Данный файл содержит в себе определения используемых для распараллеливания макросами, функциями.
Директивы препроцессора:
LS_STAG_USE_OMP
LS_STAG_USE_TBB
LS_STAG_USE_CILK
LS_STAG_USE_CUDA
*/
// ===========================================================================

namespace LS_STAG
{
#if LS_STAG_USE_CUDA
	#define LS_STAG_CUDA_BLOCKS 256
	#define LS_STAG_CUDA_BLOCKS_2D_X 32
	#define LS_STAG_CUDA_BLOCKS_2D_Y  8
	#define ELEMENTS_PER_THREAD 8

	#ifdef  __INTELLISENSE__
	#define CUDA_KERNEL(...)
	#else
	#define CUDA_KERNEL(...) <<< __VA_ARGS__ >>>
	#endif

#if defined(__CUDACC__) || !defined(__host__)
#define __host__ \
        __location__(host)
#endif /* defined(__CUDACC__) || !defined(__host__) */
#if defined(__CUDACC__) || !defined(__device__)
#define __device__ \
        __location__(device)
#endif /* defined(__CUDACC__) || !defined(__device__) */
#if defined(__CUDACC__) || !defined(__global__)
#define __global__ \
        __location__(global)
#endif /* defined(__CUDACC__) || !defined(__global__) */
#if defined(__CUDACC__) || !defined(__shared__)
#define __shared__ \
        __location__(shared)
#endif /* defined(__CUDACC__) || !defined(__shared__) */
#if defined(__CUDACC__) || !defined(__constant__)
#define __constant__ \
        __location__(constant)
#endif /* defined(__CUDACC__) || !defined(__constant__) */
#if defined(__CUDACC__) || !defined(__managed__)
#define __managed__ \
        __location__(managed)
#endif /* defined(__CUDACC__) || !defined(__managed__) */

	#define LS_STAG_SINGLE_KERNEL(func_name, ...) static __device__ void _k_s_##func_name(__VA_ARGS__) 

	#define LS_STAG_KERNEL_R(type, func_name, ...) static __device__ type _k_##func_name(const int i, __VA_ARGS__)
#define LS_STAG_KERNEL_2D_R(type, func_name, ...) static __device__ type _k_##func_name(const int i, const int j, __VA_ARGS__)

	#define EXPAND(x) x
	#define FOR_EACH_2(what, x, y) what(x, y)
	#define FOR_EACH_4(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_2(what, __VA_ARGS__))
	#define FOR_EACH_6(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_4(what, __VA_ARGS__))
	#define FOR_EACH_8(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_6(what, __VA_ARGS__))
	#define FOR_EACH_10(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_8(what, __VA_ARGS__))
	#define FOR_EACH_12(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_10(what, __VA_ARGS__))
	#define FOR_EACH_14(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_12(what, __VA_ARGS__))
	#define FOR_EACH_16(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_14(what, __VA_ARGS__))
	#define FOR_EACH_18(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_16(what, __VA_ARGS__))
	#define FOR_EACH_20(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_18(what, __VA_ARGS__))
	#define FOR_EACH_22(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_20(what, __VA_ARGS__))
	#define FOR_EACH_24(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_22(what, __VA_ARGS__))
	#define FOR_EACH_26(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_24(what, __VA_ARGS__))
	#define FOR_EACH_28(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_26(what, __VA_ARGS__))
	#define FOR_EACH_30(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_28(what, __VA_ARGS__))
	#define FOR_EACH_32(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_30(what, __VA_ARGS__))
	#define FOR_EACH_34(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_32(what, __VA_ARGS__))
	#define FOR_EACH_36(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_34(what, __VA_ARGS__))
	#define FOR_EACH_38(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_36(what, __VA_ARGS__))
	#define FOR_EACH_40(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_38(what, __VA_ARGS__))
	#define FOR_EACH_42(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_40(what, __VA_ARGS__))
	#define FOR_EACH_44(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_42(what, __VA_ARGS__))
	#define FOR_EACH_46(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_44(what, __VA_ARGS__))
	#define FOR_EACH_48(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_46(what, __VA_ARGS__))
	#define FOR_EACH_50(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_48(what, __VA_ARGS__))
	#define FOR_EACH_52(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_50(what, __VA_ARGS__))
	#define FOR_EACH_54(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_52(what, __VA_ARGS__))
	#define FOR_EACH_56(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_54(what, __VA_ARGS__))
	#define FOR_EACH_58(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_56(what, __VA_ARGS__))
	#define FOR_EACH_60(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_58(what, __VA_ARGS__))
	#define FOR_EACH_62(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_60(what, __VA_ARGS__))
	#define FOR_EACH_64(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_62(what, __VA_ARGS__))
	#define FOR_EACH_66(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_64(what, __VA_ARGS__))
	#define FOR_EACH_68(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_66(what, __VA_ARGS__))
	#define FOR_EACH_70(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_68(what, __VA_ARGS__))
	#define FOR_EACH_72(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_70(what, __VA_ARGS__))
	#define FOR_EACH_74(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_72(what, __VA_ARGS__))
	#define FOR_EACH_76(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_74(what, __VA_ARGS__))
	#define FOR_EACH_78(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_76(what, __VA_ARGS__))
	#define FOR_EACH_80(what, x, y, ...) what(x, y), EXPAND(FOR_EACH_78(what, __VA_ARGS__))

	#define FOR_EACH_NARG(...) FOR_EACH_NARG_(__VA_ARGS__, FOR_EACH_RSEQ_N())
	#define FOR_EACH_NARG_(...) EXPAND(FOR_EACH_ARG_N(__VA_ARGS__))

	#define FOR_EACH_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22,\
				_23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40, _41, _42, _43, _44,\
				_45, _46, _47, _48, _49, _50, _51, _52, _53, _54, _55, _56, _57, _58, _59, _60, _61, _62, _63, _64, _65, _66, \
				_67, _68, _69, _70, _71, _72, _73, _74, _75, _76, _77, _78, _79, _80, N, ...) N
	#define FOR_EACH_RSEQ_N() 80, 79, 78, 77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57,\
				56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32,\
				31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9,\
							8, 7, 6, 5, 4, 3, 2, 1, 0
	#define CONCATENATE(x,y) x##y
	#define FOR_EACH_(N, what, ...) EXPAND(CONCATENATE(FOR_EACH_, N)(what, __VA_ARGS__))

	#define FOR_EACH(what, ...) FOR_EACH_(FOR_EACH_NARG(__VA_ARGS__), what, __VA_ARGS__)

	#define FUNC_ARGS(type, name) type name 
	#define FUNC_PASS_ARGS(type, name) name  

	#define LS_STAG_DECLARATION_KERNEL_BASE_R(type, prefix, func_name, ...)\
		type _##prefix##_##func_name(FOR_EACH(FUNC_ARGS, __VA_ARGS__))

	#define LS_STAG_DECLARATION_KERNEL_BASE(prefix, func_name, ...)\
		LS_STAG_DECLARATION_KERNEL_BASE_R(void, prefix, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_BASE(prefix, func_name, ...)\
		__global__ LS_STAG_DECLARATION_KERNEL_BASE(prefix, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARE_KERNEL(func_name, ...)\
		LS_STAG_DECLARATION_KERNEL_BASE(c, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARE_G_KERNEL(func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL_BASE(g, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL(func_name, ...)\
		LS_STAG_DECLARE_KERNEL(func_name, const int, size, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL_2D(func_name, ...)\
		LS_STAG_DECLARE_KERNEL(func_name, const int, size, const int, size2, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL(func_name, ...)\
		LS_STAG_DECLARE_G_KERNEL(func_name, const int, size, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_2D(func_name, ...)\
		LS_STAG_DECLARE_G_KERNEL(func_name, const int, size, const int, size2, __VA_ARGS__)

	#define LS_STAG_DECLARATION_R_KERNEL(func_name, ...)\
		LS_STAG_DECLARE_KERNEL(func_name, const int, start, const int, end, __VA_ARGS__)

	#define LS_STAG_DECLARATION_R_KERNEL_2D(func_name, ...)\
		LS_STAG_DECLARE_KERNEL(func_name, const int, start, const int, end, const int, start2, const int, end2, __VA_ARGS__)

	#define LS_STAG_DECLARATION_GR_KERNEL(func_name, ...)\
		LS_STAG_DECLARE_G_KERNEL(func_name, const int, start, const int, end, __VA_ARGS__)

	#define LS_STAG_DECLARATION_GR_KERNEL_2D(func_name, ...)\
		LS_STAG_DECLARE_G_KERNEL(func_name, const int, start, const int, end, const int, start2, const int, end2, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL_R(type, func_name, ...)\
		LS_STAG_DECLARATION_KERNEL_BASE_R(type, c, func_name, const int, size, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_R(type, func_name, ...)\
		LS_STAG_DECLARE_G_KERNEL(func_name, const int, size, type*, res, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL_RT(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_KERNEL_R(type, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_RT(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_G_KERNEL_R(type, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_S_KERNEL(func_name, ...)\
		LS_STAG_DECLARATION_KERNEL_BASE(c_s, func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_S_KERNEL(func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL_BASE(g_s, func_name, __VA_ARGS__)

	#define LS_STAG_G_S_KERNEL(class_name, func_name, ...)\
		LS_STAG_DECLARATION_G_S_KERNEL(func_name, __VA_ARGS__)\
		{\
			class_name::_k_s_##func_name(FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
		}\
		LS_STAG_DECLARATION_S_KERNEL(func_name, __VA_ARGS__)\
		{\
			_g_s_##func_name<<<1,1>>>(FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
			getLastCudaError("Single kernel failed!");\
		}

	#define LS_STAG_CALL_G_KERNEL_BASE(func_name, size, ...)\
	{\
		unsigned int nt = (LS_STAG_CUDA_BLOCKS < (size)) ?  LS_STAG_CUDA_BLOCKS : (size);\
		unsigned int nb = (size + nt - 1) / nt;\
        _g_##func_name<<<nb,nt>>>(__VA_ARGS__);\
        getLastCudaError("Kernel failed!");\
	}

	#define LS_STAG_CALL_G_KERNEL_2D_BASE(func_name, size, size2, ...)\
{\
    dim3 block(LS_STAG_CUDA_BLOCKS_2D_X, LS_STAG_CUDA_BLOCKS_2D_Y);\
    dim3 grid((size + block.x - 1) / block.x, (size2 + block.y - 1) / block.y);\
    _g_##func_name<<<grid, block>>>(__VA_ARGS__);\
    getLastCudaError("2D Kernel failed!");\
}

	#define LS_STAG_CALL_G_KERNEL(func_name, ...)\
		LS_STAG_CALL_G_KERNEL_BASE(func_name, size, size, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__))

	#define LS_STAG_CALL_G_KERNEL_2D(func_name, ...)\
		LS_STAG_CALL_G_KERNEL_2D_BASE(func_name, size, size2, size, size2, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__))

	#define LS_STAG_CALL_GR_KERNEL(func_name, ...)\
		LS_STAG_CALL_G_KERNEL_BASE(func_name, end - start, start, end, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__))

	#define LS_STAG_CALL_GR_KERNEL_2D(func_name, ...)\
		LS_STAG_CALL_G_KERNEL_2D_BASE(func_name, end - start, end2 - start2, start, end, start2, end2, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__))

	#define LS_STAG_CALL_KERNEL_BASE(prefix, func_name, iEnd, iStart, ...)\
		{\
			int i = blockIdx.x * blockDim.x + threadIdx.x + iStart;\
			if (i < iEnd) prefix _k_##func_name(i, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
		}
    
	#define LS_STAG_CALL_KERNEL_2D_BASE(prefix, func_name, iEnd, iStart, iEnd2, iStart2, ...)\
    {\
        int i = blockIdx.x * blockDim.x + threadIdx.x + iStart;\
        int j = blockIdx.y * blockDim.y + threadIdx.y + iStart2;\
        if (i < iEnd && j < iEnd2) prefix _k_##func_name(i, j, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
    }

	#define LS_STAG_CALL_KERNEL(class_name, func_name, ...)\
		LS_STAG_CALL_KERNEL_BASE(class_name::, func_name, size, 0, __VA_ARGS__)

	#define LS_STAG_CALL_KERNEL_2D(class_name, func_name, ...)\
		LS_STAG_CALL_KERNEL_2D_BASE(class_name::, func_name, size, 0, size2, 0, __VA_ARGS__)

	#define LS_STAG_CALL_R_KERNEL(class_name, func_name, ...)\
		LS_STAG_CALL_KERNEL_BASE(class_name::, func_name, end, start, __VA_ARGS__)

	#define LS_STAG_CALL_R_KERNEL_2D(class_name, func_name, ...)\
		LS_STAG_CALL_KERNEL_2D_BASE(class_name::, func_name, end, start, end2, start2, __VA_ARGS__)

	#define LS_STAG_CALL_KERNEL_R(func_name, ...)\
		LS_STAG_CALL_KERNEL_BASE(_r, func_name, size, 0, __VA_ARGS__)

	#define LS_STAG_G_KERNEL(class_name, func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL(func_name, __VA_ARGS__)\
		LS_STAG_CALL_KERNEL(class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL(func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL(func_name, __VA_ARGS__)
		
	#define LS_STAG_G_KERNEL_2D(class_name, func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL_2D(func_name, __VA_ARGS__)\
		LS_STAG_CALL_KERNEL_2D(class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL_2D(func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL_2D(func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL_T(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_KERNEL(func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_T(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_G_KERNEL(func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_KERNEL_2D_T(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_KERNEL_2D(func_name, __VA_ARGS__)

	#define LS_STAG_DECLARATION_G_KERNEL_2D_T(type, func_name, ...)\
		template<typename type>\
		LS_STAG_DECLARATION_G_KERNEL_2D(func_name, __VA_ARGS__)

	#define LS_STAG_G_KERNEL_2D_T(type, class_name, func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL_2D_T(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_KERNEL_2D(class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL_2D_T(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL_2D(func_name, __VA_ARGS__)
		
	#define LS_STAG_G_KERNEL_T(type, class_name, func_name, ...)\
		LS_STAG_DECLARATION_G_KERNEL_T(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_KERNEL(class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL_T(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL(func_name, __VA_ARGS__)

	#define LS_STAG_GR_KERNEL(class_name, func_name, ...)\
		LS_STAG_DECLARATION_GR_KERNEL(func_name, __VA_ARGS__)\
		LS_STAG_CALL_R_KERNEL(class_name, func_name,  __VA_ARGS__)\
		LS_STAG_DECLARATION_R_KERNEL(func_name, __VA_ARGS__)\
		LS_STAG_CALL_GR_KERNEL(func_name, __VA_ARGS__)

	#define LS_STAG_GR_KERNEL_2D(class_name, func_name, ...)\
		LS_STAG_DECLARATION_GR_KERNEL_2D(func_name, __VA_ARGS__)\
		LS_STAG_CALL_R_KERNEL_2D(class_name, func_name,  __VA_ARGS__)\
		LS_STAG_DECLARATION_R_KERNEL_2D(func_name, __VA_ARGS__)\
		LS_STAG_CALL_GR_KERNEL_2D(func_name, __VA_ARGS__)

	#define LS_STAG_R_KERNEL(type, class_name, func_name, ...)\
		__global__ void _g_##func_name(const int size, LS_STAG_R_P(type) res, FOR_EACH(FUNC_ARGS, __VA_ARGS__))\
		{/* each block has it's own shared array filled with results of a i * b i.*/ \
			const int totalThreads = gridDim.x * blockDim.x;\
			__shared__ type sharedChache[LS_STAG_CUDA_BLOCKS]; /* ??? */ \
			const int tid = blockIdx.x * blockDim.x + threadIdx.x;\
			const int elementsPerThread = (size + totalThreads - 1) / totalThreads; \
			type threadSum = 0;\
			for (int i = tid * elementsPerThread, end = min(size, i + elementsPerThread); i < end; ++i) {\
				threadSum += class_name::_k_##func_name(i, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
			}\
			sharedChache[threadIdx.x] = threadSum; \
			__syncthreads();\
			if(threadIdx.x == 0) { threadSum = 0; for(int i = 0; i < blockDim.x; i++) threadSum += sharedChache[i]; atomicAdd(&res[0], threadSum); }\
		}

	#define LS_STAG_CALL_G_KERNEL_R(type, func_name, ...)\
	{\
		unsigned int nt = (LS_STAG_CUDA_BLOCKS < size) ?  LS_STAG_CUDA_BLOCKS : size;\
		unsigned int desiredThreads = (size + ELEMENTS_PER_THREAD - 1) / ELEMENTS_PER_THREAD;\
		unsigned int nb = (desiredThreads + nt - 1) / nt; \
		static MeshFunction<type> accGPU; static type* m_accGPU; accGPU.init(1, m_accGPU); accGPU.set(0); \
         _g_##func_name<<<nb,nt>>>(size, m_accGPU, FOR_EACH(FUNC_PASS_ARGS, __VA_ARGS__));\
        getLastCudaError("Kernel failed!");\
		return accGPU.getElem(0);\
	}

	#define LS_STAG_G_KERNEL_R(type, class_name, func_name, ...)\
		LS_STAG_R_KERNEL(type, class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL_R(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL_R(type, func_name, __VA_ARGS__)

	#define LS_STAG_G_KERNEL_RT(type, class_name, func_name, ...)\
		template<typename type>\
		LS_STAG_R_KERNEL(type, class_name, func_name, __VA_ARGS__)\
		LS_STAG_DECLARATION_KERNEL_RT(type, func_name, __VA_ARGS__)\
		LS_STAG_CALL_G_KERNEL_R(type, func_name, __VA_ARGS__)

	#define LS_STAG_PAR_SUM(class_name, func_name, type, size, ...) type acc = _c_##func_name(size, __VA_ARGS__); 

	#define LS_STAG_PAR_FOR_RANGE(class_name, func_name, start, end, ...) _c_##func_name(start, end, __VA_ARGS__);

	#define LS_STAG_PAR_FOR(class_name, func_name, size, ...) _c_##func_name(size, __VA_ARGS__);

	#define LS_STAG_PAR_FOR_2D_RANGE(class_name, func_name, start, end, start2, end2, ...) _c_##func_name(start, end, start2, end2, __VA_ARGS__);

	#define LS_STAG_PAR_FOR_2D(class_name, func_name, size, size2, ...) _c_##func_name(size, size2, __VA_ARGS__);

	#define LS_STAG_SINGLE(class_name, func_name, ...) _c_s_##func_name(__VA_ARGS__); 

	#define LS_STAG_INLINE __inline__ __host__ __device__
#else
	#define LS_STAG_INLINE inline
	#ifdef LS_STAG_USE_TBB
		#define LS_STAG_PAR_START_FOR_RANGE(start, end)\                                                           
		tbb::parallel_for(tbb::blocked_range<int>(start, end), \
			[=](const tbb::blocked_range<int>& r)\
		{\
		for (int i = r.begin(), m_end = r.end(); i < m_end; i++) {
			
		#define LS_STAG_PAR_START_FOR_2D_RANGE(start, end, start2, end2)\                                                           
		tbb::parallel_for(tbb::blocked_range<int>(start, end), \
			[=](const tbb::blocked_range<int>& r)\
		{\
		for (int i = r.begin(), m_end = r.end(); i < m_end; i++) {\
			for(int j = start2; j < end2; j++) {

		#define LS_STAG_PAR_SUM(class_name, func_name, type, size, ...)\
				type acc = tbb::parallel_reduce(tbb::blocked_range<int>(0, size), type(0),\
								[=](const tbb::blocked_range<int>& r, type acc2)->type\
								{\
									for (int i = r.begin(), m_end = r.end(); i < m_end; i++)\
										acc2 += class_name::_k_##func_name(i, __VA_ARGS__);\
									return acc2;\
								},\
								[](type x, type y)->type {return x + y; } );

		#define LS_STAG_END_PAR_FOR } } );
		#define LS_STAG_END_PAR_FOR_2D } } } );
	#else
		#define LS_STAG_END_PAR_FOR }
		#define LS_STAG_END_PAR_FOR_2D } }

		#ifdef LS_STAG_USE_CILK
			#define LS_STAG_PAR_START_FOR_RANGE(start, end)\
							cilk_for(int i = start; i < end; i++) {
								
			#define LS_STAG_PAR_START_FOR_2D_RANGE(start, end, start2, end2)\
							cilk_for(int i = start; i < end; i++) {\
							  cilk_for(int j = start2; j < end2; j++) {

			#define LS_STAG_PAR_SUM(class_name, func_name, type, size, ...)\
				cilk::reducer< cilk::op_add<type> > acc(0);\
				cilk_for(int i = 0; i < size; i++)\
					*acc += class_name::_k_##func_name(i, __VA_ARGS__);
		#else
#ifdef _MSC_VER
#define PRAGMA_OMP(x) __pragma(x)
#else
#define PRAGMA_OMP(x) _Pragma(#x)
#endif

#define LS_STAG_PAR_START_FOR_RANGE(start, end)         \
    PRAGMA_OMP(omp parallel for)                        \
    for(int i = start; i < end; i++) {

#define LS_STAG_PAR_START_FOR_2D_RANGE(start, end, start2, end2) \
    PRAGMA_OMP(omp parallel for collapse(2))                     \
    for(int i = start; i < end; i++) {                           \
        for(int j = start2; j < end2; j++) {

#define LS_STAG_PAR_SUM(class_name, func_name, type, size, ...)   \
    type acc = 0;                                                 \
    PRAGMA_OMP(omp parallel for reduction(+:acc))                 \
    for(int i = 0; i < size; i++)                                 \
        acc += class_name::_k_##func_name(i, __VA_ARGS__);
		#endif
	#endif

	#define LS_STAG_SINGLE_KERNEL(func_name, ...) static void _k_s_##func_name(__VA_ARGS__)

	#define LS_STAG_KERNEL_R(type, func_name, ...) static type _k_##func_name(const int i, __VA_ARGS__)
#define LS_STAG_KERNEL_2D_R(type, func_name, ...) static type _k_##func_name(const int i, const int j , __VA_ARGS__)

	#define LS_STAG_PAR_FOR_RANGE(class_name, func_name, start, end, ...)\
		LS_STAG_PAR_START_FOR_RANGE(start, end)\
			class_name::_k_##func_name(i, __VA_ARGS__);\
		LS_STAG_END_PAR_FOR

	#define LS_STAG_PAR_FOR_2D_RANGE(class_name, func_name, start, end, start2, end2, ...)\
		LS_STAG_PAR_START_FOR_2D_RANGE(start, end, start2, end2)\
			class_name::_k_##func_name(i, j, __VA_ARGS__);\
		LS_STAG_END_PAR_FOR_2D

	#define LS_STAG_PAR_FOR(class_name, func_name, size, ...) LS_STAG_PAR_FOR_RANGE(class_name, func_name, 0, size, __VA_ARGS__)

	#define LS_STAG_PAR_FOR_2D(class_name, func_name, size, size2, ...) LS_STAG_PAR_FOR_2D_RANGE(class_name, func_name, 0, size, 0, size2, __VA_ARGS__)

	#define LS_STAG_SINGLE(class_name, func_name, ...) class_name::_k_s_##func_name(__VA_ARGS__);

	#define LS_STAG_G_S_KERNEL(class_name, func_name, ...)
	#define LS_STAG_G_KERNEL(class_name, func_name, ...) 
	#define LS_STAG_G_KERNEL_T(type, class_name, func_name, ...) 
	#define LS_STAG_GR_KERNEL(class_name, func_name, ...) 
	#define LS_STAG_G_KERNEL_R(type, class_name, func_name, ...)
	#define LS_STAG_G_KERNEL_RT(type, class_name, func_name, ...)
	#define LS_STAG_DECLARATION_S_KERNEL(func_name, ...)
	#define LS_STAG_DECLARATION_KERNEL(func_name, ...) 
	#define LS_STAG_DECLARATION_KERNEL_T(type, func_name, ...) 
	#define LS_STAG_DECLARATION_R_KERNEL(func_name, ...) 
	#define LS_STAG_DECLARATION_KERNEL_R(type, func_name, ...)
	#define LS_STAG_DECLARATION_KERNEL_RT(type, func_name, ...)
	
	#define LS_STAG_G_KERNEL_2D(class_name, func_name, ...) 
	#define LS_STAG_G_KERNEL_2D_T(type, class_name, func_name, ...)
	#define LS_STAG_GR_KERNEL_2D(class_name, func_name, ...)  
	#define LS_STAG_DECLARATION_R_KERNEL_2D(func_name, ...) 
	#define LS_STAG_DECLARATION_KERNEL_2D(func_name, ...) 
	#define LS_STAG_DECLARATION_KERNEL_2D_T(type, func_name, ...) 

#endif

#ifdef LS_STAG_USE_CILK
	#define LS_STAG_PAR_ACC_VALUE acc.get_value()
#else
	#define LS_STAG_PAR_ACC_VALUE acc
#endif

#define LS_STAG_KERNEL_INLINE static LS_STAG_INLINE

#define LS_STAG_VIRTUAL_INLINE_CALL(class_name, func_name, ...) class_name::func_name(__VA_ARGS__)

#define LS_STAG_RETURN_PAR_SUM(class_name, func_name, type, size, ...)\
	LS_STAG_PAR_SUM(class_name, func_name, type, size, __VA_ARGS__)\
	return LS_STAG_PAR_ACC_VALUE;

#define LS_STAG_KERNEL(func_name, ...) LS_STAG_KERNEL_R(void, func_name, __VA_ARGS__)
#define LS_STAG_KERNEL_2D(func_name, ...) LS_STAG_KERNEL_2D_R(void, func_name, __VA_ARGS__)

#define LS_STAG_R_P(type) type* __restrict__
#define LS_STAG_CR_P(type) const LS_STAG_R_P(type)
#define LS_STAG_R_PD LS_STAG_R_P(double)
#define LS_STAG_CR_PD LS_STAG_CR_P(double)
#define LS_STAG_R_PI LS_STAG_R_P(int)
#define LS_STAG_CR_PI LS_STAG_CR_P(int)
#define LS_STAG_R_PCT LS_STAG_R_P(CellTypes)
#define LS_STAG_CR_PCT LS_STAG_CR_P(CellTypes)
#define LS_STAG_CR_PPI LS_STAG_CR_P(tripletInt)
#define LS_STAG_R_PC LS_STAG_R_P(coordinate)
#define LS_STAG_CR_PC LS_STAG_CR_P(coordinate)

}

#endif //_LS_STAG_LS_STAG_PARALLEL_H_
