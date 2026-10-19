cmake_minimum_required(VERSION 3.18)

project(LS-STAG)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT LS-STAG)

FILE(GLOB_RECURSE headers "src/*.h")
FILE(GLOB_RECURSE sources "src/*.cpp")
FILE(GLOB_RECURSE resources "src/*.rc")
FILE(GLOB_RECURSE cu_files "src/*.cu")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}  FILES ${sources} ${cu_files} ${headers} ${resources})
add_definitions(-DCONFIG_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/configs\") 

option(LS_STAG_USE_CUDA "Use CUDA instead of OpenMP" OFF)
option(LS_STAG_USE_TBB "Use TBB instead of OpenMP" OFF)
option(LS_STAG_USE_CILK "Use Cilk instead of OpenMP" OFF)
option(LS_STAG_USE_MKL "Use MKL (if not CUDA!)" OFF)
option(LS_STAG_USE_BICGSTAB "Use BiCGStab" OFF)
option(LS_STAG_USE_DDM "Use DDM" OFF)
option(GENERAL_SOLVER "Use general ILU instead of 5/7p" OFF)
option(LS_STAG_FULL_MG "Use MG for all equations" OFF)
option(LS_STAG_FULL_ILU "Use ILU for all equations" OFF)
option(LS_STAG_PROFILING "Build with profiling" OFF)
option(ENABLE_EXTENDED_ALIGNED_STORAGE "Enable extended aligned storage" ON)

if (ENABLE_EXTENDED_ALIGNED_STORAGE)
  add_compile_definitions("_ENABLE_EXTENDED_ALIGNED_STORAGE")
else()
  add_compile_definitions("_DISABLE_EXTENDED_ALIGNED_STORAGE")
endif()
macro(add_preprocessor_def opt_name) 
  if (${opt_name})
    add_compile_definitions(${opt_name})
  endif()
endmacro()
add_preprocessor_def(LS_STAG_USE_CUDA)
add_preprocessor_def(LS_STAG_USE_TBB)
add_preprocessor_def(LS_STAG_USE_CILK)
add_preprocessor_def(LS_STAG_USE_MKL)
add_preprocessor_def(LS_STAG_USE_BICGSTAB)
add_preprocessor_def(LS_STAG_USE_DDM)
add_preprocessor_def(GENERAL_SOLVER)
add_preprocessor_def(LS_STAG_FULL_MG)
add_preprocessor_def(LS_STAG_FULL_ILU)
add_preprocessor_def(LS_STAG_PROFILING)

# Windows: if install locations not set by user, set install prefix to "LS_STAG PATH\install".
if(WIN32 AND CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${LS_STAG_PATH}/install"
  CACHE PATH "default install path" FORCE)
endif()

# --- Архитектуры CUDA (общие) ---
if(CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0)
    set(CMAKE_CUDA_ARCHITECTURES "61;70;72;75" CACHE STRING "CUDA architectures")
else()
    set(CMAKE_CUDA_ARCHITECTURES "61;70;72;75" CACHE STRING "CUDA architectures")
endif()

# --- Создание цели и настройки в зависимости от опций ---
if (LS_STAG_USE_CUDA)
    enable_language(CUDA) 
    find_package(CUDAToolkit REQUIRED) 
    add_executable(LS-STAG ${sources} ${cu_files} ${headers} ${resources})

    set_target_properties(LS-STAG PROPERTIES
    CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}"
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(LS-STAG PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/cuda/common/inc
    ${CUDAToolkit_INCLUDE_DIRS}
)

target_link_libraries(LS-STAG PRIVATE 
    CUDA::cudart 
    CUDA::cuda_driver
)

    set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} /MAP:LS_STAG_RWDI.map")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /MAP:LS_STAG.map")

    # generate a link map so we can map crash offsets back to functions
    set_target_properties(LS-STAG PROPERTIES
        LINK_FLAGS_DEBUG ""
        LINK_FLAGS_RELEASE "/MAP:LS_STAG_release.map"
        LINK_FLAGS_RELWITHDEBINFO "/MAP:LS_STAG_rwdi.map"
        LINK_FLAGS_MINSIZEREL ""
    )

target_compile_options(LS-STAG PRIVATE 
    # CUDA математические флаги (всегда)
    $<$<COMPILE_LANGUAGE:CUDA>:-fmad=true -ftz=false -prec-div=false -prec-sqrt=false -use_fast_math>
    # Включаем стандартный препроцессор для MSVC (CUDA)
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CXX_COMPILER_ID:MSVC>>:-Xcompiler=/Zc:preprocessor>
    # GPU оптимизация по конфигурациям
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Debug>>:-G>
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:Release>>:-O3 -use_fast_math>
    $<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CONFIG:RelWithDebInfo>>:-O3 -use_fast_math -G>
    # MSVC флаги для хостового кода (C++, не CUDA)
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>>:/fp:fast>
    $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CXX_COMPILER_ID:MSVC>,$<NOT:$<CONFIG:Debug>>>:/Ox>
)

elseif(NOT LS_STAG_USE_CILK AND NOT LS_STAG_USE_TBB)
    find_package(OpenMP REQUIRED)
    add_executable(LS-STAG ${sources} ${headers} ${resources})

    set_target_properties(LS-STAG PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
    )

    target_link_libraries(LS-STAG PRIVATE OpenMP::OpenMP_CXX)

        target_compile_options(LS-STAG PRIVATE 
        $<$<CXX_COMPILER_ID:MSVC>:
            /fp:fast
            /openmp
            /wd4849                      # подавить предупреждение C4849 (collapse игнорируется)
            $<$<NOT:$<CONFIG:Debug>>:/Ox>
        >
        # ---- AddressSanitizer только в Debug ----
        $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Debug>>:/fsanitize=address>
    )

    # Необходимо также передать флаг линковщику
    target_link_options(LS-STAG PRIVATE 
        $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Debug>>:/fsanitize=address>
        $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/MAP:LS_STAG.map>
    )

else()
    # Для TBB или Cilk можно добавить соответствующие блоки позже
    message(FATAL_ERROR "TBB/Cilk not yet supported in this configuration. Disable those options.")
endif()