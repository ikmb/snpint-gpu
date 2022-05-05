#    Copyright 2019 Lars Wienbrandt, Jan Christian Kässens
#
#    This file is part of SNPInt-GPU.
#
#    SNPInt-GPU is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SNPInt-GPU is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SNPInt-GPU.  If not, see <https://www.gnu.org/licenses/>.

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED 17)

# Import local CMakeModules folder
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")


project(snpint C CXX)
cmake_minimum_required(VERSION 3.3)

execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    ${CMAKE_SOURCE_DIR}/version.h.template
    ${CMAKE_BINARY_DIR}/generated/version.h
    )

# find libraries
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)


# NVIDIA
cmake_policy(SET CMP0077 NEW)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#set(CUDA_TOOLKIT_ROOT_DIR "/opt/cuda")
#set(CUDA_BIN_PATH "/opt/cuda/bin")
#set(GPU_DEPLOYMENT_KIT_ROOT_DIR "/opt/cuda")
find_package(NVML REQUIRED)
find_package(CUDA 11.0 REQUIRED)
set(CUDA_HOST_COMPILER "/usr/bin/gcc-10")
find_package(Boost 1.65.1
    REQUIRED
    COMPONENTS
    filesystem
    program_options
)

find_package(TBB REQUIRED)

# collect linker options
set(linker_opts
  ${ZLIB_LIBRARIES}
  ${OpenMP_LIBRARIES}
  ${CMAKE_THREAD_LIBS_INIT}
  ${CUDA_LIBRARIES}
  ${NVML_LIBRARY}
#  ${TBB_LIBRARIES}
  ${Boost_LIBRARIES}
  rt
  tbb
)

# collect compiler options
set(cxx_compiler_opts
  ${OpenMP_CXX_FLAGS}
)

include_directories(
    SYSTEM "${ZLIB_INCLUDE_DIRS}"
    SYSTEM "${Boost_INCLUDE_DIRS}"
    SYSTEM "${NVML_INCLUDE_DIR}"
#    SYSTEM "${TBB_INCLUDE_DIRS}"
    "${CMAKE_BINARY_DIR}/generated"
)

set(OPT_SCORE_TYPE	double	CACHE STRING "Score precision (double or float)")

add_definitions(
    -DSCORE_TYPE=${OPT_SCORE_TYPE}
#    ${TBB_DEFINITIONS}
    -DTBB_USE_EXCEPTIONS
#    -DDEBUG
)

set(CXX_WARNINGS "-Wall -Wextra -Wshadow -Wno-deprecated-declarations -Wmissing-include-dirs -Wstrict-aliasing=2 -Wfloat-equal -Wcast-align")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_OPTIONS} ${cxx_compiler_opts} -std=c++11 ${CXX_WARNINGS} -fopenmp -march=native")

# C++ source files (DO NOT change the next line, it's a script hook)
set(cpp_sources
    main.cpp
    version.cpp
    Args.cpp
    ThreadPool.h
    MinMaxHeap.h
    Histogram.h
    IDHandler.cpp
    GPUHandler.h
    ResultHandler.h
    SNPDB.cpp
    PlinkParser.cpp
    utils.h
    Method.cpp
    ResultView.h

    # regarding host infrastructure
    hostsystem/ThreadUtils.cpp
    hostsystem/BufferAllocator.h
    hostsystem/BufferFactory.h
    hostsystem/Buffer.h
    hostsystem/DeviceCategory.cpp
    hostsystem/Device.h
    hostsystem/GPU.cpp
    hostsystem/HostSystem.cpp
    hostsystem/Spinlock.h

    # GPU-only STUFF
    GPUKernels.cu
    GPUEngine.cu
)

add_definitions (
    -D_FORCE_INLINES
    )

#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_60;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70;-gencode=arch=compute_50,code=sm_50")
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_60;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-Xcompiler;-Wall;-arch=sm_70;-gencode=arch=compute_60,code=sm_60;-gencode=arch=compute_70,code=sm_70;-gencode=arch=compute_80,code=sm_80;-gencode=arch=compute_86,code=sm_86")
set(CUDA_NVCC_FLAGS_DEBUG "${CUDA_NVCC_FLAGS_DEBUG};-O0;-g;-G")
set(CUDA_NVCC_FLAGS_RELEASE "-O3")
set(CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO};-O2;-g;-G")

cuda_add_executable(${PROJECT_NAME} ${cpp_sources})
target_link_libraries(${PROJECT_NAME} ${linker_opts})
