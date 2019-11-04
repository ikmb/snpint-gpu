//    Copyright 2019 Lars Wienbrandt, Jan Christian Kässens
//
//    This file is part of SNPInt-GPU.
//
//    SNPInt-GPU is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    SNPInt-GPU is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with SNPInt-GPU.  If not, see <https://www.gnu.org/licenses/>.

#include <sstream>
#include <ostream>
#include <algorithm>

#include "Args.h"

#include <thrust/system_error.h>
#include <thrust/system/cuda/error.h>
extern "C" {
#include <cuda_runtime.h>
}

#include "Method.h"
#include "ResultView.h"
#include "SNPDB.h"
#include "GPUHandler.h"
#include "GPUEngine.h"
#include "GPUKernels.h"

static cudaDeviceProp gpuProps;

GPUEngine::GPUEngine(const SNPDB &db, int index_, size_t idsize, size_t tablebuffersize, bool decomposed_, bool ld_, ResultView<> view, bool debug_)
    : snpdb(db), index(index_), idSize(idsize), tableBufferSize(tablebuffersize), resultView(view), decomposed(decomposed_), ld(ld_), debug(debug_)
{
    const unsigned long max_tables_per_buffer = tableBufferSize / idSize;
    resultBufferSize = max_tables_per_buffer * view.getResultSize();

    dbOffset = 0;

    devTables = NULL;
    devResults = NULL;
}

void GPUEngine::initialize() {
    cudaSetDevice(index);

    if(gpuProps.multiProcessorCount == 0) {
        checkCUDAError(cudaGetDeviceProperties(&gpuProps, 0))
    }

    if (debug) {
        std::cerr << "GPU Probs: " << index << std::endl;
        std::cerr << "MPU count: " << gpuProps.multiProcessorCount << std::endl;
        std::cerr << "Max threads per block: " << gpuProps.maxThreadsPerBlock << std::endl;
        std::cerr << "Max threads per MPU: " << gpuProps.maxThreadsPerMultiProcessor << std::endl;
        std::cerr << "Shared memory per block: " << gpuProps.sharedMemPerBlock << std::endl;
        std::cerr << "Shared memory per MPU: " << gpuProps.sharedMemPerMultiprocessor << std::endl;
        std::cerr << "Total global memory: " << gpuProps.totalGlobalMem << std::endl;
        std::cerr << "Mem bus width: " << gpuProps.memoryBusWidth << std::endl;
    }

    // allocate space for the database and copy to GPU
    size_t allocsize = tableBufferSize;
    dbOffset = tableBufferSize;
    size_t mbw = gpuProps.memoryBusWidth/8; // in bytes
    if (dbOffset % mbw) // round up to multiple of mbw
        dbOffset += mbw - (dbOffset % mbw);
    size_t dbbufsize = snpdb.getBufferSize();
    if (dbbufsize % mbw) // round up to multiple of mbw
        dbbufsize += mbw - (dbbufsize % mbw);
    allocsize = dbOffset + dbbufsize; // the DB will be located behind the table memory

    checkCUDAError(cudaMalloc(&devTables, allocsize))
    checkCUDAError(cudaMalloc(&devResults, resultBufferSize))
    // we do not need more memory, so reserve half of the total mem as heap
    checkCUDAError(cudaDeviceSetLimit(cudaLimitMallocHeapSize, gpuProps.totalGlobalMem/2))

    uploadSNPData();

    if (debug) {
        int idx = -1;
        cudaGetDevice(&idx);
        std::cout << "GPU Allocations: idx: " << index << " (" << idx << ")" << std::endl;
        std::cout << "tables: @" << std::hex << (unsigned long)devTables << ", " << std::dec << tableBufferSize << " bytes" << std::endl;
        std::cout << "results: @" << std::hex << (unsigned long)devResults << ", " << std::dec << resultBufferSize << " bytes" << std::endl;
        std::cout << "DB: @" << std::hex << (unsigned long)(devTables+dbOffset) << ", " << std::dec << snpdb.getBufferSize() << " bytes" << std::endl;
    }

    unsigned long numcases = snpdb.getCaseCountPadded();
    unsigned long numsamples = snpdb.getSampleCountPadded();
    unsigned long snpsize = snpdb.getSNPSize();

    copyConstantsToDevice(numcases, numsamples, snpsize, idSize, dbOffset, decomposed, ld);

    resultView.setBuffer(devResults, resultBufferSize);
}

void GPUEngine::free() {
    cudaFree(devTables);
    cudaFree(devResults);
}

void GPUEngine::uploadSNPData() {
    // SNP DB
    char *devDB = devTables + dbOffset;
    int idx = -1;
    cudaGetDevice(&idx);
    if (debug) {
        std::cerr << "Upload SNP data.\nIDX: " << index << "(" << idx << ")" << std::endl;
        std::cerr << "DB: @" << std::hex << (unsigned long)devDB << std::dec << std::endl;
        std::cerr << "SNP size: " << snpdb.getSNPSize() << std::endl;
    }
    checkCUDAError(cudaMemcpy(devDB, snpdb.data(), snpdb.getBufferSize(), cudaMemcpyHostToDevice))

}

void GPUEngine::runKernel(char *source, char *results, size_t idsExpected, Method method) {

    int blockSize = gpuProps.maxThreadsPerBlock/4; // empirically ;-)

    size_t numBlocks = idsExpected / blockSize;
    if(idsExpected % blockSize)
        numBlocks++;

    if (debug) {
        std::cerr << "Block size: " << blockSize << std::endl;
        std::cerr << "Tables exp: " << idsExpected << std::endl;
        std::cerr << "Num blocks: " << numBlocks << std::endl;
//    size_t heapsize;
//    cudaDeviceGetLimit(&heapsize, cudaLimitMallocHeapSize);
//    std::cerr << "Heap size: " << heapsize << "\nbs*nb*2*snpsize: " << (blockSize * numBlocks * 2 * snpdb.getSNPSize()) << "\nte*2*snpsize: " << (idsExpected * 2 * snpdb.getSNPSize()) << std::endl;
    }

    dim3 grid(numBlocks, 1);
    dim3 blocks(blockSize);

    checkCUDAError(cudaMemcpy(devTables, source, idsExpected * idSize, cudaMemcpyHostToDevice))

    cudaThreadSynchronize();

    switch (method.getType()) {
    case Method::MI3:
        MutualInformationKernel3Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                     idsExpected,
                                                     resultView);
        break;
    case Method::IG3:
        InformationGainKernel3Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                   idsExpected,
                                                   resultView);
        break;
    case Method::MI2:
        MutualInformationKernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                     idsExpected,
                                                     resultView);
        break;
    case Method::IG2:
        InformationGainKernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                   idsExpected,
                                                   resultView);
        break;
    case Method::BOOST:
        BOOSTKernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                   idsExpected,
                                                   resultView);
        break;
    case Method::LOGLIN:
        LogLinearKernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                   idsExpected,
                                                   resultView);
        break;
    case Method::LOGREG:
        LogisticRegressionKernel2Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                    idsExpected,
                                                    resultView);
        break;
    case Method::LOGREG3:
        LogisticRegressionKernel3Way<<<grid, blocks, 0>>>(reinterpret_cast<uint16_t *>(devTables),
                                                   idsExpected,
                                                   resultView);
        break;
    default:
        throw std::runtime_error("Invalid method specification");
    }

    cudaThreadSynchronize();
    checkCUDAError(cudaGetLastError())
    cudaThreadSynchronize();

    checkCUDAError(cudaMemcpy(results, devResults, resultView.getResultCount() * resultView.getResultSize(), cudaMemcpyDeviceToHost))
    cudaThreadSynchronize();
}
