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

#ifndef GPUHANDLER_H
#define GPUHANDLER_H

#include <atomic>
#include <array>
#include <cinttypes>
#include <string>
#include <vector>

#include "hostsystem/HostSystem.h"
#include "hostsystem/ThreadUtils.h"
#include "hostsystem/Buffer.h"
#include "GPUEngine.h"
#include "ThreadPool.h"

#ifdef __CUDACC_VER_MAJOR__
/*
Since GCC 5.5, NVCC fails to import C++11's to_string into the std namespace. This is a
conformant to_string version for the global namespace, just in case we're compiling with
NVCC. Otherwise, std::to_string is imported into the global namespace. Thank you, NVIDIA.
*/

#include <sstream>
template<typename T>
std::string to_string(T val) {
    std::stringstream stream;
    stream << val;
    return stream.str();
}
#else
using std::to_string;
#endif


class GPUEngine;
class SNPDB;

class GPUHandler
{
public:

    using score_type = GPUEngine::score_type;
    using id_type = std::uint32_t;

    GPUHandler(hostsystem::HostSystem &hostsystem_, const SNPDB &db, hostsystem::BufferFactory<hostsystem::CUDABuffer> &factory, size_t idSize, size_t tableBufferSize, Method method_, ResultView<> view, bool debug_)
        : hostsystem(hostsystem_),
          method(method_),
          snpdb(db),
          engines(),
          bufferFactory(factory),
		  debug(debug_)
    {
        int num_gpus = hostsystem.getGPUs().size();
        for(int i = 0; i < num_gpus; i++)
            engines.emplace_back(snpdb, hostsystem.getGPU(i).getIndex(), idSize, tableBufferSize, method_.isDetailedComputation(), method_.isLDIncluded(), view, debug);
    }

    void distributeSNPData() {
        for (auto & e : engines)
            e.uploadSNPData();
    }

    template<typename Tin, typename Tout>
    void process(typename ThreadPool<Tin, Tout>::inqueue_type &inqueue,
                 typename ThreadPool<Tin, Tout>::outqueue_type &outqueue,
                 std::atomic_bool& termination_request,
                 int threadIndex
                 ) {

        hostsystem::setThreadName("GPU:" + to_string(engines[threadIndex].getCUDAIndex()));

        engines[threadIndex].initialize();

        while(!termination_request) {
            Tin b;
            try {
                inqueue.pop(b);
            } catch (tbb::user_abort &e) {
//            	std::cerr << "GPU thread " << threadIndex << ": User abort." << std::endl;
                break;
            }

            auto target = bufferFactory.get();

            const unsigned long tablesExpected = b->getContentLength();

            target->setContentLength(tablesExpected);

            engines[threadIndex].runKernel(b->getData(), target->getData(), tablesExpected, method);

            outqueue.push(target);
        }

        engines[threadIndex].free();
    }

    void setDumpStream(std::ostream &dump) {
        for (auto &e : engines)
            e.setDumpStream(dump);
    }

private:
    hostsystem::HostSystem &hostsystem;
    Method method;
    const SNPDB &snpdb;
    std::vector<GPUEngine> engines;
    hostsystem::BufferFactory<hostsystem::CUDABuffer> &bufferFactory;
    bool debug;
};



#endif // GPUHANDLER_H
