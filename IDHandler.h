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

#ifndef IDHANDLER_H
#define IDHANDLER_H

#include <atomic>
#include <memory>
#include <chrono>
#include <iostream>
#include <utility>
#include <time.h>

#include "SNPDB.h"
#include "Method.h"
#include "ThreadPool.h"
#include "hostsystem/Buffer.h"
#include "hostsystem/HostSystem.h"
#include "hostsystem/ThreadUtils.h"

using namespace std;

class IDHandler
{
public:

    IDHandler(
            hostsystem::HostSystem &hostsystem,
            SNPDB &snpdb,
            hostsystem::BufferFactory<hostsystem::IDBuffer> &factory,
            const Method& method,
            const array<snprange,3>& sets,
            const array<snprange,3>& excludesets,
            uint64_t excluderange,
            bool quiet,
            bool debug);

    void process(tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &inqueue,
                 tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &outqueue,
                 atomic_bool& termination_request,
                 int threadIndex
                 );

    bool isFinished() {
    	return buffersRemaining == 0;
    }

    // shows the progress of the ID buffers between 0 (not started) and 1 (finished)
    double progress() {
    	return 1.0 - ((double)buffersRemaining / totalTransferBuffers);
    }

    static constexpr size_t IDSIZE_2WAY = 8;
    static constexpr size_t IDSIZE_3WAY = 12;

private:
    hostsystem::HostSystem &hostsystem;
    SNPDB &snpdb;
    unsigned order;
    hostsystem::BufferFactory<hostsystem::IDBuffer> &bufferFactory;

    unsigned long buffersRemaining;
    unsigned long residualResults;

    unsigned long totalTransferBuffers;
    unsigned long idsPerBuffer;

    void initExpectedResultCount();
    unsigned* checkBuffer(shared_ptr<hostsystem::IDBuffer> &b, size_t &items_in_buffer, tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &outqueue, unsigned *curr_ptr);

    snprange setA;
    snprange setB;
    snprange setC;
    snprange exsetA;
    snprange exsetB;
    snprange exsetC;
    uint64_t excluderange;

    bool quiet;
    bool debug;
};
#endif // IDHANDLER_H
