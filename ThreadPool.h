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

#ifndef THREADPOOL_H
#define THREADPOOL_H

#include <tbb/concurrent_queue.h>
#include <thread>
#include <vector>
#include <functional>
#include <atomic>

template<class Tin, class Tout>
class ThreadPool
{
public:
    using inqueue_type = tbb::concurrent_bounded_queue<Tin>;
    using outqueue_type = tbb::concurrent_bounded_queue<Tout>;

    using threadfunc_type = std::function<void(inqueue_type& inqueue, outqueue_type& outqueue, std::atomic_bool& terminationRequest, int threadIndex)>;

    ThreadPool(unsigned numThreads,
               inqueue_type& inqueue_,
               outqueue_type& outqueue_,
               threadfunc_type threadFunc_
               ) :
        inqueue(inqueue_),
        outqueue(outqueue_),
        threadFunc(threadFunc_),
        terminationRequest(false)
    {
        for(unsigned i = 0; i < numThreads; i++) {
            threads.emplace_back(threadFunc, std::ref(this->inqueue), std::ref(this->outqueue), std::ref(this->terminationRequest), i);
        }
    }

    ~ThreadPool() {
        terminate(false);
    }

    void terminate(bool drainInqueue = true) {
        // wait until the inqueue is drained
        while(drainInqueue && (-static_cast<int>(inqueue.size()) != static_cast<int>(threads.size()))) {
            std::this_thread::yield();
        }

        // actually cancel threads
        terminationRequest = true;
        inqueue.abort();
        for(std::thread &t : threads)
            t.join();
        threads.clear();
    }

private:
    std::vector<std::thread> threads;
    inqueue_type& inqueue;
    outqueue_type& outqueue;
    threadfunc_type threadFunc;
    std::atomic_bool terminationRequest;

};

#endif // THREADPOOL_H
