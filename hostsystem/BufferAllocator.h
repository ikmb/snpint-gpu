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

#ifndef CUDAALLOCATOR_H
#define CUDAALLOCATOR_H

#include <cassert>
#include <memory>
#include <cuda_runtime.h>
#include <iostream>
#include <cstdio>

extern "C" {
#include <unistd.h>
}

namespace hostsystem {

template<typename T>
struct CUDAAllocator {
    using value_type = T;

    // default constructors
    CUDAAllocator() noexcept {}
    template <class U> CUDAAllocator(const CUDAAllocator<U> &) {}

   /* static */ T* allocate(std::size_t n) {
        T *ptr;

       // check overflow of size_t
       size_t total_size;
       total_size = n * sizeof(T);

       cudaError_t status = cudaHostAlloc(&ptr, total_size, 0);
       if(status != cudaSuccess)
           throw std::bad_alloc();

       return ptr;
    }

    static void deallocate(T* ptr, std::size_t n) {
        (void) n;
        cudaFreeHost(ptr);
    }
};

template<typename T>
struct PageAlignedAllocator {
    using value_type = T;

    PageAlignedAllocator() noexcept {}
    template<class U> PageAlignedAllocator(const PageAlignedAllocator<U> &) {}

    /* static */ T* allocate(std::size_t n) {
        T* ptr;

        size_t total_size;
        total_size = n * sizeof(T);

        size_t page_size = sysconf(_SC_PAGESIZE);
        ptr = reinterpret_cast<T*>(aligned_alloc(page_size, total_size));
        if(!ptr) {
            throw std::bad_alloc();
        };

#ifdef LOCK_BUFFERS
        if(mlock(ptr, total_size) == -1) {
            throw std::bad_alloc();
        }
#endif
        return ptr;

    }

    static void deallocate(T* ptr, std::size_t n) {
#ifdef LOCK_BUFFERS
        munlock(ptr, n);
#else
        (void) n;
#endif
        std::free(ptr);
    }

};

} // namespace

#endif // CUDAALLOCATOR_H
