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

#ifndef BUFFER_H
#define BUFFER_H

#include <vector>
#include <list>
#include <map>
#include <mutex>
#include <memory>

#include "BufferAllocator.h"

namespace hostsystem {

template<typename T>
class BufferFactory;

template<typename T, template <typename> class Alloc>
class Buffer
{
    friend BufferFactory<Buffer<T,Alloc>>;

public:
    using value_type = T;
    using allocator = Alloc<T>;

    explicit Buffer(size_t buffer_size_)
        : buffer_size(buffer_size_), content_length(0)
    {
        buffer.resize(buffer_size_);
    }

    ~Buffer() {}

    T *getData() {
        return buffer.data();
    }

    const T *getData() const {
        return buffer.data();
    }

    void   setContentLength(size_t length) { content_length = length; }
    size_t getContentLength() const { return content_length; }

    size_t getSize() const {
        return buffer_size;
    }

private:

    const size_t buffer_size;
    std::vector<T, Alloc<T>> buffer;
    size_t content_length;
};

using IDBuffer = Buffer<char, PageAlignedAllocator>;
using CUDABuffer = Buffer<char, CUDAAllocator>;

}


#endif // BUFFER_H
