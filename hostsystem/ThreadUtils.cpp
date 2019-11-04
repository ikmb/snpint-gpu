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

// exposes the non-standard pthread_setname_np() function
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <pthread.h>

// keep this local
namespace {
    const unsigned long MAX_THREAD_NAME_LEN = 16; /* see man page PTHREAD_SETNAME_NP(3) */
}

#include <iostream>

#include <thread>
#include <cstring>

#include "ThreadUtils.h"

namespace hostsystem {

void setThreadName(const std::string& newname) {
    std::array<char, MAX_THREAD_NAME_LEN> c_newname;
    std::memcpy(c_newname.data(), newname.c_str(), std::max(newname.length(), MAX_THREAD_NAME_LEN-1));
    c_newname[MAX_THREAD_NAME_LEN-1] = 0;

    if(pthread_setname_np(pthread_self(), c_newname.data()) != 0) {

        static int messages_left = 10;
        if(messages_left-- > 0) {
            std::cerr << "Could not set thread name to " << newname << std::endl;
        } else if(messages_left-- == 0) {
            std::cerr << "This is the last time I will tell you that I could not set the thread name to " << newname << std::endl;
        }

    }
}

}
