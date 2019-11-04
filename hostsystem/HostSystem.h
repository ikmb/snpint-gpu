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

#ifndef HOSTSYSTEM_H
#define HOSTSYSTEM_H

#include <vector>

#include "GPU.h"
#include "Buffer.h"
#include "BufferFactory.h"

namespace hostsystem {

/**
 * @brief The HostSystem class represents the whole system and is meant to be used as the main entry point for all GPU activities.
 */
class HostSystem
{
public:

    /**
     * @brief HostSystem Initializes and discovers all necessary APIs.
     *
     * The default constructor performs device discovery, i.e. the GPU drivers are queried for supported
     * devices. This may take several seconds.
     */
    HostSystem(std::vector<unsigned> allowed_gpus = {});

    std::vector<GPU>& getGPUs() { return gpus; }
    GPU& getGPU(unsigned i) { return gpus.at(i); }

private:
    std::vector<GPU> gpus;
};

}
#endif // HOSTSYSTEM_H
