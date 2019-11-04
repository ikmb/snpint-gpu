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

#ifndef GPU_H
#define GPU_H

extern "C" {
#include <nvml.h>
}

#include "Device.h"

namespace hostsystem {

class GPU : Device
{
public:
    GPU(int bus, int slot);

    // Device interface
    int getBus() override;
    int getSlot() override;
    const std::string& getSerialNumber() const;

    // Convert to device indices to be used with CUDA functions
    int getIndex() const;

private:
    int bus;
    int slot;
    std::string serial;

    int cudaIndex;
    nvmlDevice_t nvmlHandle;
};

}

#endif // GPU_H
