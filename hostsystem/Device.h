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

#ifndef DEVICE_H
#define DEVICE_H

#include <string>

namespace hostsystem {

class Device {
public:
    virtual ~Device() {}

    /**
     * @brief Getter method for the logical PCI Express bus ID
     * @return the assigned bus ID
     */
    virtual int getBus() = 0;

    /**
     * @brief Getter method for the logical PCI Express slot/device ID
     * @return the assigned slot/device ID
     */
    virtual int getSlot() = 0;

};

}

#endif // DEVICE_H
