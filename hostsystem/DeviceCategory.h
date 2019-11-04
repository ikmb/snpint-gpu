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

#ifndef DEVICECATEGORY_H
#define DEVICECATEGORY_H

#include <system_error>

namespace hostsystem {

/**
 * @brief An error category to encapsulate runtime errors from the native CUDA API
 */
class GPUCategory : public std::error_category
{
public:
    GPUCategory() : std::error_category() {}

    // error_category interface
public:
    const char *name() const noexcept override;
    std::string message(int code) const override;
};

}
#endif // DEVICECATEGORY_H
