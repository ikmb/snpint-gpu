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

#include "HostSystem.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <iostream>


extern "C" {

#include <sys/types.h>
#include <dirent.h>

}


#include "BufferFactory.h"

using namespace std;


namespace hostsystem {

namespace {
struct pci_ident {
    int vendor;
    int device;
};

static bool operator==(const struct pci_ident& lhs, const struct pci_ident& rhs) {
    return (lhs.vendor == rhs.vendor) && (lhs.device == rhs.device);
}

struct pci_device {
    int bus;
    int slot;
};
}

static std::vector<struct pci_ident> pci_ident_gpu = {
    {0x10DE, 0x15F8}, // NVIDIA; Tesla P100 16G/PCIE
    {0x10DE, 0x1DB6}, // NVIDIA; Tesla V100 32G/PCIE
};

static int readHex(const std::string& file) {
    int result;
    std::ifstream f(file);
    f >> std::hex >> result;
    f.close();
    return result;
}


static std::vector<struct pci_device> findDevices(const std::vector<struct pci_ident>& idents) {
    std::vector<struct pci_device> devices;

    DIR* dir = opendir("/sys/bus/pci/devices");
    struct dirent *entry = NULL;

    while((entry = readdir(dir))) {
        if(entry->d_type == DT_LNK) {
            struct pci_ident this_device;
            std::stringstream ss;
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/vendor";

            this_device.vendor = readHex(ss.str());
            ss = std::stringstream();
            ss << "/sys/bus/pci/devices/";
            ss << entry->d_name;
            ss << "/device";
            this_device.device = readHex(ss.str());

            if(std::find(std::begin(idents), std::end(idents), this_device) != std::end(idents)) {
                int bus = 0, slot = 0;
                std::string token;
                std::istringstream name(entry->d_name);

                std::getline(name, token, ':'); // extract domain and skip it
                name >> std::hex >> bus;
                std::getline(name, token, ':'); // extract bus remainder and skip it
                name >> std::hex >> slot;


                devices.push_back( {bus, slot} );
            }
        }
    }
    return devices;
}

HostSystem::HostSystem(std::vector<unsigned> allowed_gpus)
    : gpus()
{
    std::vector<struct pci_device> devices;

    // Initialize NVML
    nvmlInit();
    devices.clear();
    devices = findDevices(pci_ident_gpu);

    // Filter GPUs
    for(const pci_device& d: devices) {
    	gpus.emplace_back(d.bus, d.slot);
    	if(find(begin(allowed_gpus), end(allowed_gpus), gpus.back().getIndex()) == end(allowed_gpus))
    		// element is not found in the allowed list -> erase
    		gpus.pop_back();
    }

}


}
