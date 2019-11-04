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

#include <cstring>
#include <cmath>
#include <iostream>


// for shared memory access:
extern "C" {
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
}

#include "utils.h"

#include "SNPDB.h"

SNPDB::SNPDB(unsigned minCases, unsigned minControls)
    : num_snps(0),
      num_cases(0),
      num_controls(0),
      num_cases_padded(0),
      num_controls_padded(0),
      word_size(1),
      snp_size(0),
      case_size(0),
      current_snp(0),
      current_case_offset {nullptr, 0},
      current_ctrl_offset {nullptr, 0},
      allow_buffer_resize(false),
      buffer_size(0),
      buffer(nullptr),
      min_cases(minCases),
      min_controls(minControls)
{

}

void SNPDB::initialize(unsigned long num_snps_, unsigned long num_cases_, unsigned long num_controls_, size_t word_size_) {

    this->num_cases = num_cases_;
    this->num_controls = num_controls_;
    this->word_size = word_size_;
    // Do NOT set the num_snps here, since the actual number of SNPs will be counted while parsing.
    // The argument is only required to reserve the correct buffer size.
    // this->num_snps = num_snps;
    this->num_snps = 0;
    current_snp = 0;

    // now that we have recorded the actual case/control counts, apply padding
    num_cases_padded = roundToMultiple(std::max(num_cases, min_cases), 16ul);
    num_controls_padded = roundToMultiple(std::max(num_controls, min_controls), 16ul);

    // count bytes required for cases and controls (4gts per byte)
    case_size = num_cases_padded / 4;
    size_t ctrl_size = num_controls_padded / 4;

    snp_size = case_size + ctrl_size;

    // round up to the next word_size
    if(snp_size % word_size)
        snp_size += word_size - (snp_size % word_size);

    buffer_size = snp_size * num_snps_;

    buffer = new unsigned char[buffer_size];
    std::memset(buffer, buffer_init_value, buffer_size);
    allow_buffer_resize = false;

    resetCurrOffsets();

}

SNPDB::~SNPDB() {
    if(buffer && buffer_fd == -1) {
        delete[] buffer;
    }

    if(buffer && buffer_fd != -1) {
        munmap(buffer, buffer_size);
        close(buffer_fd);
        buffer_fd = -1;
    }

    buffer = nullptr;
}



void SNPDB::finishSNP() {
    ++current_snp;
    ++num_snps;

    current_case_offset.bufptr = &buffer[current_snp * snp_size];
    current_case_offset.bit  = 0;

    current_ctrl_offset.bufptr = current_case_offset.bufptr + case_size;
    current_ctrl_offset.bit = 0;

}
