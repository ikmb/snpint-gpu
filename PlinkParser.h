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

#ifndef PLINKPARSER_H
#define PLINKPARSER_H

#include <vector>
#include <string>

#include "SNPDB.h"

class PlinkParser
{
public:

    explicit PlinkParser(const std::string &bim, const std::string &bed, const std::string &fam);

    void parse(SNPDB &target, size_t snp_word_size = 1);

    static const size_t input_buffer_size = 1024*1024;
private:

    std::string bimfile;
    std::string bedfile;
    std::string famfile;

    void parseFam(std::vector<SNPDB::SampleInfo> &samples);
    void parseBim(std::vector<SNPDB::SNPInfo> &snps);
    void parseBed(SNPDB &target);

};

#endif // PLINKPARSER_H
