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

#ifndef GPUENGINE_H
#define GPUENGINE_H

#include <cstdint>
#include <memory>
#include <mutex>
#include <fstream>
#include <sstream>

#include "ResultView.h"
#include "Method.h"

#define GPULOCALDBSIZE 65536 // works
//#define GPULOCALDBSIZE 131072 // out of memory error

class SNPDB;

class GPUEngine

{
public:

	// Eclipse helper
#ifndef RESULT_COMPONENTS
#define RESULT_COMPONENTS 0
#endif

    using ResultViewType = ResultView<>;
    using score_type = ResultViewType::score_type;

    GPUEngine() = delete;
    GPUEngine(const SNPDB &snpdb, int index, size_t idSize, size_t tableBufferSize, bool decomposed, bool ld, ResultView<> view, bool debug);
    ~GPUEngine(){};

    void initialize();
    void free();

    void uploadSNPData();
    void runKernel(char* source, char* results, size_t idsExpected, Method method);
    int getCUDAIndex() const { return index; }

    void setDumpStream(std::ostream &dump) {
        dumpptr = &dump;
    }

private:

    const SNPDB &snpdb;
    int index;
    size_t idSize;
    size_t tableBufferSize;
    size_t resultBufferSize;
    unsigned long dbOffset;

    char *devTables;
    char *devResults;
    ResultView<> resultView;

    bool decomposed;
    bool ld;

    bool debug;

    std::ostream *dumpptr = 0;
    static std::mutex dumpmutex;
};

#endif // GPUENGINE_H
