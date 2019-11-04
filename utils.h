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

#ifndef UTILS_H
#define UTILS_H

#include <sstream>
#include <iostream>
#include <iomanip>
#include <utility>
#include <chrono>
#include <time.h>

template<typename T>
std::string si_binary(T num, unsigned max_int_bits = 1) {
    static const int orders = 9;
    static const char* const symbols[orders] = {
        "B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"
    };

    T n = num;
    int order;
    for(order = 0; order < orders; order++) {
        if(n <= max_int_bits * 1 * 1024) {
            break;
        }
        n >>= 10;
    }

    if(order < 0)
        order = 0;
    if(order >= orders)
        order = 8;

    num >>= (10*order);
    std::stringstream ss;
    ss << num << " " << symbols[order];
    return ss.str();
}

template<typename T>
T roundToMultiple(T num, T multiple) {
    if(multiple == 0) return num;
    T remainder = num % multiple;
    if(remainder == 0) return num;

    return num + multiple - remainder;
}

// It must be: a0 <= b0 !!!
inline size_t getPairsCnt(size_t a0, size_t a1, size_t b0, size_t b1) {
    // the number of pairs is the total area minus the square part of the overlap plus the "triangle" induced by the overlap
    size_t pairscnt = (a1 - a0) * (b1 - b0); // total area
    pairscnt -= a1 <= b0 ? 0 : ((std::min(a1, b1) - b0 + 1) * (std::min(a1, b1) - b0))/2;
    return pairscnt;
}

// It must be: a0 <= b0 <= c0 !!!
inline size_t getTriplesCnt(size_t a0, size_t a1, size_t b0, size_t b1, size_t c0, size_t c1) {
    size_t triplescnt = (std::min(a1, b0) - a0) * getPairsCnt(b0, b1, c0, c1); // indepA x rest
    if (a1 <= b0)
        return triplescnt; // A has no overlaps
    if (b1 < a1) // w.l.o.g. swap rest of A and B such that a1 <= b1
        std::swap(a1, b1);
    size_t m = std::min(a1, c0);
    triplescnt += (m - b0) * getPairsCnt(m, b1, c0, c1); // (overlapAB wo. C) x ((rest of B) x C)
    triplescnt += getPairsCnt(b0, m, b0, m) * (c1 - c0); // (pairwise combs. in overlapAB wo. C) x C
    if (a1 <= c0)
        return triplescnt; // no triple overlap
    // w.l.o.g. sort (a1,b1,c1) with already a1<=b1
    if (c1 < b1) {
        std::swap(b1, c1);
        if (b1 < a1)
            std::swap(a1, b1);
    }
    triplescnt += (a1 - c0) * getPairsCnt(a1, b1, a1, c1); // overlapABC x ((rest of B) x (rest of C))
    triplescnt += getPairsCnt(c0, a1, c0, a1) * (c1 - a1); // (pairwise combs. in overlapABC) x (rest of C)
    triplescnt += ((a1 - c0) * (a1 - c0 - 1) * (a1 - c0 - 2)) / 6; // triple combs. in overlapABC
    return triplescnt;
}

inline void printProgressTime(std::ostream& out, double progress, time_t begin_timestamp) {
    using namespace std;
    out << setprecision(3) << fixed << (progress*100) << " % " << defaultfloat;

    if (progress <= 0.0) {
        out << flush;
        return;
    }

    double totsecs = (difftime(time(NULL), begin_timestamp) / progress);
    double remsecs = totsecs * (1 - progress);

    chrono::seconds ts = chrono::seconds(static_cast<long long int>(round(totsecs)));
    chrono::minutes tm(chrono::duration_cast<chrono::minutes>(ts)); // truncated
    ts -= tm;
    chrono::hours   th(chrono::duration_cast<chrono::hours>(tm)); // truncated
    tm -= th;
    chrono::duration<int,ratio<60*60*24>> td(chrono::duration_cast<chrono::duration<int,ratio<60*60*24>>>(th)); // days, truncated
    th -= td;

    chrono::seconds rs = chrono::seconds(static_cast<long long int>(round(remsecs)));
    chrono::minutes rm(chrono::duration_cast<chrono::minutes>(rs)); // truncated
    rs -= rm;
    chrono::hours   rh(chrono::duration_cast<chrono::hours>(rm)); // truncated
    rm -= rh;
    chrono::duration<int,ratio<60*60*24>> rd(chrono::duration_cast<chrono::duration<int,ratio<60*60*24>>>(rh)); // days, truncated
    rh -= rd;

    out << "\t Rem.time: ";
    if (rd.count())
        out << rd.count() << " days ";
    if (rd.count() || rh.count())
        out << rh.count() << " h ";
    if (rd.count() || rh.count() || rm.count())
        out << rm.count() << " m ";
    out << rs.count() << " s ";

    out << "\t Tot.time: ";
    if (td.count())
        out << td.count() << " days ";
    if (td.count() || th.count())
        out << th.count() << " h ";
    if (td.count() || th.count() || tm.count())
        out << tm.count() << " m ";
    out << ts.count() << " s " << flush;
}

#endif // UTILS_H
