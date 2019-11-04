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

#include <array>
#include <cmath>
#include <cstring>
#include <algorithm>

#include "hostsystem/HostSystem.h"
#include "hostsystem/ThreadUtils.h"
#include "IDHandler.h"

#include "Method.h"
#include "utils.h"


IDHandler::IDHandler(
        hostsystem::HostSystem &hostsystem_,
        SNPDB &snpdb_,
        hostsystem::BufferFactory<hostsystem::IDBuffer> &factory,
        const Method &method,
        const array<snprange,3>& sets, // sets are assumed to be sorted by their left border
        const array<snprange,3>& excludesets,
        uint64_t excluderange_,
        bool quiet_,
        bool debug_)
    : hostsystem(hostsystem_), snpdb(snpdb_), order(method.getOrder()), bufferFactory(factory),
      setA(sets[0]), setB(sets[1]), setC(sets[2]),
      exsetA(excludesets[0]), exsetB(excludesets[1]), exsetC(excludesets[2]),
      excluderange(excluderange_), quiet(quiet_), debug(debug_)
{

	totalTransferBuffers = 0; // total number of buffers to transfer: will be set later

    size_t idsize = order == 2 ? IDHandler::IDSIZE_2WAY : IDHandler::IDSIZE_3WAY;
    idsPerBuffer = bufferFactory.getBufferSize() / idsize;

    initExpectedResultCount();
}

void IDHandler::process(tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &inqueue,
             tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &outqueue,
             atomic_bool& termination_request,
             int threadIndex
             ) {

    (void) inqueue;

    // This threads produces all combinations of IDs that should be tested by the GPU. The GPU then creates the tables etc.
    hostsystem::setThreadName(string("create IDs ").append(to_string(threadIndex)));

    size_t items_in_buffer = 0;
    size_t items_should_be_in_buffer_wo_filtering = 0;
    shared_ptr<hostsystem::IDBuffer> b = bufferFactory.get(); // get a transmission buffer, wait until one is ready
    unsigned *currptr = (unsigned*) (b->getData());
    snprange locexA = exsetA;
    snprange locexB = exsetB;
    snprange locexC = exsetC;
    string chrA;
    string chrB;
    string chrC;
    unsigned long bpposA = 0;
    unsigned long bpposB = 0;
    unsigned long bpposC = 0;
    if (order == 2) {
        // the smaller right border
        unsigned rightA = setA.second < setB.second ? setA.second : setB.second;

        for (unsigned snpA = setA.first; !termination_request && snpA < rightA; snpA++) {

            // "switch" the ranges if the right borders are not in order, but only if we reached the overlap
            unsigned rightB = setB.second < setA.second && snpA >= setB.first ? setA.second : setB.second;
            if (setB.second < setA.second && snpA >= setB.first)
                swap(locexA, locexB);
            // apply exclude set A: skip all combinations with this SNP if snpA is in exclude range
            // (rudimentary implementation, could be accelerated...)
            if (snpA >= locexA.first && snpA < locexA.second) { // SNP A excluded
                // add range of B loop
                items_should_be_in_buffer_wo_filtering += rightB - max((unsigned)setB.first, snpA+1);
                // adapt buffersRemaining
                if (items_should_be_in_buffer_wo_filtering >= idsPerBuffer) {
                    buffersRemaining -= items_should_be_in_buffer_wo_filtering / idsPerBuffer;
                    items_should_be_in_buffer_wo_filtering = items_should_be_in_buffer_wo_filtering % idsPerBuffer;
                }
                continue;
            }
            if (excluderange > 0) { // if an exclude range is defined
                chrA = snpdb.getSNPInfo(snpA).chromosome;
                bpposA = snpdb.getSNPInfo(snpA).pos_bp;
            }

            for (unsigned snpB = max((unsigned)setB.first, snpA+1); !termination_request && snpB < rightB; snpB++) {

                // apply exclude set B: include this pair if snpB is not in exclude range
                if (snpB >= locexB.second || snpB < locexB.first) {
                    if (excluderange > 0) { // if an exclude range is defined
                        chrB = snpdb.getSNPInfo(snpB).chromosome;
                        bpposB = snpdb.getSNPInfo(snpB).pos_bp;
                    }
                    if (excluderange == 0 || chrA.compare(chrB) != 0 || bpposB - bpposA >= excluderange) {
                        // test snpA and snpB
                        currptr = checkBuffer(b, items_in_buffer, outqueue, currptr);
                        *currptr++ = snpA;
                        *currptr++ = snpB;
                        items_in_buffer++;
                    }
                }
                if (items_should_be_in_buffer_wo_filtering == idsPerBuffer) {
                    buffersRemaining--;
                    items_should_be_in_buffer_wo_filtering = 0;
                }
                items_should_be_in_buffer_wo_filtering++;

            } // end for snpB

        } // end for snpA

        if (debug) {
            if (termination_request)
                cout << "Termination request. " << items_in_buffer << endl;
            else
                cout << "Last buffer. " << items_in_buffer << endl;
        }
        // send the last buffer
        b->setContentLength(items_in_buffer);
        outqueue.push(b);
        buffersRemaining--; // should be zero now if not user terminated
        if (debug)
            cout << "Buffers remaining: " << buffersRemaining << endl;
    } else { // order == 3

        size_t currA = setA.second;
        size_t currB = setB.second;
        size_t currC = setC.second;

        for (size_t snpA = setA.first; !termination_request && snpA < currA; snpA++) {

            // border for b in AB overlap is the larger right limit of A and B,
            // in fact we are switching A and B here w.l.o.g. but to ensure the right limits are sorted
            if (snpA >= setB.first && currB < currA){
                swap(currA, currB);
                swap(locexA, locexB);
            }
            // borders for b and c in ABC overlap are the larger right limits from all three intervals
            // in fact, we will be switching B and C now and, if necessary, A and B afterwards,
            // such that the smallest right limit will be for a now.
            // (Note that already currA <= currB holds.)
            if (snpA >= setC.first && currC < currB) {
                swap(currB, currC);
                swap(locexB, locexC);
                if (currB < currA) {
                    swap(currA, currB);
                    swap(locexA, locexB);
                }
            }
            // apply exclude set A: simply set a flag that is recognized in the inner loops, since it is difficult to correct the progress in buffersRemainig otherwise
            bool exsnpA = false;
            if (snpA >= locexA.first && snpA < locexA.second) { // SNP A excluded
                exsnpA = true;
            }
            bool switchedBC = false;
            if (excluderange > 0) { // if an exclude range is defined
                chrA = snpdb.getSNPInfo(snpA).chromosome;
                bpposA = snpdb.getSNPInfo(snpA).pos_bp;
            }

            for (size_t snpB = max(setB.first, snpA+1); !termination_request && snpB < currB; snpB++) {

                // border for c in BC overlap is the larger right limit of B and C,
                // in fact we are switching B and C here w.l.o.g. but to ensure the right limits are sorted
                if (snpB >= setC.first && currC < currB) {
                    swap(currB, currC);
                    swap(locexB, locexC);
                    switchedBC = true;
                }
                // apply exclude set A+B: skip all combinations with this SNP if snpA or snpB is in exclude range
                // (rudimentary implementation, could be accelerated...)
                if (exsnpA || (snpB >= locexB.first && snpB < locexB.second)) { // SNP A or B excluded
                    // add range of loop C
                    items_should_be_in_buffer_wo_filtering += currC - max(setC.first, snpB+1);
                    // adapt buffersRemaining
                    if (items_should_be_in_buffer_wo_filtering >= idsPerBuffer) {
                        buffersRemaining -= items_should_be_in_buffer_wo_filtering / idsPerBuffer;
                        items_should_be_in_buffer_wo_filtering = items_should_be_in_buffer_wo_filtering % idsPerBuffer;
                    }
                    continue;
                }
                if (excluderange > 0) { // if an exclude range is defined
                    chrB = snpdb.getSNPInfo(snpB).chromosome;
                    bpposB = snpdb.getSNPInfo(snpB).pos_bp;
                    // if there is an exclude range defined we could skip the complete next loop for SNP C if A and B are already in that range
                    if (chrA.compare(chrB) == 0 && bpposB - bpposA < excluderange) {
                        // add range of loop C
                        items_should_be_in_buffer_wo_filtering += currC - max(setC.first, snpB+1);
                        // adapt buffersRemaining
                        if (items_should_be_in_buffer_wo_filtering >= idsPerBuffer) {
                            buffersRemaining -= items_should_be_in_buffer_wo_filtering / idsPerBuffer;
                            items_should_be_in_buffer_wo_filtering = items_should_be_in_buffer_wo_filtering % idsPerBuffer;
                        }
                        continue;
                    }
                }

                for (size_t snpC = max(setC.first, snpB+1); !termination_request && snpC < currC; snpC++) {

                    // apply exclude set C: include this pair if snpC is not in exclude range
                    if (snpC >= locexC.second || snpC < locexC.first) {
                        if (excluderange > 0) { // if an exclude range is defined
                            chrC = snpdb.getSNPInfo(snpC).chromosome;
                            bpposC = snpdb.getSNPInfo(snpC).pos_bp;
                        }
                        if (excluderange == 0 || ( (chrB.compare(chrC) != 0 || bpposC - bpposB >= excluderange)
//                                 && (chrA.compare(chrB) != 0 || bpposB - bpposA >= excluderange) // this was already checked before snpC loop
//                                 && (chrA.compare(chrC) != 0 || bpposC - bpposA >= excluderange) // this test is not required since SNPs A,B and C are ordered by position
                         ) ) {
                            // test snpA and snpB and snpC
                            currptr = checkBuffer(b, items_in_buffer, outqueue, currptr);
                            *currptr++ = snpA;
                            *currptr++ = snpB;
                            *currptr++ = snpC;
                            items_in_buffer++;
                        }
                    }
                    if (items_should_be_in_buffer_wo_filtering == idsPerBuffer) {
                        buffersRemaining--;
                        items_should_be_in_buffer_wo_filtering = 0;
                    }
                    items_should_be_in_buffer_wo_filtering++;

                } // end for snpC

            } // end for snpB

            if (switchedBC) { // switch back B and C if they were switched in the loop before
                swap(currB, currC);
                swap(locexB, locexC);
            }

        } // end for snpA

        if (debug) {
            if (termination_request)
                cout << "Termination request. " << items_in_buffer << endl;
            else
                cout << "Last buffer. " << items_in_buffer << endl;
        }
        // send the last buffer
        b->setContentLength(items_in_buffer);
        outqueue.push(b);
        buffersRemaining--; // should be zero now if not user terminated
        if (debug)
            cout << "Buffers remaining: " << buffersRemaining << endl;
    } // end order
    if (!quiet) {
        if (termination_request)
            cout << "ID creator terminated." << endl;
        else
            cout << "ID creator finished." << endl;
    }
}


void IDHandler::initExpectedResultCount()
{

    unsigned long long results;
    if (order == 2) {
        results = getPairsCnt(setA.first, setA.second, setB.first, setB.second);
    } else {
        results = getTriplesCnt(setA.first, setA.second, setB.first, setB.second, setC.first, setC.second);
    }
    unsigned long buffer_count = results / idsPerBuffer + (results % idsPerBuffer ? 1 : 0);
    buffersRemaining = buffer_count;
    residualResults = results % idsPerBuffer;
    if (residualResults == 0ull)
        residualResults = idsPerBuffer;
    totalTransferBuffers = buffer_count;

    if (debug) {
        using namespace std;
        cout << "Expected IDs (without range exclusion!):" << endl;
        cout << results << " ID, ";
        cout << buffersRemaining << " buffers, ";
        cout << si_binary(buffersRemaining * bufferFactory.getBufferSize());
        cout << ", " << residualResults << " IDs in last buffer." << endl;
    }
}

unsigned* IDHandler::checkBuffer(shared_ptr<hostsystem::IDBuffer> &b, size_t &items_in_buffer, tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> &outqueue, unsigned *curr_ptr) {
    if (items_in_buffer == idsPerBuffer) {
        // buffer is full -> send and get a new one
        b->setContentLength(idsPerBuffer);
        outqueue.push(b); // send
        if (debug) {
            cout << "Buffer full. " << idsPerBuffer << endl;
            cout << "Buffers remaining (without range exclusion): " << buffersRemaining << endl;
        }
        b = bufferFactory.get(); // wait for a new one
        items_in_buffer = 0;
        return (unsigned*) (b->getData());
    } else
        return curr_ptr; // no change
}
