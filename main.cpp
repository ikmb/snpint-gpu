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

/* STL */
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

/* C wrappers */
#include <csignal>

/* CUDA */
#include <cuda_runtime.h>

/* Locals */
#include "Args.h"
#include "ThreadPool.h"
#include "IDHandler.h"
#include "GPUHandler.h"
#include "ResultHandler.h"
#include "ResultView.h"
#include "SNPDB.h"
#include "PlinkParser.h"
#include "utils.h"
#include "Method.h"
#include "hostsystem/HostSystem.h"
#include "hostsystem/ThreadUtils.h"
#include "GPUKernels.h"

using namespace std;
using namespace placeholders;

using ResultViewSpec = ResultView<>;

using ctable2way = array<uint32_t, 20>; // 18 + ID
using ctable3way = array<uint32_t, 57>; // 54 + ID

static atomic<bool> user_terminate;

/**
 * @brief Signal handler registered for SIGINT reception
 *
 * This function is registered for called by the operating system if a SIGINT
 * is received (i.e. Ctrl+C). If we didn't catch it here, the program would
 * be terminated without executing the cleanup callbacks that have been registered
 * using atexit.
 *
 * @param signal The actual signal that has been received
 */
static void signalHandler(int signal) {
    if(signal == SIGINT) {
        if(!user_terminate) {
            user_terminate = true;
            cout << "\nTrying to shutdown gracefully, please be patient. Hit Ctrl+C again to terminate." << endl;
        } else {
            cout << "\nTerminating..." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

static void verifyData(const SNPDB &db, bool useGPU, bool debug) {

    stringstream errmsg;
    unsigned long cases = db.getCaseCountPadded();
    unsigned long controls = db.getControlCountPadded();
    unsigned long samples = db.getSampleCountPadded();

    if(db.getSNPSize() * 2 > GPULOCALDBSIZE && useGPU) {
        errmsg << "GPU only run supports at most " << (GPULOCALDBSIZE*2) << " samples in total (got " << db.getSampleCountPadded() << ", including padding).";
        throw runtime_error(errmsg.str());
    }

    if(debug) {
        cout << "Padded data set consists of " << cases << " cases and " << controls << " controls (" << samples << " samples in total)." << endl;
    }
}

static size_t getGPUResultSize(const Method& method) {
    size_t gpu_result_size = method.getOrder() * sizeof(ResultView<>::id_type);
    if(gpu_result_size % ResultView<>::getAlignment() != 0)
        gpu_result_size += (ResultView<>::getAlignment() - (gpu_result_size % ResultView<>::getAlignment()));

    gpu_result_size += method.getScoreFields() * sizeof(ResultView<>::score_type);
    if(gpu_result_size % ResultView<>::getAlignment() != 0)
        gpu_result_size += (ResultView<>::getAlignment() - (gpu_result_size % ResultView<>::getAlignment()));

    return gpu_result_size;
}

int main(int argc, char *argv[]) {

    checkCUDAError(cudaSuccess) // just to prevent a silly compiler warning

    Args args = Args::parseArgs(argc, argv);

    hostsystem::setThreadName("main");

    // register Ctrl+C handler for proper cleanup on interruption
    user_terminate = false;
    signal(SIGINT, signalHandler);

    {
        ifstream f("/proc/self/status");
        string line;
        int tracerPid = 0;
        static const char token[] = "TracerPid:";

        while(getline(f, line)) {
            size_t pos = line.find(token);
            if(pos != string::npos) {
                tracerPid = atoi(line.data() + pos + sizeof(token));
            }
        }
        if (tracerPid > 0)
        	cout << "WARNING: detected attached debugger!" << endl;
    }

    Method method = args.get<Method>("method");
    method.setLDAndDetails(args.count("ld"), args.count("decomposed"));

    // initialize platform
    hostsystem::HostSystem hs(args.get<vector<unsigned>>("gpu"));

    bool useGPU = hs.getGPUs().size() > 0;

    if (useGPU) {
        cout << "Using GPU(s):";
        for (const auto &gpu : hs.getGPUs())
            cout << " " << gpu.getIndex();
        cout << endl;
    }

    const unsigned num_buffers_id = useGPU ? hs.getGPUs().size()+1 : 0; // used for ID creation
    const unsigned num_buffers_per_gpu = useGPU ? args.get<unsigned>("buffers-GPU") : 0;
	if (useGPU && num_buffers_per_gpu == 0) {
		throw(runtime_error("Need to specify at least one transfer buffer for GPU->ResultProcessing!"));
	}

	string setA_str(args.get<string>("setA"));
	string setB_str(args.get<string>("setB"));
	string setC_str(args.get<string>("setC"));
	string exsetA_str(args.get<string>("excludesetA"));
	string exsetB_str(args.get<string>("excludesetB"));
	string exsetC_str(args.get<string>("excludesetC"));

	uint64_t excluderange = args.get<unsigned long>("excluderange");

    // load data
	if (!args.count("quiet"))
	    cout << "Loading data..." << endl;
    SNPDB db;

    string bim, bed, fam;
    if(args.count("snpfile")) {
        string fileset(args.get<string>("snpfile"));
        bim = fileset + string(".bim");
        bed = fileset + string(".bed");
        fam = fileset + string(".fam");
        cout << "Fileset: " << fileset << endl;
    } else {
        bim = args.get<string>("bim");
        bed = args.get<string>("bed");
        fam = args.get<string>("fam");
        cout << "Files:" << endl;
        cout << " bed: " << bed << endl;
        cout << " bim: " << bim << endl;
        cout << " fam: " << fam << endl;
    }

    PlinkParser parser(bim, bed, fam);
    if (useGPU)
        parser.parse(db, 512); // TODO the word size should be set according to GPU memory specifications (memory bus width in bytes or so)
    else // host-only
        parser.parse(db);

    if (!args.count("quiet"))
        cout << "Loading data done." << endl;

    cout << "Found " << db.getCaseCount() << " cases and " << db.getControlCount() << " controls (" << db.getSampleCount() << " samples) typed at " << db.getSNPCount() << " SNPs." << endl;
    verifyData(db, useGPU, args.count("debug"));
    if (args.count("debug"))
        cout << "SNP size: " << db.getSNPSize() << endl;

    // output filename base
    stringstream outfiless;
    if(!args.count("output")) {
        if (!args.count("snpfile"))
            outfiless << args.get<string>("bed");
        else
            outfiless << args.get<string>("snpfile");
    } else
        outfiless << args.get<string>("output");
    outfiless << "." << method.getShortName();

    GPUEngine::score_type threshold = args.get<GPUEngine::score_type>("threshold");
    unsigned long bestresults = args.get<unsigned long>("best-results");

    // User regions
    snprange setA = Args::parseSetRange(setA_str, db.getSNPCount());
    snprange setB = Args::parseSetRange(setB_str, db.getSNPCount());
    snprange setC = Args::parseSetRange(setC_str, db.getSNPCount());
    if (setA_str.compare("all") || setB_str.compare("all") || (method.getOrder() == 3 && setC_str.compare("all"))) { // using user specified SNP ranges
        cout << "User specified SNP sets:" << endl;
        cout << " Set A: " << (setA.first + 1) << " - " << (setA.second) << endl; // 1-based indices and both end included
        cout << " Set B: " << (setB.first + 1) << " - " << (setB.second) << endl; // 1-based indices and both end included
        if (method.getOrder() == 3)
            cout << " Set C: " << (setC.first + 1) << " - " << (setC.second) << endl; // 1-based indices and both end included
    }

    snprange exsetA = Args::parseSetRange(exsetA_str, db.getSNPCount());
    snprange exsetB = Args::parseSetRange(exsetB_str, db.getSNPCount());
    snprange exsetC = Args::parseSetRange(exsetC_str, db.getSNPCount());
    if (exsetA_str.compare("none") || exsetB_str.compare("none") || (method.getOrder() == 3 && exsetC_str.compare("none"))) { // using user specified SNP ranges
        cout << "User specified SNP exclude sets:" << endl;
        cout << " Exclude set A: ";
        if (exsetA.first >= exsetA.second)
            cout << "none" << endl;
        else
            cout << (exsetA.first + 1) << " - " << (exsetA.second) << endl; // 1-based indices and both end included
        cout << " Exclude set B: ";
        if (exsetB.first >= exsetB.second)
            cout << "none" << endl;
        else
            cout << (exsetB.first + 1) << " - " << (exsetB.second) << endl; // 1-based indices and both end included
        if (method.getOrder() == 3) {
            cout << " Exclude set C: ";
            if (exsetC.first >= exsetC.second)
                cout << "none" << endl;
            else
                cout << (exsetC.first + 1) << " - " << (exsetC.second) << endl; // 1-based indices and both end included
        }
    }

    if (excluderange > 0)
        cout << "User specified exclude range: " << excluderange << " bp" << endl;

    // sort sets by left border:
    if (setB.first < setA.first) {
        swap(setA, setB);
        swap(exsetA, exsetB);
    }
    if (method.getOrder() == 3) {
        if (setC.first < setA.first) {
            swap(setA, setC);
            swap(exsetA, exsetC);
        }
        if (setC.first < setB.first) {
            swap(setB, setC);
            swap(exsetB, exsetC);
        }
    }

    array<snprange,3> sets = {setA, setB, setC};
    array<snprange,3> exsets = {exsetA, exsetB, exsetC};


    // buffer sizes and result size
    size_t gpuresultsize = useGPU ? getGPUResultSize(method) : 0;
    size_t idsize = useGPU ? (method.getOrder() == 2 ? IDHandler::IDSIZE_2WAY : IDHandler::IDSIZE_3WAY) : 0;
    size_t buffersize = args.get<size_t>("buffer-size");
    size_t gpubuffersize = useGPU ?
            gpuresultsize * (buffersize / gpuresultsize) // take the provided buffer size as GPU buffer size, but floored to exactly fit a multiple of the results
            : 0; // not specified for CPU-only run
    if (useGPU)
        buffersize = (buffersize / gpuresultsize) * idsize; // adjusted to keep only the IDs for each result

    // define a view on the results
    ResultViewSpec view(method.getOrder(), method.getScoreFields(), gpubuffersize); // gpubuffersize == 0 for CPU-only runs

    ResultHandler<typename ResultViewSpec::id_type, typename ResultViewSpec::score_type> resultHandler(
    			db,
                outfiless.str(),
                useGPU ? num_buffers_per_gpu : 1, // provide one thread for each available buffer or a single one for CPU-only runs
                threshold,
                bestresults,
                view,
                method,
                args
                );

    // main() return value
    int retval = EXIT_FAILURE;
    // progress
    time_t procBegin = time(NULL);
    double progress = 0.0;

    if (useGPU) {

        ///////////////////////////////
        //  hardware accelerated run //
        ///////////////////////////////

        hostsystem::BufferFactory<hostsystem::IDBuffer> idBufferFactory(buffersize);
        hostsystem::BufferFactory<hostsystem::CUDABuffer> gpuBufferFactory(gpubuffersize);

        if (args.count("debug")) {
            cout << "ID buffer size: " << buffersize << " CUDA buffer size: " << gpubuffersize
                << " Items per buffer: " << (gpubuffersize / gpuresultsize) << endl;
        }

        for(unsigned j = 0; j < num_buffers_id; ++j) {
            idBufferFactory.preallocateBuffer(); // allocate buffers for the IDs
        }

        for(unsigned i = 0; i < hs.getGPUs().size(); ++i) {
            cudaSetDevice(hs.getGPU(i).getIndex());
            for(unsigned j = 0; j < num_buffers_per_gpu; ++j) {
                gpuBufferFactory.preallocateBuffer();
            }
        }

        if (!args.count("quiet"))
            cout << "Preloaded buffers." << endl;


        // set up transmission queues
        tbb::concurrent_bounded_queue<shared_ptr<hostsystem::IDBuffer>> dummyQueue, gpuQueue;
        tbb::concurrent_bounded_queue<shared_ptr<hostsystem::CUDABuffer>> resultProcessorQueue;

        // set up queue handlers
        GPUHandler gpuHandler(hs, db, gpuBufferFactory, idsize, buffersize, method, view, args.count("debug"));
        if (!args.count("quiet"))
            cout << "Prepared GPU(s)." << endl;

        IDHandler idHandler(hs, db, idBufferFactory, method, sets, exsets, excluderange, args.count("quiet"), args.count("debug")); // will generate buffers containing IDs

        // starting the ID creator thread
        ThreadPool<shared_ptr<hostsystem::IDBuffer>, shared_ptr<hostsystem::IDBuffer>> idProcessor(1, // one ID creator
                                                                                           dummyQueue,
                                                                                           gpuQueue,
                                                                                           bind(&IDHandler::process,
                                                                                                     &idHandler, _1, _2, _3, _4));

        // starting GPU handler threads (one thread per GPU)
        ThreadPool<shared_ptr<hostsystem::IDBuffer>, shared_ptr<hostsystem::CUDABuffer>> gpuProcessor(hs.getGPUs().size(),
                                                                                          gpuQueue,
                                                                                          resultProcessorQueue,
                                                                                          bind(&GPUHandler::process<shared_ptr<hostsystem::IDBuffer>, shared_ptr<hostsystem::CUDABuffer>>,
                                                                                                    &gpuHandler, _1, _2, _3, _4));

        // starting result handler threads (one thread per GPU buffer)
        ThreadPool<shared_ptr<hostsystem::CUDABuffer>, shared_ptr<hostsystem::IDBuffer>> resultProcessor(num_buffers_per_gpu,
                                                                                             resultProcessorQueue,
                                                                                             dummyQueue,
                                                                                             bind(&ResultHandler<typename ResultViewSpec::id_type, typename ResultViewSpec::score_type>::process<shared_ptr<hostsystem::CUDABuffer>, shared_ptr<hostsystem::IDBuffer>>,
                                                                                                       &resultHandler, _1, _2, _3, _4));

        if (!args.count("quiet"))
            cout << "Let's go..." << endl;
        procBegin = time(NULL); // set timestamp of beginning of process
        bool terminated = false;
        while(!terminated) {
            this_thread::sleep_for(chrono::milliseconds(500)); // poll twice a second

            if(user_terminate || idHandler.isFinished()) {
                idProcessor.terminate(false);
                gpuProcessor.terminate(idHandler.isFinished()); // true if normal finish, false if user interrupt
                if (!args.count("quiet"))
                    cout << "\nGPUs have terminated. " << endl;
                resultProcessor.terminate(idHandler.isFinished()); // true if normal finish, false if user interrupt
                if (!args.count("quiet"))
                    cout << "Result processing has terminated. " << endl;
                terminated = true;
            }

            progress = idHandler.progress();

            if (!args.count("quiet") && !user_terminate) {

                // print progress
                if(isatty(STDOUT_FILENO) && !args.count("debug"))
                    cout << "\r                                                                                  \r";
                if(args.count("debug")) {
                // Note that the actual numbers are calculated as actual queue size + putters - getters,
                // so if the queue is empty and 5 threads are waiting for a buffer, it will display as "-5".
                    cout << "ID -> [" << gpuQueue.size()
                              << "] -> GPU -> [" << resultProcessorQueue.size()
                              << "] -> Result Processing ";
                }

                printProgressTime(cout, progress, procBegin);

                // If we're on a terminal, just update the line instead of making a new line to not spam the output
                // This won't work when the output is piped (i.e. not a terminal), then insert a newline instead.
                if(!isatty(STDOUT_FILENO) || args.count("debug"))
                    cout << endl;
            }
        }

    } else { // !useGPU

        ////////////////////////////
        //      CPU-only run      //
        ////////////////////////////

        // determine the correct result size
        size_t resultsize = view.getResultSize();

        mutex m;

        bool decomp = args.count("decomposed");
        bool ld = args.count("ld");

        snprange locexA = exsetA;
        snprange locexB = exsetB;
        snprange locexC = exsetC;
        string chrA;
        string chrB;
        string chrC;
        unsigned long bpposA = 0;
        unsigned long bpposB = 0;
        unsigned long bpposC = 0;

        if (!args.count("quiet"))
            cout << "Let's go..." << endl;
        procBegin = time(NULL); // set timestamp of beginning of process
        cout << setprecision(4) << fixed << ((double) 0) << " %" << flush;

        if (method.getOrder() == 2) { // 2way

            double pairs = (double) getPairsCnt(setA.first, setA.second, setB.first, setB.second); // cast for progress print

            unsigned long currpair = 0;

            // the smaller right border
            unsigned rightA = setA.second < setB.second ? setA.second : setB.second;

            for(unsigned snpA = setA.first; snpA < rightA && !user_terminate; ++snpA) {

                // "switch" the ranges if the right borders are not in order, but only if we reached the overlap
                unsigned rightB = setB.second < setA.second && snpA >= setB.first ? setA.second : setB.second;
                if (setB.second < setA.second && snpA >= setB.first)
                    swap(locexA, locexB);

                // apply exclude set A: skip all combinations with this SNP if snpA is in exclude range
                if (snpA >= locexA.second || snpA < locexA.first) { // SNP A is not excluded

                    if (excluderange > 0) { // if an exclude range is defined
                        chrA = db.getSNPInfo(snpA).chromosome;
                        bpposA = db.getSNPInfo(snpA).pos_bp;
                    }

#ifndef DEBUG_SINGLETHREAD
#pragma omp parallel for
#endif
                    for(unsigned snpB = max((unsigned)setB.first, snpA+1); snpB < rightB; ++snpB) {
                        // apply exclude set B: include this pair if snpB is not in exclude range
                        if (snpB >= locexB.second || snpB < locexB.first) {
                            if (excluderange > 0) { // if an exclude range is defined
                                chrB = db.getSNPInfo(snpB).chromosome;
                                bpposB = db.getSNPInfo(snpB).pos_bp;
                            }
                            if (excluderange == 0 || chrA.compare(chrB) || bpposB - bpposA >= excluderange) {

                                // test snpA and snpB
                                ctable2way t;
                                uint32_t numcases, numctrls;
                                generateContingencyTable2Way(t.data(), &numcases, &numctrls, db.data(), db.getSNPSize(), db.getCaseCountPadded(), snpA, snpB);

                                char* resultbuf = (char*)malloc(resultsize*sizeof(char)); // TODO could be improved
                                ResultViewSpec v(view); // local copy
                                v.setBuffer(resultbuf, resultsize);

                                switch (method.getType()) {
                                case Method::MI2:
                                    MutualInformationKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                    break;
                                case Method::IG2:
                                    InformationGainKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                    break;
                                case Method::BOOST:
                                    BOOSTKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                    break;
                                case Method::LOGLIN:
                                    LogLinearKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                    break;
                                case Method::LOGREG:
                                    LogisticRegressionKernel2WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                    break;
                                default: // 3way or unknown kernel (should never get here)
                                    throw runtime_error("Undefined kernel!");
                                }

                                m.lock();
                                resultHandler.processResultView(v, 0);
                                m.unlock();

                                free(resultbuf);
                            }
                        } // end if snpB is not excluded
                    } // end for snpB
                } // end if snpA is not in exclude range

                // progress
                currpair += rightB - max((unsigned)setB.first, snpA+1);
                progress = currpair/pairs;

                if (!args.count("quiet")) {
                    if(isatty(STDOUT_FILENO) && !args.count("debug"))
                        cout << "\r                                                                                  \r";
                    if (args.count("debug"))
                        cout << currpair << "/" << (unsigned long) pairs << " ";
                    printProgressTime(cout, progress, procBegin);
                    if(!isatty(STDOUT_FILENO) || args.count("debug"))
                        cout << endl;
                }
            } // end for snpA

        } else { // 3way

            double triples = (double) getTriplesCnt(setA.first, setA.second, setB.first, setB.second, setC.first, setC.second); // cast for progress print

            unsigned long currtriple = 0;

            size_t currA = setA.second;
            size_t currB = setB.second;
            size_t currC = setC.second;

            for(size_t snpA = setA.first; snpA < currA && !user_terminate; ++snpA) {
                // border for b in AB overlap is the larger right limit of A and B,
                // in fact we are switching A and B here w.l.o.g. but to ensure the right limits are sorted
                if (snpA >= setB.first && currB < currA) {
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
                // apply exclude set A: simply set a flag that is recognized in the inner loops, since it is difficult to correct the progress otherwise,
                // anyway this could delay the main process significantly only if the exclusion area is large
                // TODO if necessary, I will fix this in the future
                bool exsnpA = false;
                if (snpA >= locexA.first && snpA < locexA.second) { // SNP A excluded
                    exsnpA = true;
                }
                bool switchedBC = false;
                if (excluderange > 0) { // if an exclude range is defined
                    chrA = db.getSNPInfo(snpA).chromosome;
                    bpposA = db.getSNPInfo(snpA).pos_bp;
                }

                for(size_t snpB = max(setB.first, snpA+1); snpB < currB && !user_terminate; ++snpB) {
                    // border for c in BC overlap is the larger right limit of B and C,
                    // in fact we are switching B and C here w.l.o.g. but to ensure the right limits are sorted
                    if (snpB >= setC.first && currC < currB) {
                        swap(currB, currC);
                        swap(locexB, locexC);
                        switchedBC = true;
                    }
                    // apply exclude set A+B: skip all combinations with this SNP if snpA or snpB is in exclude range
                    if (!exsnpA && (snpB >= locexB.second || snpB < locexB.first)) { // SNP A and B not excluded
                        if (excluderange > 0) { // if an exclude range is defined
                            chrB = db.getSNPInfo(snpB).chromosome;
                            bpposB = db.getSNPInfo(snpB).pos_bp;
                        }
                        // if there is an exclude range defined we could skip the complete next loop for SNP C if A and B are already in that range
                        if (excluderange == 0 || chrA.compare(chrB) || bpposB - bpposA >= excluderange) {

#ifndef DEBUG_SINGLETHREAD
#pragma omp parallel for
#endif
                            for(size_t snpC = max(setC.first, snpB+1); snpC < currC; ++snpC) {
                                // apply exclude set C: include this pair if snpC is not in exclude range
                                if (snpC >= locexC.second || snpC < locexC.first) {
                                    if (excluderange > 0) { // if an exclude range is defined
                                        chrC = db.getSNPInfo(snpC).chromosome;
                                        bpposC = db.getSNPInfo(snpC).pos_bp;
                                    }
                                    if (excluderange == 0 || ( (chrB.compare(chrC) != 0 || bpposC - bpposB >= excluderange)
                                    // && (chrA.compare(chrB) != 0 || bpposB - bpposA >= excluderange) // this was already checked before snpC loop
                                    // && (chrA.compare(chrC) != 0 || bpposC - bpposA >= excluderange) // this test is not required since SNPs A,B and C are ordered by position
                                     ) ) {
                                        // test snpA and snpB and snpC
                                        ctable3way t;
                                        uint32_t numcases, numctrls;
                                        generateContingencyTable3Way(t.data(), &numcases, &numctrls, db.data(), db.getSNPSize(), db.getCaseCountPadded(), snpA, snpB, snpC);

                                        char* resultbuf = (char*)malloc(resultsize*sizeof(char)); // TODO could be improved
                                        ResultViewSpec v(view); // local copy
                                        v.setBuffer(resultbuf, resultsize);

                                        switch (method.getType()) {
                                        case Method::MI3:
                                            MutualInformationKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                            break;
                                        case Method::IG3:
                                            InformationGainKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                            break;
                                        case Method::LOGREG3:
                                            LogisticRegressionKernel3WayCore(t.data(), numcases, numctrls, v.getDeviceResult(0), ld, decomp);
                                            break;
                                        default: // 2way or unknown kernel (should never get here)
                                            throw runtime_error("Undefined kernel!");
                                        }

                                        m.lock();
                                        resultHandler.processResultView(v, 0);
                                        m.unlock();

                                        free(resultbuf);
                                    }
                                }
                            } // end for snpC
                        }
                    }

                    // progress
                    currtriple += currC - max(setC.first, snpB+1);
                    progress = currtriple/triples;

                    if (!args.count("quiet")) {
                        if(isatty(STDOUT_FILENO) && !args.count("debug"))
                            cout << "\r                                                                                  \r";
                        if (args.count("debug"))
                            cout << currtriple << "/" << (unsigned long) triples << " ";
                        printProgressTime(cout, progress, procBegin);
                        if(!isatty(STDOUT_FILENO) || args.count("debug"))
                            cout << endl;
                    }
                } // end for snpB
                if (switchedBC) { // switch back B and C if they were switched in the loop before
                    swap(currB, currC);
                    swap(locexB, locexC);
                }
            } // end for snpA
        } // end 3way
    } // end if(useGPU) ... else ...

    if (!user_terminate) { // finished
        cout << "\nResult file: " << resultHandler.getResultFileName() << endl;
        if (!args.count("quiet"))
            cout << "\nWriting results..." << endl;
        string cmdline;
        for (int i = 1; i < argc; i++) // all args but the binary
            cmdline += string(" ") + string(argv[i]);
        resultHandler.flush(cmdline);
        retval = EXIT_SUCCESS;
    } else { // User terminate or timeout
        cout << "Terminating at:" << endl;
        printProgressTime(cout, progress, procBegin);
        cout << endl;
    }

    return retval;
}
