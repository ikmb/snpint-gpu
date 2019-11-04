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

#ifndef RESULTHANDLER_H
#define RESULTHANDLER_H

#include <atomic>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <fstream>

#include "MinMaxHeap.h"
#include "Histogram.h"

#include "SNPDB.h"
#include "GPUEngine.h"
#include "ResultView.h"
#include "ThreadPool.h"
#include "Method.h"
#include "Args.h"
#include "hostsystem/Buffer.h"
#include "hostsystem/ThreadUtils.h"

template<typename Id, typename Score>
class ResultHandler
{
public:
    explicit ResultHandler(const SNPDB &db, const std::string &output_file_, int num_threads_, Score threshold_, unsigned long bestResults_, ResultView<Id,Score> &v, Method &method_, Args &args)
        : snpdb(db),
		  output_file(output_file_),
          view(v),
          heaps(num_threads_),
          histograms(),
          num_threads(num_threads_),
          threshold(threshold_),
          bestResults(bestResults_),
          method(method_),
          debug(args.count("debug"))
    {
        useThreshold = !std::isnan(threshold);
        useResults = (bestResults > 0) || useThreshold;

        heaps.resize(num_threads);

        // histograms
        initHistograms(args, v);
    }

    void processResultView(ResultView<Id,Score> &v, int threadIndex) {

		// scan through all results
        if (v.getScoreCount()) { // ensure that this method has at least one result
            for(auto it = v.begin(); it != v.end(); ++it) {

                // add all scores of the current result to the histogram set, if enabled
                if (!histograms.empty()) {
                    for(size_t i = 0; i < v.getScoreCount(); i++) {
                        if (useHistogram[i])
                            histograms[threadIndex][i].add(it.getScore(i));
                    }
                }

                // don't need results if we're working on a permutation; only histograms
                if(!useResults)
                    continue;

                // sorting out by threshold is enabled
                if(useThreshold) {
                    if(it.getScore(0) >= threshold)
                        heaps[threadIndex].push(*it);

                // push into heap, maintaining its capacity
                } else {
                    if(heaps[threadIndex].getSize() >= bestResults) {
                        if(heaps[threadIndex].getMin().getScore(0) < it.getScore(0)) {
                            heaps[threadIndex].popMin();
                            heaps[threadIndex].push(*it);
                        }
                    } else {
                        heaps[threadIndex].push(*it);
                    }
                }
            }
        }
    }

    template<typename Tin, typename Tout>
    void process(typename ThreadPool<Tin, Tout>::inqueue_type &inqueue,
                 typename ThreadPool<Tin, Tout>::outqueue_type &outqueue,
                 std::atomic_bool& termination_request,
                 int threadIndex
                 ) {

        ResultView<Id,Score> v(view); // create a local copy of the view
        (void) outqueue;

        hostsystem::setThreadName("ResultH:" + std::to_string(threadIndex));


        while(!termination_request) {
            Tin srcbuf;
            try {
                inqueue.pop(srcbuf);
            } catch(tbb::user_abort &e) {
//            	std::cerr << "Result thread " << threadIndex << ": Abort." << std::endl;
                break;
            }

			// set view to current buffer
            v.setBuffer(srcbuf->getData(), (srcbuf->getContentLength())*v.getResultSize());

            processResultView(v, threadIndex);
        }
    }

    void flush(const string &cmdline) {

        if(useResults)
            flushResults(cmdline);

        if(!histograms.empty())
            flushHistograms();

    }

    string getResultFileName() const {
        return output_file + ".scores";
    }


private:
    const SNPDB &snpdb;
    std::string output_file;
    ResultView<Id,Score> view;
    std::vector<jcommon::MinMaxHeap<Result<Id,Score>>> heaps;
    std::vector<std::vector<jcommon::Histogram<Score>>> histograms;

    int num_threads;
    Score threshold;
    unsigned long bestResults;

    Method method;

    bool useThreshold;
    bool useResults;
    bool useHistogram[Method::maxScoreComponents];

    bool debug;

    void initHistograms(Args &args, ResultView<Id,Score> &v) {
        // histogram resolutions and boundaries
        unsigned histogram_resolutions[Method::maxScoreComponents];
        Score histogram_minima[Method::maxScoreComponents];
        Score histogram_maxima[Method::maxScoreComponents];

        Score last_min = std::nan("");
        Score last_max = std::nan("");
        if (args.count("hmin"))
            last_min = args.get<Score>("hmin");
        if (args.count("hmax"))
            last_max = args.get<Score>("hmax");
        for (int i = 0; i < Method::maxScoreComponents; i++) {
            std::stringstream sshres; sshres << "hres" << i;
            std::stringstream sshmin; sshmin << "hmin" << i;
            std::stringstream sshmax; sshmax << "hmax" << i;
            histogram_resolutions[i] = args.get<unsigned>(sshres.str());
            if (args.count(sshmin.str())) {
                histogram_minima[i] = args.get<Score>(sshmin.str());
                last_min = histogram_minima[i];
            } else {
                histogram_minima[i] = last_min;
            }
            if (args.count(sshmax.str())) {
                histogram_maxima[i] = args.get<Score>(sshmax.str());
                last_max = histogram_maxima[i];
            } else {
                histogram_maxima[i] = last_max;
            }
            // if histogram is activated, check if the user has provided the mandatory boundaries. otherwise, deactivate histogram
            if (histogram_resolutions[i] && (std::isnan(histogram_minima[i]) || std::isnan(histogram_maxima[i]))) {
                std::cerr << "WARNING: No valid boundaries found for histogram " << i << ". Deactivating histogram." << std::endl;
                histogram_resolutions[i] = 0;
            }
        }
        // --hres0 has precedence before --hres, but if the former has not been applied at all, the value is set to zero now. So, take the latter then.
        if (histogram_resolutions[0] == 0)
            histogram_resolutions[0] = args.get<unsigned>("hres");

        if (args.count("debug")) {
            bool usehist = false;
            for (unsigned int i = 0; i < v.getScoreCount(); i++)
                usehist |= histogram_resolutions[i] != 0;
            if (usehist) {
                std::cout << "Histograms:" << std::endl;
                for (unsigned int i = 0; i < v.getScoreCount(); i++)
                    if (histogram_resolutions[i])
                        std::cout << i << ":\thres = " << histogram_resolutions[i] << "\thmin = " << histogram_minima[i] << "\thmax = " << histogram_maxima[i] << std::endl;
            }
        }

        // init histograms
        bool uh = false;
        for (int i = 0; i < Method::maxScoreComponents; i++) {
            useHistogram[i] = histogram_resolutions[i] != 0;
            if (useHistogram[i])
                uh = true;
        }

        if (uh) {
            histograms.resize(num_threads);
            for(auto &h: histograms)
                for (size_t i = 0; i < v.getScoreCount(); i++)
                    h.emplace_back(histogram_resolutions[i], histogram_minima[i], histogram_maxima[i]);
        }
    }

    void flushResults(const string &cmdline) {
        // build filename
        std::stringstream ss;
        ss << output_file << ".scores";
        std::ofstream results {ss.str()};

        // how many results did we collect?
        unsigned long collected_results = 0ul;
        if (debug)
            std::cout << "\nResults from " << heaps.size() << " heaps: ";
        for (const jcommon::MinMaxHeap<Result<Id,Score>> &h : heaps) {
        	collected_results += h.getSize();
        	if (debug)
        	    std::cout << h.getSize() << " ";
        }

        // how many do we want?
        unsigned long actual_results = collected_results;
        if (!useThreshold && collected_results > bestResults)
            actual_results = bestResults;
        if (debug)
            std::cout << "\nTotal from heaps: " << collected_results << "\nTotal to be saved: " << actual_results << std::endl;

        // print command line
        results << "#" << cmdline << endl;

        // print header
        method.printHeader(results);

        // collect <actual_results> top elements from all our heaps
        for(unsigned long i = 0; i < actual_results; i++) {
            auto it = std::max_element(std::begin(heaps), std::end(heaps), [](const jcommon::MinMaxHeap<Result<Id,Score>> &lhs, const jcommon::MinMaxHeap<Result<Id,Score>> &rhs) {
                if(lhs.isEmpty()) return true;
                if(rhs.isEmpty()) return false;
                return lhs.getMax() < rhs.getMax();
            });

            // print SNP names
            auto &r = it->getMax();
            for(unsigned id = 0; id < r.getID().size(); id++) {
				results << snpdb.getSNPInfo(r.getID(id)).variant_id << "\t";
			}
            // print rest
            results << it->popMax();
        }

        results.close();

        std::cout << "Saved " << actual_results << " results." << std::endl;
    }

    void flushHistograms() {
        if (debug)
            std::cout << "Merging " << histograms.size() << " histograms." << std::endl;
        for(size_t i = 1; i < histograms.size(); i++) {
            for(size_t j = 0; j < view.getScoreCount(); j++)
            	if (useHistogram[j])
            		histograms[0][j] += histograms[i][j];
        }

        if (debug)
            std::cout << "Writing histogram." << std::endl;
        for(size_t hist = 0; hist < view.getScoreCount(); hist++) {
        	if (useHistogram[hist]) {
				const jcommon::Histogram<Score> &current = histograms[0][hist];

				// build filename
				std::stringstream ss;
				ss << output_file;
				if (view.getScoreCount() > 1)
					ss << "." << hist;
				ss << ".hist";

				std::ofstream histfile {ss.str()};

				current.writeToFile(histfile);

				histfile.close();
        	}
        }
    }
};

#endif // RESULTHANDLER_H
