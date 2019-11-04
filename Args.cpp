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

#include <iostream>
#include <iomanip>
#include <cmath>
#include <thread>

#include <boost/program_options.hpp>

#include "Args.h"
#include <version.h>
#include "GPUEngine.h"

#include "Method.h"

namespace bpo = boost::program_options;

using namespace std;
using namespace bpo;

/**
 * @brief Constructs, parses and verifies all command-line arguments
 *
 * This function constructs the Args object and parses all given command-line
 * arguments. If unknown arguments and/or missing obligatory arguments are
 * detected, this function does not return and instead prints an appropriate
 * help message and calls exit(), executing all exit handlers registered
 * up to the parseArgs call.
 *
 * @param argc argument count
 * @param argv argument vector
 * @return an rvalue reference of a new Args object containing all defined or defaulted CLI arguments
 */
Args Args::parseArgs(int argc, char *argv[]) {
    Args args {argc, argv};

    if(args.count("help")) {
        args.printHelp(argv[0], cout);
        exit(EXIT_SUCCESS);
    }

    if(args.count("version")) {
        args.printVersion();
        exit(EXIT_SUCCESS);
    }

    if(!(args.count("snpfile") || (args.count("bim") && args.count("bed") && args.count("fam")))) {
        args.printHelp(argv[0], cout);
        throw runtime_error("You need to specify an input file basename OR individual file names with --{bim,bed,fam}. You may not specify both or neither.");
    }

    if(!args.count("method")) {
        args.printHelp(argv[0], cout);
        throw runtime_error("Missing required option --method/-m");
    }

    return args;
}

ostream &operator<<(ostream &out, const Args &args) {
    variables_map::const_iterator it;

    long name_width = 0;
    long value_width = 0;
    long orig_width = out.width();

    // collect field widths
    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        long this_width = static_cast<long>(it->first.length());

        if(this_width > name_width) {
            name_width = this_width;
        }

        this_width = 0;

        if(it->second.value().type() == typeid(string)) {
            this_width = static_cast<long>(it->second.as<string>().length());
        }

        if(it->second.value().type() == typeid(int)) {
            this_width = static_cast<long>(log10(it->second.as<int>()));
        }

        if(this_width > value_width) {
            value_width = this_width;
        }
    }

    // dump option values
    out.setf(ios::adjustfield, ios::left);

    for(it = args.vars.begin(); it != args.vars.end(); ++it) {
        out.width(name_width + 2); // more space for the colon
        out << (it->first + ":") << endl;

        out.width(value_width);

        if(it->second.value().type() == typeid(string)) {
            out << it->second.as<string>() << endl;
        } else if(it->second.value().type() == typeid(int)) {
            out << it->second.as<int>() << endl;
        } else {
            out << "(unknown)" << endl;
        }
    }

    resetiosflags(ios::adjustfield);
    out.width(orig_width);

    return out;
}

Args::Args(int argc, char *argv[]) :
    opts_regular("Program options"),
	opts_regular_GPU("Additional program options for GPU architecture"),
    opts_hidden("Hidden options"),
	opts_hidden_GPU("Additional hidden options for GPU architecture"),
	opts_hres("Hidden options for histograms, never printed")
	{

    opts_regular.add_options()
    ("help,h", "produce this help message and terminates")
    ("version", "prints version information and terminates")
    ("quiet", "reduces output to a minimum")
    ("output,o", value<string>(), "output file(s prefix, if applicable)")
    ("method,m", value<Method>(), "select method (see below)")
    ("best-results,n", value<unsigned long>()->default_value(1000UL), "select only the N best results, 0 for no results (e.g. histogram only)")
	("decomposed", "decompose result values (e.g. information gain is decomposed into mutual information subvalues)")
    ("threshold,t", value<GPUEngine::score_type>()->default_value(numeric_limits<GPUEngine::score_type>::quiet_NaN()), "threshold for result selection")
    ("ld", "additionally calculate LD (r-square). 2way: From all possible solutions the highest and the lowest value are provided. 3way: all pairwise (highest) r-squares are provided.")
    ("hres", value<unsigned>()->default_value(0), "histogram resolution in bits for first score component, others: --hres# (e.g. --hres2) (default = 0 for no histogram)")
    ("hmin", value<GPUEngine::score_type>(), "left border of histogram (inclusive) for first score component, others: --hmin# (e.g. --hmin2); for the leftmost component, this option is mandatory, others use the value of the previous one, if activated with --hres#")
    ("hmax", value<GPUEngine::score_type>(), "right border of histogram (exclusive) for first score component, others: --hmax# (e.g. --hmax2); for the leftmost component, this option is mandatory, others use the value of the previous one, if activated with --hres#")
    ("bim", value<string>(), "BIM file input (can be used instead of the <SNP input file> basename")
    ("bed", value<string>(), "BED file input (can be used instead of the <SNP input file> basename")
    ("fam", value<string>(), "FAM file input (can be used instead of the <SNP input file> basename")
    ("setA", value<string>()->default_value("all"), "Subset of SNPs for interaction tests at first position in pair or triple. 1-based SNP index range or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("setB", value<string>()->default_value("all"), "Subset of SNPs for interaction tests at second position in pair or triple. 1-based SNP index range or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("setC", value<string>()->default_value("all"), "Subset of SNPs for interaction tests at third position in triple. 1-based SNP index range or \"all\" (default). Range may as well be a single SNP or have an open end.")
    ("excludesetA", value<string>()->default_value("none"), "Subset of SNPs for interaction tests to be excluded at first position in pair or triple. 1-based SNP index range or \"none\" (default). Range may as well be a single SNP.")
    ("excludesetB", value<string>()->default_value("none"), "Subset of SNPs for interaction tests to be excluded at second position in pair or triple. 1-based SNP index range or \"none\" (default). Range may as well be a single SNP.")
    ("excludesetC", value<string>()->default_value("none"), "Subset of SNPs for interaction tests to be excluded at third position in triple. 1-based SNP index range or \"none\" (default). Range may as well be a single SNP.")
    ("excluderange", value<unsigned long>()->default_value(0), "exclude test if a pair (pairwise tests) or one pair of a triple (3way tests) is within this range (in bp)")
    ;

    opts_regular_GPU.add_options()
    ("gpu", value<vector<unsigned>>()->multitoken()->default_value(vector<unsigned>(),""), "restrict to these GPUs (separate indices by whitespace, defaults to no GPU acceleration)")
    ;

    opts_hidden.add_options()
    ("debug", "produce lots of debug output")
    ("snpfile", value<string>(), "SNP input file (filename base for .bim/.bed/.fam (automatic parameter for positional argument #1)")
	;

    opts_hidden_GPU.add_options()
    ("buffer-size", value<size_t>()->default_value(256*1024*1024), "Size for transmission buffers (GPU->host) in bytes.")
    ("buffers-GPU", value<unsigned>()->default_value(16), "Number of transmission buffers (GPU->host) to keep around.")
    ;

	// need one for each possible result component
    opts_hres.add_options()
	("hres0", value<unsigned>()->default_value(0), "see --hres")
	("hmin0", value<GPUEngine::score_type>(), "see --hmin")
	("hmax0", value<GPUEngine::score_type>(), "see --hmax")
	("hres1", value<unsigned>()->default_value(0), "see --hres")
	("hmin1", value<GPUEngine::score_type>(), "see --hmin")
	("hmax1", value<GPUEngine::score_type>(), "see --hmax")
	("hres2", value<unsigned>()->default_value(0), "see --hres")
	("hmin2", value<GPUEngine::score_type>(), "see --hmin")
	("hmax2", value<GPUEngine::score_type>(), "see --hmax")
	("hres3", value<unsigned>()->default_value(0), "see --hres")
	("hmin3", value<GPUEngine::score_type>(), "see --hmin")
	("hmax3", value<GPUEngine::score_type>(), "see --hmax")
	("hres4", value<unsigned>()->default_value(0), "see --hres")
	("hmin4", value<GPUEngine::score_type>(), "see --hmin")
	("hmax4", value<GPUEngine::score_type>(), "see --hmax")
	("hres5", value<unsigned>()->default_value(0), "see --hres")
	("hmin5", value<GPUEngine::score_type>(), "see --hmin")
	("hmax5", value<GPUEngine::score_type>(), "see --hmax")
	("hres6", value<unsigned>()->default_value(0), "see --hres")
	("hmin6", value<GPUEngine::score_type>(), "see --hmin")
	("hmax6", value<GPUEngine::score_type>(), "see --hmax")
	("hres7", value<unsigned>()->default_value(0), "see --hres")
	("hmin7", value<GPUEngine::score_type>(), "see --hmin")
	("hmax7", value<GPUEngine::score_type>(), "see --hmax")
	("hres8", value<unsigned>()->default_value(0), "see --hres")
	("hmin8", value<GPUEngine::score_type>(), "see --hmin")
	("hmax8", value<GPUEngine::score_type>(), "see --hmax")
	("hres9", value<unsigned>()->default_value(0), "see --hres")
	("hmin9", value<GPUEngine::score_type>(), "see --hmin")
	("hmax9", value<GPUEngine::score_type>(), "see --hmax")
	("hres10", value<unsigned>()->default_value(0), "see --hres")
	("hmin10", value<GPUEngine::score_type>(), "see --hmin")
	("hmax10", value<GPUEngine::score_type>(), "see --hmax")
	;
    assert(Method::maxScoreComponents == 11); // this is the maximum number of score components for all currently available methods

    opts_positional.add("snpfile", 1);

    parse(argc, argv);
}

void Args::parse(int argc, char *argv[]) {
    bpo::options_description all_options;

    // combine all options
    all_options.add(opts_regular).add(opts_hidden).add(opts_hres);

    all_options.add(opts_regular_GPU).add(opts_hidden_GPU);

    // do the actual parsing
    store(command_line_parser(argc, argv).options(all_options).positional(opts_positional).run(), vars);
    notify(vars);
}

bool Args::isDefined(const string &optname) const {
    bool found = false;
    found = !!this->opts_regular.find_nothrow(optname, false); // return null (-> false) if option has not been found
    found |= !!this->opts_hidden.find_nothrow(optname, false);
    found |= !!this->opts_hres.find_nothrow(optname, false);
    found |= !!this->opts_regular_GPU.find_nothrow(optname, false);
    found |= !!this->opts_hidden_GPU.find_nothrow(optname, false);
    return found;
}

void Args::printHelp(const string &progname, ostream &out) const {
    out << "SNPInt-GPU  Copyright (C) 2019  Lars Wienbrandt, Jan Christian Kässens" << endl;
    out << "This program comes with ABSOLUTELY NO WARRANTY;\nThis is free software, and you are welcome to redistribute it under certain conditions;\nSee file COPYING for details.\n" << endl;
    out << "Usage: " << progname << " <SNP input file> [options]" << endl << endl;
    out << opts_regular << endl;
    out << opts_regular_GPU << endl;
//    out << opts_hidden << endl;
//    out << opts_hidden_GPU << endl;
    out << endl;

    out << "The following methods are available for the -m/--method option:" << endl;
    for(const auto& m: Method::getMethods()) {
        if(m.type != Method::Type::INVALID)
            out << "  " << left << setw(12) << m.shortName << "" << m.descriptiveName << endl;
    }
    out << endl;

    printVersion();
}

/* static */
snprange Args::parseSetRange(string rset, size_t maximum) {
    snprange range;
    istringstream rs(rset);
    string r;

    if (!rset.compare("all")) {
        return make_pair(0, maximum); // maximum range for string "all"
    }

    if (!rset.compare("none")) {
        return make_pair(0,0); // empty range for string "none"
    }

    size_t pos = rset.find('-');
    if (pos == string::npos) { // no range, single element
        size_t idx = stoull(rset);
        if (idx > 0 && idx <= maximum) {
            range.first = idx-1; // indices are provided 1-based but stored zero-based
            range.second = idx;  // second is excluding
        } else
            throw runtime_error("No valid range argument or range exceeds maximum!");
    } else { // range
        size_t idx1 = stoull(rset.substr(0,pos)); // parsing will stop at the dash
        size_t idx2;
        if (pos == rset.size()-1)
            idx2 = maximum; // open range "x-"
        else
            idx2 = stoull(rset.substr(pos+1)); // parsing after the dash
        if (idx1 > idx2 || idx1 <= 0 || idx2 <= 0 || idx1 > maximum || idx2 > maximum)
            throw runtime_error("No valid range argument or range exceeds maximum!");
        range.first = idx1-1; // indices provided 1-based
        range.second = idx2;   // second is excluding
    }
    return range;
}

/* static */
void Args::printVersion() {
    cout << "This is version " << prog_version << ", compiled on " << prog_timestamp << endl;
    cout << "Send bugs to " << prog_bugaddress << endl;
}
