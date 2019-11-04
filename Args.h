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

#ifndef ARGS_H_
#define ARGS_H_

#include <string>
#include <utility>

#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

using namespace std;

/**
 * Class for storing and retrieving command-line arguments.
 */
class Args {
public:

	static Args parseArgs(int argc, char *argv[]);

    /**
     * Returns the argument with the given (long) name. The template argument
     * specifies the return type you want the parameter casted to.
     * @param name argument name
     * @return the variable value
     */
    template<typename T> T get(const string &name) const {
        auto where = vars.find(name);
        if(where == end(vars)) {
            if(!isDefined(name))
                throw invalid_argument("Option undefined: " + name + " (This is a bug)");
            else
                throw out_of_range("Option has not been specified and does not have a default value associated: " + name);
        }
        return where->second.as<T>();
    }

    /**
     * Counts the number of argument occurences. Mainly useful for boolean switches
     * and/or counting flags, such as verbosity or debug mode.
     * @param name argument name
     * @return argument value
     */
    unsigned int count(const string &name) const {
        if(!isDefined(name))
            throw invalid_argument("Option undefined: " + name + " (This is a bug)");
        return vars.count(name);
    }

    bool operator()(const string &name) const {
        return count(name) > 0;
    }

    /**
     * Prints a help message.
     * @param progname the program name, usually argv[0]
     * @param out output stream, e.g. cout
     */
    void printHelp(const string &progname, ostream &out) const;

    // Parses a range string of the form "x-y" where x and y are 1-based SNP indices, both including!
    // Single elements and open ranges of the form "x-" are allowed as well as the string "all" which
    // means "from 1 to max".
    // The provided maximum is the maximum right value of a provided range, i.e. 1-based and including
    // (which is the same as 0-based and excluding... ;-))
    // The returned pair is developer-friendly 0-based and forms an interval where its left limit is included
    // and its right limit is excluded.
    // Throws an error if the range exceeds in any kind the provided maximum.
    static pair<size_t,size_t> parseSetRange(string rset, size_t maximum);

    static void printVersion();

    Args(Args&& other) = default;

protected:
    /** Constructs the arguments list and adds all defined options */
    Args();
    Args(Args const &);
    void operator=(Args const &);

    void parse(int argc, char *argv[]);
    bool isDefined(const string &optname) const;

    bpo::options_description opts_regular;        /**< regular options, shown on help output */
    bpo::options_description opts_regular_GPU; /**< regular options, shown on help output */
    bpo::options_description opts_hidden;         /**< hidden options */
    bpo::options_description opts_hidden_GPU;  /**< hidden options */
    bpo::options_description opts_hres;           /**< hidden options for histograms, never printed */
    bpo::positional_options_description opts_positional;    /**< positional options (without name) */

    bpo::variables_map vars;    /**< this is where the option values go */

    /** dump all options to the given stream */
    friend ostream &operator<<(ostream &out, const Args &args);

private:
    /** parses the main() options into the variable map */
    Args(int argc, char *argv[]);
};

#endif /* ARGS_H_ */
