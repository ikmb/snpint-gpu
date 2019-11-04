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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <vector>
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <limits>
#include <stdexcept>
#include <fstream>

namespace jcommon {

template<class T>
class Histogram {

public:

	// min is inclusive, max is exclusive; resolution is in bits
    // DO NOT USE std::numeric_limits<T>::min() for floating point types since it returns only the minimum positive normalized value!
	explicit Histogram(int resolution_, T min = std::numeric_limits<T>::lowest(), T max = std::numeric_limits<T>::max())
	: resolution(resolution_), mind((double)min), maxd((double)max), mini((uint64_t)min), maxi((uint64_t)max),
	  histogram(1 << resolution_)
	{
		// this is the multiplication factor for the interval mapping
		map_factor = std::numeric_limits<uint64_t>::max();
		map_factor /= maxd - mind;
		//std::cout << map_factor << " " << std::hex << min << " " << max << std::endl;

		// reserve space for the corrsponding limits array
		limitsi.clear();
		limitsd.clear();
		if (std::is_integral<T>()) {
			limitsi.resize(1 << resolution);
			//std::cout << "creating " << (1 << resolution) << " fields for integer histogram." << std::endl;
		} else {
			limitsd.resize(1 << resolution);
			//std::cout << "creating " << (1 << resolution) << " fields for fp histogram." << std::endl;
		}

		// generate the limits array
        for (std::size_t i = 0; i < histogram.size(); ++i) {
			double limit_tmp = (i * (maxd - mind)) / (1ull << resolution);
			//std::cout << std::dec << i << " " << limit_tmp << " ";
			if (std::is_integral<T>()) {
				limit_tmp = ceil(limit_tmp); // have to take ceiling if T is an integer type!
				//std::cout << limit_tmp << " ";
				limitsi[i] = limit_tmp + mini;
				//std::cout << std::hex << limitsi[i] << " " << std::dec << limitsi[i] << std::endl;
			} else {
				//std::cout << limit_tmp << " ";
				limitsd[i] = limit_tmp + mind;
				//std::cout << limitsd[i] << std::endl;
			}
		}
	}

	void add(T value_p, uint64_t increment = 1) {

		if (std::is_integral<T>()) {
			uint64_t value = (uint64_t) value_p;
			if (value < mini)
				histogram_n_ov+=increment;
			else if (value >= maxi)
				histogram_p_ov+=increment;
			else {
				// map the value to a uint64_t
				uint64_t bucket = (value - mini) * map_factor;
				//std::cout << std::dec << value << " " << std::hex << value << " " << bucket << " ";
				// the bucket is directly determined by the first 'resolution' bits
				bucket >>= 64-resolution;
				//std::cout << std::dec << bucket << std::endl;
				histogram[bucket]+=increment;
			}
		} else {
			double value = (double) value_p;
			if (value < mind)
				histogram_n_ov+=increment;
			else if (value >= maxd)
				histogram_p_ov+=increment;
			else {
				// map the value to a uint64_t
				uint64_t bucket = (value - mind) * map_factor;
				//std::cout << std::dec << value << " " << std::hex << value << " " << bucket << " ";
				// the bucket is directly determined by the first 'resolution' bits
				bucket >>= 64-resolution;
				//std::cout << std::dec << bucket << std::endl;
				histogram[bucket]+=increment;
			}
		}

	}

	Histogram& operator+=(const Histogram& h) {
		if (resolution != h.resolution
				|| mind > h.mind || mind < h.mind  // exact equality is desired here,
				|| maxd > h.maxd || maxd < h.maxd  // but != results in a compiler warning for floating point types
				|| mini != h.mini || maxi != h.maxi)
			throw std::invalid_argument("Unable to accumulate histograms of different size or of different ranges!");

        for (std::size_t i = 0; i < histogram.size(); ++i) {
			histogram[i] += h.histogram[i];
		}
		histogram_p_ov += h.histogram_p_ov;
		histogram_n_ov += h.histogram_n_ov;
		return *this;
	}

	// returns the histogram data
	const std::vector<uint64_t>& data() const {
		return histogram;
	}

	// returns the values of left and right overflows from the histogram borders
	std::pair<uint64_t,uint64_t> dataOV() const {
		return std::make_pair(histogram_n_ov, histogram_p_ov);
	}

	// sets the contents of the histogram to zero
	void clear() {
		histogram_n_ov = 0;
		histogram_p_ov = 0;
		histogram.clear();
		histogram.resize(1 << resolution);
	}

	// return the left (inclusive) interval border for each bucket
	// the vector is empty if the internal type is no integral type!
	const std::vector<uint64_t>& dataLimitsInt() const {
			return limitsi;
	}

	// return the left (inclusive) interval border for each bucket
	// the vector is empty if the internal type is no floating point type!
	const std::vector<double>& dataLimitsDouble() const {
			return limitsd;
	}

	// return a pair of the internally used minimum and maximum values
	std::pair<uint64_t,uint64_t> minMaxInt() const {
		return std::make_pair(mini,maxi);
	}

	// return a pair of the internally used minimum and maximum values
	std::pair<double,double> minMaxDouble() const {
		return std::make_pair(mind,maxd);
	}

	void writeToFile(std::ofstream &histfile) const {
		histfile << std::setprecision(16) << std::fixed;

		// write header
		histfile << "# ";
		if (std::is_integral<T>())
			histfile << "int\t" << resolution << "\t" << mini << "\t" << maxi << std::endl;
		else
			histfile << "float\t" << resolution << "\t" << mind << "\t" << maxd << std::endl;

		// write underflow if not zero
		bool leading_zero = true;
		if(histogram_n_ov != 0) {
			leading_zero = false;
			if (std::is_integral<T>())
				histfile << "<" << mini << "\t";
			else
				histfile << "<" << mind << "\t";
			histfile << histogram_n_ov << std::endl;
		}

		// find last non-zero bucket
		std::size_t upper_limit;
		if (std::is_integral<T>())
			upper_limit = limitsi.size();
		else
			upper_limit = limitsd.size();
		if(histogram_p_ov == 0) {
			std::size_t i = upper_limit;
			while (i > 0 && histogram[i-1] == 0) {
				i--;
			}
			upper_limit = i;
		}

		// write buckets, starting with first non-zero one and stopping with last non-zero one
		for(std::size_t i = 0; i < upper_limit; i++) {
			if(!leading_zero || histogram[i] != 0) {
				leading_zero = false;
				if (std::is_integral<T>())
					histfile << limitsi[i];
				else
					histfile << limitsd[i];
				histfile << "\t";
				histfile << histogram[i] << std::endl;
			}
		}

		// write overflow if not zero
		if(histogram_p_ov != 0) {
			if (std::is_integral<T>())
				histfile << ">=" << maxi << "\t";
			else
				histfile << ">=" << maxd << "\t";
			histfile << histogram_p_ov << std::endl;
		}
	}

private:
	int resolution;
	double mind;
	double maxd;
	uint64_t mini;
	uint64_t maxi;
	double map_factor;

	std::vector<uint64_t> histogram;

	// to prevent precision problems limits are always stored in double or uint64_t depending on type of T
	std::vector<double> limitsd;
	std::vector<uint64_t> limitsi;

	uint64_t histogram_p_ov = 0;
	uint64_t histogram_n_ov = 0;

};

} /* namespace jcommon */
#endif /* HISTOGRAM_H */
