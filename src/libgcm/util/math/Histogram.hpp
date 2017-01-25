#ifndef LIBGCM_HISTOGRAM_HPP
#define LIBGCM_HISTOGRAM_HPP

#include <algorithm>

#include <libgcm/util/infrastructure/infrastructure.hpp>

namespace gcm {

class Histogram {
public:
	template<typename ForwIter>
	Histogram(const ForwIter begin, const ForwIter end, const size_t numberOfBins) {
		assert_true(begin != end);
		_min = (real)(*std::min_element(begin, end));
		_max = (real)(*std::max_element(begin, end));
		assert_gt(numberOfBins, 1);
		const real _binSize = (_max - _min) / real(numberOfBins);
		
		if (max() == min()) {
			_bins.resize(numberOfBins, 0);
			for (ForwIter it = begin; it != end; ++it) {
				++(_bins[0]);
			}
		} else {
			_bins.resize(numberOfBins + 1, 0);
			for (ForwIter it = begin; it != end; ++it) {
				++(_bins[(size_t)((*it - min()) / _binSize)]);
			}
			*std::next(_bins.rbegin()) += _bins.back();
			_bins.pop_back();
		}
	}
	
	real min() const { return _min; }
	
	real max() const { return _max; }
	
	size_t binsNumber() const { return _bins.size(); }
	
	const std::vector<size_t>& binCounts() const { return _bins; }
	
	real binSize() const { return (max() - min()) / (real)binsNumber(); }
	
	std::vector<real> binCenters() const {
		std::vector<real> ans(binsNumber());
		for (size_t i = 0; i < binsNumber(); i++) {
			ans[i] = min() + (real(i) + 0.5) * binSize();
		}
		return ans;
	}
	
	real mean() const {
		return std::inner_product(
			_bins.begin(), _bins.end(), binCenters().begin(), real(0)) /
				std::accumulate(_bins.begin(), _bins.end(), real(0));
	}
	
	friend std::ostream& operator<<(std::ostream& os, const Histogram& hist) {
		os << "binCenter" << "\t" << "binCount" << std::endl;
		const std::vector<real>  values = hist.binCenters();
		const std::vector<size_t>& counts = hist.binCounts();
		assert_eq(values.size(), counts.size());
		for (size_t i = 0; i < values.size(); i++) {
			os << values[i] << "\t" << counts[i] << std::endl;
		}
		return os;
	}
	
	
private:
	real _min, _max;
	std::vector<size_t> _bins;
	
};

}

#endif // LIBGCM_HISTOGRAM_HPP
