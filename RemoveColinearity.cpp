#include <ranges>
#include <algorithm>
#include <map>
#include <numeric>
#include <unordered_set>
#include <unordered_map>

using Triplet = std::tuple<size_t, size_t, size_t>;
constexpr std::vector<Triplet> ComputeColinearTriplets(const std::vector<Point>& pointset) {
	std::vector<Triplet> colinearTriplets;
	int counter{ 0 };
	const size_t numPoints{ pointset.size() };
	// tetrahedal number n*n+1*n+2/6 -> number of triples we consider
	for (size_t i = 0; i < numPoints - 2; i++) {
		for (size_t j = i + 1; j < numPoints - 1; j++) {
			for (size_t k = j + 1; k < numPoints; k++) {
				if (Orientation(pointset[i], pointset[j], pointset[k]) == 0) {
					colinearTriplets.emplace_back(i, j, k);
					++counter;
				}
			}
		}
	}
	return colinearTriplets;
}

std::vector<size_t> RemoveCollinearityFromPointset(const std::vector<Point>& pointset) {
	const size_t numPoints{pointset.size()};
	std::unordered_set<size_t> deletedIndices;
	auto colinearTriplets = ComputeColinearTriplets(pointset);
	std::vector<std::tuple<size_t, size_t, size_t>> remainingTriplets;
	while (!colinearTriplets.empty()) {
		remainingTriplets.clear();
		
		std::vector<int> occurences(numPoints, 0);
		// Determine max occurence
		for (auto const& [indexOne, indexTwo, indexThree] : colinearTriplets) {
			occurences[indexOne] += 1;
			occurences[indexTwo] += 1;
			occurences[indexThree] += 1;
		}
		int max_val = 0;
		int max_idx{ -1 };
		for (int i = 0; i < numPoints; i++) {
			if (occurences[i] > max_val) {
				max_val = occurences[i];
				max_idx = i;
			}
		}
		// Delte max occured point.
		if (max_idx != -1) {
			deletedIndices.insert(max_idx);
			for (auto const& [p1, p2, p3] : colinearTriplets) {
				if (!((p1 == max_idx) || (p2 == max_idx) || (p3 == max_idx))) {
					remainingTriplets.emplace_back(p1, p2, p3);
				}
			}
		}
		colinearTriplets = remainingTriplets;
	}
	std::vector<size_t> remainingIndices;
	for (size_t i = 0; i < pointset.size(); ++i) {
		if (!deletedIndices.contains(i)) {
			remainingIndices.push_back(i);
		}
	}
	return remainingIndices;
}


