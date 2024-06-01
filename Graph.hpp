#pragma once
#include <utility>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include "Geometry.hpp"
#include <cassert>

using size_t = std::size_t;
using Edge = std::pair<size_t,size_t>;
namespace {
	const std::vector<std::vector<size_t>> InitAdjList(const size_t numNodes, const std::vector<Edge>& edges) {
		//Generate adj list.
		std::vector<std::vector<size_t>>tmpAdjList(numNodes);
		for (const auto& [u, v] : edges) {
			tmpAdjList[u].push_back(v);
			tmpAdjList[v].push_back(u);
		}
		return tmpAdjList;
	}
}
struct Graph {
	Graph(const size_t numNodes, const std::vector<Edge>& edges, const std::vector<Point>& pointSet)
		: numNodes(numNodes), numEdges(edges.size()), numPoints(pointSet.size()), points(pointSet), edges(edges), adjList(InitAdjList(numNodes, edges)) {
		// Simple assignment to first n points
		for (size_t i = 0; i < numNodes;++i) {
			mapVerticesToPoints[i] = i;
			usedPoints.insert(i);
		
		};
		for (size_t i = numNodes; i < numPoints;++i) {
			freePoints.insert(i);
		}
		}

	// adjacencylist representation
	const size_t numNodes;
	const size_t numEdges;
	const size_t numPoints;
	const std::vector<Point>points;
	const std::vector<Edge>edges;
	const std::vector<std::vector<size_t>>adjList;

	std::unordered_set<size_t>freePoints;
	std::unordered_set<size_t>usedPoints;
	std::unordered_map<size_t, size_t> mapVerticesToPoints;

	// Point query.
	[[nodiscard]] inline const Point& GetConstPosOfNode(const size_t nodeId) const {
		return points[mapVerticesToPoints.at(nodeId)];
	}
	
	// Edge relations
	[[nodiscard]] bool inline AreAdjacent(const Edge& lhs, const Edge& rhs) const {
		return lhs.first == rhs.first || lhs.first == rhs.second || rhs.first == lhs.second || rhs.second == lhs.second;
	}
	// TODO: Include weighted crossing in case of overlaps.
	[[nodiscard]] int inline DoEdgesIntersect(const Edge& lhs, const Edge& rhs) const {
		if (AreAdjacent(lhs, rhs)) {
			return 0;
		}
		return DoIntersect(GetConstPosOfNode(lhs.first), GetConstPosOfNode(lhs.second), GetConstPosOfNode(rhs.first), GetConstPosOfNode(rhs.second), numNodes);
	}
	
	[[nodiscard]] inline int ComputeCrossings() const {
		int crossings{ 0 };
		// #pragma omp parallel for reduction (+:crossings) <- use this if you want to parallelize the code
		for (size_t i = 0; i < numEdges - 1; ++i) {
			const auto& edge = edges[i];
			for (size_t j = i + 1; j < numEdges; ++j) {
				crossings += DoEdgesIntersect(edge,edges[j]);
		}
		
	}
	return crossings;
	}
	
	};
