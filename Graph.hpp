#pragma once
#include <utility>
#include <tuple>
#include <map>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include "Geometry.hpp"
#include "Node.hpp"
#include <cassert>
#include <iostream>
#include <list>
#include <queue>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/basic/SList.h>

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
	Graph(const size_t numNodes, const std::vector<Node>& nodes, const std::vector<Edge>& edges, const std::vector<Point>& pointSet)
		: numNodes(nodes.size()), numEdges(edges.size()), numPoints(pointSet.size()), points(pointSet), edges(edges), nodes(nodes), adjList(InitAdjList(numNodes, edges)) {

		

		queueNodes.push(0);
		queuePoints.push(0);
		

		// Assign all other points as free 
		for (size_t i = 0; i < numPoints;++i) {
			freePoints.insert(i);
		}


		while (!queueNodes.empty()) {


			size_t nodeID = queueNodes.front();
			size_t pointID = queuePoints.front();
			queueNodes.pop();
			queuePoints.pop();
			mapVerticesToPoints[nodeID] = pointID ;

			std::vector<size_t> adjNodes = GetAdjacentNodes(nodeID,usedNodes);
			std::tuple<std::queue<size_t>,std::unordered_set<size_t>> tmpNodes = AppendQueue(queueNodes, adjNodes, usedNodes);
			queueNodes = std::get<0>(tmpNodes);
			usedNodes = std::get<1>(tmpNodes);
			std::tuple<std::queue<size_t>,std::unordered_set<size_t>> tmp = AppendQueue( queuePoints, GetClosestPoints(adjNodes.size(), points[pointID], points) ,usedPoints);
			queuePoints = std::get<0>(tmp); 
			usedPoints =  std::get<1>(tmp);
		}
		
		}


		// TODO  O(n²)
			// find closest point to center 
			// sort Points from central point
			// use starting node 
			// map adjacent nodes to closest points

	// adjacencylist representation
	const size_t numNodes;
	const size_t numEdges;
	const size_t numPoints;
	std::vector<Node> nodes;
	std::vector<Point>points;
	const std::vector<Edge>edges;
	const std::vector<std::vector<size_t>>adjList;


	std::unordered_set<size_t>freePoints;
	std::unordered_set<size_t>usedPoints;
	std::unordered_set<size_t>usedNodes;
	std::unordered_map<size_t, size_t> mapVerticesToPoints;
	std::vector<Point> orderedPoints;

	std::queue<size_t> queuePoints ;
	std::queue<size_t> queueNodes ; 

	std::vector<std::vector<Point>> pointClusters;
    std::vector<std::vector<Node>> NodeClusters;

	[[nodiscard]] inline std::tuple<std::queue<size_t>,std::unordered_set<size_t>> AppendQueue( std::queue<size_t> queue,  std::vector<size_t> vector, std::unordered_set<size_t> used){
		for (const auto& element : vector) {
			if (used.find(element) == used.end() ) {
				queue.push(element);
				used.insert(element);
			}
    	}
		return std::make_tuple(queue, used) ;
	}


	[[nodiscard]] inline std::vector<size_t> GetAdjacentNodes( size_t nodeID,  std::unordered_set<size_t> usedNodes){
		std::vector<size_t> adjNodes;
		// std::cout << " GetAdjacentNodes " << std::endl;
		for (const auto& idx : adjList[nodeID]){
			if (usedNodes.find(idx) == usedNodes.end() ) {
				adjNodes.push_back(idx) ;
			}
		}
		return adjNodes;
    }
	


	[[nodiscard]] inline std::vector<size_t> GetClosestPoints( size_t adjSize, const Point& p1, const std::vector<Point>points ){
		std::vector<std::tuple<float, Point , size_t>> closestPoints;
	    for (const auto& q : freePoints) {
			Point p2 = points[q];
			float dis = euclideanDistance(p1,  p2);
			closestPoints.push_back(std::make_tuple(dis, p2, q));
		}
		std::sort(closestPoints.begin(), closestPoints.end());

		// only printing
		for(const auto& tuple : closestPoints) {
        	// std::cout << "Tuple: " << std::get<0>(tuple) << ", " << std::get<2>(tuple) << std::endl;
			
    	}
		std::vector<size_t> indices(adjSize);
		int idx = 0;
		int i = 0;
		while(i < adjSize){
			int element = std::get<2>(closestPoints[i]);
			if (usedPoints.find(element) == usedPoints.end() )
			{
				indices[idx] = element;
				usedPoints.insert(element);
				freePoints.erase(element);
				idx++;
			}
			i++;
		}
		return indices;
	}


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

	void toOgdfGraph(ogdf::Graph& ogdfGraph) const {
        std::unordered_map<size_t, ogdf::node> nodeMap;

        // Add nodes
        for (size_t i = 0; i < numNodes; ++i) {
            ogdf::node n = ogdfGraph.newNode();
            nodeMap[i] = n;
        }

        // Add edges
        for (const auto& edge : edges) {
            ogdfGraph.newEdge(nodeMap[edge.first], nodeMap[edge.second]);
        }
    }


	[[nodiscard]] inline std::vector<std::vector<Point>> getPointClusters(std::vector<Point>points, int numClusters){
		std::vector<std::vector<Point>> pointClusters(numClusters);
		for (auto point : points){
			pointClusters[point.GetCluster()].push_back(point);
		}
		return pointClusters;
	}


    // [[nodiscard]] inline std::vector<std::vector<Point>> getNodeClusters(ogdf::SList<ogdf::SimpleCluster *> clusters){

	// }


	void initializeCentroids(vector<Point>& centroids, const vector<Point>& points, int k) {
		for (int i = 0; i < k; ++i) {
			centroids[i] = points[i];
		}
	}

	
	void balancedAssignClusters(vector<Point>& points, const vector<Point>& centroids, vector<int> targetClusterSizes) {
	int k = targetClusterSizes.size();
    vector<vector<pair<double, int>>> distances(points.size(), vector<pair<double, int>>(k));

    // Calculate distances from each point to each centroid
    for (int i = 0; i < points.size(); ++i) {
        for (int j = 0; j < k; ++j) {
            distances[i][j] = {euclideanDistance(points[i], centroids[j]), j};
        }
        sort(distances[i].begin(), distances[i].end());
    }

    // Create a balanced assignment
    vector<int> clusterSizes(k, 0);

    for (int i = 0; i < points.size(); ++i) {
        for (int j = 0; j < k; ++j) {
            int cluster = distances[i][j].second;
            if (clusterSizes[cluster] < targetClusterSizes[cluster]) {
                points[i].SetCluster(cluster);
                clusterSizes[cluster]++;
                break;
            }
        }
    }
}





	void updateCentroids(vector<Point>& centroids, const vector<Point>& points, int k) {
		vector<int> count(k, 0);
		vector<Point> newCentroids(k, {0, 0, 0});
		
		for (const auto& point : points) {
			newCentroids[point._cluster]._x += point._x;
			newCentroids[point._cluster]._y += point._y;
			count[point._cluster]++;
		}
		
		for (int i = 0; i < k; ++i) {
			if (count[i] != 0) {
				newCentroids[i]._x /= count[i];
				newCentroids[i]._y /= count[i];
			}
		}
		
		centroids = newCentroids;
	}

	bool hasConverged(const vector<Point>& centroids, const vector<Point>& oldCentroids) {
		for (int i = 0; i < centroids.size(); ++i) {
			if (euclideanDistance(centroids[i], oldCentroids[i]) > 1) {
				return false;
			}
		}
		return true;
	}

	void kMeans( vector<Point>& points, int k, vector<int> clusterSizes) {
		vector<Point> centroids(k);
		initializeCentroids(centroids, points, k);

		vector<Point> oldCentroids(k);

		for(int tmp = 0; tmp < 10000; tmp++ ) {
			oldCentroids = centroids;

			balancedAssignClusters(points, centroids, clusterSizes);
			updateCentroids(centroids, points, k);

			if (hasConverged(centroids, oldCentroids)) {
				break;
			}
		}
	}

	/**
	 * @brief Reduces the number of points in a cluster.
	 *
	 * This function reduces the number of points in a cluster by deleting collinear points.
	 * @param points the points of a certain cluster
	 */
	void reductCluster(std::vector<Point> points, int newNumPoints){
		vector<vector<double>> radialGradient(points.size(), vector<double> (points.size())) ; 
		for(int i = 0; i < points.size() ; i++){ 
			for(int j = 0; j < points.size(); j++) { 
				if (i == j){
					radialGradient[i][j] = - std::numeric_limits<double>::infinity();
				}
				else{ 
				radialGradient[i][j] = computeGradient(points[i],points[j]);
				}
			}
		}
		std::unordered_map<double, int> countMap;
		for (int i = 0; i < points.size(); ++i){
			for (int j = 0; j < points.size(); ++j) {
				countMap[radialGradient[i][j]]++;
			}
			for (const auto& pair : countMap) {
				if(pair.second != (double)1){
					std::cout << i << " Value: " << pair.first << ", Count: " << pair.second << " ;" << std::endl;
				}
			
			}
			countMap.clear();
		}
	}

	std::vector<vector<Node>> assignClustersToNodes(ogdf::SList<ogdf::SimpleCluster *> clusters){
		int clusterIndex = 0;
		std::vector<vector<Node>> nodeClusters(clusters.size());
		std::vector<int> clusterSizes; 
		for (auto cluster : clusters) {
			clusterSizes.push_back(cluster->m_size);

			for (ogdf::node v : cluster->nodes()) {
				nodes[v->index()].setCluster(clusterIndex);
				nodeClusters[clusterIndex].push_back(nodes[v->index()]);
			}
			
			clusterIndex++;
		}
		return nodeClusters;


    }


};

