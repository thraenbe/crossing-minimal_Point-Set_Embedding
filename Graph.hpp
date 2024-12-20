#pragma once
#include "Geometry.hpp"
#include "Node.hpp"

#include <utility>
#include <tuple>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <ranges>
#include <queue>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/basic/SList.h>


using size_t = std::size_t;
struct Edge : public std::pair<size_t, size_t> {
    size_t id;

    // Constructor with id
    Edge(size_t u, size_t v, size_t id) : std::pair<size_t, size_t>(u, v), id(id) {}

    // Constructor without id (default id = 0)
    Edge(size_t u, size_t v) : std::pair<size_t, size_t>(u, v), id(0) {}

};

namespace {
	const std::vector<std::vector<size_t>> InitAdjList(const size_t numNodes, const std::vector<Edge>& edges) {
		//Generate adj list.
		std::cout << "I\n";
		std::vector<std::vector<size_t>>tmpAdjList(numNodes);
		for (const auto& e : edges) {
			auto& u = e.first;
			auto& v = e.second;
			tmpAdjList[u].push_back(v);
			tmpAdjList[v].push_back(u);
		}
	// 	for (size_t i = 0; i < tmpAdjList.size(); ++i) {
	// 		std::cout << "Node " << i << ": ";
	// 		for (const auto& neighbor : tmpAdjList[i]) {
	// 			std::cout << neighbor << " ";
	// 		}
    // std::cout << "\n";
	// 	}
		return tmpAdjList;
	}
}
struct Graph {
	Graph(const size_t numNodes, const std::vector<Node>& nodes, const std::vector<Edge>& edges, const std::vector<Point>& pointSet, const int org_width , const int org_height)
		: numNodes(nodes.size()), numEdges(edges.size()), numPoints(pointSet.size()), points(pointSet),  edges(edges), nodes(nodes), adjList(InitAdjList(numNodes, edges)), org_width(org_width), org_height(org_height) {
		

		width = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._x < b._x;
	    })->_x;

	    height = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._y < b._y;
	    })->_y;


		queueNodes.push(0);
		queuePoints.push(0);
		

		// // Assign all other points as free 
		// for (size_t i = 0; i < numPoints;++i) {
		// 	freePoints.insert(i);
		// }
		std::cout << "m\n";
		for (size_t i = 0; i < numNodes; i++) {
			mapVerticesToPoints[i] = i;
		}
		std::cout << "m\n";
}


		// TODO  O(n²)
			// find closest point to center 
			// sort Points from central point
			// use starting node 
			// map adjacent nodes to closest points

	// adjacencylist representation
	int width; 
	int height;
	const int org_width;
	const int org_height;
	const size_t numNodes;
	const size_t numEdges;
	const size_t numPoints;
	std::vector<Node> nodes;
	std::vector<Point> points;
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
	

	// [[nodiscard]] inline std::tuple<std::queue<size_t>,std::unordered_set<size_t>> AppendQueue( std::queue<size_t> queue,  std::vector<size_t> vector, std::unordered_set<size_t> used){
	// 	for (const auto& element : vector) {
	// 		if (used.contains(element)) {
	// 			queue.push(element);
	// 			used.insert(element);
	// 		}
    // 	}
	// 	return {queue, used};
	// }

	// [[nodiscard]] inline std::vector<size_t> GetAdjacentNodes(const size_t nodeID,  std::unordered_set<size_t> usedNodes){
	// 	std::vector<size_t> adjNodes;
	// 	for (const auto idx : adjList[nodeID]){
	// 		if (usedNodes.contains(idx)) {
	// 			adjNodes.push_back(idx) ;
	// 		}
	// 	}
	// 	return adjNodes;
    // }
	
	// [[nodiscard]] inline std::vector<size_t> GetClosestPoints( const size_t adjSize, const Point& p1, const std::vector<Point>&points ){
	// 	std::vector<std::tuple<double, Point , size_t>> closestPoints;
	//     for (const auto q : freePoints) {
	// 		const Point& p2 = points[q];
	// 		const double dis = L2DistSquared(p1,  p2);
	// 		closestPoints.emplace_back(dis, p2, q);
	// 	}
	// 	std::ranges::sort(closestPoints);

	// 	std::vector<size_t> indices(adjSize);
	// 	int idx = 0;
	// 	int i = 0;
	// 	while(i < adjSize){
	// 		const auto element = std::get<2>(closestPoints[i]);
	// 		if (usedPoints.contains(element))
	// 		{
	// 			indices[idx] = element;
	// 			usedPoints.insert(element);
	// 			freePoints.erase(element);
	// 			idx++;
	// 		}
	// 		i++;
	// 	}
	// 	return indices;
	// }

	[[nodiscard]] inline std::vector<Edge> getAdjEdges(const Node &n) const {
        std::vector<Edge> adjEdges;
		for (auto &e : edges){
			if (e.first == n.GetId() || e.second == n.GetId()){
				adjEdges.push_back(e);
			}
		}
		assert(adjEdges.size() == adjList[n.GetId()].size());
        return adjEdges;
    }


	// Point query.
	[[nodiscard]] inline const Point& GetConstPosOfNode(const size_t nodeId) const {
		return points[mapVerticesToPoints.at(nodeId)];
	}

[[nodiscard]] inline Point NodeToPoint(const size_t nodeId) const {
    const auto& node = nodes[nodeId];  // Annahme: nodes ist ein Container von `Node`-Objekten
    return Point(node.GetId(), node.GetX(), node.GetY());  // Rückgabe des Point-Objekts, ohne eine lokale Referenz zu erzeugen
}

	
	// Edge relations
	[[nodiscard]] int inline DoEdgesIntersect(const Edge& lhs, const Edge& rhs) const {
		const auto& p1 = GetConstPosOfNode(lhs.first);
       	const auto& p2 = GetConstPosOfNode(lhs.second);
		const auto& q1 = GetConstPosOfNode(rhs.first);
       	const auto& q2 = GetConstPosOfNode(rhs.second);
       	auto DoOverlap = [&](const Point& commonPoint, const Point& p, const Point& q) {
          	if (Orientation(commonPoint, p, q) == 0 && !InBoundingBox(p, commonPoint, q)) {
               return numNodes;
           }
           return (size_t) 0;
       };

       if (lhs.first== rhs.first) {
           return DoOverlap(p1,p2,q2);
       }
       if (lhs.first== rhs.second) {
           return DoOverlap(p1,p2,q1);
       }
       if (lhs.second== rhs.first) {
           return DoOverlap(p2,p1,q2);
       }
       if (lhs.second== rhs.second) {
           return DoOverlap(p2,p1,q1);
       }
	return DoIntersect(p1,p2,q1,q2, numNodes);
	}

	[[nodiscard]] int inline DoEdgesIntersectVar(const Edge& lhs, const Edge& rhs) const {
		
		const auto& p1 = NodeToPoint(lhs.first);
		
       	const auto& p2 = NodeToPoint(lhs.second);
		const auto& q1 = NodeToPoint(rhs.first);
       	const auto& q2 = NodeToPoint(rhs.second);
       	auto DoOverlap = [&](const Point& commonPoint, const Point& p, const Point& q) {
          	if (Orientation(commonPoint, p, q) == 0 && !InBoundingBox(p, commonPoint, q)) {
               return numNodes;
           }
           return (size_t) 0;
       };

       if (lhs.first== rhs.first) {
           return DoOverlap(p1,p2,q2);
       }
       if (lhs.first== rhs.second) {
           return DoOverlap(p1,p2,q1);
       }
       if (lhs.second== rhs.first) {
           return DoOverlap(p2,p1,q2);
       }
       if (lhs.second== rhs.second) {
           return DoOverlap(p2,p1,q1);
       }
	return DoIntersect(p1,p2,q1,q2, numNodes);
	}

	
	[[nodiscard]] inline int ComputeCrossings(int& counter) const {
		int crossings{ 0 };
		int tmp = 0;
		counter = 0;
		#pragma omp parallel for reduction (+:crossings) //<- use this if you want to parallelize the code
		for (size_t i = 0; i < numEdges - 1; ++i) {
			const auto& edge = edges[i];
			for (size_t j = i + 1; j < numEdges; ++j) {
				tmp = DoEdgesIntersect(edge,edges[j]);
				crossings += tmp;
				if (tmp > 1){
					counter += 1;
				}
			}	
		}
		// std::cout << " Collinear Crossings: " << counter << std::endl;
		return crossings;
	}

	[[nodiscard]] inline int ComputeEffectiveCrossings(int& counter) const {
		int crossings{ 0 };
		int tmp = 0;
		#pragma omp parallel for reduction (+:crossings) //<- use this if you want to parallelize the code
		for (size_t i = 0; i < numEdges - 1; ++i) {
			const auto& edge = edges[i];
			for (size_t j = i + 1; j < numEdges; ++j) {
				tmp = DoEdgesIntersectVar(edge,edges[j]);
				crossings += tmp;
				if (tmp > 1){
					counter += 1;
				}
			}	
		}
		// std::cout << " Collinear Crossings: " << counter << std::endl;
		return crossings;
	}


	[[nodiscard]] inline int ComputeEffectiveCrossings(const std::vector<size_t>& edgeWeights) const {
		int crossings{ 0 };
		int tmp = 0;
		#pragma omp parallel for reduction (+:crossings) //<- use this if you want to parallelize the code
		for (size_t i = 0; i < numEdges - 1; ++i) {
			const auto& edge = edges[i];
			for (size_t j = i + 1; j < numEdges; ++j) {
				tmp = DoEdgesIntersectVar(edge,edges[j]);
				if (tmp == 1) {
					tmp = edgeWeights[i]*edgeWeights[j]; 
				}
				crossings += tmp;
			}
		}
		return crossings;
	}


	void toOgdfGraph(ogdf::Graph& ogdfGraph, ogdf::GraphAttributes& GA) const {
        std::unordered_map<size_t, ogdf::node> nodeMap;

        // Add nodes
        for (size_t i = 0; i < numNodes; ++i) {
            ogdf::node n = ogdfGraph.newNode();
            nodeMap[i] = n;
			GA.x(n) = nodes[i].GetX();
			GA.y(n) = nodes[i].GetY();
			GA.label(n) = std::to_string(nodes[i].GetId());
			
        }

        // Add edges
        for (const auto& edge : edges) {
            ogdfGraph.newEdge(nodeMap[edge.first], nodeMap[edge.second]);
        }
    }


	[[nodiscard]] inline std::vector<std::vector<Point>> getPointClusters(const std::vector<Point>& points, size_t numClusters){
		std::vector<std::vector<Point>> pointClusters(numClusters);
		for (const auto& point : points){
			// std::cout << point.GetId() << " Cluster: " << point.GetCluster() << "     Postion: (" << point.GetX() << " , " << point.GetY() << ")" << std::endl;
			pointClusters[point.GetCluster()].push_back(point);
		}
		return pointClusters;
	}


	void initializeCentroids(vector<Point>& centroids, const vector<Point>& points, const size_t k) {
		for (int i = 0; i < k; ++i) {
			centroids[i] = points[i];
		}
	}

	
	void balancedAssignClusters(vector<Point>& points, const vector<Point>& centroids, vector<size_t> targetClusterSizes) {
	const size_t k = targetClusterSizes.size();
    vector<vector<pair<double, size_t>>> squaredDistances(points.size(), vector<pair<double, size_t>>(k));
    // Calculate distances from each point to each centroid
    for (int i = 0; i < points.size(); ++i) {
        for (int j = 0; j < k; ++j) {
            squaredDistances[i][j] = {L2DistSquared(points[i], centroids[j]), j};
        }
        std::ranges::sort(squaredDistances);
    }

    // Create a balanced assignment
    vector<size_t> clusterSizes(k, 0);

    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = 0; j < k; ++j) {
            const auto cluster = squaredDistances[i][j].second;
            if (clusterSizes[cluster] < targetClusterSizes[cluster]) {
                points[i].SetCluster(cluster);
                clusterSizes[cluster]++;
                break;
            }
        }
    }
}


	void updateCentroids(vector<Point>& centroids, const vector<Point>& points, const size_t k) {
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
			if (L2DistSquared(centroids[i], oldCentroids[i]) > 1) {
				return false;
			}
		}
		return true;
	}

	void kMeans( vector<Point>& points, const int k, const vector<size_t>& clusterSizes) {
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

int randomCrossingNumber(int& counter){
	for(int i = 0; i < nodes.size(); i++){
		mapVerticesToPoints[i] = i;
	}
	return ComputeCrossings(counter);
}

/*
	assign Clusters for final Local Clustering from MatchCluster Values
*/
std::vector<size_t> assignClusters(Graph& G, const std::vector<size_t>& clusterSizes){
	std::vector<size_t> newClusterSizes(clusterSizes.size());
	std::vector<vector<Node>> newNodeClusters(G.NodeClusters.size()) ;
    std::cout << mapVerticesToPoints.size() << " " << NodeClusters.size() << " " << G.NodeClusters.size()<< " "  << clusterSizes.size() << std::endl;

 	for (int i = 0; i < mapVerticesToPoints.size(); i++){
		newNodeClusters[i] = G.NodeClusters[mapVerticesToPoints[i]];
		newClusterSizes[i] = clusterSizes[mapVerticesToPoints[i]];
	}
	G.NodeClusters = newNodeClusters;
	NodeClusters = newNodeClusters;
	return newClusterSizes;
}


void sortDirectionVectors(std::vector<vector<Point>>& directions){

	int xmin = 0;
	int ymin = 0;
	int ymax = height;
	int xmax = width;

	//sort points of vector N by distance of Y to Ymax
	std::ranges::sort(directions[4], [ymax](const Point& a, const Point& b) {
        return a.distanceToLineY(ymax) < b.distanceToLineY(ymax);
    });
	Point reference(0, xmax, ymax);
    // Sort points of vector NE by distance of each point to (10, 10)
    std::ranges::sort(directions[0], [&reference](const Point& a, const Point& b) {
        return a.squaredDistanceTo(reference) < b.squaredDistanceTo(reference);
    });
	std::ranges::sort(directions[5], [xmax](const Point& a, const Point& b) {
        return a.distanceToLineX(xmax) < b.distanceToLineX(xmax);
    });

    // Sort points of vector NE by distance of each point to (10, 10)
	reference._y = ymin;
    std::ranges::sort(directions[1], [&reference](const Point& a, const Point& b) {
        return a.squaredDistanceTo(reference) < b.squaredDistanceTo(reference);
    });

	std::ranges::sort(directions[6], [ymin](const Point& a, const Point& b) {
        return a.distanceToLineY(ymin) < b.distanceToLineY(ymin);
    });

    // Sort points of vector NE by distance of each point to (10, 10)
	reference._x = xmin;
    std::ranges::sort(directions[2], [&reference](const Point& a, const Point& b) {
        return a.squaredDistanceTo(reference) < b.squaredDistanceTo(reference);
    });
  
	std::ranges::sort(directions[7], [xmin](const Point& a, const Point& b) {
    
        return a.distanceToLineX(xmin) < b.distanceToLineX(xmin);
    });

	reference._y = ymax;
    std::ranges::sort(directions[3], [&reference](const Point& a, const Point& b) {
        return a.squaredDistanceTo(reference) < b.squaredDistanceTo(reference);
    });
}
  
  
	/**
	 * @brief Clusters all points locally into the given clusterSizes.
	 *
	 * @param clusterSizes Vector with size of number of Clusters.
	 * @param xmax width.
	 * @param ymax height.
	 */
void manClustering(vector<size_t> clusterSizes, int xmax, int ymax){
	// sort clustersizes ascending
	std::cout << xmax << " " << ymax << "\n";


	
	int xmin = 0;
	int ymin = 0;

	std::vector<vector<Point>> directions(8);


	directions[0] = points;

	directions[1] = points;
	directions[2] = points; 
	directions[3] = points;
	directions[4] = points;
	directions[5] = points;
	directions[6] = points; 
	directions[7] = points;

	std::cout << "Sort Direction Vectors \n";
	sortDirectionVectors(directions);
	std::cout << "Finished Sorting Direction Vectors \n";

	// Print the entire matrix
    // for (int i = 0; i < 8; ++i) {
    //     std::cout << "Direction " << i << ":\n";
    //     for (const auto& point : directions[i]) {
    //         std::cout << point.GetId() << " ";
    //     }
    //     std::cout << "\n";
    // }
	// std::cout << "Points: \n";
    //     for (const auto& point : points) {
    //         std::cout << point.GetId() << " ";
    //     }
    //     std::cout << "\n";


	assert(freePoints.empty());

	assert(points.size() == numPoints);

	for (size_t i = 0; i < points.size();++i) {
		freePoints.insert(i).second;
	}

	assert(freePoints.size() == numPoints);

	std::cout << "FreePoints reassigned \n";

    const size_t count = std::accumulate(clusterSizes.begin(),clusterSizes.end(),(size_t)0);
		std::cout << count << points.size() <<  " ----------------------------- \n";
	assert(count == points.size());

	int direction = 0;
	std::vector<Point> pointsInCurrentCluster;
	for(auto i = 0; i < clusterSizes.size(); ++i){
		direction = i % 8;
		const auto& currentDirection = directions[direction];
		assert(directions[direction].size() == points.size());
		size_t s = 0;
		// std::cout << "Starting While Loop: " << i << std::endl;
		for(auto idx = 0; idx < points.size(); ++idx){
			assert(currentDirection.size() > idx);
			assert(s < clusterSizes[i]);
			const auto& currentPointId = currentDirection[idx].GetId();
			auto& currentPoint = points[currentPointId];
			//assert(currentPointId == points[currentPointId].GetId());
			if (freePoints.contains(currentPointId)){
				// if (clusterSizes[i] - s < freePoints.size()) {
				// 	if (isCollinear(currentPoint, pointsInCurrentCluster )){
				// 		freePoints.erase(currentPointId);
				// 		currentPoint.SetCluster(i);
				// 		pointsInCurrentCluster.push_back(currentPoint);
				// 		s++;
				// 	}
				// }else{ 
				freePoints.erase(currentPointId);
				currentPoint.SetCluster(i);
				pointsInCurrentCluster.push_back(currentPoint);
				s++;
				//}
			}
			if (s == clusterSizes[i]){
				break;
			}
		}
	}
	std::cout << "Clustering Complete" << freePoints.size() << " " << clusterSizes[0] << " " << clusterSizes.size() << std::endl;
	assert(freePoints.empty());

}



	/**
	 * @brief Reduces the number of points in a cluster.
	 *
	 * This function reduces the number of points in a cluster by deleting collinear points.
	 * @param points the points of a certain cluster
	 */
	void reductCluster(const std::vector<Point>& points, const int newNumPoints){
		vector<vector<double>> radialGradient(points.size(), vector<double> (points.size(),std::numeric_limits<double>::infinity())) ;
		const size_t numPoints = points.size();
		for(int i = 0; i < numPoints-1; ++i){
			for(int j = i+1; j < numPoints; ++j) {
				const auto grad = computeGradient(points[i],points[j]);
				radialGradient[i][j] = grad;
				radialGradient[j][i] = -grad;
			}
		}
		// TODO: missing functionality
		std::unordered_map<double, int> countMap;
		for (int i = 0; i < numPoints; ++i){
			for (int j = 0; j < numPoints; ++j) {
				countMap[radialGradient[i][j]]++;
			}
			for (const auto& pair : countMap) {
				if(pair.second != (double)1){
					// std::cout << i << " Value: " << pair.first << ", Count: " << pair.second << " ;" << std::endl;
				}
			
			}
			countMap.clear();
		}
	}

	std::vector<vector<Node>> assignClustersToNodes(const ogdf::SList<ogdf::SimpleCluster *>& clusters){
		int clusterIndex = 0;
		std::vector<vector<Node>> nodeClusters(clusters.size());
		for (auto* cluster : clusters) {
			for (const ogdf::node& v : cluster->nodes()) {
				Node &node = nodes[v->index()];
				node.SetCluster(clusterIndex);
				nodeClusters[clusterIndex].push_back(node);
			}
			
			clusterIndex++;
		}

		return nodeClusters;


    }
};




