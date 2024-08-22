#pragma once
#include "Graph.hpp"

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/geometry.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterPlanarity.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/ClustererModule.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/PlanarizationGridLayout.h>
#include <ogdf/planarity/PlanarSubgraphModule.h>
#include <ogdf/planarity/PlanarizerMixedInsertion.h>
#include <ogdf/planarity/CrossingMinimizationModule.h>
#include <ogdf/energybased/SpringEmbedderKK.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>

#include <algorithm> // For std::find
#include <limits>
#include <thread>
#include <future>
#include <chrono>



// Define the Boost graph type

typedef boost::property< boost::edge_weight_t, float, boost::property< boost::edge_index_t, int > >EdgeProperty;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty > BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor VertexDescriptor;


class Minimizer : public ogdf::PlanarizerMixedInsertion {
public:
    void callDoCall(ogdf::PlanRep& pr, int cc, const ogdf::EdgeArray<int>* pCostOrig,
                          const ogdf::EdgeArray<bool>* pForbiddenOrig, const ogdf::EdgeArray<uint32_t>* pEdgeSubGraphs,
                          int& crossingNumber) {
        doCall(pr, cc, pCostOrig, pForbiddenOrig, pEdgeSubGraphs, crossingNumber);
    }

    std::vector<ogdf::DPoint> getNodePositions(const ogdf::GraphAttributes &GA) const {
        std::vector<ogdf::DPoint> positions;
        std::cout << GA.constGraph().nodes.size() << std::endl;
        for (ogdf::node v : GA.constGraph().nodes) {
            positions.push_back(GA.point(v));
            std::cout << "X: " << GA.point(v).m_x << "   Y: " << GA.point(v).m_y << std::endl;
        }
        return positions;
    }
};




void createClusterGraph(ogdf::Graph& G, ogdf::GraphAttributes& GA, Graph& myGraph, std::vector<Node> nodes ){
    std::vector<size_t> tmp;
    std::unordered_map<int, ogdf::node> nodeMap;
    for (auto node : nodes){
        ogdf::node newNode = G.newNode();
        GA.x(newNode) = node.getX();
        GA.y(newNode) = node.getY();
        GA.label(newNode) = std::to_string(node.getId());
        tmp.push_back(node.getId());
        nodeMap[node.getId()] = newNode;
        
    }
    for (int i = 0; i < nodes.size();i++){
        for (auto target : myGraph.adjList[nodes[i].getId()]){
            auto it = std::find(tmp.begin(), tmp.end(), target);
            if (it != tmp.end()){
                if (nodes[i].getId() < target){ // Every Edge is only appended once istead of twice
                    G.newEdge(nodeMap[nodes[i].getId()], nodeMap[target]);
                }
            }
        }
    }
}

void springEmbedding(ogdf::Graph& G, ogdf::GraphAttributes& GA){
    // Instantiate the Kamada-Kawai spring embedder
    ogdf::SpringEmbedderKK KK;

    // Set parameters (optional)
    // KK.setMaxGlobalIterations(10);
    // KK.setStopTolerance(1e-4);
    // KK.setDesLength(30.0); // Set desired edge length

    std::cout << "Positions before Spring Embedder:" << std::endl;
    for (ogdf::node n : G.nodes) {
        std::cout << "Node " << GA.label(n) << " position: (" << GA.x(n) << ", " << GA.y(n) << ");    " ;
    }
    std::cout << std::endl;
    // Apply the layout algorithm
    if (ogdf::isConnected(G)){        
        KK.call(GA);
    }

    std::cout << std::endl;

    // Print node positions
    std::cout << "Positions after Spring Embedder:" << std::endl;
    for (ogdf::node n : G.nodes) {
        std::cout << "Node " << GA.label(n) << " position: (" << GA.x(n) << ", " << GA.y(n) << ");    " ;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    
}
// void crossMin(ogdf::Graph& G, ogdf::GraphAttributes& GA, std::vector<Node> nodes){
//     double maxXNodes = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
//         return a._x < b._x;
//     })->_x;

//     double maxYNodes = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
//         return a._y < b._y;
//     })->_y;

//     iterativeCrossMin(G, GA, maxXNodes, maxYNodes, 100, 100 );
// }


long getWeight(Point point, Node node, int num){
    long weight;
    weight =  (long) num * (1 / (euclideanDistance(point, node.nodeToPoint()) + 1));
    return weight;
}

std::vector<Node> allignBoundingBoxes(std::vector<Point> points, std::vector<Node> nodes) {
    if (points.empty() || nodes.empty()) return nodes;
    // Get min x and y of nodes
    double minXNodes = std::min_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._x < b._x;
    })->_x;

    double minYNodes = std::min_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._y < b._y;
    })->_y;

    // Get max x and y of nodes
    double maxXNodes = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._x < b._x;
    })->_x;

    double maxYNodes = std::max_element(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) {
        return a._y < b._y;
    })->_y;

    // Get min x and y of points
    double minXPoints = std::min_element(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a._x < b._x;
    })->_x;

    double minYPoints = std::min_element(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a._y < b._y;
    })->_y;

    // Get max x and y of points
    double maxXPoints = std::max_element(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a._x < b._x;
    })->_x;

    double maxYPoints = std::max_element(points.begin(), points.end(), [](const Point& a, const Point& b) {
        return a._y < b._y;
    })->_y;

    // Align and stretch points
    for (auto& point : nodes) {
        // Align points to minimal x and y of nodes
        point._x = point._x - minXPoints + minXNodes;
        point._y = point._y - minYPoints + minYNodes;

        // Stretch points horizontally
        double widthNodes = maxXNodes - minXNodes;
        double widthPoints = maxXPoints - minXPoints;
        if (widthPoints != 0) {
            point._x = minXNodes + (point._x - minXNodes) * (widthNodes / widthPoints);
        }

        // Stretch points vertically
        double heightNodes = maxYNodes - minYNodes;
        double heightPoints = maxYPoints - minYPoints;
        if (heightPoints != 0) {
            point._y = minYNodes + (point._y - minYNodes) * (heightNodes / heightPoints);
        }
        point._x = (int) point._x;
        point._y = (int) point._y;
    }

    return nodes;
}



BoostGraph convertClusterToBoostGraph(std::vector<Point> points, std::vector<Node> nodes, int num){
    
    std::vector<Node> allignedNodes;
    // std::cout << "Point " << points[0].GetX() << points[0].GetY() << "Point " << points[1].GetX() << points[1].GetY() << std::endl ; 
    allignedNodes = allignBoundingBoxes(points, nodes);
    // std::cout << "Point " << allignedPoints[0].GetX() << allignedPoints[0].GetY() << "Point " << allignedPoints[1].GetX() << allignedPoints[1].GetY() << std::endl ; 
    std::vector<VertexDescriptor> tmp;
    BoostGraph BoostGraph;
    for (auto point : points) {
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        tmp.push_back(v);
        // std::cout << "Added vertex " << v << " for Point(" << point._x << ", " << point._y << ")\n";
    }
    for (int j = points.size(); j < allignedNodes.size() + points.size() ; j++){
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        for (int i = 0; i < tmp.size(); i++){
            long weight = getWeight(points[i], allignedNodes[j-points.size()], num);
            boost::add_edge(j , i , EdgeProperty(weight), BoostGraph);
            // std::cout << "Match: " << j << " " << i << "   Weight: " << weight << "  Node: " << nodes[j-points.size()].getX() << " " << nodes[j-points.size()].getY() << " " << "Point:  "<< allignedNodes[i].getX() << " " << allignedNodes[i].getY() << std::endl;
        }
    }
    return BoostGraph;
}

std::vector<size_t> minimumWeightMatching(BoostGraph& BoostGraph, int clusterSize){
    std::vector<VertexDescriptor> mate(boost::num_vertices(BoostGraph));
    std::atomic<bool> cancelFlag(false);  // Shared flag for cancellation
    std::promise<void> promise;
    std::future<void> future = promise.get_future();
    std::thread t([&promise, &BoostGraph, &mate]() {
    try {
        boost::maximum_weighted_matching(BoostGraph, &mate[0]);
        promise.set_value();  // Signal that the function finished successfully
        std::cout << "Max weight matching completed." << std::endl;
    } catch (const std::exception& e) {
        promise.set_exception(std::current_exception());
        std::cerr << "Exception during matching: " << e.what() << std::endl;

    }
    });

    // Wait for the function to finish or timeout after 10 seconds
    if (future.wait_for(std::chrono::seconds(10)) == std::future_status::timeout) {
        std::cout << "Function didn't finish in 10 seconds, aborting..." << std::endl;
        cancelFlag.store(true);

        // Join the thread after cancellation

        // Set mate to fallback solution

        for (int i = 0; i < clusterSize; i++) {
            mate[i] = i;
        }
        for (int i = 0; i < clusterSize; i++) {
            mate[i + clusterSize] = i;
        }

        std::cout << "Fallback solution applied." << std::endl;
    } else {
        std::cout << "Function finished in time." << std::endl;
        t.join();  // Wait for the thread to finish
        std::cout << "Found a weighted matching:" << std::endl;
        std::cout << "Matching size is " << boost::matching_size(BoostGraph, &mate[0]) << ", total weight is " << boost::matching_weight_sum(BoostGraph, &mate[0]) << std::endl;
        std::cout << std::endl;

    }

    return mate;
}

void mapMatching(std::vector<size_t> mate, Graph& myGraph, int cluster){
//     int i = 0;
//     for (auto m : mate){
//         std::cout << i++ << " " << m << std::endl; 
//     }
    for (int i = 0; i < myGraph.NodeClusters[cluster].size(); i++){
        int nodeId = myGraph.NodeClusters[cluster][i].getId();
        int pointId = myGraph.pointClusters[cluster][mate[i+myGraph.pointClusters[cluster].size()]].GetId() ;
        myGraph.mapVerticesToPoints[nodeId] = pointId; 
        // std::cout << "Matched Node: " << nodeId << " with Point: " << pointId << std::endl; 
    }
}

void ogdfToNodes(ogdf::Graph& G, ogdf::GraphAttributes& GA, Graph& myGraph, int cluster){
    int i = 0;
    std::cout << "NodeClusters " << myGraph.NodeClusters.size() << "\n" << " " << myGraph.NodeClusters[0].size() << std::endl;
    for (ogdf::node v : G.nodes){
        myGraph.NodeClusters[cluster][i]._x = GA.x(v);
        myGraph.NodeClusters[cluster][i]._y = GA.y(v);
        if (myGraph.NodeClusters[cluster][i]._id != std::stoi(GA.label(v))){
            std::cout << "IDs Do Not match \n";
        }
        i++;
    }
}

std::vector<int> computeVirtualClusterSizes(int numberOfClusters, int numberOfNodes){
    assert(numberOfClusters < numberOfNodes);

    int rest = numberOfNodes % numberOfClusters;
    int q = (numberOfNodes - rest) / numberOfClusters;

    std::vector<int> virtualClusterSizes(numberOfClusters);
    for (int i = 0; i < numberOfClusters; i++){
        if (rest > 0){
            virtualClusterSizes[i] = q+1;
            rest--;
        }
        else {
            virtualClusterSizes[i] = q;
        }
    }
    return virtualClusterSizes;
}

std::vector<Point> computeCentroids(std::vector<vector<Point>> pointClusters){
    std::vector<Point> centroids(pointClusters.size());
    double sumX = 0; 
    double sumY = 0;
    for(int i = 0; i < pointClusters.size(); i++){
        // compute center of all points in cluster
        for(auto point : pointClusters[i]){
            sumX += point.GetX();
            sumY += point.GetY();
        }
        sumX = sumX / pointClusters[i].size();
        sumY = sumY / pointClusters[i].size();
        centroids[i].SetPosition(sumX, sumY);
    }
    return centroids;
}

std::vector<Node> addNodesToClusterGraph(std::vector<Point> centroids){
    std::vector<Node> clusterNodes(centroids.size());
    for(int i = 0; i < centroids.size(); i++){
        clusterNodes[i]._x = centroids[i].GetX();
        clusterNodes[i]._y = centroids[i].GetY();
        clusterNodes[i]._id = i;
        clusterNodes[i].setCluster(i);
    }
    return clusterNodes;
}
std::pair<std::vector<Edge>, std::vector<int>> addEdgesToClusterGraph(std::vector<vector<int>> edgeWeightsArray){
    //Add all possible edges to the Graph
    std::vector<Edge> edges;
    std::vector<int> edgeWeights;
    
    for (size_t i = 0; i < edgeWeightsArray.size(); ++i) {
        for (size_t j = 0; j + i < edgeWeightsArray.size(); ++j) {
            Edge edge;
            edge.first = i;
            edge.second = j;
            edges.push_back(edge);
            edgeWeights.push_back(edgeWeightsArray[i][j]);
        }
    }
    return std::make_pair(edges, edgeWeights);
}

std::vector<vector<int>> compEdgeWeightToClusterGraph(Graph myGraph, int numOfClusters){
    //The weight represents the number of Connecitons between the clusters. 
    std::vector<vector<int>> edgeWeights(numOfClusters, std::vector<int>(numOfClusters, 0));

    for(int i = 0; i < myGraph.NodeClusters.size(); i++){
        for(auto node : myGraph.NodeClusters[i]){
            for(auto target : myGraph.adjList[node._id]){
                if (node.getCluster() != myGraph.nodes[target].getCluster()){
                    edgeWeights[node.getCluster()][myGraph.nodes[target].getCluster()] += 1;
                }
            }
        }
    }
    for (int i = 0; i < numOfClusters; ++i) {
        for (int j = 0; j < numOfClusters; ++j) {
            std::cout << edgeWeights[i][j] << " ";
        }
        std::cout << std::endl;
    }


    return edgeWeights;
}

Graph createClusterGraph(Graph& myGraph, int numberOfClusters){
    std::cout << "compute virtual sizes\n";
    std::vector<int> virtualClusterSizes = computeVirtualClusterSizes(numberOfClusters, myGraph.nodes.size());
    std::cout << "virtual local clustering \n";

    myGraph.manClustering(virtualClusterSizes, myGraph.width, myGraph.height);

    std::cout << "compute centroids \n";

    myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, numberOfClusters);
    std::vector<Point> centroids = computeCentroids(myGraph.pointClusters);
    std::vector<Node> clusterNodes = addNodesToClusterGraph(centroids);
    std::pair<std::vector<Edge>, std::vector<int>> edgePair = addEdgesToClusterGraph(compEdgeWeightToClusterGraph(myGraph, numberOfClusters));
    std::vector<Edge> clusterEdges = edgePair.first;
    std::vector<int> edgeWeights = edgePair.second;
    Graph ClusterGraph(numberOfClusters, clusterNodes, clusterEdges, centroids, myGraph.width, myGraph.height);
    ClusterGraph.NodeClusters.push_back(ClusterGraph.nodes);

    return ClusterGraph;

}


void matchClusters(Graph& ClusterGraph){
    
    std::vector<size_t> mate;
    ogdf::Graph G;
    ogdf::GraphAttributes GA = ogdf::GraphAttributes(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
    std::vector<Point> centroids;
    BoostGraph BoostGraph = BoostGraph;
        std::cout << "creating cluster Graph \n";
    createClusterGraph(G, GA, ClusterGraph, ClusterGraph.nodes );
    springEmbedding(G, GA); //TODO Add weights of edges
    std::cout << "finished Spring Embedding \n";
    ogdfToNodes(G, GA, ClusterGraph, 0);

    std::cout << "convert cluster to Boost Graph \n";

    BoostGraph = convertClusterToBoostGraph(ClusterGraph.points, ClusterGraph.NodeClusters[0], ClusterGraph.NodeClusters[0].size()*ClusterGraph.points.size());
    mate = minimumWeightMatching(BoostGraph, ClusterGraph.NodeClusters[0].size());

    for(int i = 0; i < mate.size(); i++){
        std::cout << mate[i] << std::endl;
    }
    mapMatching(mate, ClusterGraph, 0);

}
 


std::unordered_map<size_t, size_t> matching(Graph& myGraph, std::vector<vector<Node>> nodeClusters, std::unordered_set<size_t>freePoints, std::unordered_set<size_t>usedPoints, std::vector<Point>points ){

    std::vector<ogdf::Graph> graphClusters(nodeClusters.size());
    std::vector<ogdf::GraphAttributes> graphAttributesClusters(nodeClusters.size());
    for (int i = 0 ; i < graphClusters.size(); i++){
        BoostGraph BoostGraph ;
        std::vector<size_t> mate;

        graphAttributesClusters[i] = ogdf::GraphAttributes(graphClusters[i],ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
        createClusterGraph(graphClusters[i], graphAttributesClusters[i], myGraph, nodeClusters[i]);
                // for (auto p : myGraph.points){
                //     std::cout << "ID" << p.GetId() << " C" << p.GetCluster() << std::endl; 
                // }
        springEmbedding(graphClusters[i], graphAttributesClusters[i]);
                // for (auto p : myGraph.points){
                //     std::cout << "ID" << p.GetId() << " C" << p.GetCluster() << std::endl; 
                // }
        // crossMin(graphClusters[i], graphAttributesClusters[i], myGraph.nodes);
        ogdfToNodes(graphClusters[i], graphAttributesClusters[i], myGraph, i);
        BoostGraph = convertClusterToBoostGraph(myGraph.pointClusters[i], myGraph.NodeClusters[i], myGraph.nodes.size() * myGraph.edges.size() * myGraph.points.size());
        mate = minimumWeightMatching(BoostGraph, myGraph.NodeClusters[i].size());
                // for (auto l : myGraph.pointClusters){
                //     for (auto p : l){ 
                //     std::cout << "ID" << p.GetId() << " C" << p.GetCluster() << std::endl; 
                //     }
                // }
        mapMatching(mate, myGraph, i);
        
    }
    // Print Results    
    // for (int i = 0; i < myGraph.mapVerticesToPoints.size(); i++){
    //     std::cout << "NodeId: " << i <<" PointId: " << myGraph.mapVerticesToPoints[i] << std::endl;
    // }
    return myGraph.mapVerticesToPoints;
}
