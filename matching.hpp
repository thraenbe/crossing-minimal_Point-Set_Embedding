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
#include <ogdf/basic/EdgeArray.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/graph/graphviz.hpp>

#include <algorithm> // For std::find
#include <limits>
#include <thread>
#include <future>
#include <chrono>



// Define the Boost graph type

using Weight = // boost::multiprecision::cpp_dec_float_50;
    boost::multiprecision::number<
        boost::multiprecision::cpp_dec_float<50>,
        boost::multiprecision::et_off >;

typedef boost::property< boost::edge_weight_t, Weight, boost::property< boost::edge_index_t, int > >EdgeProperty;
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

void createClusterGraph(ogdf::Graph& G, ogdf::GraphAttributes& GA, Graph& myGraph, std::vector<int> edgeWeights){
    std::vector<Node> nodes = myGraph.nodes;
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
    for (auto edge : myGraph.edges){
        G.newEdge(nodeMap[edge.first], nodeMap[edge.second]);
        
    }
}

void springEmbedding(ogdf::Graph& G, ogdf::GraphAttributes& GA, std::vector<int> edgeWeights){
    // Instantiate the Kamada-Kawai spring embedder
    ogdf::SpringEmbedderKK KK;
    ogdf::EdgeArray<double>eLength(G);
    int idx =0 ;
    for(auto e : eLength){
        e = 100 * (1/ (1 + edgeWeights[idx] ));
        idx++;
    }
    

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
        KK.call(GA, eLength);
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

void springEmbedding(ogdf::Graph& G, ogdf::GraphAttributes& GA){
    // Instantiate the Kamada-Kawai spring embedder
    ogdf::SpringEmbedderKK KK;
    

    // Set parameters (optional)
    // KK.setMaxGlobalIterations(10);
    // KK.setStopTolerance(1e-4);
    // KK.setDesLength(30.0); // Set desired edge length

    // std::cout << "Positions before Spring Embedder:" << std::endl;
    // for (ogdf::node n : G.nodes) {
    //     std::cout << "Node " << GA.label(n) << " position: (" << GA.x(n) << ", " << GA.y(n) << ");    " ;
    // }
    std::cout << std::endl;
    // Apply the layout algorithm
    if (ogdf::isConnected(G)){        
        KK.call(GA);
    }

    std::cout << std::endl;

    // Print node positions
    std::cout << "Spring Embedder completed" << std::endl;
    // for (ogdf::node n : G.nodes) {
    //     std::cout << "Node " << GA.label(n) << " position: (" << GA.x(n) << ", " << GA.y(n) << ");    " ;
    // }
    // std::cout << std::endl;
    // std::cout << std::endl;
    
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


double getWeight(Point point, Node node, int num){
    double weight;
    double distance = euclideanDistance(point, node.nodeToPoint());
    std::cout << distance << std::endl;
    assert(distance >= 0);
    weight = (double) num * (1 / (distance + 1));
    if (weight == 0){
        weight += 1;
    }

    assert(weight > 0);
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
    std::cout << points.size() << " " << nodes.size() << " \n\n\n\n" ;
    for (auto point : points) {
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        tmp.push_back(v);
        // std::cout << "Added vertex " << v << " for Point(" << point._x << ", " << point._y << ")\n";
    }
    for (int j = points.size(); j < allignedNodes.size() + points.size() ; j++){
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        for (int i = 0; i < tmp.size(); i++){
            double weight = getWeight(points[i], allignedNodes[j-points.size()], num);
            
            assert(weight > 0);
            boost::add_edge(j , i , EdgeProperty(weight), BoostGraph);
            std::cout << "Match: " << j << " " << i <<"  NUM: " << num << "   Weight: " << weight << "  Node: " << nodes[j-points.size()].getX() << " " << nodes[j-points.size()].getY() << " " << "Point:  "<< allignedNodes[i].getX() << " " << allignedNodes[i].getY() << std::endl;
        }
    }
    return BoostGraph;
}

std::vector<size_t> minimumWeightMatching(BoostGraph& BoostGraph, int clusterSize){
    std::vector<VertexDescriptor> mate(boost::num_vertices(BoostGraph));
    try {
        boost::maximum_weighted_matching(BoostGraph, &mate[0]);
        
        std::cout << "Max weight matching completed. Size: " << boost::matching_size(BoostGraph, &mate[0]) << "   Weight: "<< boost::matching_weight_sum(BoostGraph, &mate[0])  << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception during matching: " << e.what() << std::endl;

    }

    return mate;
}

void mapMatching(std::vector<size_t> mate, Graph& myGraph, int cluster){
    int i = 0;
    for (auto m : mate){
        std::cout << i++ << " " << m << std::endl; 
    }
    for (int i = 0; i < myGraph.NodeClusters[cluster].size(); i++){
        int nodeId = myGraph.NodeClusters[cluster][i].getId();
        assert(i+myGraph.pointClusters[cluster].size()< mate.size());
        int pointId = myGraph.pointClusters[cluster][mate[i+myGraph.pointClusters[cluster].size()]].GetId() ;
        std::cout << pointId << std::endl;
        assert(pointId < myGraph.points.size());
        assert(nodeId < myGraph.nodes.size());
        myGraph.mapVerticesToPoints[nodeId] = pointId; 
    }
}

void mapMatching(std::vector<size_t> mate, Graph& myGraph){
    std::cout << "Starting to match vertices" << std::endl;
    for (int i = 0; i < myGraph.nodes.size(); i++){
        int nodeId = myGraph.nodes[i].getId();
        assert(i+myGraph.nodes.size() < mate.size());
        int pointId = mate[i+myGraph.nodes.size()] ;
        std::cout << pointId << " " << myGraph.points.size() << std::endl;
        assert(pointId < myGraph.points.size());
        assert(nodeId < myGraph.nodes.size());
        myGraph.mapVerticesToPoints[nodeId] = pointId; 
        // std::cout << "Matched Node: " << nodeId << " with Point: " << pointId << std::endl; 
    }
}

void ogdfToNodes(ogdf::Graph& G, ogdf::GraphAttributes& GA, Graph& myGraph, int cluster){
    int i = 0;
    std::cout << "NodeCluster: " << cluster << " of " << " " << myGraph.NodeClusters.size() << std::endl;
    for (ogdf::node v : G.nodes){
        myGraph.NodeClusters[cluster][i]._x = GA.x(v);
        myGraph.NodeClusters[cluster][i]._y = GA.y(v);
        if (myGraph.NodeClusters[cluster][i]._id != std::stoi(GA.label(v))){
            std::cout << "IDs Do Not match \n";
        }
        i++;
    }
}

std::vector<int> computeVirtualClusterSizes(int numberOfClusters, int numberOfPoints){
    assert(numberOfClusters < numberOfPoints);

    int rest = numberOfPoints % numberOfClusters;
    int q = (numberOfPoints - rest) / numberOfClusters;

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
        for (size_t j = i + 1; j < edgeWeightsArray.size(); ++j) {
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

std::pair<Graph, std::vector<int>> createClusteringGraph(Graph& myGraph, int numberOfClusters){
    std::cout << "compute virtual sizes\n";
    std::vector<int> virtualClusterSizes = computeVirtualClusterSizes(numberOfClusters, myGraph.points.size());
    std::cout << "virtual local clustering \n";

    myGraph.manClustering(virtualClusterSizes, myGraph.width, myGraph.height);

    myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, numberOfClusters);
    std::vector<Point> centroids = computeCentroids(myGraph.pointClusters);
    std::vector<Node> clusterNodes = addNodesToClusterGraph(centroids);
    std::pair<std::vector<Edge>, std::vector<int>> edgePair = addEdgesToClusterGraph(compEdgeWeightToClusterGraph(myGraph, numberOfClusters));
    std::vector<Edge> clusterEdges = edgePair.first;
    std::vector<int> edgeWeights = edgePair.second;
    Graph ClusterGraph(numberOfClusters, clusterNodes, clusterEdges, centroids, myGraph.width, myGraph.height);
    ClusterGraph.NodeClusters.push_back(ClusterGraph.nodes);
    std::cout << "Initial Cluster Layout Graph:      Edges: " << clusterEdges.size() << "     edgeWeights: " << edgeWeights.size() << "   Nodes: " << centroids.size() << std::endl;

    return std::make_pair(ClusterGraph, edgeWeights);
}


void matchClusters(Graph& ClusterGraph, std::vector<int> edgeWeights){
    
    std::vector<size_t> mate;
    ogdf::Graph G;
    ogdf::GraphAttributes GA = ogdf::GraphAttributes(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
    
    BoostGraph BoostGraph;
        std::cout << "creating cluster Graph \n";
    createClusterGraph(G, GA, ClusterGraph, edgeWeights );
    std::cout << " Number of edges: " << edgeWeights.size() << std::endl;


    springEmbedding(G, GA, edgeWeights); //TODO Add weights of edges
    std::cout << "finished Spring Embedding \n";
    ogdfToNodes(G, GA, ClusterGraph, 0);

    std::cout << "convert cluster to Boost Graph \n";

    BoostGraph = convertClusterToBoostGraph(ClusterGraph.points, ClusterGraph.NodeClusters[0], ClusterGraph.NodeClusters[0].size()*ClusterGraph.points.size());
    mate = minimumWeightMatching(BoostGraph, ClusterGraph.NodeClusters[0].size());

    // for(int i = 0; i < mate.size(); i++){
    //     std::cout << mate[i] << std::endl;
    // }
    mapMatching(mate, ClusterGraph);

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
        std::cout << "Matching successful! " << std::endl;
        
    }
    // Print Results    
    // for (int i = 0; i < myGraph.mapVerticesToPoints.size(); i++){
    //     std::cout << "NodeId: " << i <<" PointId: " << myGraph.mapVerticesToPoints[i] << std::endl;
    // }
    return myGraph.mapVerticesToPoints;
}
