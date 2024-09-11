#pragma once
#include "Graph.hpp"

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/geometry.h>

#include <ogdf/cluster/ClusterPlanarity.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/ClustererModule.h>
#include <ogdf/planarity/PlanarizationLayout.h>
#include <ogdf/planarity/PlanarSubgraphModule.h>
#include <ogdf/planarity/PlanarizerMixedInsertion.h>
#include <ogdf/energybased/SpringEmbedderKK.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/graph/graphviz.hpp>

#include <algorithm> // For std::find
#include <iostream>
#include <limits>
#include <thread>
#include <future>
#include <chrono>
#include <ranges>



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




void createClusterGraph(ogdf::Graph& G, ogdf::GraphAttributes& GA, const Graph& myGraph, const std::vector<Node>& nodes ){
    std::vector<size_t> tmp;
    std::unordered_map<size_t, ogdf::node> nodeMap;
    for (const auto& node : nodes){
        ogdf::node newNode = G.newNode();
        GA.x(newNode) = node.GetX();
        GA.y(newNode) = node.GetY();
        GA.label(newNode) = std::to_string(node.GetId());
        tmp.push_back(node.GetId());
        nodeMap[node.GetId()] = newNode;
        
    }
    for (const auto& node : nodes) {
            for (const auto target : myGraph.adjList[node.GetId()]){
            if (std::ranges::find(tmp,target)!= tmp.end()){
                if (node.GetId() < target){ // Every Edge is only appended once istead of twice
                    G.newEdge(nodeMap[node.GetId()], nodeMap[target]);
                }
            }
        }
    }
}

void createClusterGraph(ogdf::Graph& G, ogdf::GraphAttributes& GA, const Graph& myGraph, const std::vector<size_t>& edgeWeights){
    std::vector<size_t> tmp;
    std::unordered_map<size_t, ogdf::node> nodeMap;
    for (const auto& node : myGraph.nodes){
        ogdf::node newNode = G.newNode();
        GA.x(newNode) = node.GetX();
        GA.y(newNode) = node.GetY();
        GA.label(newNode) = std::to_string(node.GetId());
        tmp.push_back(node.GetId());
        nodeMap[node.GetId()] = newNode;
    }
    for (const auto& edge : myGraph.edges){
        G.newEdge(nodeMap[edge.first], nodeMap[edge.second]);
        
    }
}

void springEmbedding(const ogdf::Graph& G, ogdf::GraphAttributes& GA, const std::vector<size_t>& edgeWeights){

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

void springEmbedding(const ogdf::Graph& G, ogdf::GraphAttributes& GA){
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


inline float getWeight(const Point& point, const Node& node, const size_t num){
    float weight =  num * (1 / (L2Dist(point, node.NodeToPoint()) + 1));
    if (weight == 0){
        weight += 1;
    }
    assert(weight != 0);
    return weight;
}

namespace
{
    template <typename T>
    std::array<double,4>GetMinMaxXY(const std::vector<T>& nodes)
    {
        // TODO:: Assert that input is really in positive halfplane
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        double maxX = -1;
        double maxY = -1;
        for (const auto& node : nodes){
            const double x = node.GetX();
            const double y = node.GetY();
            minX = std::min(minX,x);
            maxX = std::max(maxX,x);
            minY = std::min(minY,y);
            maxY = std::max(maxY,y);
        }
        return {minX, minY, maxX, maxY};
    }
}
    std::vector<Node> allignBoundingBoxes(std::vector<Point>& points,std::vector<Node>& nodes) {
        if (points.empty() || nodes.empty())
        {
            return nodes;
        }
        // Get min x and y of nodes
        const auto minMaxPoints = GetMinMaxXY(points);
        const auto minMaxNodes = GetMinMaxXY(nodes);
        assert(minMaxPoints.size() == 4);
        assert (minMaxNodes.size() == 4);
        const auto& [minXPoints,minYPoints,maxXPoints,maxYPoints] = minMaxPoints;
        const auto& [minXNodes,minYNodes,maxXNodes,maxYNodes] = minMaxNodes; // Changed Points To Nodes
        assert(minXPoints <= maxXPoints);
        assert(minYPoints<= maxYPoints);
        assert(minXNodes <= maxXNodes);
        assert(minYNodes <= maxYNodes);

        // Align and stretch points
        const double widthNodes = maxXNodes - minXNodes;
        const double heightNodes = maxYNodes - minYNodes;
        const double widthPoints = maxXPoints - minXPoints;
        const double heightPoints = maxYPoints - minYPoints;
        const bool stretchWidth{widthPoints > 0};
        const bool stretchHeight{heightPoints > 0};
        for (auto& node : nodes) {
            // Align points to minimal x and y of nodes
            node._x = node._x - minXPoints + minXNodes;
            node._y = node._y - minYPoints + minYNodes;
            if (stretchWidth) {
                node._x = minXNodes + (node._x - minXNodes) * (widthNodes / widthPoints);
            }

            // Stretch points vertically
            if (stretchHeight) {
                node._y = minYNodes + (node._y - minYNodes) * (heightNodes / heightPoints);
            }
            node._x = (int) node._x;
            node._y = (int) node._y;
        }

        return nodes;
    }


    BoostGraph convertClusterToBoostGraph(std::vector<Point>& points,std::vector<Node>& nodes, const size_t num){

        // std::cout << "Point " << points[0].GetX() << points[0].GetY() << "Point " << points[1].GetX() << points[1].GetY() << std::endl ;
        auto allignedNodes = allignBoundingBoxes(points, nodes);
        const auto minMaxNodes = GetMinMaxXY(allignedNodes);
        const auto& [minXNodes,minYNodes,maxXNodes,maxYNodes] = minMaxNodes;
        
        // std::cout << "Point " << allignedPoints[0].GetX() << allignedPoints[0].GetY() << "Point " << allignedPoints[1].GetX() << allignedPoints[1].GetY() << std::endl ;
        std::vector<VertexDescriptor> tmp;
        BoostGraph BoostGraph;
        std::cout  << points.size() << " " << allignedNodes.size() << std::endl;
        assert(points.size() >=  allignedNodes.size());
        for (const auto& _ : points) {
            VertexDescriptor v = boost::add_vertex(BoostGraph);
            tmp.push_back(v);
            // std::cout << "Added vertex " << v << " for Point(" << point._x << ", " << point._y << ")\n";
        }
        for (auto j = points.size(); j < allignedNodes.size() + points.size() ; j++){
            VertexDescriptor v = boost::add_vertex(BoostGraph);
            for (auto i = 0; i < tmp.size(); i++){
                float weight = getWeight(points[i], allignedNodes[j-points.size()], num);
                if (weight == 0){
                    weight += 1;
                }

                assert(weight != 0);
                boost::add_edge(j , i , EdgeProperty(static_cast<int>(weight)), BoostGraph);
                // std::cout << "Match: " << j << " " << i <<"  NUM: " << num << "   Weight: " << weight << "  Node: " << nodes[j-points.size()].GetX() << " " << nodes[j-points.size()].GetY() << " " << "Point:  "<< allignedNodes[i].GetX() << " " << allignedNodes[i].GetY() << std::endl;
            }

        }
        assert(boost::num_vertices(BoostGraph) == points.size() + nodes.size());
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

    void simpleAssign(Graph& myGraph, int size){
        std::cout << " ------------------------------------------ Simple Assignment for matching \n \n";
        for (int i = 0; i < size ; i++){
            myGraph.mapVerticesToPoints[i] = i; 
        }
    }

    void mapMatching(std::vector<size_t> mate, Graph& myGraph, int cluster){
        std::cout << std::endl;
        for (int i = 0; i < myGraph.NodeClusters[cluster].size(); i++){
            int nodeId = myGraph.NodeClusters[cluster][i].GetId();
            assert(i+myGraph.pointClusters[cluster].size()< mate.size());
            int pointId = myGraph.pointClusters[cluster][mate[i+myGraph.pointClusters[cluster].size()]].GetId() ;
            if (pointId >= myGraph.points.size()){
                simpleAssign(myGraph,myGraph.NodeClusters[cluster].size() );
                break;
            }
            assert(pointId < myGraph.points.size());
            assert(nodeId < myGraph.nodes.size());
            myGraph.mapVerticesToPoints[nodeId] = pointId; 
        }
    }



    void mapMatching(std::vector<size_t> mate, Graph& myGraph){
        std::cout << "Starting to match vertices of 1. Clustering" << std::endl;
        for (int i = 0; i < myGraph.nodes.size(); i++){
            int nodeId = myGraph.nodes[i].GetId();
            assert(i+myGraph.nodes.size() < mate.size());
            int pointId = mate[i+myGraph.nodes.size()];
            std::cout << pointId << " " << myGraph.points.size() << std::endl;
            if (pointId >= myGraph.points.size()){
                simpleAssign(myGraph,myGraph.nodes.size() );
                break;
            }
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
            assert(myGraph.NodeClusters[cluster][i]._id == std::stoi(GA.label(v)));
            i++;
        }
    }

    void ogdfToNodes(ogdf::Graph& G, ogdf::GraphAttributes& GA, Graph& myGraph){
        int i = 0;
        std::cout << "Clustering Graph " << myGraph.NodeClusters.size() << std::endl;
        for (ogdf::node v : G.nodes){
            myGraph.nodes[i]._x = GA.x(v);
            myGraph.nodes[i]._y = GA.y(v);
            if (myGraph.nodes[i]._id != std::stoi(GA.label(v))){
                std::cout << "IDs Do Not match \n";
            }
            assert(myGraph.nodes[i]._id == std::stoi(GA.label(v)));
            i++;
        }
    }

    std::vector<size_t> computeVirtualClusterSizes(const size_t numberOfClusters,const size_t numberOfPoints){
        assert(numberOfClusters < numberOfPoints);

        auto rest = numberOfPoints % numberOfClusters;
        auto q = (numberOfPoints - rest) / numberOfClusters;

        std::vector<size_t> virtualClusterSizes(numberOfClusters);
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

    std::vector<Point> computeCentroids(const std::vector<vector<Point>>& pointClusters){
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

    std::vector<Node> addNodesToClusterGraph(const std::vector<Point>& centroids){
        std::vector<Node> clusterNodes(centroids.size());
        for(int i = 0; i < centroids.size(); i++){
            clusterNodes[i]._x = centroids[i].GetX();
            clusterNodes[i]._y = centroids[i].GetY();
            clusterNodes[i]._id = i;
            clusterNodes[i].SetCluster(i);
        }
        return clusterNodes;
    }

    std::pair<std::vector<Edge>, std::vector<size_t>> addEdgesToClusterGraph(const std::vector<vector<size_t>>& edgeWeightsArray){
        //Add all possible edges to the Graph
        std::vector<Edge> edges;
        std::vector<size_t> edgeWeights;

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

    std::vector<vector<size_t>> compEdgeWeightToClusterGraph(Graph& myGraph, const size_t numOfClusters){
        //The weight represents the number of Connecitons between the clusters.
        std::vector<vector<size_t>> edgeWeights(numOfClusters, std::vector<size_t>(numOfClusters, 0));

        for(auto i = 0; i < myGraph.NodeClusters.size(); ++i){
            for(const auto& node : myGraph.NodeClusters[i]){
                for(auto target : myGraph.adjList[node._id]){
                    if (node.GetCluster() != myGraph.nodes[target].GetCluster()){
                        edgeWeights[node.GetCluster()][myGraph.nodes[target].GetCluster()] += 1;
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

    std::pair<Graph, std::vector<size_t>> createClusteringGraph(Graph& myGraph, size_t numberOfClusters){
        std::cout << "compute virtual sizes\n";
        std::vector<size_t> virtualClusterSizes = computeVirtualClusterSizes(numberOfClusters, myGraph.points.size());
        std::cout << "virtual local clustering \n";

        myGraph.manClustering(virtualClusterSizes, myGraph.width, myGraph.height);

        std::cout << "compute centroids \n";

        myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, numberOfClusters);
        std::vector<Point> centroids = computeCentroids(myGraph.pointClusters);
        std::vector<Node> clusterNodes = addNodesToClusterGraph(centroids);
        std::pair<std::vector<Edge>, std::vector<size_t>> edgePair = addEdgesToClusterGraph(compEdgeWeightToClusterGraph(myGraph, numberOfClusters));
        std::vector<Edge> clusterEdges = edgePair.first;
        std::vector<size_t> edgeWeights = edgePair.second;
        Graph ClusterGraph(numberOfClusters, clusterNodes, clusterEdges, centroids, myGraph.width, myGraph.height);
        ClusterGraph.NodeClusters.push_back(ClusterGraph.nodes);
        std::cout << "Initial Cluster Layout Graph:      Edges: " << clusterEdges.size() << "     edgeWeights: " << edgeWeights.size() << "   Nodes: " << centroids.size() << std::endl;

        return std::make_pair(ClusterGraph, edgeWeights);
    }

    void matchClusters(Graph& ClusterGraph, const std::vector<size_t>& edgeWeights){

        std::vector<size_t> mate;
        ogdf::Graph G;
        ogdf::GraphAttributes GA = ogdf::GraphAttributes(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);

        BoostGraph BoostGraph;
        std::cout << "creating cluster Graph \n";
        createClusterGraph(G, GA, ClusterGraph, edgeWeights );
        std::cout << " Number of edges: " << edgeWeights.size() << std::endl;


        springEmbedding(G, GA, edgeWeights); //TODO Add weights of edges
        std::cout << "finished Spring Embedding \n";
        ogdfToNodes(G, GA, ClusterGraph);

        std::cout << "convert cluster to Boost Graph \n";

        BoostGraph = convertClusterToBoostGraph(ClusterGraph.points, ClusterGraph.NodeClusters[0], ClusterGraph.NodeClusters[0].size()*ClusterGraph.points.size());
        mate = minimumWeightMatching(BoostGraph, ClusterGraph.NodeClusters[0].size());

        mapMatching(mate, ClusterGraph);
    }



    std::unordered_map<size_t, size_t> matching(Graph& myGraph, const std::vector<vector<Node>>& nodeClusters, std::unordered_set<size_t>&freePoints, std::unordered_set<size_t>&usedPoints, std::vector<Point>&points ){
        assert(myGraph.pointClusters.size() == myGraph.NodeClusters.size());
        std::vector<ogdf::Graph> graphClusters(nodeClusters.size());
        std::vector<ogdf::GraphAttributes> graphAttributesClusters(nodeClusters.size());
        for (int i = 0 ; i < graphClusters.size(); i++){
            BoostGraph BoostGraph ;
            std::vector<size_t> mate;


            graphAttributesClusters[i] = ogdf::GraphAttributes(graphClusters[i],ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
            createClusterGraph(graphClusters[i], graphAttributesClusters[i], myGraph, nodeClusters[i]);
            springEmbedding(graphClusters[i], graphAttributesClusters[i]);
            ogdfToNodes(graphClusters[i], graphAttributesClusters[i], myGraph, i);
            std::cout << "print Sizes" << myGraph.points.size() << " " << myGraph.nodes.size() << " " << myGraph.NodeClusters.size() <<  " \n" ;

            for(int idx = 0; idx < myGraph.NodeClusters.size(); idx++){
                std::cout << myGraph.pointClusters[idx].size() << " " << myGraph.NodeClusters[idx].size() << std::endl;
            }
            assert(myGraph.pointClusters[i].size() >= myGraph.NodeClusters[i].size());
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
        for (int i = 0; i < myGraph.mapVerticesToPoints.size(); i++){
            std::cout << "(" << i <<", " << myGraph.mapVerticesToPoints[i] << ");   " ;
        }
        std::cout << std::endl;

        return myGraph.mapVerticesToPoints;
    }