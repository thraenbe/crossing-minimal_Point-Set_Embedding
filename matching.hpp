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

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>

#include <algorithm> // For std::find
#include <limits>



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




void createClusterGraph(ogdf::Graph& G,Graph& myGraph, std::vector<Node> nodes ){
    std::vector<size_t> tmp;
    std::unordered_map<int, ogdf::node> nodeMap;
    for (auto node : nodes){
        ogdf::node newNode = G.newNode(node.getId());
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

void springEmbedding(ogdf::Graph& G){

}

void crossingMinimization(ogdf::Graph& G){
    // Prepare the PlanRep and other required parameters
    ogdf::PlanRep pr(G);
    int cc = 0;
    int crossingNumber;
    ogdf::EdgeArray<int> *pCostOrig = nullptr;
    ogdf::EdgeArray<bool> *pForbiddenOrig = nullptr;
    ogdf::EdgeArray<uint32_t> *pEdgeSubGraphs = nullptr;

    Minimizer minimizer;

    ogdf::GraphAttributes TMP(G,  ogdf::GraphAttributes::nodeGraphics);

    ogdf::PlanarizerMixedInsertion mixedInsertion;
    mixedInsertion.call(pr, cc, crossingNumber, pCostOrig, pForbiddenOrig, pEdgeSubGraphs);
    // minimizer.callDoCall(pr, cc, pCostOrig, pForbiddenOrig, pEdgeSubGraphs, crossingNumber);
    
    ogdf::GraphAttributes GA(pr,  ogdf::GraphAttributes::nodeGraphics);
    // minimizer.getNodePositions(GA);

    ogdf::PlanarizationLayout pl;
    // std::cout << "SEGME" << std::endl ;
    // Retrieve and store the node positions in graph attributes

    // Print positions of the nodes
    for (auto v : pr.nodes) {
        // std::cout << "SEGME" << std::endl ;
        // std::cout << "Node " << v->index() << ": ("
        //           << GA.x(v) << ", " << GA.y(v) << ")" << std::endl;
    }
}

int getWeight(Point point, Node node, int num){
    int weight;
    weight =  (int) num * (1 / (euclideanDistance(point, node.nodeToPoint()) + 1));
    return weight;
}

std::vector<Point> allignBoundingBoxes(std::vector<Point> points, std::vector<Node> nodes) {
    if (points.empty() || nodes.empty()) return points;

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
    for (auto& point : points) {
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
    }

    return points;
}



BoostGraph convertClusterToBoostGraph(std::vector<Point> points, std::vector<Node> nodes, int num){
    // Add a node for each point
    std::cout << points.size() << " " << nodes.size() << std::endl;
    points = allignBoundingBoxes(points, nodes);
    std::vector<VertexDescriptor> tmp;
    BoostGraph BoostGraph;
    for (auto point : points) {
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        tmp.push_back(v);
        // std::cout << "Added vertex " << v << " for Point(" << point._x << ", " << point._y << ")\n";
    }
    for (int j = points.size(); j < points.size() + nodes.size() ; j++){
        VertexDescriptor v = boost::add_vertex(BoostGraph);
        for (int i = 0; i < tmp.size(); i++){
            int weight = getWeight(points[i], nodes[j-points.size()], num);
            boost::add_edge(j , i , EdgeProperty(weight), BoostGraph);
            // std::cout << "Match: " << j << " " << i << "   Weight: " << weight << "  Node: " << nodes[j-points.size()].getX() << " " << nodes[j-points.size()].getY() << " " << "Point:  "<< points[i].GetX() << " " << points[i].GetY() << std::endl;
        }
    }
    return BoostGraph;
}

std::vector<size_t> minimumgWeightMatching(BoostGraph& BoostGraph){
    std::vector<VertexDescriptor> mate(boost::num_vertices(BoostGraph));
    // Find the maximum cardinality matching
    std::cout << "Starting max weight matching..." << std::endl;
    try {
        boost::maximum_weighted_matching(BoostGraph, &mate[0]);
        std::cout << "Max weight matching completed." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception during matching: " << e.what() << std::endl;
    }

    std::cout << "Found a weighted matching:" << std::endl;
    std::cout << "Matching size is " << boost::matching_size(BoostGraph, &mate[0])
              << ", total weight is " << boost::matching_weight_sum(BoostGraph, &mate[0])
              << std::endl;
    std::cout << std::endl;

    // Print the matching result
    for (std::size_t i = 0; i < mate.size(); ++i) {
        if (mate[i] != BoostGraph::null_vertex() && i < mate[i]) {
            std::cout << "Matched: " << i << " - " << mate[i] << std::endl;
        }
    }



}
 
// Uncomment following

std::unordered_map<size_t, size_t> matching(Graph& myGraph, std::vector<vector<Node>> nodeClusters, std::unordered_set<size_t>freePoints, std::unordered_set<size_t>usedPoints, std::vector<Point>points ){

    std::vector<ogdf::Graph> graphClusters(nodeClusters.size());
    for (int i = 0 ; i < graphClusters.size(); i++){
        std::cout << "Graph: " << i << std::endl;
        std::cout << "A" << std::endl;
        createClusterGraph(graphClusters[i], myGraph, nodeClusters[i]);
        std::cout << "B" << std::endl;
        crossingMinimization(graphClusters[i] );
        BoostGraph BoostGraph ;
        std::cout << myGraph.pointClusters[i].size() << " " << myGraph.NodeClusters[i].size() << std::endl;
        BoostGraph = convertClusterToBoostGraph(myGraph.pointClusters[i],myGraph.NodeClusters[i], myGraph.nodes.size() + myGraph.edges.size() + myGraph.points.size());
        minimumgWeightMatching(BoostGraph);
        
    }

    for (auto graph : graphClusters){
        if (graph.edges.size() > 0) {
            std::cout << graph.edges.size() << ": " ;
            std::cout << graph.firstEdge()->target()->index() << " - " << graph.firstEdge()->source()->index() << std::endl;

        }
    }
    return myGraph.mapVerticesToPoints;
    // for (auto cluster : clusters){
    //     ogdf::Graph optimizedClusterGraph = crossingMinimization(cluster, points, edges);
    //     BoostGraph boostGraph = convertClusterToBoostGraph(cluster, points);
    //     allignBoundingBoxes(boostGraph, points);
    //     minimumgWeightMatching(boostGraph);


    // }
}

// Function to convert OGDF SimpleCluster to Boost Graph
BoostGraph convertOGDFToBoostGraph(const ogdf::ClusterGraph &cg) {
    BoostGraph bg;

    // Mapping OGDF nodes to Boost graph vertices
    std::unordered_map<ogdf::node, VertexDescriptor> node_map;

    // // Add all nodes to the Boost graph
    // for (auto v : cg.nodes) {
    //     VertexDescriptor vd = boost::add_vertex(bg);
    //     node_map[v] = vd;
    // }

    // // Add all edges to the Boost graph
    // for (auto e : cg.edges) {
    //     auto source = e->source();
    //     auto target = e->target();
    //     boost::add_edge(node_map[source], node_map[target], bg);
    // }

    return bg;
}
