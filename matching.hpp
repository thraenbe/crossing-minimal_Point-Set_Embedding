#pragma once
#include <ogdf/basic/Graph.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/cluster/ClusterPlanarity.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/ClustererModule.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include "Graph.hpp"
#include <algorithm> // For std::find




// Define the Boost graph type
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> BoostGraph;
typedef boost::graph_traits<BoostGraph>::vertex_descriptor VertexDescriptor;



void createClusterGraph(ogdf::Graph& G,Graph myGraph, std::vector<Node> nodes ){
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

// ogdf::Graph crossingMinimization(cluster, points, edges){
//     ogdf::Graph optimizedClusterGraph;

// }


BoostGraph convertClusterToBoostGraph(){

}
void allignBoundingBoxes(){

}
void minimumgWeightMatching(){

}
 
// Uncomment following

// std::unordered_map<size_t, size_t> matching(Graph myGraph, std::vector<vector<Node>> nodeClusters, std::unordered_set<size_t>freePoints, std::unordered_set<size_t>usedPoints, std::vector<Point>points ){
    
    
//     for (auto cluster : clusters){
//         ogdf::Graph optimizedClusterGraph = crossingMinimization(cluster, points, edges);
//         BoostGraph boostGraph = convertClusterToBoostGraph(cluster, points);
//         allignBoundingBoxes(boostGraph, points);
//         minimumgWeightMatching(boostGraph);


//     }
// }

// // Function to convert OGDF SimpleCluster to Boost Graph
// BoostGraph convertOGDFToBoostGraph(const ogdf::ClusterGraph &cg) {
//     BoostGraph bg;

//     // Mapping OGDF nodes to Boost graph vertices
//     std::unordered_map<ogdf::node, VertexDescriptor> node_map;

//     // Add all nodes to the Boost graph
//     for (auto v : cg.nodes) {
//         VertexDescriptor vd = boost::add_vertex(bg);
//         node_map[v] = vd;
//     }

//     // Add all edges to the Boost graph
//     for (auto e : cg.edges) {
//         auto source = e->source();
//         auto target = e->target();
//         boost::add_edge(node_map[source], node_map[target], bg);
//     }

//     return bg;
// }
