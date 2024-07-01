#include"JSONParser.hpp"
#include <random>
#include <iostream>
#include <cassert>
#include <iterator>
#include <limits>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/cluster/ClusterGraph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/SList.h>
#include <ogdf/planarity/PlanarizationLayout.h>


 

using size_t = std::size_t;

Graph read(size_t numNodes, std::vector<std::pair<size_t, size_t>>  edges, std::vector<Point> positions, int width, int height, JSONParser jp, std::string path) {
    // Width and height are kind of useless but have to be written back.
    jp.jsonToGraph(numNodes, edges, positions, width, height, path);
    const int n = (int) numNodes;
    const int m = edges.size();
    const int f = positions.size();

	std::cout << "#Nodes : " << n << "   #Edges: " << m << "    #Positions: " << f << std::endl ;
    assert(f >= n);

    Graph G(n, edges, positions);


    return G;
}

void write(Graph G, std::string path, std::string name, int width, int height, JSONParser jp ){
    auto createOutputFileName = [&](){
    	return name + "_"+std::to_string(G.ComputeCrossings());
    };
    auto WriteGraph = [&]() {
        jp.graphToJson(G, path, createOutputFileName(), width, height); };

    size_t crossNum = G.ComputeCrossings();
    std::cout << "Starting crossings" << std::endl;
    std::cout << crossNum << std::endl;
    WriteGraph();
}

void printVector(const std::vector<std::vector<size_t>>& vec) {
    for (const auto& innerVec : vec) {
        for (size_t elem : innerVec) {
            std::cout << elem << " " << std::endl;
        }
    }
}



int main(int argc, char* argv[]) {


    std::vector<Point> positions;
    std::vector<std::pair<size_t, size_t>> edges;
    size_t numNodes{0};

    int width;
    int height;

    JSONParser jp;

    const std::string filepath = "../graphs/";
    const std::string outPath = "../results";
    const std::string filename = "g1";

    Graph myGraph = read(numNodes, edges, positions, width, height, jp, filepath + filename);
    std::cout << " Graph was read succesfully " << std::endl ;

    int k = 3;
    // std::vector<Point>points = myGraph.points;
    myGraph.kMeans(myGraph.points, k);



    ogdf::Graph G; 
    ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);



    myGraph.toOgdfGraph(G);
    

    ogdf::SList<ogdf::SimpleCluster *> clusters; 
    ogdf::ClusterGraph ClusterGraph;

    ogdf::Clusterer CG(G);

    double avgCIdx = CG.averageCIndex();
    for (ogdf::node v : G.nodes) {
        double CIdx = CG.computeCIndex(G, v); 
    }
    std::cout << avgCIdx << std::endl;

    CG.setRecursive(false);
    CG.setStopIndex(0.2);

    

    CG.setAutomaticThresholds(8);


    CG.computeClustering(clusters);
    CG.createClusterGraph(ClusterGraph); // not really working 


    int numberOfClusters = ClusterGraph.numberOfClusters();


    std::cout << " Number of Clusters: " << clusters.size() << std::endl ;



    // Print the clusters
    int clusterIndex = 0;
    for (auto cluster : clusters) {
        std::cout << "Cluster " << clusterIndex++ << ": ";
        for (ogdf::node v : cluster->nodes()) {
            std::cout << v->index() << " ";
        }
        std::cout << std::endl;
    }

    // Apply a layout algorithm
    ogdf::PlanarizationLayout layout;
    layout.call(GA);

    // Export the graph to a SVG file
    ogdf::GraphIO::write(GA, "../results/output.svg", ogdf::GraphIO::drawSVG);


    ogdf::GraphIO::write(G, "../results/clustered_graph.gml");
    std::cout << " Graph was written back succesfully " << std::endl ;
}

    
    

