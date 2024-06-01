#include"JSONParser.hpp"
#include <random>
#include <iostream>
#include <cassert>
#include <iterator>
#include <limits>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/fileformats/GraphIO.h>

 

using size_t = std::size_t;

Graph read(size_t numNodes, std::vector<std::pair<size_t, size_t>>  edges, std::vector<Point> positions, int width, int height, JSONParser jp, std::string path) {
    // Width and height are kind of useless but have to be written back.
    jp.jsonToGraph(numNodes, edges, positions, width, height, path);
    const int n = (int) numNodes;
    const int m = edges.size();
    const int f = positions.size();
	std::cout << "#Nodes : " << n << "   #Edges: " << m << "    #Positions: " << f << std::endl;
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


int main(int argc, char* argv[]) {



    std::vector<Point> positions;
    std::vector<std::pair<size_t, size_t>> edges;

    size_t numNodes{0};
    int width;
    int height;
    JSONParser jp;

    const std::string filepath = "../graphs/";
    const std::string outPath = "results";
    const std::string filename = "g5";

    Graph G = read(numNodes, edges, positions, width, height, jp, filepath + filename);



    write(G, outPath, filename, width, height, jp);
}

    
    

