#include"JSONParser.hpp"
#include"matching.hpp"
#include <iostream>
#include <cassert>
#include <chrono>
#include <limits>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/SList.h>
#include <boost/graph/adjacency_list.hpp>


 

using size_t = std::size_t;

Graph read(size_t numNodes,std::vector<Node> nodes, std::vector<std::pair<size_t, size_t>>  edges, std::vector<Point> positions, int width, int height, JSONParser jp, const std::string& path) {
    // Width and height are kind of useless but have to be written back.
    jp.jsonToGraph(numNodes, nodes, edges, positions, width, height, path);
    const auto n = nodes.size();
    const auto m = edges.size();
    const auto f = positions.size();

	std::cout << "#Nodes : " << n << "   #Edges: " << m << "    #Positions: " << f << "     #Width: " << width << "    height: " << height << std::endl ;
    assert(f >= n);

    Graph G(n, nodes, edges, positions, width, height);


    return G;
}

void write(Graph G, const std::string& path, const std::string& name, int width, int height, JSONParser jp ){
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

void printVector(const std::vector<std::vector<Node>>& vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "Cluster " << i << "  : " ;
        for (auto elem : vec[i]) {
            std::cout << elem.GetId() << " " ;
        }
        std::cout << std::endl;
    }
}

std::vector<size_t> computeClusterSizes(const Graph& myGraph, const std::vector<size_t>& initialSizes, size_t numClusters){
    std::vector<size_t>newClusterSizes(initialSizes.size(),0);
    auto tmpPoints = (myGraph.numPoints - myGraph.numNodes);
    const size_t quotient = tmpPoints / numClusters  ;
    std::cout << tmpPoints << std::endl;

    for (int i = 0; i < initialSizes.size(); i++) {
        newClusterSizes[i] = initialSizes[i]+quotient;
        if (const size_t remainder = tmpPoints % numClusters;remainder > 0){
            ++newClusterSizes[i];
            tmpPoints-- ;
        }
    }
    return newClusterSizes;
}

std::vector<int> computation( int clustering_automatic_threshold, double clustering_stop_idx, int kk_des_edge_len, int kk_global_iterations, double kk_stop_tolerance, const string& graphfile ){


    /*
     *1.READ GRAPH
     */
    std::vector<Point> positions;
    std::vector<Node> nodes;
    std::vector<std::pair<size_t, size_t>> edges;
    size_t numNodes{0};
    int width;
    int height;
    JSONParser jp;
    const std::string filepath = "../graphs/";
    const std::string outPath = "../results";
    Graph myGraph = read(numNodes, nodes, edges, positions, width, height, jp, filepath + graphfile);
    std::cout << " Graph was read succesfully " << std::endl ;

    /*
     *1.1 Convert graph to OGDF graph
     */

    ogdf::Graph G; 
    ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
    myGraph.toOgdfGraph(G, GA);

    /*
     *2. Compute cluster graph C(G)
     */
    ogdf::SList<ogdf::SimpleCluster *> clusters;
    ogdf::Clusterer CG(G);
    /*
    double avgCIdx = CG.averageCIndex();

    for (ogdf::node v : G.nodes) {
        double CIdx = CG.computeCIndex(G, v); 
    }
    std::cout << avgCIdx << std::endl;
    */
    CG.setRecursive(false);
    CG.setStopIndex(clustering_stop_idx); // 0.6
    CG.setAutomaticThresholds(clustering_automatic_threshold); // 10
    CG.computeClustering(clusters);
    // std::cout << " Number of Clusters: " << clusters.size() << " Cluster Graph Clusters: " << ClusterGraph.numberOfClusters()<<  std::endl ;
    std::vector<size_t> clusterSizes;
    myGraph.NodeClusters = myGraph.assignClustersToNodes(clusters);

    for (const auto& cluster : myGraph.NodeClusters){
        clusterSizes.push_back(cluster.size());
    }
    auto augmentedClusterSizes = computeClusterSizes(myGraph, clusterSizes, clusters.size());
    std::cout << "Creating Cluster Graph\n";
    std::pair<Graph, std::vector<size_t>> pair = createClusteringGraph(myGraph, myGraph.NodeClusters.size());
    std::cout << "Matching Cluster Graph\n";

    matchClusters(pair.first, pair.second);
    

    std::cout << "Old Cluster Sizes: " << std::endl;

    for (auto cluster : clusterSizes) {
        std::cout << ", " << cluster ;
    }
    std::cout << std::endl;

    const auto newClusterSizes = pair.first.assignClusters(myGraph, clusterSizes);

    std::cout << "New Cluster Sizes: " << std::endl;

    for (const auto c: newClusterSizes)
    {
        std::cout << ", " << c ;
    }
    std::cout << std::endl;
    for (int i = 0; i < newClusterSizes.size(); i++) {
        std::cout << ", " << myGraph.NodeClusters[i].size() ;
    }
    std::cout << std::endl;



    assert(pair.first.NodeClusters[0].size() == myGraph.NodeClusters[0].size());

    std::cout << "Clusters Succesfully assigned" << std::endl;

    // myGraph.kMeans(myGraph.points,clusters.size() , clusterSizes);
    myGraph.manClustering(newClusterSizes, myGraph.width, myGraph.height);
    std::cout << "Local Clustering Succesful" << std::endl;

    myGraph.reductCluster(myGraph.points,clusterSizes[0]);
    printVector(myGraph.NodeClusters);

    myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, clusters.size());

    matching(myGraph, myGraph.NodeClusters, myGraph.freePoints, myGraph.usedPoints, myGraph.points);

    int crossings = myGraph.ComputeCrossings();    

    std::cout << "The Graph has " << crossings << " crossings \n";


    // Apply a layout algorithm
    // ogdf::PlanarizationLayout layout;
    
    // layout.call(GA);

    // Export the graph to a SVG file
    // ogdf::GraphIO::write(GA, "../results/output.svg", ogdf::GraphIO::drawSVG);


    // ogdf::GraphIO::write(G, "../results/clustered_graph.gml");
    // std::cout << " Graph was written back succesfully " << std::endl ;

    std::vector<int> results;
    results.push_back( crossings);
    std::cout << " Compute Random Cossing Number: \n";   
    results.push_back(  myGraph.randomCrossingNumber());
    results.push_back( clusterSizes.size());

    return results;


}


int main(int argc, char* argv[]) {
    const string graphFile{"g9"};     // n = 1200.
    int clustering_automatic_threshold = 10;
    double clustering_stop_idx = 0.6;
    // TODO:: the following are currently unused - remove in case they are not required.
    int kk_des_edge_len{1};
    int kk_global_iterations{100};
    double kk_stop_tolerance{.1};

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const auto crossingVector = computation(clustering_automatic_threshold, clustering_stop_idx, kk_des_edge_len, kk_global_iterations, kk_stop_tolerance, graphFile);
    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
    const auto crossings = crossingVector[0];
    std::cout << crossings << std::endl;
    const auto randomCrossings = crossingVector[1];
    std::cout << randomCrossings << std::endl;

    return 0;

}

    
    

