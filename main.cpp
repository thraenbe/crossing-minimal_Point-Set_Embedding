#include"JSONParser.hpp"
#include"matching.hpp"
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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

#include <xlsxwriter.h>

 

using size_t = std::size_t;

Graph read(size_t numNodes,std::vector<Node> nodes, std::vector<std::pair<size_t, size_t>>  edges, std::vector<Point> positions, int width, int height, JSONParser jp, std::string path) {
    // Width and height are kind of useless but have to be written back.
    jp.jsonToGraph(numNodes, nodes, edges, positions, width, height, path);
    const int numN = (int) numNodes;
    const int n = nodes.size();
    const int m = edges.size();
    const int f = positions.size();

	std::cout << "#Nodes : " << n << "   #Edges: " << m << "    #Positions: " << f << "     #Width: " << width << "    height: " << height << std::endl ;
    assert(f >= n);

    Graph G(n, nodes, edges, positions, width, height);


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

void printVector(const std::vector<std::vector<Node>> vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << "Cluster " << i << "  : " ;
        for (auto elem : vec[i]) {
            std::cout << elem.getId() << " " ;
        }
        std::cout << std::endl;
    }
}

std::vector<int> computeClusterSizes(Graph myGraph, std::vector<int> clusterSizes, int numClusters){
    int tmpPoints = myGraph.numPoints - myGraph.numNodes;
    int quotient = tmpPoints / numClusters  ;
    std::cout << tmpPoints << std::endl;

    for (int i = 0; i < clusterSizes.size(); i++) {
        int rest = tmpPoints % numClusters;
        if (rest > 0){
            clusterSizes[i] += quotient + 1;
            tmpPoints-- ;
        }
        else{
            clusterSizes[i] += quotient;
        }
    }
    return clusterSizes;
}

std::vector<int> computation( int clustering_automatic_threshold, double clustering_stop_idx, int kk_des_edge_len, int kk_global_iterations, double kk_stop_tolerance, string graphfile ){

    std::vector<Point> positions;
    std::vector<Node> nodes;
    std::vector<std::pair<size_t, size_t>> edges;
    size_t numNodes{0};

    int width;
    int height;

    JSONParser jp;

    const std::string filepath = "../graphs/";
    const std::string outPath = "../results";
    const std::string filename = graphfile;

    Graph myGraph = read(numNodes, nodes, edges, positions, width, height, jp, filepath + filename);

    std::cout << " Graph was read succesfully " << std::endl ;



    ogdf::Graph G; 
    ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);



    myGraph.toOgdfGraph(G, GA);
    

    ogdf::SList<ogdf::SimpleCluster *> clusters; 
    // ogdf::ClusterGraph ClusterGraph(G);
    


    ogdf::Clusterer CG(G);

    double avgCIdx = CG.averageCIndex();
    for (ogdf::node v : G.nodes) {
        double CIdx = CG.computeCIndex(G, v); 
    }
    std::cout << avgCIdx << std::endl;

    CG.setRecursive(false);
    CG.setStopIndex(clustering_stop_idx); // 0.6

    CG.setAutomaticThresholds(clustering_automatic_threshold); // 10

    CG.computeClustering(clusters);
    // CG.createClusterGraph(ClusterGraph); // not really working 
    
    
    // std::cout << " Number of Clusters: " << clusters.size() << " Cluster Graph Clusters: " << ClusterGraph.numberOfClusters()<<  std::endl ;




    // Print the clusters
    std::vector<int> clusterSizes;

    
    myGraph.NodeClusters = myGraph.assignClustersToNodes(clusters);



    for (auto cluster : myGraph.NodeClusters){
        clusterSizes.push_back(cluster.size());
    }

    clusterSizes = computeClusterSizes(myGraph, clusterSizes, clusters.size());

    for (auto cluster : clusterSizes) {
        std::cout << "Cluster size " << cluster  << std::endl;
    }
    
    std::cout << "Creating Cluster Graph\n";
    Graph ClusterGraph = createClusterGraph(myGraph, myGraph.NodeClusters.size());
    std::cout << "Matching Cluster Graph\n";

    matchClusters(ClusterGraph);




    // myGraph.kMeans(myGraph.points,clusters.size() , clusterSizes);
    myGraph.manClustering(clusterSizes, myGraph.width, myGraph.height);
    std::cout << "Local Clustering Succesful" << std::endl;

    myGraph.reductCluster(myGraph.points,clusterSizes[0]);
    printVector(myGraph.NodeClusters);

    myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, clusters.size());

    matching(myGraph, myGraph.NodeClusters, myGraph.freePoints, myGraph.usedPoints, myGraph.points);

    int crossings = myGraph.ComputeCrossings();    

    std::cout << "The Graph has " << crossings << " crossings \n";


    // Apply a layout algorithm
    ogdf::PlanarizationLayout layout;
    
    layout.call(GA);

    // Export the graph to a SVG file
    ogdf::GraphIO::write(GA, "../results/output.svg", ogdf::GraphIO::drawSVG);


    ogdf::GraphIO::write(G, "../results/clustered_graph.gml");
    std::cout << " Graph was written back succesfully " << std::endl ;

    std::vector<int> crossingTypes;
    crossingTypes.push_back( crossings);
    crossingTypes.push_back(  myGraph.randomCrossingNumber());

    return crossingTypes;


}


int main(int argc, char* argv[]) {
    string graphfile;
    int clustering_automatic_threshold = 10;
    double clustering_stop_idx = 0.6;
    int kk_des_edge_len;
    int kk_global_iterations;
    double kk_stop_tolerance;

    int crossings;
    int randomCrossings;

    lxw_workbook  *workbook  = workbook_new("../results/results.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);

    // Write some data.
    worksheet_write_string(worksheet, 0, 0, "Graph", NULL);
    worksheet_write_string(worksheet, 1, 0, "Crossings", NULL);
    worksheet_write_string(worksheet, 2, 0, "RandomCrossings", NULL);
    for(int clustering_automatic_threshold = 5; clustering_automatic_threshold < 50; clustering_automatic_threshold = clustering_automatic_threshold+5 )
    for(double j = 0.1; j < 1; j = j + 0.1){
        clustering_stop_idx = j;
        for( int i = 0; i < 6; i++){
            graphfile = "g" + std::to_string(i+1);
            crossings = computation(clustering_automatic_threshold, clustering_stop_idx, kk_des_edge_len, kk_global_iterations, kk_stop_tolerance, graphfile).front() ;
            randomCrossings = computation(clustering_automatic_threshold, clustering_stop_idx, kk_des_edge_len, kk_global_iterations, kk_stop_tolerance, graphfile)[1] ;

            worksheet_write_number(worksheet, j*10*4 , i + 1, i+1, NULL);
            worksheet_write_number(worksheet, j*10*4 +1 , i + 1, crossings, NULL);
            worksheet_write_number(worksheet, j*10*4 +2 , i + 1, randomCrossings, NULL);
        }
    }

    workbook_close(workbook);

    return 0;

}

    
    

