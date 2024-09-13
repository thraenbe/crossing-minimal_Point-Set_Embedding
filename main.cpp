#include "JSONParser.hpp"
#include "matching.hpp"
#include "switch-crossmin.hpp"
#include <iostream>
#include <cassert>
#include <chrono>
#include <limits>
#include <ogdf/basic/Graph.h>
#include <ogdf/graphalg/Clusterer.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/SList.h>
#include <ogdf/graphalg/ModifiedNibbleClusterer.h>
#include <boost/graph/adjacency_list.hpp>

#include <xlsxwriter.h>
 

using size_t = std::size_t;

Graph read(size_t numNodes,std::vector<Node> nodes, std::vector<std::pair<size_t, size_t>>  edges, std::vector<Point> positions, int width, int height, JSONParser jp, const std::string& path) {
    // Width and height are kind of useless but have to be written back.
    jp.jsonToGraph(numNodes, nodes, edges, positions, width, height, path);
    const auto n = nodes.size();
    const auto m = edges.size();
    const auto f = positions.size();

	std::cout << "\n #Nodes : " << n << "   #Edges: " << m << "    #Positions: " << f << "     #Width: " << width << "    height: " << height << std::endl ;
    assert(f >= n);

    Graph G(n, nodes, edges, positions, width, height);


    return G;
}

void write(Graph G, const std::string& path, const std::string& name, int width, int height, JSONParser jp ){
    auto createOutputFileName = [&](){
        int collinear = 0;
    	return name + "_"+std::to_string(G.ComputeCrossings(collinear));
    };
    auto WriteGraph = [&]() {
        jp.graphToJson(G, path, createOutputFileName(), width, height); };
        int collinear = 0;
    size_t crossNum = G.ComputeCrossings(collinear);
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


void calculateBetweennessCentrality(ogdf::Graph &G, ogdf::EdgeArray<double> &edgeBetweenness) {
    ogdf::NodeArray<double> nodeBetweenness(G, 0.0);

    for (ogdf::node v : G.nodes) {
        // Single-source shortest paths (SSSP) and dependencies
        std::stack<ogdf::node> S;
        ogdf::NodeArray<std::vector<ogdf::node>> P(G);
        ogdf::NodeArray<int> sigma(G, 0);
        ogdf::NodeArray<int> d(G, -1);
        ogdf::NodeArray<double> delta(G, 0.0);
        sigma[v] = 1;
        d[v] = 0;

        std::queue<ogdf::node> Q;
        Q.push(v);

        while (!Q.empty()) {
            ogdf::node w = Q.front();
            Q.pop();
            S.push(w);

            for (ogdf::adjEntry adj : w->adjEntries) {
                ogdf::node u = adj->twinNode();

                if (d[u] < 0) {
                    Q.push(u);
                    d[u] = d[w] + 1;
                }

                if (d[u] == d[w] + 1) {
                    sigma[u] += sigma[w];
                    P[u].push_back(w);
                }
            }
        }

        while (!S.empty()) {
            ogdf::node w = S.top();
            S.pop();
            for (ogdf::node u : P[w]) {
                double c = ((double)sigma[u] / sigma[w]) * (1.0 + delta[w]);
                delta[u] += c;

                ogdf::edge e = u->firstAdj()->theEdge();
                edgeBetweenness[e] += c;
            }
        }
    }

    for (ogdf::edge e : G.edges) {
        edgeBetweenness[e] /= 2.0;  // Betweenness is counted twice for undirected graphs
    }
}

// Girvan-Newman algorithm implementation
void girvanNewman(ogdf::Graph &G, int desiredClusters, std::vector<ogdf::List<ogdf::node>> &clusters) {
    while (true) {
        // Step 1: Compute edge betweenness centrality
        ogdf::EdgeArray<double> edgeBetweenness(G, 0.0);
        calculateBetweennessCentrality(G, edgeBetweenness);

        // Step 2: Find and remove the edge with the highest betweenness centrality
        double maxBetweenness = -1.0;
        ogdf::edge maxEdge = nullptr;

        for (ogdf::edge e : G.edges) {
            if (edgeBetweenness[e] > maxBetweenness) {
                maxBetweenness = edgeBetweenness[e];
                maxEdge = e;
            }
        }

        if (maxEdge != nullptr) {
            G.delEdge(maxEdge);
            std::cout << "Removed edge with max betweenness: (" << maxEdge->source()->index() 
                      << ", " << maxEdge->target()->index() << ") with centrality " 
                      << maxBetweenness << std::endl;
        }

        // Step 3: Check connected components and store clusters
        ogdf::NodeArray<int> component(G);  // Component array for each node
        int numComponents = ogdf::connectedComponents(G, component);  // Get the number of components

        // Clear and store clusters
        clusters.clear();
        clusters.resize(numComponents);

        for (ogdf::node v : G.nodes) {
            int compID = component[v];  // Get the component ID for each node
            clusters[compID].pushBack(v);
        }

        // Print the clusters
        std::cout << "Number of clusters: " << numComponents << std::endl;
        for (int i = 0; i < numComponents; ++i) {
            std::cout << "Cluster " << i << ":";
            for (ogdf::node v : clusters[i]) {
                std::cout << " " << v->index();
            }
            std::cout << std::endl;
        }

        if (numComponents >= desiredClusters) {
            std::cout << "Reached desired number of clusters!" << std::endl;
            break;
        }

        if (G.numberOfEdges() == 0) {
            std::cout << "No more edges to remove." << std::endl;
            break;
        }
    }
}


std::vector<size_t> computeClusterSizes( Graph& myGraph, const std::vector<size_t>& initialSizes, size_t numClusters){
    std::vector<size_t>newClusterSizes(initialSizes.size(),0);
    auto tmpPoints = (myGraph.numPoints - myGraph.numNodes);
    const size_t quotient = tmpPoints / numClusters  ;
    std::cout << "cluster Sizes  and (Inital Sizes) for final Matching" << std::endl;
    

    for (int i = 0; i < initialSizes.size(); i++) {
        newClusterSizes[i] = initialSizes[i]+quotient;
        if (const size_t remainder = tmpPoints % numClusters;remainder > 0){
            ++newClusterSizes[i];
            tmpPoints-- ;
        }
        std::cout << newClusterSizes[i] <<" (" << initialSizes[i] << ") ; " ;
    }
    std::cout << std::endl;
    return newClusterSizes;
}

std::vector<ogdf::List<ogdf::node>> createOgdfGraph(ogdf::Graph& G, ogdf::GraphAttributes& GA, const Graph& myGraph ){
    std::vector<size_t> tmp;
    std::unordered_map<size_t, ogdf::node> nodeMap;
    std::vector<ogdf::List<ogdf::node>> ogdfClusters(myGraph.NodeClusters.size());
    for (const auto& node : myGraph.nodes){
        const ogdf::node newNode = G.newNode();
        assert(myGraph.mapVerticesToPoints.at(node.GetId()) < myGraph.points.size());
        Point point = myGraph.points[myGraph.mapVerticesToPoints.at(node.GetId())];
        assert(node.GetCluster() < ogdfClusters.size());
        ogdfClusters[node.GetCluster()].pushBack(newNode) ; 
        assert(point.GetX() >= 0 && point.GetY()  >= 0);
        assert(point.GetX() == static_cast<int>(point.GetX()) && point.GetY() == static_cast<int>(point.GetY()));
        GA.x(newNode) = (int) point.GetX();
        GA.y(newNode) = (int) point.GetY();
        assert(newNode->index() == node.GetId());
        GA.label(newNode) = std::to_string(point.GetId());
        tmp.push_back(node.GetId());
        nodeMap[node.GetId()] = newNode;
        
    }
    for (const auto& node : myGraph.nodes) {
        for (const auto target : myGraph.adjList[node.GetId()]){
            if (std::ranges::find(tmp,target)!= tmp.end()){
                if (node.GetId() < target){ // Every Edge is only appended once istead of twice
                    G.newEdge(nodeMap[node.GetId()], nodeMap[target]);
                }
            }
        }
    }
    return ogdfClusters;
}


ogdf::List<ogdf::node> getFreeOgdfPositions( Graph & myGraph,ogdf::Graph& G, ogdf::GraphAttributes& GA){
ogdf::List<ogdf::node> freePositions;
    std::vector<int> usedPoints;
    std::vector<int> freePoints;
    if(myGraph.points.size() > myGraph.nodes.size()){ 
        for(int i = 0; i < myGraph.nodes.size(); i++){
            usedPoints.push_back(myGraph.mapVerticesToPoints.at(i));
        }
        assert(usedPoints.size() == myGraph.numNodes );
        std::sort(usedPoints.begin(), usedPoints.end());
        int idx = 0;
        for(int i = 0; i < myGraph.points.size(); i++){
           // std::cout << usedPoints[idx] << "   " << i << "  " << idx << std::endl;
            if(usedPoints[idx] == i ){
                idx += 1;
            }else{
                freePoints.push_back(i);
            }
        }

        // std::cout << usedPoints.size() << "  ----------------------   " << freePoints.size() << std::endl ;
        assert(areValuesUnique(myGraph.mapVerticesToPoints));
        assert(usedPoints.size() + freePoints.size() == myGraph.points.size());
        for(auto i : freePoints){
            ogdf::node newNode = G.newNode();
            GA.x(newNode) = myGraph.points[i].GetX();
            GA.y(newNode) = myGraph.points[i].GetY();
            GA.label(newNode) = std::to_string(i);
            assert(std::stoi(GA.label(newNode)) >= 0 );
            assert(std::stoi(GA.label(newNode)) < myGraph.points.size());
            freePositions.pushBack(newNode) ;
        }
    }
    return freePositions;

}

void nodeClusteringFallback( Graph& myGraph, int& numberOfClusters){
    std::cout << "start Fallback Clustering \n";
    int fallBackClusterSize = (int) myGraph.numNodes/10;
    std::vector<std::vector<Node>> fallBackClusters(fallBackClusterSize);
    int count = 0;
    int cluster = 0;
    for(int i = 0; i < myGraph.numNodes; i++){
        if (count == 10){
            count = 0;
            cluster++;
        }
        Node& node = myGraph.nodes[i];
        node.SetCluster(cluster);
        fallBackClusters[cluster].push_back(node);
        count++;
    }
    std::cout << " finished Fallback Clustering \n";
    myGraph.NodeClusters = fallBackClusters;
    numberOfClusters = fallBackClusters.size();
}
void nodeClustering( Graph& myGraph, ogdf::Graph& G, int& numberOfClusters, const int clustering_stop_idx, const int clustering_automatic_threshold){
    ogdf::SList<ogdf::SimpleCluster *> clusters;
    ogdf::Clusterer CG(G);
    CG.setRecursive(false);
    CG.setStopIndex(clustering_stop_idx); // 0.6
    CG.setAutomaticThresholds(clustering_automatic_threshold); // 10

    CG.computeClustering(clusters);
    std::cout << " Number of Clusters: " << clusters.size()<<  std::endl ;
    myGraph.NodeClusters = myGraph.assignClustersToNodes(clusters);
    
    numberOfClusters = clusters.size();

}

std::vector<size_t> computation( const int numberOfOuterLoopsGlobal , const int numberOfSamples, const int numberOfOuterLoopsMove, const int numberOfOuterLoopsSwitch, int clustering_automatic_threshold, const double clustering_stop_idx, const int kk_des_edge_len, const int kk_global_iterations, const double kk_stop_tolerance, const string& graphfile ){
    std::vector<size_t> results(15,0);
    ogdf::Graph secondG; 
    ogdf::GraphAttributes secondGA(secondG, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate | ogdf::GraphAttributes::nodeWeight);


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
    ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate | ogdf::GraphAttributes::nodeWeight);
    myGraph.toOgdfGraph(G, GA);

    /*
     *2. Compute cluster graph C(G)
     */

    std::cout << "start clustering \n";
    int numberOfClusters;

    if(myGraph.edges.size() > 1000000){
        nodeClusteringFallback(myGraph, numberOfClusters);
    }
    else{
        nodeClustering(myGraph, G, numberOfClusters, clustering_stop_idx, clustering_automatic_threshold);
    }
    std::vector<size_t> clusterSizes;

    for (const auto& cluster : myGraph.NodeClusters){
        clusterSizes.push_back(cluster.size());
    }

    /*
     *3. Match ClusterGraph
     */

    std::cout << "Creating Cluster Graph\n";

    
    std::pair<Graph, std::vector<size_t>> pair = createClusteringGraph(myGraph, myGraph.NodeClusters.size());
    std::cout << "Matching Cluster Graph\n";

    matchClusters(pair.first, pair.second);

    const auto newClusterSizes = pair.first.assignClusters(myGraph, clusterSizes);
    auto augmentedClusterSizes = computeClusterSizes(myGraph, newClusterSizes, numberOfClusters);

    assert(pair.first.NodeClusters[0].size() == myGraph.NodeClusters[0].size());

    std::cout << "Clusters Succesfully assigned" << std::endl;

    /*
     *4. Cluster Points and Match
     */

    myGraph.manClustering(augmentedClusterSizes, myGraph.width, myGraph.height);
    std::cout << "Local Clustering Succesful" << std::endl;

    myGraph.reductCluster(myGraph.points, clusterSizes[0]);

    std::cout << "Get Point Clusters" << std::endl;

    myGraph.pointClusters = myGraph.getPointClusters(myGraph.points, numberOfClusters);

    std::cout << "Start Matching" << std::endl;

    myGraph.mapVerticesToPoints = matching(myGraph, myGraph.NodeClusters, myGraph.freePoints, myGraph.usedPoints, myGraph.points);
    if (!areValuesUnique(myGraph.mapVerticesToPoints)){
        simpleAssign(myGraph, myGraph.nodes.size());
    }
    assert(areValuesUnique(myGraph.mapVerticesToPoints));

    int collinear = 0;
    
    /*
     *5. Swaping and Moving Nodes to reduce Crossings
     */


    results[4] =  myGraph.ComputeCrossings(collinear);    
    results[8] = collinear;
    //IterativeCrossMin for Clusters
    ogdf::Graph crossminG;
    ogdf::GraphAttributes crossminGA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate);
    auto ogdfCluster = createOgdfGraph(crossminG, crossminGA, myGraph);
    auto ogdfFreePoints = getFreeOgdfPositions(myGraph, secondG, secondGA);
    std::cout << "Starting Iterative Crossing Minimization (SWITCH)\n";
    for (int i = 0 ; i < myGraph.NodeClusters.size(); i++){
        iterativeCrossMinSwitch(crossminG, crossminGA, numberOfOuterLoopsSwitch, numberOfSamples, ogdfCluster[i], myGraph.numNodes);
        std::cout << "    Completed Local Crossmin for Cluster: " << i  << std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << "Completed Iterative Crossing Minimization (Local)\n";
    for (const auto &n : crossminG.nodes){
        auto nodeId = n->index();
        auto pointId = std::stoi(crossminGA.label(n));
        assert(nodeId < myGraph.nodes.size());
        assert(pointId < myGraph.points.size());
        myGraph.mapVerticesToPoints.at(nodeId) = pointId;
    }

    results[5] = myGraph.ComputeCrossings(collinear);    
    results[9] = collinear;
    ogdf::List<ogdf::node> allNodes ;
    crossminG.allNodes(allNodes);
    
    if (ogdfFreePoints.size() > 0 ){ 
        std::cout << "Starting Iterative Crossing Minimization (Global)\n";
        for (int i = 0; i < numberOfOuterLoopsGlobal  ; i++){
            iterativeCrossMinMove(crossminG, crossminGA, secondGA, numberOfOuterLoopsMove, numberOfSamples, ogdfFreePoints, myGraph.numNodes );
            iterativeCrossMinSwitch(crossminG, crossminGA, numberOfOuterLoopsSwitch, numberOfSamples, allNodes, numNodes);
            std::cout << "     Completed Global Crossmin for Iteration: " << i << " of " << numberOfOuterLoopsGlobal << std::endl; ;
        }
        std::cout << std::endl << std::endl;
    }
    else{
        std::cout << "Skipped Iterative Crossing Minimization (Global)\n";
        for (int i = 0; i < numberOfOuterLoopsGlobal ; i++){
            iterativeCrossMinSwitch(crossminG, crossminGA, numberOfOuterLoopsSwitch, numberOfSamples, allNodes, numNodes);
        }
    }


    for (const auto &n : crossminG.nodes){
        auto nodeId = n->index();
        auto pointId = std::stoi(crossminGA.label(n));
        assert(nodeId < myGraph.nodes.size());
        assert(pointId < myGraph.points.size());
        // std::cout << "NODE: " << nodeId << "   PointId: " << pointId << std::endl;
        myGraph.mapVerticesToPoints.at(nodeId) = pointId;
    }

    /*
        6. Write Back
    */

    std::cout << "Finished Iterative Crossing Minimization (MOVE)\n";
    results[6] = myGraph.ComputeCrossings(collinear);    
    results[10] = collinear;

    // Apply a layout algorithm
    // ogdf::PlanarizationLayout layout;
    
    // layout.call(GA);

    // Export the graph to a SVG file
    // ogdf::GraphIO::write(GA, "../results/output.svg", ogdf::GraphIO::drawSVG);


    // ogdf::GraphIO::write(G, "../results/clustered_graph.gml");
    // std::cout << " Graph was written back succesfully " << std::endl ;

    results[0] = myGraph.nodes.size();
    results[1] = myGraph.NodeClusters.size();
    results[2] = (int) myGraph.nodes.size() / myGraph.NodeClusters.size();
    results[3] = myGraph.randomCrossingNumber(collinear);
    results[7] = collinear;

    return results;


}


int main(int argc, char* argv[]) {
    string graphFile{"g9"};     // n = 1200.
    int clustering_automatic_threshold = 15;
    double clustering_stop_idx = 0.6;
    int numberOfOuterLoopsMove = 1;
    int numberOfOuterLoopsSwitch = 1;
    int numberOfOuterLoopsGlobal = 1;

  // TODO:: the following are currently unused - remove in case they are not required.
    int kk_des_edge_len{1};
    int kk_global_iterations{100};
    double kk_stop_tolerance{.1};

    const std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    const std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;


    std::cout << "opening Workbook" << std::endl;

    lxw_workbook  *workbook  = workbook_new("../results/results-real-clustering.xlsx");
    lxw_worksheet *worksheet = workbook_add_worksheet(workbook, NULL);

    worksheet_write_string(worksheet, 0, 0, "Graph", NULL);
    worksheet_write_string(worksheet, 0, 1, "Parameter", NULL);
    worksheet_write_string(worksheet, 0, 2, "Number of Nodes", NULL);
    worksheet_write_string(worksheet, 0, 3, "Number Of Clusters", NULL);
    worksheet_write_string(worksheet, 0, 4, "Average Cluster Size", NULL);
    worksheet_write_string(worksheet, 0, 5, "Average Cluster Size", NULL);
    worksheet_write_string(worksheet, 0, 6, "Random Crossings", NULL);
    worksheet_write_string(worksheet, 0, 7, "Crossings before Iterative Crossing Minimization ", NULL);
    worksheet_write_string(worksheet, 0, 8, "Final Crossings", NULL);


    // for(int clustering_automatic_threshold = 5; clustering_automatic_threshold < 50; clustering_automatic_threshold = clustering_automatic_threshold+5 )

    int idx = 1;
    for( int i = 45; i < 150; i++){
        worksheet_write_number(worksheet, i*10 + 1, 0, i, NULL);
        idx++;
        for (int k = 1 ; k < 1.5; k++){ 
            for(double j = 1; j < 51; j += 5){ 
                for(int p = 1; p < 1.3; p += 1 ){ 
                    clustering_automatic_threshold = j;
                    int numberOfSamples = 1;


                    string graphfile = "rand_" + std::to_string(i+1);

                    auto crossingVector = computation(numberOfOuterLoopsGlobal, numberOfSamples, numberOfOuterLoopsMove, numberOfOuterLoopsSwitch , 
                                                        clustering_automatic_threshold, 
                                                        clustering_stop_idx, 
                                                        kk_des_edge_len, 
                                                        kk_global_iterations, 
                                                        kk_stop_tolerance, 
                                                        graphfile);


                    worksheet_write_number(worksheet,  idx, 1 , j, NULL);
                    worksheet_write_number(worksheet,  idx, 2, p, NULL);
                    worksheet_write_number(worksheet,  idx, 3, crossingVector[0], NULL);
                    worksheet_write_number(worksheet,  idx, 4, crossingVector[1], NULL);
                    worksheet_write_number(worksheet,  idx, 5, crossingVector[2], NULL);
                    worksheet_write_number(worksheet,  idx, 6, crossingVector[3], NULL);
                    worksheet_write_number(worksheet,  idx, 7, crossingVector[4], NULL);
                    worksheet_write_number(worksheet,  idx, 8, crossingVector[5], NULL);
                    worksheet_write_number(worksheet,  idx, 9, crossingVector[6], NULL);
                    worksheet_write_number(worksheet,  idx, 10, crossingVector[7], NULL);
                    worksheet_write_number(worksheet,  idx, 11, crossingVector[8], NULL);
                    worksheet_write_number(worksheet,  idx, 12, crossingVector[9], NULL);
                    worksheet_write_number(worksheet,  idx, 13, crossingVector[10], NULL);

                    std::cout << "Graph " << i+1 << " computed successful" << std::endl;
                    idx++;
                }
            }
        }
    }

    workbook_close(workbook);

    return 0;

}

    
    

