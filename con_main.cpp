#include "JSONParser.hpp"
#include "matching.hpp"
#include "graphGrowing.hpp"
#include "light_crossmin.hpp"
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

Graph read(size_t numNodes,std::vector<Node> nodes, std::vector<Edge>  edges, std::vector<Point> positions, int& width, int& height, JSONParser jp, const std::string& path) {
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
    std::cout << numClusters << " Check for Float \n";
    const size_t quotient = (int) tmpPoints / numClusters  ;
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

void assignNodePosToPointPos( Graph& myGraph ){
    for (auto& node : myGraph.nodes){
        assert(myGraph.mapVerticesToPoints.at(node.GetId()) < myGraph.points.size());
        assert(areValuesUnique(myGraph.mapVerticesToPoints));

        Point &point = myGraph.points[myGraph.mapVerticesToPoints.at(node.GetId())];
        assert(point.GetX() >= 0 && point.GetY()  >= 0);
        assert(point.GetX() == static_cast<int>(point.GetX()) && point.GetY() == static_cast<int>(point.GetY()));
        node._x = (int) point.GetX();
        node._y = (int) point.GetY();
        node._pointId = point.GetId();
    }
}


std::vector<size_t> getFreePositions( Graph &myGraph){
    std::vector<Point> freePoints;
    std::vector<int> usedPointIds;
    std::vector<size_t> freePointIds;
    assert(areValuesUnique(myGraph.mapVerticesToPoints));
    if(myGraph.points.size() > myGraph.nodes.size()){ 
        for(int i = 0; i < myGraph.nodes.size(); i++){
            usedPointIds.push_back(myGraph.mapVerticesToPoints.at(i));
        }
        assert(usedPointIds.size() == myGraph.numNodes );
        std::sort(usedPointIds.begin(), usedPointIds.end());
        int idx = 0;
        for(int i = 0; i < myGraph.points.size(); i++){
           // std::cout << usedPoints[idx] << "   " << i << "  " << idx << std::endl;
            if(usedPointIds[idx] == i ){
                idx += 1;
            }else{
                freePointIds.push_back(i);
                freePoints.push_back(myGraph.points[i]);
            }
        }

        // std::cout << usedPoints.size() << "  ----------------------   " << freePoints.size() << std::endl ;
        assert(areValuesUnique(myGraph.mapVerticesToPoints));
        assert(usedPointIds.size() + freePointIds.size() == myGraph.points.size());
    }
    return freePointIds;
}

void mapVerticesToPoints(Graph &myGraph){
    assert(areValuesUnique(myGraph.mapVerticesToPoints));
    for (const auto &n : myGraph.nodes){
        auto nodeId = n.GetId();
        auto pointId = n._pointId;
        assert(nodeId < myGraph.nodes.size());
        assert(pointId < myGraph.points.size());
        assert(pointId == myGraph.mapVerticesToPoints.at(nodeId));
        myGraph.mapVerticesToPoints.at(nodeId) = pointId;
    }
    assert(areValuesUnique(myGraph.mapVerticesToPoints));
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
void nodeClustering( Graph& myGraph, ogdf::Graph& G, int& numberOfClusters, const int clustering_automatic_threshold){
    ogdf::SList<ogdf::SimpleCluster *> clusters;
    ogdf::Clusterer CG(G);
    CG.setRecursive(false);
    CG.setAutomaticThresholds(clustering_automatic_threshold); // 10

    CG.computeClustering(clusters);
    std::cout << " Number of Clusters: " << clusters.size()<<  std::endl ;
    myGraph.NodeClusters = myGraph.assignClustersToNodes(clusters);
    
    numberOfClusters = clusters.size();

}

int setNumberOfOuterLoopsGlobal(size_t numEdges, int numNodes){
    if (numEdges > 100000){
        return 0;
    }
    else if (numEdges >= 10000 || numNodes > 1000){
        return 1;
    }
    else if (numEdges >= 1000){
        return 10;
    }
    else if (1000 > numEdges){
        return 100;
    }
    
    return 0;
    
}
int setNumberOfSamples(size_t numEdges, int numNodes){
    if (numEdges >= 100000){
        return 1;
    }
    else if (numEdges > 10000){
        return 3;
    }
    else if (numEdges > 1000){
        return 10;
    }
    else if (1000 > numEdges){
        return static_cast<int>(round(numNodes/2));
    }
    return 0;
    
}
int setClusteringAutomaticThreshold(size_t numEdges, int numNodes){
    int roundedSqrt = static_cast<int>(round(sqrt(numNodes)));
    return roundedSqrt;
}

std::vector<size_t> computation(  int numberOfOuterLoopsGlobal , int numberOfSamples,  int numberOfOuterLoopsMove,  int numberOfOuterLoopsLocal, int clustering_automatic_threshold, const string& graphfile ){
    std::vector<size_t> results(15,0);

    /*
     *1.READ GRAPH
     */
    std::vector<Point> positions;
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    size_t numNodes{0};
    int width;
    int height;
    JSONParser jp;
    const std::string filepath = "../graphs/";
    const std::string outPath = "../results";
    Graph myGraph = read(numNodes, nodes, edges, positions, width, height, jp, filepath + graphfile);
    std::cout << " Graph was read succesfully " << std::endl ;

    /*
     * 1.1 Set Parameters
    */



    /*
     *1.2 Convert graph to OGDF graph
     */

    ogdf::Graph G; 
    ogdf::GraphAttributes GA(G, ogdf::GraphAttributes::nodeGraphics | ogdf::GraphAttributes::edgeGraphics | ogdf::GraphAttributes::nodeLabel | ogdf::GraphAttributes::edgeStyle | ogdf::GraphAttributes::nodeStyle | ogdf::GraphAttributes::nodeTemplate | ogdf::GraphAttributes::nodeWeight);
    myGraph.toOgdfGraph(G, GA);

    /*
     *2. Compute cluster graph C(G)
     */

    std::cout << "start clustering \n";
    int numberOfClusters;

    if(myGraph.edges.size() > 130000){
        nodeClusteringFallback(myGraph, numberOfClusters);
    }
    else if (false){
        std::cout << "Try Graph growing Clustering \n" ;
        numberOfClusters = setClusteringAutomaticThreshold(myGraph.edges.size(), myGraph.numNodes);
        std::vector<std::vector<Node>> BFSClusters(numberOfClusters);
        vector<int> membership = graphGrowingPartition(myGraph.adjList, numberOfClusters);

        cout << membership.size() << "Cluster membership:\n";
        assert(membership.size() == myGraph.nodes.size());
        for (int i = 0; i < membership.size(); ++i) {
            assert(i < myGraph.nodes.size());
            assert(i >= 0);
            int cluster = membership[i];
            Node& node = myGraph.nodes[i];
            // std::cout << "Node " << i << " -> Cluster " << cluster << endl;
            node.SetCluster(cluster);
            BFSClusters[cluster].push_back(node);
        }
        myGraph.NodeClusters = BFSClusters;
    }
    else{
        nodeClustering(myGraph, G, numberOfClusters, clustering_automatic_threshold);
    }
    std::vector<size_t> clusterSizes;

    for (const auto& cluster : myGraph.NodeClusters){
        std::cout << cluster.size() << std::endl;
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
    std::cout << "assigned Clusters \n ";
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
    int crossings = myGraph.ComputeCrossings(collinear);
    
    /*
     *5. Swaping and Moving Nodes to reduce Crossings
     */

    //IterativeCrossMin for Clusters
    assignNodePosToPointPos(myGraph);
    std::vector<size_t> freePoints = getFreePositions(myGraph);

    std::cout << "Crossings: " <<     myGraph.ComputeCrossings(collinear) << std::endl;
    std::cout << "Starting Iterative Local Crossing Minimization \n";



    iterativeCrossMinLocal(myGraph, numberOfOuterLoopsLocal, freePoints, numberOfSamples, numberOfOuterLoopsMove, myGraph.numNodes);
    assert(areValuesUnique(myGraph.mapVerticesToPoints));

    std::cout << "Crossings: " <<   myGraph.ComputeCrossings(collinear) << std::endl;
    std::cout << "Starting Global Iterative Crossing Minimization \n";
    

    
    //iterativeCrossMinGlobal(myGraph, numberOfOuterLoopsGlobal, numberOfOuterLoopsMove, freePoints, numberOfSamples, myGraph.numNodes, crossings);
    assert(areValuesUnique(myGraph.mapVerticesToPoints));

    std::cout << "Crossings: " <<     myGraph.ComputeCrossings(collinear) << std::endl;

    /*
        6. Write Back
    */

    for (auto m : myGraph.mapVerticesToPoints){
        std::cout << m.second << std::endl;
    }
    assert(areValuesUnique(myGraph.mapVerticesToPoints) );

    write(myGraph, "../results" , graphfile, width, height, jp);


    return results;


}


int main(int argc, char* argv[]) {
    
    string graphfile{"g1"};     // n = 1200.
    int clustering_automatic_threshold = 6;

    int numberOfOuterLoopsGlobal = 1000;
    int numberOfSamples = 1000 ;
    int numberOfOuterLoopsMove = 4;
    int numberOfOuterLoopsLocal = 10000;

    for (int i = 0; i < 10; i++){
    auto crossingVector = computation(numberOfOuterLoopsGlobal, numberOfSamples, numberOfOuterLoopsMove, numberOfOuterLoopsLocal, 
                                        clustering_automatic_threshold, 
                                        graphfile);
    }

    return 0;

}

    
    

