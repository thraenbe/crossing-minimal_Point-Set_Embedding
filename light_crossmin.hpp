#include "Geometry.hpp"
#include "Node.hpp"
#include "matching.hpp"
#include "rtree.h"
#include "Graph.hpp"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/geometry.h> // For DPoint and DSegment
#include <omp.h>
#include <ostream>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>


void getEdgeCrossingsWithRTree(const RTree<int, int, 2, float> &RTree, const Graph &G, const Edge &e,
                                                       unsigned int &cr, 
                                                       const int numNodes) {
    int px = (int) G.nodes[e.first].GetX();  
    int py = (int) G.nodes[e.first].GetY();
    int qx = (int) G.nodes[e.second].GetX();
    int qy = (int) G.nodes[e.second].GetY();

    

   Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));

    std::vector<int> treeRes;
    int hits = RTree.Search(rect.min, rect.max, &treeRes);
    for (int i = 0; i < hits; ++i) {
        // std::cout << G.nodes[e.first].GetX() << " " << py << " " << qx << " " << qy << " ############ " << e.id << "\n";
        if(e.id != G.edges[treeRes[i]].id){
        cr = cr + G.DoEdgesIntersectVar(e, G.edges[treeRes[i]]);
        assert(e.id != G.edges[treeRes[i]].id);
        }
    }
}


void removeFromRTree(RTree<int, int, 2, float> &RTree, Graph &G, const Node &n) {
    for (auto const &e : G.getAdjEdges(n)) {
        int px = (int) G.nodes[e.first].GetX();  
        int py = (int) G.nodes[e.first].GetY();
        int qx = (int) G.nodes[e.second].GetX();
        int qy = (int) G.nodes[e.second].GetY();
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Remove(rect.min, rect.max, e.id);
    }
}

void insertIntoRTree(RTree<int, int, 2, float> &RTree, Graph &G, const Node &n) {
    for (auto const &e : G.getAdjEdges(n)) {
        int px = (int) G.nodes[e.first].GetX();  
        int py = (int) G.nodes[e.first].GetY();
        int qx = (int) G.nodes[e.second].GetX();
        int qy = (int) G.nodes[e.second].GetY();
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, e.id);
    }
}

bool isNumeric(const std::string &str) {
    return !str.empty() && std::all_of(str.begin(), str.end(), ::isdigit);
}

std::vector<int> generate_random_indices(int k, int n) {
    // Create a vector with values from 0 to n-1
    if (k> n){
        k = n;
    }
    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) {
        indices[i] = i;
    }
    
    // Use random_device and mt19937 for randomness
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Shuffle the indices vector
    std::shuffle(indices.begin(), indices.end(), gen);
    
    // Resize the vector to contain only the first k elements
    indices.resize(k);
    
    return indices;
}

void swapNodePositions(Node &n, Node &m){
	auto tmp_x = n.GetX();
	auto tmp_y = n.GetY();
	n._x = m._x;
	n._y = m._y;
	m._x = tmp_x;
	m._y = tmp_y;
}
void swapPositions(Node &n, Point &m){
	auto tmp_x = n.GetX();
	auto tmp_y = n.GetY();
	n._x = m._x;
	n._y = m._y;
	m._x = tmp_x;
	m._y = tmp_y;
}

std::vector<Point> freePointsToMove(std::vector<Point> points, std::vector<size_t> freePoints){
    std::vector<Point> freePointsInCluster;
    for (auto &p : points){
        std::cout << p.GetId() << " " << p.GetCluster() << std::endl;
        if (std::ranges::find(freePoints, p.GetId()) != std::ranges::end(freePoints)){
            freePointsInCluster.push_back(p);
            std::cout << p.GetId() << std::endl;
        }
    }
    for (auto f : freePoints){
        std::cout << "TEST " << f << std::endl;
    }
    return freePointsInCluster;
}

std::vector<std::size_t> freePointsToMove(const Graph &G, const size_t cluster) {
    std::vector<std::size_t> freePoints;
    std::vector<std::size_t> usedPoints;
    std::vector<std::size_t> allPoints; 

    // Sammeln der IDs der genutzten Punkte
    for (auto n : G.nodes) { // G.NodeClusters[cluster]
        usedPoints.push_back(G.mapVerticesToPoints.at(n._id));
    }

    // Sammeln der IDs aller Punkte im Cluster
    for (auto p : G.pointClusters[cluster]) {
        allPoints.push_back(p._id);
    }

    // Sortieren der Vektoren für symmetrische Differenz
    std::sort(allPoints.begin(), allPoints.end());
    std::sort(usedPoints.begin(), usedPoints.end());

    // Berechnung der freien Punkte
    std::set_difference(
        allPoints.begin(), allPoints.end(), 
        usedPoints.begin(), usedPoints.end(), 
        std::back_inserter(freePoints)
    );
    for (auto f : freePoints){
        std::cout << "Point " << f << std::endl;
    }
        // Print allPoints
    std::cout << "All points: ";
    for (int point : allPoints) {
        std::cout << point << " ";
    }
    std::cout << std::endl;

    // Print usedPoints
    std::cout << "Used points: ";
    for (int point : usedPoints) {
        std::cout << point << " ";
    }
    std::cout << std::endl;



    return freePoints;
}

std::vector<std::size_t> freePointsToMove(const Graph &G) {
    std::vector<std::size_t> freePoints;
    std::vector<std::size_t> usedPoints;
    std::vector<std::size_t> allPoints; 

    // Sammeln der IDs der genutzten Punkte
    for (auto n : G.nodes) { // G.NodeClusters[cluster]
        usedPoints.push_back(G.mapVerticesToPoints.at(n._id));
    }

    // Sammeln der IDs aller Punkte im Cluster
    for (auto p : G.points) {
        allPoints.push_back(p._id);
    }

    // Sortieren der Vektoren für symmetrische Differenz
    std::sort(allPoints.begin(), allPoints.end());
    std::sort(usedPoints.begin(), usedPoints.end());

    // Berechnung der freien Punkte
    std::set_difference(
        allPoints.begin(), allPoints.end(), 
        usedPoints.begin(), usedPoints.end(), 
        std::back_inserter(freePoints)
    );
    for (auto f : freePoints){
        std::cout << "Point " << f << std::endl;
    }
        // Print allPoints
    std::cout << "All points: ";
    for (int point : allPoints) {
        std::cout << point << " ";
    }
    std::cout << std::endl;

    // Print usedPoints
    std::cout << "Used points: ";
    for (int point : usedPoints) {
        std::cout << point << " ";
    }
    std::cout << std::endl;



    return freePoints;
}



void move(Graph &G, RTree<int, int, 2, float> &RTree,std::vector<size_t> &freePoints, std::vector<Node> &nodes,const int numberOfSamples, const int numNodes, Graph &C){
    auto randomSamples = generate_random_indices(numberOfSamples, nodes.size());
    // std::cout << "MOVE \n";
    assert(areValuesUnique(C.mapVerticesToPoints));
    for (auto i : randomSamples) {
        assert( i < nodes.size());
        auto &n = G.nodes[nodes[i].GetId()];
        removeFromRTree(RTree, C, n);
        int n_x = (int) n._x ;
        int n_y = (int) n._y;
        size_t n_id = n.GetId();
        int n_pointId;
        bool swap = false;
        
        unsigned int cr_before = 0;
        unsigned int cr_after = 0;
        auto adjEdges = C.getAdjEdges(n);
        for (auto const &a : adjEdges) {
            getEdgeCrossingsWithRTree(RTree, C, a, cr_before, numNodes);
        }
        auto randomPoints = generate_random_indices(freePoints.size(), freePoints.size());
        // std::cout << "2.\n";

        for (auto &p_idx : randomPoints) {
            
            auto &point = G.points[freePoints[p_idx]];
            cr_after = 0;

            assert(areValuesUnique(C.mapVerticesToPoints)); 

            // for (auto point : C.mapVerticesToPoints) {
            //     std::cout << point.first << "," << point.second << " # ";
            // }
            // std::cout << "\n";

            swapPositions(n, point);
            C.mapVerticesToPoints.at(n.GetId()) = point.GetId();
            // std::cout << n.GetId() << " " << n._pointId << " " << point.GetId() << "   p   \n";

            // for (auto point : C.mapVerticesToPoints) {
            //     std::cout << point.first << "," << point.second << " # ";
            // }
            // std::cout << "\n";


            assert(areValuesUnique(C.mapVerticesToPoints));
            // #pragma omp parallel for reduction(+:cr_after)
            for (int i = 0; i < adjEdges.size(); ++i) {
                getEdgeCrossingsWithRTree(RTree, C, adjEdges[i], cr_after, numNodes);
            }
            
            assert(areValuesUnique(C.mapVerticesToPoints));
            int tmp = 0;
            // std::cout << cr_before << " " << cr_after << "------------- " << G.ComputeCrossings(tmp) << " " << C.ComputeCrossings(tmp) << " " << G.ComputeEffectiveCrossings(tmp) << " " << C.ComputeEffectiveCrossings(tmp) << std::endl;
            
            if (cr_after <= cr_before) {
                // std::cout << cr_before << " " << cr_after << "-----------------------------" << std::endl;
                swap = true;
                cr_before = cr_after;
                n_x = n._x;
                n_y = n._y;
                n_id = n.GetId();
                n_pointId = point.GetId();
                
            }
            // Swap
            swapPositions(n, point);
            C.mapVerticesToPoints.at(n.GetId()) = n._pointId;

        }

        if (swap){
            // std::cout << C.mapVerticesToPoints.at(n_id) << " " << n_id << " " <<n_pointId << "   fffffff   \n";
            freePoints.push_back(C.mapVerticesToPoints.at(n_id)); 
            n._x = n_x;
            n._y = n_y;
            n._pointId = n_pointId;
            C.mapVerticesToPoints.at(n_id) = n_pointId;
            auto it = std::ranges::find(freePoints, n_pointId);
            freePoints.erase(it);
        }
        insertIntoRTree(RTree, C, n);
        // std::cout << "3\n";
    }
    assert(areValuesUnique(C.mapVerticesToPoints));
}

void switchNodes(Graph &G, RTree<int, int, 2, float> RTree, std::vector<Node> &nodesToSwitch, const int numberOfSamples, const int numNodes, Graph &C ){
    for (auto &n_id : nodesToSwitch) {
        Node &n = G.nodes[n_id.GetId()];
        auto randomSamples = generate_random_indices(numberOfSamples, nodesToSwitch.size());
        removeFromRTree(RTree, C, n);
        int n_x = (int) n.GetX();
        int n_y = (int) n.GetY();
        int n_pointID = n._pointId;
        int n_idx = n.GetId();
        int m_x = -1;
        int m_y = -1;
        int m_pointID;
        int m_idx ;
        bool swap = false;
        int tmp_x;
        unsigned int cr_before = 0;
        unsigned int cr_after = 0;
        std::vector<Edge> n_adjEdges = C.getAdjEdges(n);

        for (auto i : randomSamples) {
            cr_after = 0; 
            cr_before = 0;
            assert(i < nodesToSwitch.size());
            auto &m = G.nodes[nodesToSwitch[i].GetId()] ;
            if (n.GetId() == m.GetId() ) {
                break;
            }

            std::vector<Edge> m_adjEdges = C.getAdjEdges(m);

            #pragma omp parallel for reduction(+:cr_before)
            for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, C, n_adjEdges[idx], cr_before, numNodes);
            }

            insertIntoRTree(RTree, C, n);
            removeFromRTree(RTree, C, m);

            assert(areValuesUnique(C.mapVerticesToPoints));

            #pragma omp parallel for reduction(+:cr_before)
            for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, C, m_adjEdges[idx], cr_before, numNodes);
            }

            swapNodePositions(n, m);
            assert(m._pointId == C.mapVerticesToPoints.at(m._id));
            assert(n._pointId == C.mapVerticesToPoints.at(n._id));
            // m._pointId = C.mapVerticesToPoints.at(m._id);
            // n._pointId = C.mapVerticesToPoints.at(n._id);

            C.mapVerticesToPoints[n.GetId()]= m._pointId;
            C.mapVerticesToPoints[m.GetId()]= n._pointId;

            insertIntoRTree(RTree, C, m);
            removeFromRTree(RTree, C, n);


            #pragma omp parallel for reduction(+:cr_after)
            for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, C, n_adjEdges[idx], cr_after, numNodes);
            }

            insertIntoRTree(RTree, C, n);
            removeFromRTree(RTree, C, m);
            
            #pragma omp parallel for reduction(+:cr_after)
            for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, C, m_adjEdges[idx], cr_after, numNodes);
            }
            insertIntoRTree(RTree, C, m);
            removeFromRTree(RTree, C, n);

            if (cr_after <= cr_before) {
                swap = true;
                m_pointID = m._pointId;
                n_pointID = n._pointId;
                m_idx = m.GetId();
                n_idx = n.GetId();
                n_x = n._x ;
                n_y = n._y ; 
                m_x = m._x; 
                m_y = m._y;

                                    std::cout << "m_pointID = " << m._pointId << " " << C.mapVerticesToPoints[n.GetId()] << "\n";
                                    std::cout << "n_pointID = " << n._pointId << " " <<  C.mapVerticesToPoints[m.GetId()] << "\n";
                                    std::cout << "m_idx = " << m.GetId() << "\n";
                                    std::cout << "n_id = " << n.GetId() << "\n";
                                    std::cout << "n_x = " << n._x << "\n";
                                    std::cout << "n_y = " << n._y << "\n";
                                    std::cout << "m_x = " << m._x << "\n";
                                    std::cout << "m_y = " << m._y << "\n";


            }
            swapNodePositions(n,m);
            C.mapVerticesToPoints[n.GetId()]= n._pointId;
            C.mapVerticesToPoints[m.GetId()]= m._pointId;
        }
        if (swap) {
            auto& m = C.nodes[m_idx];
                                    std::cout << "FINal:" "\n";
                                    std::cout << "m_pointID = " << m._pointId << " " << C.mapVerticesToPoints[n.GetId()] << "\n";
                                    std::cout << "n_pointID = " << n._pointId << " " <<  C.mapVerticesToPoints[m.GetId()] << "\n";
                                    std::cout << "m_idx = " << m.GetId() <<  "\n";
                                    std::cout << "n_id = " << n.GetId() << "\n";
                                    std::cout << "n_x = " << n._x << "\n";
                                    std::cout << "n_y = " << n._y << "\n";
                                    std::cout << "m_x = " << m._x << "\n";
                                    std::cout << "m_y = " << m._y << "\n\n";
            
            C.mapVerticesToPoints[n.GetId()]= m_pointID;
            C.mapVerticesToPoints[m.GetId()]= n_pointID;
            m._pointId = n_pointID;
            n._pointId = m_pointID;
            n._x = n_x;
            n._y = n_y;
            m._x = m_x; 
            m._y = m_y;

            assert(m._pointId == C.mapVerticesToPoints.at(m._id));
            assert(n._pointId == C.mapVerticesToPoints.at(n._id));

        }
        
        insertIntoRTree(RTree, C, n);
    }
}
std::vector<Node> createCNodes(const Graph &G, const Graph &WCG,  size_t cluster){
    assert(cluster < WCG.numNodes );
    size_t id_Node_overhead = G.numNodes;

    size_t numNodes = WCG.numNodes + G.nodes.size() ; 
    std::vector<Node> nodes(numNodes);
    std::vector<Node> nodesG = G.NodeClusters[cluster];

    for (size_t i = 0; i < G.NodeClusters.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < G.NodeClusters[i].size(); ++j) {
            Node node = G.nodes[G.NodeClusters[i][j]._id];
            nodes[id_Node_overhead+i]._x += node._x;
            nodes[id_Node_overhead+i]._y += node._y;
        }
        nodes[id_Node_overhead+i]._x = nodes[id_Node_overhead+i]._x / G.NodeClusters[i].size();
        nodes[id_Node_overhead+i]._y = nodes[id_Node_overhead+i]._y / G.NodeClusters[i].size();
        nodes[id_Node_overhead+i]._id = id_Node_overhead + i;
        std::cout << nodes[id_Node_overhead+i]._id << " " << nodes[id_Node_overhead+i]._pointId << " " << nodes[id_Node_overhead+i]._x << " " << nodes[id_Node_overhead+i]._y <<"  \n";

    }

    for (auto n :  G.nodes){
        nodes[n._id] = n;
    }
    return nodes;
}

std::vector<Edge> createCEdges(const Graph &G,const  Graph &WCG, const size_t cluster){
    assert(cluster < WCG.numNodes );
    size_t id_Node_overhead = G.nodes.size();
    size_t id_Point_overhead = G.points.size();
    std::vector<Edge> edges = {}; 

    int idx = 0;
    for (auto e : G.edges){
        assert(e.first < G.nodes.size());
        assert(e.second < G.nodes.size());
        if (G.nodes[e.first]._cluster == cluster && G.nodes[e.second]._cluster == cluster ){
            // std::cout << "1\n";
            edges.push_back(Edge(e.first, e.second));
            edges.back().id = idx;
            idx++;
        }
        else if (G.nodes[e.first]._cluster == cluster ){
            // std::cout << "2\n";
            // assert(e.first != G.nodes[e.second]._cluster);
            assert(G.nodes[e.first]._cluster != G.nodes[e.second]._cluster);
            if (cluster < G.nodes[e.second]._cluster ){
                edges.push_back(Edge(e.first, G.nodes[e.second]._cluster + id_Node_overhead )) ;
            }
            else {
                edges.push_back(Edge(e.first, G.nodes[e.second]._cluster + id_Node_overhead )) ;
            }
            assert(e.second < G.numNodes);
            edges.back().id = idx;
            idx++;
        }
        else if (G.nodes[e.second]._cluster == cluster ){
            // std::cout << "3\n";
            // assert(e.second != G.nodes[e.first]._cluster);
            assert(G.nodes[e.first]._cluster != G.nodes[e.second]._cluster);
            if (cluster < G.nodes[e.first]._cluster ){
                edges.push_back(Edge(e.second, G.nodes[e.first]._cluster + id_Node_overhead )) ;
            }
            else {
                edges.push_back(Edge(e.second, G.nodes[e.first]._cluster + id_Node_overhead )) ;
            }
            edges.back().id = idx;
            idx++;
        }
    }

    std::vector<Point> points(G.points.size()+ WCG.points.size());
    for (auto e :  edges){
        // std::cout << e.first << " " << e.second << "\n";
    }

    return edges;
}

void assignCtoG(std::unordered_map<size_t, size_t> &G, const std::unordered_map<size_t, size_t> &C, const size_t max_elements){


    for (const auto& [key, value] : C) {
        if (key < max_elements) {
            G.at(key) = value;
        }
    }
}


void iterativeCrossMinLocal( Graph &G, const int numberOfOuterLoops, std::vector<size_t> &freePoints, const int numberOfSamples, const int numberOfOuterLoopsMove,
                                const int numNodes, Graph &WCG ){
    // initialize RTree data struct to speed up crossing computation
    areValuesUnique(G.mapVerticesToPoints);


    for (int cl = 0; cl < G.NodeClusters.size(); cl++) {
        std::cout  << "Nodes \n";
        auto nodes = createCNodes(G, WCG, cl);
        std::cout  << "Edges \n";
        auto edges = createCEdges(G, WCG, cl);
        std::cout  << "Points \n";
        std::vector<Point> points = G.points;
        for (auto c : G.NodeClusters){
            points.push_back(Point());
        };


        freePoints = freePointsToMove( G, cl);
        std::cout  << "Test\n";
        Graph C(nodes.size(), nodes, edges, points, G.width, G.height );
        for (int i = 0; i < G.nodes.size(); i++){
            C.mapVerticesToPoints.at(i) = G.mapVerticesToPoints.at(i);
        }
        for (int i = 0; i < C.nodes.size() - G.nodes.size(); i++){
            int nodeId = i + G.nodes.size();
            C.mapVerticesToPoints.at(nodeId) = G.points.size()+ i;
        }
        areValuesUnique(C.mapVerticesToPoints);


        RTree<int, int, 2, float> RTree;
        for (auto const &e : edges) {
            int px = (int) nodes[e.first].GetX();  
            int py = (int) nodes[e.first].GetY();
            int qx = (int) nodes[e.second].GetX();
            int qy = (int) nodes[e.second].GetY();
            Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
            RTree.Insert(rect.min, rect.max, e.id);
        }
        std::cout  << "b\n";
        for (int outerLoop = 0 ; outerLoop < numberOfOuterLoops; outerLoop++){
            std::cout << outerLoop << " of " << numberOfOuterLoops << std::endl;
            std::vector<Node> nodesToSwitch ;
            for (auto nTS : G.NodeClusters[cl]){
                Node &node = G.nodes[nTS._id];
                nodesToSwitch.push_back(node);
            }
            std::vector<Point> &PointsToSwitch = G.pointClusters[cl];
            if (nodesToSwitch.size() < PointsToSwitch.size()){
                for(int outerLoopMove = 0; outerLoopMove < numberOfOuterLoopsMove; outerLoopMove++){
        
                    move(G, RTree, freePoints, nodesToSwitch ,numberOfSamples, numNodes, C);

                }
            }
            assert(G.NodeClusters.size() == G.pointClusters.size());
            std::cout << G.NodeClusters.size() << " Switch \n";
            assert(areValuesUnique(C.mapVerticesToPoints));
            // switchNodes(G, RTree, nodesToSwitch, numberOfSamples, numNodes, C);
        }
        assignCtoG(G.mapVerticesToPoints, C.mapVerticesToPoints, G.nodes.size());
        int counter = 0;
        std::cout << G.ComputeCrossings(counter) << " cross \n";

    }
    
}



void iterativeCrossMinGlobal( Graph &G, const int numberOfOuterLoops, const int numberOfOuterLoopsMove, std::vector<size_t> &freePoints, const int numberOfSamples,
                                const int numNodes) {
    // initialize RTree data struct to speed up crossing computation
    RTree<int, int, 2, float> RTree;
    for (auto const &e : G.edges) {
        int px = (int) G.nodes[e.first].GetX();  
        int py = (int) G.nodes[e.first].GetY();
        int qx = (int) G.nodes[e.second].GetX();
        int qy = (int) G.nodes[e.second].GetY();
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, e.id);
    }

    freePoints = freePointsToMove(G);

    for (int outerLoop = 0 ; outerLoop < numberOfOuterLoops; outerLoop++){
        std::cout << outerLoop << " of " << numberOfOuterLoops << std::endl;
        std::vector<Node> &nodesToSwitch = G.nodes;
        std::vector<Point> &pointsToSwitch = G.points;

        if (pointsToSwitch.size() > nodesToSwitch.size()){ 
            for(int outerLoopMove = 0; outerLoopMove < numberOfOuterLoopsMove; outerLoopMove++){
                
                move(G, RTree, freePoints , G.nodes ,numberOfSamples, numNodes, G);

            }
        }
        assert(areValuesUnique(G.mapVerticesToPoints));

        switchNodes(G, RTree, nodesToSwitch, numberOfSamples, numNodes, G);
        assert(areValuesUnique(G.mapVerticesToPoints));
    }
}

