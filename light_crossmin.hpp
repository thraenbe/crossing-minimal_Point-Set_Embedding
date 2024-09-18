#include "rtree.h"
#include "Graph.hpp"
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/geometry.h> // For DPoint and DSegment
#include <omp.h>
#include <ranges>


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
        if(e.id != G.edges[treeRes[i]].id){
        cr = cr + G.DoEdgesIntersect(e, G.edges[treeRes[i]]);
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

void move(Graph &G, RTree<int, int, 2, float> &RTree,std::vector<size_t> &freePoints, std::vector<Node> &nodes,const int numberOfSamples, const int numNodes ){
    auto randomSamples = generate_random_indices(numberOfSamples, nodes.size());
    for (auto i : randomSamples) {
        assert( i < nodes.size());
        auto &n = G.nodes[nodes[i].GetId()];
        removeFromRTree(RTree, G, n);
        int n_x = (int) n._x ;
        int n_y = (int) n._y;
        int n_pointId;
        bool swap = false;
        
        unsigned int cr_before = 0;
        unsigned int cr_after = 0;
        auto adjEdges = G.getAdjEdges(n);
        for (auto const &a : adjEdges) {
            getEdgeCrossingsWithRTree(RTree, G, a, cr_before, numNodes);
        }
        auto randomPoints = generate_random_indices(freePoints.size(), freePoints.size());
        for (auto &p_idx : randomPoints) {
            auto &point = G.points[freePoints[p_idx]];
            cr_after = 0;

            // Swap 
            swapPositions(n, point);
            G.mapVerticesToPoints[n.GetId()] = point.GetId();

            #pragma omp parallel for reduction(+:cr_after)
            for (int i = 0; i < adjEdges.size(); ++i) {
                getEdgeCrossingsWithRTree(RTree, G, adjEdges[i], cr_after, numNodes);
            }

            
            if (cr_after <= cr_before) {
                swap = true;
                cr_before = cr_after;
                n_x = n._x;
                n_y = n._y;
                n_pointId = point.GetId();
            }
            // Swap
            swapPositions(n, point);
            G.mapVerticesToPoints[n.GetId()] = n._pointId;

        }

        if (swap){
            freePoints.push_back(n._pointId); 
            n._x = n_x;
            n._y = n_y;
            n._pointId = n_pointId;
            G.mapVerticesToPoints[n.GetId()] = n_pointId;
            auto it = std::ranges::find(freePoints, n_pointId);
            freePoints.erase(it);
        }
        insertIntoRTree(RTree, G, n);
    }
}

void switchNodes(Graph &G, RTree<int, int, 2, float> RTree, std::vector<Node> &nodesToSwitch, const int numberOfSamples, const int numNodes ){
    for (auto &n_id : nodesToSwitch) {
        Node &n = G.nodes[n_id.GetId()];
        auto randomSamples = generate_random_indices(numberOfSamples, nodesToSwitch.size());
        removeFromRTree(RTree, G, n);
        int n_x = (int) n.GetX();
        int n_y = (int) n.GetY();
        int n_pointID = n._pointId;
        int m_x = -1;
        int m_y = -1;
        int m_pointID;
        int m_idx ;
        bool swap = false;
        int tmp_x;
        unsigned int cr_before = 0;
        unsigned int cr_after = 0;
        std::vector<Edge> n_adjEdges = G.getAdjEdges(n);


        for (auto i : randomSamples) {
            cr_after = 0; 
            cr_before = 0;
            assert(i < nodesToSwitch.size());
            auto &m = G.nodes[nodesToSwitch[i].GetId()] ;
            if (n.GetId() == m.GetId() ) {
                break;
            }

            std::vector<Edge> m_adjEdges = G.getAdjEdges(m);

            #pragma omp parallel for reduction(+:cr_before)
            for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, G, n_adjEdges[idx], cr_before, numNodes);
            }

            insertIntoRTree(RTree, G, n);
            removeFromRTree(RTree, G, m);
            

            #pragma omp parallel for reduction(+:cr_before)
            for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, G, m_adjEdges[idx], cr_before, numNodes);
            }

            swapNodePositions(n, m);
            G.mapVerticesToPoints[n.GetId()]= m._pointId;
            G.mapVerticesToPoints[m.GetId()]= n._pointId;

            insertIntoRTree(RTree, G, m);
            removeFromRTree(RTree, G, n);


            #pragma omp parallel for reduction(+:cr_after)
            for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, G, n_adjEdges[idx], cr_after, numNodes);
            }

            insertIntoRTree(RTree, G, n);
            removeFromRTree(RTree, G, m);
            
            #pragma omp parallel for reduction(+:cr_after)
            for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, G, m_adjEdges[idx], cr_after, numNodes);
            }
            insertIntoRTree(RTree, G, m);
            removeFromRTree(RTree, G, n);

            

            if (cr_after <= cr_before) {
                swap = true;
                m_pointID = n._pointId;
                n_pointID = m._pointId;
                m_idx = m.GetId();
            }
            swapNodePositions(n,m);
            G.mapVerticesToPoints[n.GetId()]= n._pointId;
            G.mapVerticesToPoints[m.GetId()]= m._pointId;
            
        }
        if (swap) {
            auto& m = G.nodes[m_idx];
            G.mapVerticesToPoints[n.GetId()]= n_pointID;
            G.mapVerticesToPoints[m.GetId()]= m_pointID;
            m._pointId = m_pointID;
            n._pointId = n_pointID;

            swapNodePositions(m,n);
        }
        
        insertIntoRTree(RTree, G, n);
    }
}


void iterativeCrossMinLocal( Graph &G, const int numberOfOuterLoops, std::vector<size_t> &freePoints, const int numberOfSamples, const int numberOfOuterLoopsMove,
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

    for (int outerLoop = 0 ; outerLoop < numberOfOuterLoops; outerLoop++){
        // std::cout << outerLoop << " of " << numberOfOuterLoops << std::endl;
        for (int l = 0; l < G.NodeClusters.size(); l++) {
            std::vector<Node> nodesToSwitch ;
            for (auto nTS : G.NodeClusters[l]){
                Node &node = G.nodes[nTS._id];
                nodesToSwitch.push_back(node);
            }
            std::vector<Point> &PointsToSwitch = G.pointClusters[l];
            if (nodesToSwitch.size() < PointsToSwitch.size()){
                for(int outerLoopMove = 0; outerLoopMove < numberOfOuterLoopsMove; outerLoopMove++){
        
                    move(G, RTree, freePoints, nodesToSwitch ,numberOfSamples, numNodes);

                }
            }
            assert(G.NodeClusters.size() == G.pointClusters.size());

            switchNodes(G, RTree, nodesToSwitch, numberOfSamples, numNodes);
        }
    }
}



void iterativeCrossMinGlobal( Graph &G, const int numberOfOuterLoops, const int numberOfOuterLoopsMove, std::vector<size_t> &freePoints, const int numberOfSamples,
                                const int numNodes, int crossings) {
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


    for (int outerLoop = 0 ; outerLoop < numberOfOuterLoops; outerLoop++){
        // std::cout << outerLoop << " of " << numberOfOuterLoops << std::endl;
        std::vector<Node> &nodesToSwitch = G.nodes;
        std::vector<Point> &pointsToSwitch = G.points;

        if (pointsToSwitch.size() > nodesToSwitch.size()){ 
            for(int outerLoopMove = 0; outerLoopMove < numberOfOuterLoopsMove; outerLoopMove++){
                
               move(G, RTree, freePoints , G.nodes ,numberOfSamples, numNodes);

            }
        }

        //switchNodes(G, RTree, nodesToSwitch, numberOfSamples, numNodes);
        
    }
}








