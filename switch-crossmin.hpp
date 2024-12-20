#include "rtree.h"
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/geometry.h> // For DPoint and DSegment
#include <omp.h>

//using namespace ogdf;

int doIntersect(const ogdf::edge &e1, const ogdf::edge &e2, const ogdf::GraphAttributes &GA, const int numNodes);
int orientation(ogdf::node const &p, ogdf::node const &q, ogdf::node const &r, const ogdf::GraphAttributes &GA);
int onSegment(ogdf::node const &p, ogdf::node const &q, ogdf::node const &r, const ogdf::GraphAttributes &GA);

void getEdgeCrossingsWithRTree(const RTree<int, int, 2, float> &RTree, const ogdf::Graph &G,
                                                       const ogdf::GraphAttributes &GA, const ogdf::edge &e,
                                                       unsigned int &cr, 
                                                       const int numNodes) {
    int px = (int) GA.x(e->source());
    int py = (int) GA.y(e->source());
    int qx = (int) GA.x(e->target());
    int qy = (int) GA.y(e->target());
    

   Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));

    ogdf::Array<ogdf::edge, int> edges;
    G.allEdges(edges);
    std::vector<int> treeRes;
    int hits = RTree.Search(rect.min, rect.max, &treeRes);

    for (int i = 0; i < hits; ++i) {
        cr = cr + doIntersect(e, edges[treeRes[i]], GA, numNodes);
    }
}
int doIntersect(const ogdf::edge &e1, const ogdf::edge &e2,
                                        const ogdf::GraphAttributes &GA, const int numNodes) {
    int o1 = orientation(e1->source(), e1->target(), e2->source(), GA);
    int o2 = orientation(e1->source(), e1->target(), e2->target(), GA);
    int o3 = orientation(e2->source(), e2->target(), e1->source(), GA);
    int o4 = orientation(e2->source(), e2->target(), e1->target(), GA);
    if (o1 != o2 && o3 != o4) return 1;
    if (o1 == 0 && onSegment(e1->source(), e2->source(), e1->target(), GA)) return numNodes;
    if (o2 == 0 && onSegment(e1->source(), e2->target(), e1->target(), GA)) return numNodes;
    if (o3 == 0 && onSegment(e2->source(), e1->source(), e2->target(), GA)) return numNodes;
    if (o4 == 0 && onSegment(e2->source(), e1->target(), e2->target(), GA)) return numNodes;
    return 0;
}
int onSegment(ogdf::node const &p, ogdf::node const &q, ogdf::node const &r,
                                      const ogdf::GraphAttributes &GA) {
    return ((int) GA.x(q) <= std::max((int) GA.x(p), (int) GA.x(r)) &&
            (int) GA.x(q) >= std::min((int) GA.x(p), (int) GA.x(r)) &&
            (int) GA.y(q) <= std::max((int) GA.y(p), (int) GA.y(r)) &&
            (int) GA.y(q) >= std::min((int) GA.y(p), (int) GA.y(r)));
}
int orientation(ogdf::node const &p, ogdf::node const &q, ogdf::node const &r,
                                        const ogdf::GraphAttributes &GA) {
    long long val = ((long long) ((int) GA.y(q) - (int) GA.y(p))) * ((long long) ((int) GA.x(r) - (int) GA.x(q))) -
                    ((long long) ((int) GA.x(q) - (int) GA.x(p))) * ((long long) ((int) GA.y(r) - (int) GA.y(q)));
    if (val == 0) return 0;
    return (val > 0) ? 1 : 2;
}



void removeFromRTree(RTree<int, int, 2, float> &RTree, const ogdf::Graph &G, const ogdf::GraphAttributes &GA,
                                const ogdf::node &n) {
    for (auto const &a : n->adjEntries) {
        int px = (int) GA.x(a->theEdge()->source());
        int py = (int) GA.y(a->theEdge()->source());
        int qx = (int) GA.x(a->theEdge()->target());
        int qy = (int) GA.y(a->theEdge()->target());
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Remove(rect.min, rect.max, a->theEdge()->index());
    }
}

void insertIntoRTree(RTree<int, int, 2, float> &RTree, const ogdf::Graph &G, const ogdf::GraphAttributes &GA,
                                const ogdf::node &n) {
    for (auto const &a : n->adjEntries) {
        int px = (int) GA.x(a->theEdge()->source());
        int py = (int) GA.y(a->theEdge()->source());
        int qx = (int) GA.x(a->theEdge()->target());
        int qy = (int) GA.y(a->theEdge()->target());
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, a->theEdge()->index());
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



void iterativeCrossMinSwitch(const ogdf::Graph &G, ogdf::GraphAttributes &GA, const int numberOfOuterLoops, const int numberOfSamples,
                                const ogdf::List<ogdf::node> nodesToSwitch, 
                                const int numNodes) {
    // initialize RTree data struct to speed up crossing computation
    RTree<int, int, 2, float> RTree;
    for (auto const &e : G.edges) {
        int px = (int) GA.x(e->source());
        int py = (int) GA.y(e->source());
        int qx = (int) GA.x(e->target());
        int qy = (int) GA.y(e->target());
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, e->index());
    }
    // std::random_device rd;
    // std::mt19937 gen(rd());
    //
    // Outer loop
    //
    ogdf::Array<ogdf::node, int> nodes;
    G.allNodes(nodes);
    
    for (int l = 0; l < numberOfOuterLoops; l++) {
        for (auto const &n : nodesToSwitch) {
            auto randomSamples = generate_random_indices(numberOfSamples, nodesToSwitch.size());
            removeFromRTree(RTree, G, GA, n);
            int n_x = (int)GA.x(n);
            int n_y = (int)GA.y(n);
            int n_pointID = std::stoi(GA.label(n));
            int m_x = -1;
            int m_y = -1;
            int m_pointID;
            ogdf::node m_node;
            unsigned int cr_before = 0;
            unsigned int cr_after = 0;
            ogdf::List<ogdf::edge> n_adjEdges;
            ogdf::List<ogdf::edge> m_adjEdges;
            n->adjEdges(n_adjEdges);
            #pragma omp parallel for reduction(+:cr_before)
            for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                getEdgeCrossingsWithRTree(RTree, G, GA, *n_adjEdges.get(idx), cr_before, numNodes);
            }

            for (int i = 0; i < randomSamples.size(); ++i) {
                const auto &m = *nodesToSwitch.get(i) ;
                if (n == m ) {
                    break;
                }
                m->adjEdges(m_adjEdges);

                removeFromRTree(RTree, G, GA, m);
                insertIntoRTree(RTree, G, GA, n);

                #pragma omp parallel for reduction(+:cr_before)
                for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                    getEdgeCrossingsWithRTree(RTree, G, GA, *m_adjEdges.get(idx), cr_before, numNodes);
                }

                insertIntoRTree(RTree, G, GA, m);
                removeFromRTree(RTree, G, GA, n);

                
                auto tmp_x = GA.x(n);
                auto tmp_y = GA.y(n);
                GA.x(n) = GA.x(m);
                GA.y(n) = GA.y(m);
                GA.x(m) = tmp_x;
                GA.y(m) = tmp_y;
                #pragma omp parallel for reduction(+:cr_after)
                for (int idx = 0; idx < n_adjEdges.size(); ++idx) {
                    getEdgeCrossingsWithRTree(RTree, G, GA, *n_adjEdges.get(idx), cr_after, numNodes);
                }
                removeFromRTree(RTree, G, GA, m);
                insertIntoRTree(RTree, G, GA, n);
                #pragma omp parallel for reduction(+:cr_after)
                for (int idx = 0; idx < m_adjEdges.size(); ++idx) {
                    getEdgeCrossingsWithRTree(RTree, G, GA, *m_adjEdges.get(idx), cr_after, numNodes);
                }
                insertIntoRTree(RTree, G, GA, m);
                removeFromRTree(RTree, G, GA, n);

                if (cr_after <= cr_before) {
                    cr_before = cr_after;
                    n_x = (int)GA.x(n);
                    n_y = (int)GA.y(n);
                    n_pointID = std::stoi(GA.label(m));
                    m_x = (int)GA.x(m);
                    m_y = (int)GA.y(m);
                    m_pointID = std::stoi(GA.label(n));
                    m_node = m;
                    assert(m_x >= 0 && m_y >= 0);
                }else{ // check better results with or without else
                tmp_x = GA.x(n);
                tmp_y = GA.y(n);
                GA.x(n) = GA.x(m);
                GA.y(n) = GA.y(m);
                GA.x(m) = tmp_x;
                GA.y(m) = tmp_y;

                }
            }
            if (m_x >= 0 && m_y >= 0){
                GA.x(n) = n_x;
                GA.y(n) = n_y;
                GA.label(n) = std::to_string(n_pointID);
                GA.x(m_node) = m_x;
                GA.y(m_node) = m_y;
                GA.label(m_node) = std::to_string(m_pointID);
            }
            insertIntoRTree(RTree, G, GA, n);
            
        }
    }
}



// Ensure that nodes lie in range [0,xdim] x [0,ydim] -> will only sample positions in this box
void iterativeCrossMinMove(const ogdf::Graph &G, ogdf::GraphAttributes &GA, ogdf::GraphAttributes &secondGA, const int numberOfOuterLoops, const int numberOfSamples, ogdf::List<ogdf::node> positions, const int numNodes) {
    // initialize RTree data struct to speed up crossing computation
    RTree<int, int, 2, float> RTree;
    for (auto const &e : G.edges) {
        int px = (int) GA.x(e->source());
        int py = (int) GA.y(e->source());
        int qx = (int) GA.x(e->target());
        int qy = (int) GA.y(e->target());
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, e->index());
    }
    
    //
    // Outer loop
    //
    ogdf::Array<ogdf::node, int> nodes;
    G.allNodes(nodes);
    for (int l = 0; l < numberOfOuterLoops; l++) {
        for (auto const &n : G.nodes) {
            auto randomSamples = generate_random_indices(numberOfSamples, positions.size());
            removeFromRTree(RTree, G, GA, n);
            int n_x = (int)GA.x(n);
            int n_y = (int)GA.y(n);
            int n_pointId;
            int m_x = -1;
            int m_y = -1;
            int m_pointId;
            ogdf::node& point = positions.front() ;
            if (isNumeric(GA.label(n)) ){
                int pointId =  std::stoi(GA.label(n));
                } 
                else { 
                    break;
                }
            unsigned int cr_before = 0;
            unsigned int cr_after = 0;
            for (auto const &a : n->adjEntries) {
                getEdgeCrossingsWithRTree(RTree, G, GA, a->theEdge(), cr_before, numNodes);
            }
            ogdf::List<ogdf::edge> adjEdges;
            n->adjEdges(adjEdges);
            for (int idx : randomSamples) {
                assert(idx < positions.size());
                point = *positions.get(idx);

                // Swap 
                std::swap(GA.x(n), secondGA.x(point));
                std::swap(GA.y(n), secondGA.y(point));

                
                #pragma omp parallel for reduction(+:cr_after)
                for (int i = 0; i < adjEdges.size(); ++i) {
                    getEdgeCrossingsWithRTree(RTree, G, GA, *adjEdges.get(i), cr_after, numNodes);
                }
                if (cr_after <= cr_before) {
                    
                    cr_before = cr_after;
                    n_x = (int)GA.x(n);
                    n_y = (int)GA.y(n);

                    if (isNumeric(GA.label(n)) && isNumeric(secondGA.label(point))){
                        n_pointId =  std::stoi(secondGA.label(point));
                    } 
                    else { 
                        break;
                    }

                    m_x = (int) secondGA.x(point);
                    m_y = (int) secondGA.y(point);

                    if (isNumeric(GA.label(n)) && isNumeric(secondGA.label(point))){
                        m_pointId =  std::stoi(GA.label(n));
                    } 
                    else { 
                        break;
                    }
                }
                // Swap
                std::swap(GA.x(n), secondGA.x(point));
                std::swap(GA.y(n), secondGA.y(point));
            }

            if (m_x >= 0 && m_y >= 0){
                GA.x(n) = n_x;
                GA.y(n) = n_y;
                GA.label(n) = std::to_string(n_pointId);
                secondGA.x(point) = m_x;
                secondGA.y(point) = m_y;
                secondGA.label(point) = std::to_string(m_pointId);
            }
            insertIntoRTree(RTree, G, GA, n);
        }
    }
}



void iterativeCrossMin(const ogdf::Graph &G, ogdf::GraphAttributes &GA, const int xdim, const int ydim, const int numberOfOuterLoops, const int numberOfPositionSamples, const int numNodes) {
    // initialize RTree data struct to speed up crossing computation
    RTree<int, int, 2, float> RTree;
    for (auto const &e : G.edges) {
        int px = (int) GA.x(e->source());
        int py = (int) GA.y(e->source());
        int qx = (int) GA.x(e->target());
        int qy = (int) GA.y(e->target());
        Rect rect(std::min(px, qx), std::min(py, qy), std::max(px, qx), std::max(py, qy));
        RTree.Insert(rect.min, rect.max, e->index());
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    //
    // Outer loop
    //
    ogdf::Array<ogdf::node, int> nodes;
    G.allNodes(nodes);
    std::uniform_int_distribution<> noderange(0, G.numberOfNodes() - 1);
    // we only sample integer points here
    std::uniform_int_distribution<> xrange(0, xdim);
    std::uniform_int_distribution<> yrange(0, ydim);
    for (int l = 0; l < numberOfOuterLoops; l++) {
        for (auto const &n : G.nodes) {
            removeFromRTree(RTree, G, GA, n);
            int x = (int)GA.x(n);
            int y = (int)GA.y(n);
            unsigned int cr_before = 0;
            unsigned int cr_after = 0;
            for (auto const &a : n->adjEntries) {
                getEdgeCrossingsWithRTree(RTree, G, GA, a->theEdge(), cr_before, numNodes);
            }
            ogdf::List<ogdf::edge> adjEdges;
            n->adjEdges(adjEdges);
            for (int i = 0; i < numberOfPositionSamples; ++i) {
                GA.x(n) = xrange(gen);
                GA.y(n) = yrange(gen);
                #pragma omp parallel for reduction(+:cr_after)
                for (int i = 0; i < adjEdges.size(); ++i) {
                    getEdgeCrossingsWithRTree(RTree, G, GA, *adjEdges.get(i), cr_after, numNodes);
                }
                if (cr_after <= cr_before) {
                    cr_before = cr_after;
                    x = (int)GA.x(n);
                    y = (int)GA.y(n);
                }
            }
            GA.x(n) = x;
            GA.y(n) = y;
            insertIntoRTree(RTree, G, GA, n);
        }
    }
}
