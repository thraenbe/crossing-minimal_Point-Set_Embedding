#include <iostream>
#include <vector>
#include <queue>
#include <cstdlib>
#include <ctime>
#include <algorithm>


using namespace std;

// BFS-based Graph Growing Algorithm
vector<int> graphGrowingPartition(const vector<vector<size_t>> &graph, int num_clusters) {
    int num_nodes = graph.size();
    int cluster_size = num_nodes / num_clusters;
    
    // To store the cluster membership of each node
    vector<int> membership(num_nodes, -1);

    // Track the size of each cluster
    vector<int> cluster_sizes(num_clusters, 0);

    // Randomly select seed nodes for each cluster
    srand(time(0));
    vector<int> seeds;
    for (int i = 0; i < num_clusters; ++i) {
        int seed;
        do {
            seed = rand() % num_nodes;
        } while (membership[seed] != -1); // Ensure unique seed
        seeds.push_back(seed);
        membership[seed] = i;
        cluster_sizes[i] += 1;
    }

    // BFS queues to grow each cluster
    queue<pair<int, int>> bfs_queue;
    for (int i = 0; i < num_clusters; ++i) {
        bfs_queue.push({seeds[i], i});
    }

    // Grow clusters until all nodes are assigned
    while (!bfs_queue.empty()) {
        auto [node, cluster_id] = bfs_queue.front();
        bfs_queue.pop();

        for (int neighbor : graph[node]) {
            if (membership[neighbor] == -1 && cluster_sizes[cluster_id] < cluster_size) {
                membership[neighbor] = cluster_id;
                cluster_sizes[cluster_id] += 1;
                bfs_queue.push({neighbor, cluster_id});
            }
        }
    }

    // If some nodes are unassigned, assign them randomly
    for (int i = 0; i < num_nodes; ++i) {
        if (membership[i] == -1) {
            int smallest_cluster = std::min_element(cluster_sizes.begin(), cluster_sizes.end()) - cluster_sizes.begin();
            membership[i] = smallest_cluster;
            cluster_sizes[smallest_cluster] += 1;
        }
    }

    return membership;
}
