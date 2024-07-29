#include <iostream>

class Node {
// private:
//     int id;
//     int x;
//     int y;
//     int cluster;

public:
    // Constructor
    Node() = default;
        Node(int id, int x, int y, int cluster = -1) : _id(id), _x(x), _y(y) , _cluster(cluster) {}


    int getCluster() {
        return _cluster;
    }
    // Getters
    int getId() const {
        return _id;
    }

    int getX() const {
        return _x;
    }

    int getY() const {
        return _y;
    }

    void setCluster(int cluster){
        _cluster = cluster;
    }


    // Method to display the node's details
    void display() const {
        std::cout << "Node ID: " << _id << ", x: " << _x << ", y: " << _y << "cluster: " << _cluster << std::endl;
    }

    int _x{ 0 };
	int _y{ 0 };
	int _id{ 0 };
	int _cluster{0};
	
};