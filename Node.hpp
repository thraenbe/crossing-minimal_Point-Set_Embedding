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
        Node(int id, double x, double y, int cluster = -1) : _id(id), _x(x), _y(y) , _cluster(cluster) {}


    int getCluster() {
        return _cluster;
    }
    // Getters
    int getId() const {
        return _id;
    }

    double getX() const {
        return _x;
    }

    double getY() const {
        return _y;
    }

    void setCluster(int cluster){
        _cluster = cluster;
    }

    Point nodeToPoint() const {
        return Point(0, _x, _y);
    }


    // Method to display the node's details
    void display() const {
        std::cout << "Node ID: " << _id << ", x: " << _x << ", y: " << _y << "cluster: " << _cluster << std::endl;
    }

    double _x{ 0 };
	double _y{ 0 };
	int _id{ 0 };
	int _cluster{0};
	
};