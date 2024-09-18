#pragma once
#include <iostream>


class Node {
public:
    // Constructor
    Node() = default;
        Node(const size_t id, const double x, const double y,const size_t cluster) : _id(id), _x(x), _y(y) , _cluster(cluster) {}
        Node(const size_t id, const double x, const double y) : _id(id), _x(x), _y(y){}

    // Getters
    [[nodiscard]] size_t GetId() const {
        return _id;
    }

    [[nodiscard]] double GetX() const {
        return _x;
    }

    [[nodiscard]] double GetY() const {
        return _y;
    }

    [[nodiscard]] size_t GetCluster()const  {
            return _cluster;
        }

    void SetCluster(const size_t cluster){
        _cluster = cluster;
    }

    [[nodiscard]] Point NodeToPoint() const {
        return{0, _x, _y};
    }


    // Method to display the node's details
    void Display() const {
        std::cout << "Node ID: " << _id << ", x: " << _x << ", y: " << _y << "cluster: " << _cluster << std::endl;
    }

    size_t _id{ 0 };
    double _x{ 0 };
	double _y{ 0 };
	size_t _cluster{0};
    size_t _pointId{0};
	
};