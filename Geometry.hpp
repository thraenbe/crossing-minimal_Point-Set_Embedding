#pragma once
#include <cmath>
#include <map>
#include <vector>

// Required geometry for this contest -> points are on integer coordinates
using size_t = std::size_t;
using namespace std;


class Point {

public:
	Point() = default;
	    Point(int id, double x, double y, int cluster = 0) :  _id(id), _x(x), _y(y), _cluster(cluster) {}

    // SetPosition methods
    void SetPosition(const std::pair<double, double>& newPos) noexcept {
        _x = newPos.first;
        _y = newPos.second;
    }
    void SetPosition(double x, double y) noexcept {
        _x = x;
        _y = y;
	}

    // SetCluster method
    void SetCluster(int cluster) noexcept {
        _cluster = cluster;
    }

    // GetCluster method
    int GetCluster() const noexcept {
        return _cluster;
    }
	Point operator+ (const Point& other) const
	{
		return Point(_x + other._x, _y + other._y , _id);
	}

	Point operator- (const Point& other) const
	{
		return Point(_x - other._x, _y - other._y, _id);
	}

	void operator+=(const Point& other) {
		_x += other._x;
		_y += other._y;
	}

	void operator*=(const int scalar) {
		_x *= scalar;
		_y *= scalar;
	}

	bool operator== (const Point& other) const {
		return _x == other._x && _y == other._y;
	}
	bool operator<( const Point& other) const {
		if (_x != other._x) {
			return _x < other._x;
		}
		return _y < other._y;
	}

	double GetX() const {
		return _x;
	}
	double GetY() const {
		return _y;
	}

	int GetId() const  {
		return _id;
	}
	
	double _x{ 0 };
	double _y{ 0 };
	int _id{ 0 };
	int _cluster{0};
	};

using Segment = std::pair<Point, Point>;

	// Products
[[nodiscard]] inline double DotProd(const Point& lhs, const Point& rhs) {
	return lhs._x * rhs._x + lhs._y * rhs._y;
}
[[nodiscard]] inline double CrossProd(const Point& lhs, const Point& rhs) {
	return lhs._x * rhs._y - lhs._y * rhs._x;
}
	// Distance.
[[nodiscard]] inline int L1Dist(const Point& lhs, const Point& rhs) {
	return std::abs(lhs._x - rhs._x) + std::abs(lhs._y + rhs._y);
}
[[nodiscard]] inline double L2DistSquared(const Point& lhs, const Point& rhs) {
	return DotProd(lhs - rhs, lhs - rhs);
}

[[nodiscard]] inline double L2Dist(const Point& lhs, const Point& rhs) {
	return std::sqrt(L2DistSquared(lhs, rhs));
}

// Geometric functions
// check if query lies in range [a,b]
[[nodiscard]] inline bool InRange(const int left, const int query, const int right) {
	return query >= left && query <= right;
}

[[nodiscard]] inline bool InInteriorOfRange(const int left, const int query, const int right) {
	return query > left && query < right;
}

[[nodiscard]] inline bool InRangeUnordered(const int a, const int query, const int b) {
	if (a < b) {
		return InRange(a, query, b);
	}
	return InRange(b, query, a);
}

[[nodiscard]] inline bool InInteriorOfRangeUnordered(const int a, const int query, const int b) {
	if (a < b) {
		return InRangeUnordered(a, query, b);
	}
	return InRangeUnordered(b, query, a);
}

// Check if point q lies in bounding box defined by p and q
[[nodiscard]] inline bool InBoundingBox(const Point& p, const Point& query, const Point& q)
{
	return InRangeUnordered(p._x, query._x, q._x) && InRangeUnordered(p._y, query._y, q._y);
}

[[nodiscard]] inline bool PointInInteriorOfSegment(const Point& p, const Point& query, const Point& q)
{
	return InInteriorOfRangeUnordered(p._x, query._x, q._x) && InInteriorOfRangeUnordered(p._y, query._y, q._y);
}


// Returns orientation for tripple of points
[[nodiscard]] inline int Orientation(const Point& p1, const Point& p2, const Point& p3)
{
	int val = (p2._y - p1._y) * (p3._x - p2._x)
		- (p2._x - p1._x) * (p3._y - p2._y);
	if (val == 0)
		return 0;
	return (val > 0) ? 1 : 2;
}

// The main function that returns true if line segment 'p1q1'
// and 'p2q2' intersect.
[[nodiscard]] inline int DoIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2, const int numNodes)
{
	// Find the four orientations needed for general and
	// special cases
	int o1 = Orientation(p1, q1, p2);
	int o2 = Orientation(p1, q1, q2);
	int o3 = Orientation(p2, q2, p1);
	int o4 = Orientation(p2, q2, q1);

	// General case
	if (o1 != o2 && o3 != o4)
		return 1;

	// Special Cases
	// p1, q1 and p2 are collinear and p2 lies on segment p1q1
	if (o1 == 0 && InBoundingBox(p1, p2, q1)) return numNodes;

	// p1, q1 and q2 are collinear and q2 lies on segment p1q1
	if (o2 == 0 && InBoundingBox(p1, q2, q1)) return numNodes;

	// p2, q2 and p1 are collinear and p1 lies on segment p2q2
	if (o3 == 0 && InBoundingBox(p2, p1, q2)) return numNodes;

	// p2, q2 and q1 are collinear and q1 lies on segment p2q2
	if (o4 == 0 && InBoundingBox(p2, q1, q2)) return numNodes;

	return 0; // Doesn't fall in any of the above cases
}
[[nodiscard]] inline bool DoIntersect(const Segment& lhs, const Segment& rhs, int numNodes) {
	return DoIntersect(lhs.first, lhs.second, rhs.first, rhs.second, numNodes);
}
[[nodiscard]] inline int sqr(int number) {
    return number * number;
}

[[nodiscard]] inline float euclideanDistance(const Point& p, const Point& q ){
	return sqrt(sqr(p._x - q._x) + sqr(p._y - q._y));
}

[[nodiscard]] inline float computeGradient(const Point& p, const Point& q ){
	double gradient;
	if (p._x != q._x){ 
		gradient =  (double)(p._y - q._y) / (double)(p._x - q._x);
	}
	else{
		gradient = std::numeric_limits<double>::infinity();
	}
	return gradient;
}

