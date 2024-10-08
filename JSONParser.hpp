#pragma once
#include <string>
#include "json.hpp"
#include "Graph.hpp"



class JSONParser {
public:
    JSONParser();
    void jsonToGraph(size_t& numNodes,std::vector<Node>& nodes, std::vector<Edge>& edges, std::vector<Point>& points, int& width, int& height, const std::string& filepath);
    void graphToJson(Graph& G, const std::string& filepath, const std::string& fname, int width, int height);

private:
    nlohmann::json j;
};

