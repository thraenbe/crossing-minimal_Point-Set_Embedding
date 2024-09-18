#pragma once
#include <fstream>
#include "JSONParser.hpp"
#include <list>
#include <vector>
#include <iostream>

using namespace nlohmann;

JSONParser::JSONParser() = default;


void JSONParser::jsonToGraph(size_t& numNodes,std::vector<Node>& nodes ,std::vector<Edge>& edges, std::vector<Point>& points, int& width, int& height, const std::string& filepath) {
    std::ifstream is(filepath+".json");
    is >> j;
    for (auto const& p : j["points"]) {
        points.emplace_back(p["id"],p["x"],p["y"]);
    }
    int idx = 0;
    for (auto const& e : j["edges"]) {
        edges.push_back({ e["source"],e["target"],idx });
        idx++;
    }
    for (auto const& n : j["nodes"]){
        nodes.emplace_back(n["id"],n["x"],n["y"]);
    }
    numNodes = j["nodes"].size();
    width = j["width"];
    height = j["height"];
}

void JSONParser::graphToJson(Graph& G, const std::string& filepath, const std::string& fname, int width, int height) {
    std::cout << "\nOutputting JSON to " << filepath << std::endl;

    j.clear();
    
    
    for (size_t i = 0; i < G.numNodes; i++){
        const auto& node = G.nodes[i]; 
        const auto& pos = G.GetConstPosOfNode(i);
        j["nodes"].push_back({ {"id", i},
                              {"x", pos.GetX()},
                              {"y", pos.GetY()}});
    }

    for (auto const& e : G.edges) {
        j["edges"].push_back({ {"source", e.first},
                              {"target", e.second} });
    }

    std::vector<std::tuple<int, int>>::iterator it;
    int i = 0;
    for (size_t i = 0; i < G.numPoints;++i) {
        const auto& point = G.points[i];
        j["points"].push_back({ {"id",i},
        {"x", point.GetX()},
        {"y", point.GetY()}});
    }
    j["width"] = width;
    j["height"] = height;


    // Write file
    std::cout << filepath + "/" + fname + ".json" << std::endl;
    std::ofstream os(filepath + "/" + fname + ".json");
    os << j;
}
