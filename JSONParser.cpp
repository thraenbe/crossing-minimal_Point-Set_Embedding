#pragma once
#include <fstream>
#include "JSONParser.hpp"
#include <list>
#include <vector>
#include <iostream>

using namespace nlohmann;

JSONParser::JSONParser() = default;


void JSONParser::jsonToGraph(size_t& numNodes ,std::vector<std::pair<size_t,size_t>>& edges, std::vector<Point>& points, int& width, int& height, const std::string& filepath) {
    std::ifstream is(filepath+".json");
    is >> j;
    for (auto const& p : j["points"]) {
        points.push_back({p["x"], p["y"]});
    }
    for (auto const& e : j["edges"]) {
        edges.push_back({ e["source"],e["target"] });
    }
    numNodes = j["nodes"].size();
    width = j["width"];
    height = j["height"];
}

void JSONParser::graphToJson(Graph& G, const std::string& filepath, const std::string& fname, int width, int height) {
    std::cout << "\nOutputting JSON to " << filepath << std::endl;

    j.clear();
    //push node information at the back of the "nodes" list
    for (auto const& [nodeId,pointId] : G.mapVerticesToPoints) {
        const auto& point = G.points[pointId];
        j["nodes"].push_back({ {"id", int(nodeId)},
                              {"x", point.GetX()},
                              {"y", point.GetY()}});
    }

    for (auto const& [src,target] : G.edges) {
        j["edges"].push_back({ {"source", src},
                              {"target", target} });
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
