// Copyright 2017-2019, Nicholas Sharp and the Polyscope contributors.
// http://polyscope.run.
#include "polyscope/polyscope.h"

namespace polyscope {

namespace state {

bool initialized = false;
double lengthScale = 1.0;
std::tuple<glm::vec3, glm::vec3> boundingBox;
glm::vec3 center{0, 0, 0};
std::map<std::string, std::map<std::string, Structure*>> structures;
std::function<void()> userCallback;
std::set<Widget*> widgets;

// CHANGED: for DDG
size_t facePickIndStart;
size_t edgePickIndStart;
size_t halfedgePickIndStart;
MeshSubset subset;
int currVertexIndex = -1;
int currFaceIndex = -1;
int currEdgeIndex = -1;
int deleteVertexIndex = -1;
int deleteFaceIndex = -1;
int deleteEdgeIndex = -1;

} // namespace state
} // namespace polyscope
