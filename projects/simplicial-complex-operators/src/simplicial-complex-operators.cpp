// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }
    
    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    std::vector<Eigen::Triplet<size_t>> tripletList;
    SparseMatrix<size_t> result(mesh->nVertices(), mesh->nEdges());

    for (Vertex v : mesh->vertices()) {
        size_t vid = v.getIndex();
        for (Edge e : v.adjacentEdges()) {
            tripletList.push_back(Eigen::Triplet<size_t>(vid, e.getIndex(), 1));
        }
    }
    result.setFromTriplets(tripletList.begin(), tripletList.end());

    return result;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    std::vector<Eigen::Triplet<size_t>> tripletList;
    SparseMatrix<size_t> result(mesh->nFaces(), mesh->nEdges());

    for (Face f : mesh->faces()) {
        size_t fid = f.getIndex();
        for (Edge e : f.adjacentEdges()) {
            tripletList.push_back(Eigen::Triplet<size_t>(fid, e.getIndex(), 1));
        }
    }
    result.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return result; // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector<size_t> encodingVector = Vector<size_t>::Zero(mesh->nVertices());
    for (size_t vid : subset.vertices) {
        encodingVector(vid) = 1;
    }
    return encodingVector;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector<size_t> encodingVector = Vector<size_t>::Zero(mesh->nEdges());
    for (size_t eid : subset.edges) {
        encodingVector(eid) = 1;
    }
    return encodingVector;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    Vector<size_t> encodingVector = Vector<size_t>::Zero(mesh->nFaces());
    for (size_t fid : subset.faces) {
        encodingVector(fid) = 1;
    }
    return encodingVector;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    MeshSubset star = subset.deepCopy();

    // get all the adjacent elements of the vertices in the given subset
    for (size_t vid : subset.vertices) {
        for (Edge e : mesh->vertex(vid).adjacentEdges()) {
            star.edges.insert(e.getIndex());
        }

        for (Face f : mesh->vertex(vid).adjacentFaces()) {
            star.faces.insert(f.getIndex());
        }
    }

    // get all the adjacent elements of the vertices in the given subset
    for (size_t eid : subset.edges) {
        for (Vertex v : mesh->edge(eid).adjacentVertices()) {
            star.vertices.insert(v.getIndex());
        }

        for (Face f : mesh->edge(eid).adjacentFaces()) {
            star.faces.insert(f.getIndex());
        }
    }

    // get all the adjacent elements of the faces in the given subset
    for (size_t fid : subset.faces) {
        for (Edge e : mesh->face(fid).adjacentEdges()) {
            star.edges.insert(e.getIndex());
        }

        for (Vertex v : mesh->face(fid).adjacentVertices()) {
            star.vertices.insert(v.getIndex());
        }
    }

    return star; 
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    MeshSubset closure = subset.deepCopy();

    for (size_t fid : subset.faces) {
        Face f = mesh->face(fid);

        for (Vertex v : f.adjacentVertices()) {
            closure.vertices.insert(v.getIndex());
        }
        for (Edge e : f.adjacentEdges()) {
            closure.edges.insert(e.getIndex());
        }
    }

    for (size_t eid : subset.edges) {
        Edge e = mesh->edge(eid);

        for (Vertex v : e.adjacentVertices()) {
            closure.vertices.insert(v.getIndex());
        }
    }

    return closure; 
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    MeshSubset closureOfStar = closure(star(subset));
    MeshSubset starOfClosure = star(closure(subset));
    closureOfStar.deleteSubset(starOfClosure);
    return closureOfStar;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    for (size_t fid : subset.faces) {
        Face f = mesh->face(fid);
        for (auto e : f.adjacentEdges()) {
            if (subset.edges.find(e.getIndex()) == subset.edges.end())
                return false;
        }

        for (auto v : f.adjacentVertices()) {
            if (subset.vertices.find(v.getIndex()) == subset.vertices.end())
                return false;
        }
    }

    for (size_t eid : subset.edges) {
        Edge e = mesh->edge(eid);
        for (auto v : e.adjacentVertices()) {
            if (subset.vertices.find(v.getIndex()) == subset.vertices.end())
                return false;
        }
    }

    return true; 
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    if (!isComplex(subset)) return -1;

    int degree = 2;
    if (subset.faces.empty()) {
        // if the subset has no faces, the maximum dgree is 1
        degree--;
    }
    else if (subset.edges.empty()) {
        // if the subset has no faces & edges, the maximum dgree is 0
        degree--;
        return degree;
    }

    if (degree == 2) {
        for (size_t eid : subset.edges) {
            bool isInSubset = false;
            for (auto f : mesh->edge(eid).adjacentFaces()) {
                if (subset.faces.find(f.getIndex()) != subset.faces.end())
                    isInSubset = true;
            }
            if (!isInSubset) return -1;
        }

        for (size_t vid : subset.vertices) {
            bool isInSubset = false;
            for (auto e : mesh->vertex(vid).adjacentEdges()) {
                if (subset.edges.find(e.getIndex()) != subset.edges.end())
                    isInSubset = true;
            }
            if (!isInSubset) return -1;
        }
    }
    else if (degree == 1) {
        for (size_t vid : subset.vertices) {
            bool isInSubset = false;
            for (auto e : mesh->vertex(vid).adjacentEdges()) {
                if (subset.edges.find(e.getIndex()) != subset.edges.end())
                    isInSubset = true;
            }
            if (!isInSubset) return -1;
        }
    }

    return degree;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    MeshSubset boundary;

    int degree = isPureComplex(subset);
   
    if (degree == -1) return boundary;

    if (degree == 2) {
        for (size_t eid : subset.edges) {
            Edge e = mesh->edge(eid);

            bool onlyOne = false;
            for (auto f : e.adjacentFaces()) {
                if (subset.faces.find(f.getIndex()) != subset.faces.end()) {
                    if (onlyOne) {
                        onlyOne = false;
                        break;
                    }
                    else {
                        onlyOne = true;
                    }    
                }
            }

            if (onlyOne)
                boundary.edges.insert(eid);
        }
    }
    else if (degree == 1) {
        for (size_t vid : subset.vertices) {
            Vertex v = mesh->vertex(vid);

            bool onlyOne = false;
            for (auto e : v.adjacentEdges()) {
                if ((subset.edges.find(e.getIndex()) != subset.edges.end())) {
                    if (onlyOne) {
                        onlyOne = false;
                        break;
                    }
                    else {
                        onlyOne = true;
                    }
                }
            }

            if (onlyOne)
                boundary.vertices.insert(vid);
        }
    }
    
    return closure(boundary); // placeholder
}