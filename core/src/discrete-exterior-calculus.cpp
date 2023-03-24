// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    
    std::vector<Eigen::Triplet<double>> tripList;
    for (Vertex v : mesh.vertices()) {
        auto vid = v.getIndex();
        double A = this->barycentricDualArea(v);
        tripList.push_back(Eigen::Triplet<double>(vid, vid, A));
    }

    SparseMatrix<double> result(mesh.nVertices(), mesh.nVertices());
    result.setFromTriplets(tripList.begin(), tripList.end());

    return result; 
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    std::vector<Eigen::Triplet<double>> tripList;
    for (Edge e : mesh.edges()) {
        size_t eid = e.getIndex();
        double cot = 0.0;
        Halfedge he1 = e.halfedge();
        Halfedge he2 = he1.twin();
        cot += cotan(he1);
        cot += cotan(he2);
        tripList.push_back(Eigen::Triplet<double>(eid, eid, 0.5 * cot));
    }
    Eigen::SparseMatrix<double> result(mesh.nEdges(), mesh.nEdges());
    result.setFromTriplets(tripList.begin(), tripList.end());
    return result; 
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    std::vector<Eigen::Triplet<double>> tripList;
    for (Face f : mesh.faces()) {
        auto fid = f.getIndex();
        std::vector<double> edgeList;
        for (Edge e : f.adjacentEdges()) {
            edgeList.push_back(this->edgeLength(e));
        }
        double s = 0.0;
        for (auto l : edgeList) {
            s += l;
        }
        double A = 0.0;
        for (auto l : edgeList) {
            A *= s - l;
        }
        A = 1.0 / sqrt(s * A);
        tripList.push_back(Eigen::Triplet<double>(fid, fid, A));
    }
    Eigen::SparseMatrix<double> result(mesh.nFaces(), mesh.nFaces());
    result.setFromTriplets(tripList.begin(), tripList.end());

    return result; 
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    std::vector<Eigen::Triplet<double>> tripList;

    for (Edge e : mesh.edges()) {
        auto eid = e.getIndex();
        Halfedge he = e.halfedge();
        tripList.push_back(Eigen::Triplet<double>(eid, e.firstVertex().getIndex(), -1.0));
        tripList.push_back(Eigen::Triplet<double>(eid, e.secondVertex().getIndex(), 1.0));
    }

    Eigen::SparseMatrix<double> result(mesh.nEdges(), mesh.nVertices());
    result.setFromTriplets(tripList.begin(), tripList.end());

    return result;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    std::vector<Eigen::Triplet<double>> tripList;

    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            Edge e = he.edge();
            double coeff = he.orientation()? 1.0 : -1.0;
            tripList.push_back(Eigen::Triplet<double>(f.getIndex(), e.getIndex(), coeff));
        }
    }

    Eigen::SparseMatrix<double> result(mesh.nFaces(), mesh.nEdges());
    result.setFromTriplets(tripList.begin(), tripList.end());
    return result;
}

} // namespace surface
} // namespace geometrycentral