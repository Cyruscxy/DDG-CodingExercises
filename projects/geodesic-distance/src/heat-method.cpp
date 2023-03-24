// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    this->A = this->geometry->laplaceMatrix();

    const double h = this->geometry->meanEdgeLength();

    this->F = this->geometry->massMatrix() + h * h * this->A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    
    FaceData<Vector3> vectorField(*this->mesh);

    for ( auto f : this->mesh->faces() )
    {
        Vector3 gradient{ 0, 0, 0 };
        for ( auto he : f.adjacentHalfedges() )
        {
            gradient += this->geometry->halfedgeVector(he) * u[he.next().tipVertex().getIndex()];
        }
        gradient = cross(this->geometry->faceNormal(f), gradient);
        vectorField[f] = gradient.normalize();
    }

    return vectorField;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    // TODO

    Vector<double> integrated_divergence(this->mesh->nVertices());
    for ( auto v : this->mesh->vertices() )
    {
    	double grad = 0.0;
		for ( auto he : v.outgoingHalfedges() )
		{
            double cot1 = this->geometry->cotan(he);
            auto X1 = X[he.face()];
            grad += 0.5 * cot1 * dot(X1, this->geometry->halfedgeVector(he));

            double cot2 = this->geometry->cotan(he.twin());
            auto X2 = X[he.twin().face()];
            grad += 0.5 * cot2 * dot(X2, this->geometry->halfedgeVector(he));

		}
        integrated_divergence[v.getIndex()] = grad;
    }

    return integrated_divergence;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // TODO
    Vector<double> phi = Vector<double>::Zero(delta.rows());

    auto flow = this->F;
    auto u = geometrycentral::solvePositiveDefinite<double>(flow, delta);
    auto b = this->computeDivergence(this->computeVectorField(u));
    auto Lc = this->A;
    phi = solvePositiveDefinite<double>(Lc, b);
    phi = -phi;

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}