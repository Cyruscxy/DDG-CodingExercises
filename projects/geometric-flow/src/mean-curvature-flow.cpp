// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    auto C = geometry->laplaceMatrix();
    auto result = M + h * C;
    
    return result;
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions

    DenseMatrix<double> f0(mesh->nVertices(), 3);
    auto M = geometry->massMatrix();
    auto A = this->buildFlowOperator(M, h);

    for (Vertex v : mesh->vertices()) {
        auto& p = geometry->inputVertexPositions[v];
        f0.row(v.getIndex()) << p.x, p.y, p.z;
    }

    Vector<double> x0 = M * f0.col(0);
    Vector<double> y0 = M * f0.col(1);
    Vector<double> z0 = M * f0.col(2);

    auto xh = geometrycentral::solvePositiveDefinite(A, x0);
    auto yh = geometrycentral::solvePositiveDefinite(A, y0);
    auto zh = geometrycentral::solvePositiveDefinite(A, z0);

    for (auto v : mesh->vertices()) {
        auto index = v.getIndex();
        geometry->inputVertexPositions[v].x = xh[index];
        geometry->inputVertexPositions[v].y = yh[index];
        geometry->inputVertexPositions[v].z = zh[index];
    }
}