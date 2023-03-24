// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    // building area matrix in complex plane
    std::vector<Eigen::Triplet<std::complex<double>>> areaTripList;

    for (auto bl : this->mesh->boundaryLoops()) {
        for (auto he : bl.adjacentHalfedges()) {
            const auto i = he.tailVertex().getIndex();
            const auto j = he.tipVertex().getIndex();

            areaTripList.emplace_back(Eigen::Triplet<std::complex<double>>(i, j, std::complex<double>(0.0, -0.5)));
            areaTripList.emplace_back(Eigen::Triplet<std::complex<double>>(j, i, std::complex<double>(0.0, 0.5)));
        }
    }

    SparseMatrix<std::complex<double>> areaMatrix(this->mesh->nVertices(), this->mesh->nVertices());
    areaMatrix.setFromTriplets(areaTripList.begin(), areaTripList.end());

    // EC = ED - A
    return this->geometry->complexLaplaceMatrix() - areaMatrix;

}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    SurfaceMesh& m = this->geometry->mesh;
	VertexData<Vector2> mapping(m);

    auto complexMapping = solveInversePowerMethod(buildConformalEnergy());

    for ( auto v : this->mesh->vertices() ) {
        mapping[v] = Vector2{ complexMapping[v.getIndex()].real(), complexMapping[v.getIndex()].imag() };
    }

    return mapping;
}