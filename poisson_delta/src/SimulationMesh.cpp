#include "SimulationMesh.hpp"

SimulationMesh::SimulationMesh(int nx, int ny, int order)
    : mesh(mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::TRIANGLE)),
      fec(order, mesh.Dimension()),
      fespace(&mesh, &fec)
{
    // Aquí podrías añadir refinamiento de malla si quisieras
}