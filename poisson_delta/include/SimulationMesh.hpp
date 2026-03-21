#ifndef SIMULATION_MESH_HPP
#define SIMULATION_MESH_HPP

#include "mfem.hpp"

class SimulationMesh {
private:
    mfem::Mesh mesh;
    mfem::H1_FECollection fec;
    mfem::FiniteElementSpace fespace;

public:
    SimulationMesh(int nx, int ny, int order = 1);
    
    mfem::FiniteElementSpace* GetSpace() { return &fespace; }
    mfem::Mesh* GetMesh() { return &mesh; }
};

#endif