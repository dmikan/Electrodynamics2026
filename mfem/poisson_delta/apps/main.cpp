#include "SimulationMesh.hpp"
#include "PoissonProblem.hpp"
#include "Visualizer.hpp"

int main() {
    // 1. Geometría
    SimulationMesh sim(100, 100);

    // 2. Física
    PoissonProblem physics(sim.GetSpace());
    physics.SetSourceDelta(0.5, 0.5, 1.0);
    physics.SetBoundaryConditions();
    physics.Assemble();

    // 3. Cálculo
    mfem::GridFunction sol(sim.GetSpace());
    physics.Solve(sol);

    // 4. Visualización (Llamada limpia al método estático)
    Visualizer::SendToGLVis(*sim.GetMesh(), sol, 19916);
    std::cout << "\n--- Simulation Complete ---" << std::endl;
    std::cout << "The data has been sent to GLVis." << std::endl;
    std::cout << "Press Enter to close the program and the connection..." << std::endl;
    std::cin.get();

    return 0;
}