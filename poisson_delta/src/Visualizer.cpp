#include "Visualizer.hpp"
#include <iostream>

void Visualizer::SendToGLVis(mfem::Mesh& mesh, mfem::GridFunction& x, int port) {
    char vishost[] = "localhost";
    mfem::socketstream sol_sock(vishost, port);

    if (sol_sock.is_open()) {
        sol_sock.precision(8);
        // Enviamos la malla, la solución y comandos de cámara de GLVis
        sol_sock << "solution\n" << mesh << x << "keys R2mc\n" << std::flush;
        std::cout << "Data sent to GLVis on port " << port << std::endl;
    } else {
        std::cerr << "Unable to connect to GLVis at " << vishost << ":" << port << std::endl;
    }
}