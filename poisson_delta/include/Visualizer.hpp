#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "mfem.hpp"

class Visualizer {
public:
    // Un método estático para no tener que instanciar la clase
    static void SendToGLVis(mfem::Mesh& mesh, mfem::GridFunction& x, int port = 19916);
};

#endif