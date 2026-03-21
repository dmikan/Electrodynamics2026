#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main()
{
   // 1. Create mesh square 50x50 of triangles.
   Mesh mesh = Mesh::MakeCartesian2D(100, 100, Element::TRIANGLE);
   
   // 2. Finite element space (H1, order 1)
   H1_FECollection fec(1, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fec);

   // 3. Boundary conditions (Dirichlet = 0 on the edges)
   Array<int> ess_tdof_list;
   if (mesh.bdr_attributes.Size())
   {
      Array<int> ess_bdr(mesh.bdr_attributes.Max());
      ess_bdr = 1; 
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 4. Right-hand side: Delta of Dirac in the center (0.5, 0.5)
   LinearForm b(&fespace);
   
   double x_centro = 0.5;
   double y_centro = 0.5;
   double amplitud = 1.0;
   
   DeltaCoefficient delta(x_centro, y_centro, amplitud);
   
   b.AddDomainIntegrator(new DomainLFIntegrator(delta));
   b.Assemble();

   // 5. Operator: Laplacian (-Delta)
   BilinearForm a(&fespace);
   ConstantCoefficient one(1.0);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));
   a.Assemble();

   // 6. Solve Ax = B
   GridFunction x(&fespace);
   x = 0.0;

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   CG(*A, B, X, 1, 200, 1e-12, 0.0);

   // 7. Recover solution and save files
   a.RecoverFEMSolution(X, b, x);
   
   ofstream sol_ofs("solution.gf");
   x.Save(sol_ofs);
   ofstream mesh_ofs("mesh.mesh");
   mesh.Print(mesh_ofs);

   // 8. Visualization in GLVis (Port 19916 for automatic redirection)
   char vishost[] = "localhost";
   int  visport   = 19916; 
   socketstream sol_sock(vishost, visport);
   
   if (sol_sock.is_open())
   {
      sol_sock.precision(8);
      // "keys R2mc" to reset camera, 2D mode, show mesh and center
      sol_sock << "solution\n" << mesh << x << "keys R2mc\n" << flush;
      
      std::cout << "Point load simulation sent to GLVis (Port " << visport << ")." << std::endl;
      std::cout << "Press Enter to finish..." << std::endl;
      std::cin.get();
   }
   else
   {
      std::cout << "Could not connect to GLVis at " << vishost << ":" << visport << std::endl;
   }

   return 0;
}