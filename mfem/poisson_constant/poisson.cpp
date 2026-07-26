#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main()
{
   // 1. Create a square mesh of 10x10 elements directly.
   mfem::Mesh mesh = mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::TRIANGLE);
   
   // 2. Create a finite element space. We use polynomials of order 1 (linear) in 2D.
   mfem::H1_FECollection fec(1, mesh.Dimension());
   mfem::FiniteElementSpace fespace(&mesh, &fec);

   // 3. Boundary conditions (Dirichlet homogeneous). We mark all edges to have a solution of 0 there.
   mfem::Array<int> ess_tdof_list;
   if (mesh.bdr_attributes.Size())
   {
      mfem::Array<int> ess_bdr(mesh.bdr_attributes.Max());
      ess_bdr = 1; 
      fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   }

   // 4. The right-hand side (the source / Linear Form). It represents the integral of (1 * test function).
   mfem::LinearForm b(&fespace);
   mfem::ConstantCoefficient one(1.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   b.Assemble();

   // 5. The operator (the stiffness matrix / Bilinear Form). It represents the Laplacian (-Delta).
   mfem::BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));
   a.Assemble();

   // 6. Solve the system Ax = B. We create the solution vector x (initialized to 0).
   mfem::GridFunction x(&fespace);
   x = 0.0;

   mfem::OperatorPtr A;
   mfem::Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   // Use a simple solver (Conjugate Gradient).
   mfem::CG(*A, B, X, 1, 200, 1e-12, 0.0);

   // 7. Recover and save the solution.
   a.RecoverFEMSolution(X, b, x);
   
   std::ofstream sol_ofs("solution.gf");
   x.Save(sol_ofs);
   std::ofstream mesh_ofs("mesh.mesh");
   mesh.Print(mesh_ofs);

   // 8. Send to GLVIS (To see it in the browser)
   char vishost[] = "localhost";
   int  visport   = 19916; 
   mfem::socketstream sol_sock(vishost, visport);
   
   if (sol_sock.is_open())
   {
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << "keys R2mc\n" << std::flush;
      
      std::cout << "Visualization sent to GLVis server on port " << visport << "." << std::endl;
      std::cout << "Press Enter to close the program and the connection..." << std::endl;
      std::cin.get();
   }
   else
   {
      std::cout << "Could not connect to GLVis server at " << vishost << ":" << visport << std::endl;
      std::cout << "Make sure the server is running." << std::endl;
   }

   std::cout << "Simulation finished successfully!" << std::endl;

   return 0;
}