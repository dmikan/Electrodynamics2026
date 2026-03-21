#ifndef POISSON_PROBLEM_HPP
#define POISSON_PROBLEM_HPP

#include "mfem.hpp"

class PoissonProblem {
private:
    mfem::FiniteElementSpace* fespace;
    mfem::BilinearForm a;
    mfem::LinearForm b;
    mfem::DeltaCoefficient* delta_coeff;
    mfem::ConstantCoefficient diff_coeff;
    mfem::Array<int> ess_tdof_list;

public:
    PoissonProblem(mfem::FiniteElementSpace* fes);
    ~PoissonProblem();
    void SetSourceDelta(double x, double y, double amplitude);
    void SetBoundaryConditions();
    void Assemble();
    void Solve(mfem::GridFunction& x);
};

#endif