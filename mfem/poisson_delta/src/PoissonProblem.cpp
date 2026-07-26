#include "PoissonProblem.hpp"

PoissonProblem::PoissonProblem(mfem::FiniteElementSpace* fes)
    : fespace(fes), a(fes), b(fes), delta_coeff(nullptr) { }

void PoissonProblem::SetSourceDelta(double x, double y, double amplitude) {
    if (delta_coeff) delete delta_coeff;
    delta_coeff = new mfem::DeltaCoefficient(x, y, amplitude);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(*delta_coeff));
}

PoissonProblem::~PoissonProblem() {
    delete delta_coeff;
}

void PoissonProblem::SetBoundaryConditions() {
    if (fespace->GetMesh()->bdr_attributes.Size()) {
        mfem::Array<int> ess_bdr(fespace->GetMesh()->bdr_attributes.Max());
        ess_bdr = 1;
        fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }
}

void PoissonProblem::Assemble() {
    diff_coeff.constant = 1.0;    
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(diff_coeff));
    a.Assemble();
    b.Assemble();
}

void PoissonProblem::Solve(mfem::GridFunction& x) {
    x = 0.0;
    mfem::OperatorPtr A;
    mfem::Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
    
    mfem::CG(*A, B, X, 1, 200, 1e-12, 0.0);
    a.RecoverFEMSolution(X, b, x);
}