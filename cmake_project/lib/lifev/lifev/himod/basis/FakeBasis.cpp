#include <lifev/himod/basis/FakeBasis.hpp>

namespace LifeV
{

void
FakeBasis::
evaluateBasis (     Basis1DAbstract::MBMatrix_type& phi,
                    Basis1DAbstract::MBMatrix_type& dphi,
                    const Basis1DAbstract::MBVector_type& eigenvalues,
                    const QuadratureRule* quadrule) const
{
    UInt dim = eigenvalues.size();
    phi.resize (dim);
    dphi.resize (dim);

    for (UInt i = 0; i < dim; ++i)
    {
        phi[i].resize (quadrule->nbQuadPt() );
        dphi[i].resize (quadrule->nbQuadPt() );
    }

    for (UInt p = 0; p < dim; ++p)
        for (UInt n = 0; n < quadrule->nbQuadPt(); ++n)
        {
            phi[p][n]  = 1.;
            dphi[p][n] = 0.;

        }
}

Real
FakeBasis::
evalSinglePoint ( const Real& /*eigen*/, const Real& /*yh*/) const
{
    return 1.;
}

}//End lifev namespace
