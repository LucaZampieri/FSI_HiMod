#include <lifev/himod/basis/EducatedBasisDD.hpp>

namespace LifeV
{

void
EducatedBasisDD::
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

    // These functions should be normalized in respect to L^2 ( 0 , 1 )
    Real norm2 = std::sqrt (2.0);

    for (UInt p = 0; p < dim; ++p)
        for (UInt n = 0; n < quadrule->nbQuadPt(); ++n)
        {
            phi[p][n]  = norm2 * std::sin ( (p + 1) * M_PI * quadrule->quadPointCoor (n, 0) );
            dphi[p][n] = norm2 * std::cos ( (p + 1) * M_PI * quadrule->quadPointCoor (n, 0) ) * M_PI * (p + 1);

        }
}

Real
EducatedBasisDD::
evalSinglePoint ( const Real& eigen, const Real& yh) const
{
    return std::sqrt (2.0) * std::sin ( eigen * M_L * yh );
}

}//End lifev namespace
