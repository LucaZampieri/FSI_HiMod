#include <lifev/himod/basis/EducatedBasisRR.hpp>
#include <lifev/himod/basis/EducatedBasisFunctorImplemented.hpp>

// This is for the root finding algorithm, we should find something better.
#include <lifev/core/algorithm/NonLinearBrent.hpp>

namespace LifeV
{

Real
EducatedBasisDR::
Next()
{
    Real eps1 = 1e-12;

    if (M_icurrent == 0)
    {
        M_ptrFunctor = static_cast<functor_ptrType> ( new  EducatedBasisFunctorDR (M_mu, M_chi, M_L) );
    }


    Real Homega;
    Real M_a, M_b;
    M_icurrent++;


    M_a = M_icurrent * M_PI + M_PI / 2;
    M_b = M_PI * (M_icurrent + 1);
    Homega = NonLinearBrent<functor_type> (*M_ptrFunctor, M_a , M_b - eps1, eps1 / 10., 100);

    return Homega / M_L;
}


void
EducatedBasisDR::
evaluateBasis (          Basis1DAbstract::MBMatrix_type& phi,
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
    Real norm;
    Real omega;

    for (UInt p = 0; p < dim; ++p)
    {
        //Chi c'è in eigenvalues di p?
        omega = eigenvalues[p] * M_L;

        //B=1
        //A=1
        norm   = 1.0 / std::sqrt ( (2 * M_L - std::sin (2 * M_L) ) / 4. );
        /*
                std::cout<<" p = "<<p;
                std::cout<<" A = "<<A;
                std::cout<<" M_chi = "<<M_chi;
                std::cout<<" M_L = "<<M_L;
                std::cout<<" M_mu = "<<M_mu;
                std::cout<<" norm = "<<norm;
                std::cout<<" omega = "<<omega;
                std::cout<<" gamma = "<<gamma<<std::endl;
        */
        for (UInt n = 0; n < quadrule->nbQuadPt(); ++n)
        {

            phi[p][n]  =  norm * (  std::sin ( omega * quadrule->quadPointCoor (n, 0) ) );

            dphi[p][n] = norm * omega * ( std::cos ( omega * quadrule->quadPointCoor (n, 0) ) );

        }
    }

}



//gli viene passato un valore nell'intervallo di riferimento [0,1]
Real
EducatedBasisDR::
evalSinglePoint ( const Real& eigen, const Real& yh) const
{

    //B=1
    //A = 1
    Real omega  = eigen * M_L;
    Real norm   = 1.0 / std::sqrt ( (2 * M_L - std::sin (2 * M_L) ) / 4. );

    return      norm * std::sin ( omega * yh );
}


} // End LifeV namespace
