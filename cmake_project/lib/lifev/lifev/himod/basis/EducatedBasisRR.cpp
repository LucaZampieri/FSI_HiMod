#include <lifev/himod/basis/EducatedBasisRR.hpp>
#include <lifev/himod/basis/EducatedBasisFunctorImplemented.hpp>

// This is for the root finding algorithm, we should find something better.
#include <lifev/core/algorithm/NonLinearBrent.hpp>

namespace LifeV
{

Real
EducatedBasisRR::
Next()
{
    Real eps1 = 1e-12;

    Real Homega;
    Real M_a, M_b;

    // We had analysed the function and the zero are located into the intervals ((i-1)*M_PI + M_PI/2,  i*M_PI + M_PI/2 )
    // but there are some exception depending on where the quadratic factor changes its sign.
    if ( M_icurrent == 0 )
    {
        M_x = M_chi * M_L / M_mu; //At this x the quadratic term changes its sign.
        M_star = false;
        M_istar  = static_cast<UInt> (std::floor ( M_x / M_PI - 0.5 ) + 1); //i = M_istar  (i-1)*M_PI + M_PI/2 < M_x < i*M_PI + M_PI/2
        M_ptrFunctor = static_cast<functor_ptrType> ( new  EducatedBasisFunctorRR (M_mu, M_chi, M_L) );

        if ( M_istar == 0)
        {
            // In this case the first eigenvalue is located into the interval (0, pi/2)
            // so we have to pay attention that the algorithm does not converge on the wrong point which is 0, that is a root of the function, but it is not an eigen value.
            // We made the math and found out that in M_a the first derivative of f is greater than zero so the root we were looking for is on its right because the function in this case in concave.
            M_a = std::min (M_x, M_PI / 4.);
            M_b = M_PI / 2.0;
            ++M_icurrent;
            Homega = NonLinearBrent<functor_type> (*M_ptrFunctor, M_a , M_b - eps1, eps1 / 10., 100);
            return Homega / M_L;
        }
        else
        {
            ++M_icurrent;
            M_a = (M_icurrent * M_PI - M_PI / 2.0);
            M_b = (M_icurrent * M_PI + M_PI / 2.0);

            if (M_a < M_x && M_b > M_x)
            {
                if (!M_star)
                {
                    M_b = (M_a + M_b) / 2.0;
                    M_star = true;
                }
                else
                {
                    M_a = (M_a + M_b) / 2.0;
                    M_star = false;
                    ++M_icurrent;
                }
            }
            else
            {
                ++M_icurrent;
            }
        }
    }
    else
    {
        M_a = (M_icurrent * M_PI - M_PI / 2.0);
        M_b = (M_icurrent * M_PI + M_PI / 2.0);

        if (M_a < M_x && M_b > M_x)
        {
            if (!M_star)
            {
                M_b = (M_a + M_b) / 2.0;
                M_star = true;
            }
            else
            {
                M_a = (M_a + M_b) / 2.0;
                M_star = false;
                ++M_icurrent;
            }
        }
        else
        {
            ++M_icurrent;
        }
    }

    Homega = NonLinearBrent<functor_type> (*M_ptrFunctor, M_a + eps1, M_b - eps1, 1e-6, 100); //TO be changed so that user can change toll and maxiter
    return Homega / M_L; // We return omega without hat.
}


void
EducatedBasisRR::
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
    Real A;
    Real norm;
    Real omega;
    Real gamma;

    for (UInt p = 0; p < dim; ++p)
    {
        //Chi c'è in eigenvalues di p?
        omega = eigenvalues[p] * M_L;

        //B=1
        A     = M_chi * M_L / (M_mu * omega);
        gamma = std::sin (2.0 * omega) / (4.0 * omega);
        norm  = 1.0 / std::sqrt ( ( A * A * (0.5 - gamma) + (0.5 + gamma) + A * (1.0 - std::cos (2.0 * omega) ) / (2.0 * omega) ) );
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

            phi[p][n]  =  norm * ( A * std::sin ( omega * quadrule->quadPointCoor (n, 0) )
                                   +  std::cos ( omega * quadrule->quadPointCoor (n, 0) ) );

            dphi[p][n] = norm * omega * ( A * std::cos ( omega * quadrule->quadPointCoor (n, 0) )
                                          - std::sin ( omega * quadrule->quadPointCoor (n, 0) ) );

        }
    }

}




Real
EducatedBasisRR::
evalSinglePoint ( const Real& eigen, const Real& yh) const
{

    //B=1
    Real omega  = eigen * M_L;
    Real A      = M_chi * M_L / (M_mu * omega);
    Real gamma  = std::sin (2.0 * omega) / (4.0 * omega);
    Real norm   = 1.0 / std::sqrt ( ( A * A * (0.5 - gamma) + (0.5 + gamma) + A * (1.0 - std::cos (2.0 * omega) ) / (2.0 * omega) ) );

    return      norm * ( A * std::sin ( omega * yh ) +  std::cos ( omega * yh ) );
}


} // End LifeV namespace
