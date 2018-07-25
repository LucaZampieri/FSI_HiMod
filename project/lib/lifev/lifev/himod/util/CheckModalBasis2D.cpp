#include <lifev/himod/util/CheckModalBasis2D.hpp>

namespace LifeV
{

CheckModalBasis2D::
CheckModalBasis2D (const modalbasis_ptrType& modalbasis ) : M_modalbasis ( modalbasis )
{
}

void CheckModalBasis2D::
VerifyOrthonormality() const
{
    std::cout << std::endl;
    std::cout << "Checking orthonormality... " << std::endl;
    std::cout << std::endl;
    std::cout << "In theory, the basis should be orthonormal. So when you try to compute " << std::endl;
    std::cout << "the integral on the slice of Phi_j * Phi_i. You should obtain a Kronecker Delta as a result." << std::endl;
    std::cout << "The result depends a lot on the order of the quadrature rule, if the number of modes increases you will have problems." << std::endl;
    std::cout << std::endl;

    UInt m = M_modalbasis->mtot();
    Real err( 0 );
    std::vector<Real> delta( m * m, 0.0 );
    for ( UInt j( 0 ); j != m; ++j )
    {
        for ( UInt i( 0 ); i != m; ++i )
        {

            delta[i + m * j] = M_modalbasis->compute1_PhiPhi( i, j );

            if( i == j )
            {
                err += std::abs( delta[i + m * j] - 1 );
            }
            else
            {
                err += std::abs(delta[i + m * j] );
            }
        }
    }
    std::cout << "Now we compute the integral and compare it with the kronecker delta, obtaining a matrix that should be " << std::endl;
    std::cout << "equal to an identity matrix, we compute the difference and, then, the sum of the absolute value of the result" << std::endl;
    std::cout << "The result should be close to zero." << std::endl;
    std::cout << std::endl;
    std::cout << "The sum of the absolute value is: " << err << std::endl;
    std::cout << std::endl;

}

void CheckModalBasis2D::
VerifyBC( const Real& mu, const Real& chi ) const
{
    std::cout << "We want to verify that the basis satisfies the required boundary conditions so we report the value on the boundary." << std::endl;
    UInt ntheta = M_modalbasis->qrTheta().nbQuadPt() - 1;
    UInt nrho = M_modalbasis->qrRho().nbQuadPt() - 1;

    std::cout << "Here we put the coefficient, be aware that they may not make sense for a Dirichlet or Neumann BC." << std::endl;
    std::cout << "mu <- " << mu << std::endl;
    std::cout << "chi <- " << chi << std::endl;

    std::cout << "If you use a quadrature rule that doesn't have a node in the boundary" << std::endl;
    std::cout << "the results may be inaccurate, so we report here the position of the point where the basis is actually evaluated" << std::endl;

    std::cout << "First quadrature node in rho <- ";
    std::cout << M_modalbasis->qrRho().quadPointCoor( 0, 0 ) << std::endl;
    std::cout << "First quadrature node in theta <- ";
    std::cout << M_modalbasis->qrTheta().quadPointCoor( 0, 0 ) << std::endl;
    std::cout << "Last quadrature node in rho <- ";
    std::cout << M_modalbasis->qrRho().quadPointCoor( nrho, 0 ) << std::endl;
    std::cout << "Last quadrature node in theta <- ";
    std::cout << M_modalbasis->qrTheta().quadPointCoor( ntheta, 0 ) << std::endl;

    for ( UInt k( 0 ); k != M_modalbasis->mtot(); ++k )
    {

        std::cout << "###########################" << std::endl;
        std::cout << " m <- " << k + 1 << std::endl;
        std::cout << " Dirichlet  " << "<- \t";
        std::cout << chi * M_modalbasis->phirho( k, nrho ) << std::endl;

        std::cout << " Neumann " << "<- \t";
        std::cout << mu * M_modalbasis->dphirho( k, nrho )<< std::endl;

        std::cout << " Robin " << "<- \t";
        std::cout << chi* M_modalbasis->phirho( k, nrho ) + mu* M_modalbasis->dphirho( k, nrho ) << std::endl;
    }

}
} //End lifev namespace
