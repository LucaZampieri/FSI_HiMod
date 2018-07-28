#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/himod/modalbasis/ModalSpaceCircular.hpp>

namespace LifeV
{

void ModalSpaceCircular::
addSliceBC( const std::string& BC, const Real& mu, const Real& Chi )
{
   M_genbasisRhoTheta = EducatedBasisFactory::instance().createObject( BC );
   M_genbasisRhoTheta -> setRho( M_Rho );
   M_genbasisRhoTheta -> setTheta( M_Theta );
   M_genbasisRhoTheta -> setMu( mu );
   M_genbasisRhoTheta -> setChi( Chi );
   M_genbasisRhoTheta -> setNumberModes(M_mtot);
   eigensProvider();
}

void ModalSpaceCircular::
evaluateBasis()
{
    M_genbasisRhoTheta -> evaluateBasis( M_phirho, M_dphirho, M_phitheta, M_dphitheta, M_eigenvalues, M_quadruleRho, M_quadruleTheta );
}

void ModalSpaceCircular::
showMe() const
{
    std::cout << "---- MODAL SPACE CIRCULAR SHOWME ---" << std::endl;
    std::cout << "Rho = "	<<	M_Rho					<< std::endl;
    std::cout << "Theta = "	<<	M_Theta					<< std::endl;
    std::cout << "M = "		<<	M_mtot					<< std::endl;
    std::cout << "------------------------------------" << std::endl;
}

void ModalSpaceCircular::
eigensProvider()
{

    M_genbasisRhoTheta -> computeEigenvalues();
    
    for(UInt i = 0; i != ( M_genbasisRhoTheta->eigenValues() ).size(); ++i)
    {
    	M_eigenvalues.push_back( ( M_genbasisRhoTheta->eigenValues() )[i] );
    }  
	return;
}


//Fourier Coefficients in the case g is x-independent
std::vector<Real> ModalSpaceCircular::
fourierCoefficients( const function_Type& g ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mtot, 0.0 );
    
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >		evaluate_g;

    // first resize
    evaluate_g.resize( M_quadruleRho->nbQuadPt() );

    // second resize
    for ( UInt k = 0; k != evaluate_g.size(); ++k )
    {
        evaluate_g[k].resize(M_quadruleTheta->nbQuadPt() );
    }
    // evaluation of g
    for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
        {
            evaluate_g[n][m] = g( 0 , 0 , M_quadruleRho -> quadPointCoor (n, 0) * M_Rho , M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mtot; ++k )
    {
        // extraction of the sub-indeces
//        UInt p_k = ( M_eigenvalues[k] - 1 ); 

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt(2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )							//y
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )						//z
                FourCoeff[k] += evaluate_g[n][m]*                                           // function evaluated in the right nodes
                                M_phirho[k][ n] * normrho *
                                M_phitheta[k][m] * normtheta *     			            // \f$ \phi_k\f$ evaluated on the nodes
                                M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho*                // Jacobian
                                M_quadruleRho->weight( n ) * M_Rho *
                                M_quadruleTheta->weight( m ) * M_Theta;    // weights

    }
    return FourCoeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real ModalSpaceCircular::
fourierCoeffPointWise ( const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt(2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
            coeff += f ( 0 , x , M_quadruleRho->quadPointCoor ( n , 0 ) * M_Rho , M_quadruleTheta->quadPointCoor ( m , 0 ) * M_Theta , 0 ) *
                     M_phirho[k][n] * normrho * M_phitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) * M_Rho *
                     M_Theta * M_Rho * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );

    return coeff;
}

//Compute the integral on modal basis
//PHI PHI
Real ModalSpaceCircular::
compute1_PhiPhi( const UInt& j, const UInt& k ) const
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += normrho * M_phirho [j][n] *
                   normrho * M_phirho [k][n] * 
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
                   M_quadruleRho->weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n < M_quadruleTheta -> nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute2_PhiPhi( const UInt& j, const UInt& k ) const 
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += normrho * M_phirho [j][n] *
                   normrho * M_phirho [k][n] * 
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho * M_Rho *
                   M_quadruleRho->weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n != M_quadruleTheta -> nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute1_DrhoPhiPhi( const UInt& j, const UInt& k ) const 
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;
    
    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_dphirho[j][n] * normrho * normrho * // the former normrho is due to normalization, the latter to derivative
                   M_phirho[k][n] * normrho *
                   M_quadruleRho->quadPointCoor( n, 0 ) *M_Rho * M_quadruleRho -> weight(n) * M_Rho;
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute2_DrhoPhiPhi( const UInt& j, const UInt& k ) const 
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_dphirho[j][n] * normrho * normrho *
                   M_phirho[k][n] * normrho *
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho * M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
                   M_quadruleRho -> weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute2_DrhoPhiDrhoPhi ( const UInt& j, const UInt& k ) const
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for (UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n)
    {
        coeff_rho += M_dphirho [j][n] * normrho * normrho * // the former normrho is due to normalization, the latter to derivative
                   M_dphirho [k][n] * normrho * normrho *
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *  M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
                   M_quadruleRho->weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute1_DrhoPhiDrhoPhi ( const UInt& j, const UInt& k ) const
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for (UInt n = 0; n < M_quadruleRho->nbQuadPt(); ++n)
    {
        coeff_rho += M_dphirho [j][n] * normrho *
                   M_dphirho [k][n] * normrho *
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n < M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute1_DthetaPhiPhi (const UInt& j, const UInt& k) const 
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_phirho[j][n] * normrho *
                   M_phirho[k][n] * normrho *
                   M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
                   M_Rho * M_quadruleRho->weight( n );
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += ( 1. / ( 2. * M_PI ) ) * M_dphitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_quadruleTheta->weight( n ) * M_Theta;
    }

    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute0_DthetaPhiDthetaPhi( const UInt& j, const UInt& k ) const 
{
    Real coeff_rho = 0.0;
    Real coeff_theta = 0.0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_phirho[j][n] * normrho *
                   M_phirho[k][n] * normrho *
                   M_Rho * M_quadruleRho->weight( n );
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += ( 1. / ( 2 * M_PI ) ) * M_dphitheta[j][n] * normtheta *
                   (1. / ( 2 * M_PI ) ) * M_dphitheta[k][n] * normtheta *
                   M_Theta * M_quadruleTheta->weight( n );
    }

    return coeff_rho * coeff_theta;
}

//PHI
Real ModalSpaceCircular::
compute_Phi( const UInt& k ) const
{
    Real coeff_rho = 0;
    Real coeff_theta = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_phirho[k][n] * normrho * M_Rho * 
        		     M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho * 
        			 M_quadruleRho->weight( n ) * M_Rho;
    }

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[k][n] * normtheta * M_Theta * M_quadruleTheta->weight( n );
    }
    return coeff_rho * coeff_theta;
}

Real ModalSpaceCircular::
compute_rho_PhiPhi( const UInt& j, const UInt& k ) const 
{
    Real coeff_rho = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        coeff_rho += M_phirho[j][n] * normrho *
                   M_phirho [k][n] * normrho *
                   M_quadruleRho->weight( n ) * M_Rho;
    }

    return coeff_rho;
}

Real ModalSpaceCircular::
compute_theta_PhiPhi( const UInt& j, const UInt& k ) const
{
    Real coeff_theta = 0.0;

	Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleTheta->nbQuadPt(); ++n )
    {
        coeff_theta += M_phitheta[j][n] * normtheta *
                   M_phitheta[k][n] * normtheta *
                   M_quadruleTheta->weight( n ) * M_Theta;
    }

    return coeff_theta;
}

Real
ModalSpaceCircular::
compute_R11( const UInt& j, const UInt& k, const function_Type& mu, const Real& x ) const
{
    Real coeff = 0;

	Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {

            coeff += mu( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                         M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 0 ) *
                     M_Rho * M_phirho[j][n] * M_phirho[k][n] * M_quadruleRho->quadPointCoor( n, 0 ) *
                     M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->weight( n ) *
                     M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m );
        }
    }
    return coeff;
}

Real
ModalSpaceCircular::
compute_R10( const UInt& j, const UInt& k, const function_Type& beta, const Real& x ) const
{
    Real coeff = 0;

	Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
            coeff += beta( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                           M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 0 ) *
                     M_Rho * M_phirho[j][n] * M_phirho[k][n] * 
                     M_quadruleRho->quadPointCoor( n, 0 ) *M_quadruleRho->quadPointCoor( n, 0 ) *
                     M_quadruleRho->weight( n ) *
                     M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m );
        }
    }

    return coeff;
}


Real
ModalSpaceCircular::
compute_R00( const UInt& j, const UInt& k, const function_Type& mu, const function_Type& beta, const function_Type& sigma, const Real& x ) const
{
    Real coeff = 0;

   	Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
            coeff += mu( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                         M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 0 ) *
                     (
                     	 normrho * normrho * normrho * normrho * M_Rho * M_Rho * M_Rho *
                         M_dphirho[j][n] * M_dphirho[k][n] 
                         * M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->quadPointCoor( n, 0 ) * 
                         M_quadruleRho->weight( n ) *  
                         normtheta * normtheta * M_Theta *
                         M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m ) +
                         
                         normrho * normrho * normrho * M_Rho * M_Rho *
                         M_dphirho[j][n] * M_phirho[k][n] * M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->weight( n ) *
                         M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m ) *
                         normtheta * normtheta * M_Theta +
                         
                         M_phirho[j][n] * M_phirho[k][n] * normrho * normrho * M_Rho * M_quadruleRho->weight( n ) *   
                         M_dphitheta[j][m] * M_dphitheta[k][m] * normtheta * normtheta * M_Theta *
                         ( 1. / (2 * M_PI) ) * ( 1. / (2 * M_PI) ) *
                         M_quadruleTheta->weight( m )
                     );
            coeff += beta( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                           M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 1 ) *
                     M_dphirho[j][n] * M_phirho[k][n] * normrho * normrho * normrho * M_Rho * M_Rho * M_Rho *
                     M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->quadPointCoor( n, 0 ) *  
                     M_quadruleRho->weight( n ) *
                     normtheta * normtheta * M_Theta *
                     M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m );

            coeff += beta( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                           M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 2 ) *
                     normrho * normrho * M_Rho * M_Rho *
                     M_phirho[j][n] * M_phirho[k][n] * M_quadruleTheta->quadPointCoor( m, 0 ) * M_quadruleRho->weight( n ) *
                     normtheta * normtheta * M_Theta * ( 1. / ( 2 * M_PI ) ) *
                     M_dphitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m );

            coeff += sigma( 0, x, M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho,
                            M_quadruleTheta->quadPointCoor( m, 0 ) * M_Theta, 0 ) *
                     normrho * normrho * M_Rho * M_Rho * M_Rho *
                     M_phirho[j][n] * M_phirho[k][n] * 
                     M_quadruleRho->quadPointCoor( n, 0 ) * M_quadruleRho->quadPointCoor( n, 0 ) *  
                     M_quadruleRho->weight( n ) *
                     normtheta * normtheta * M_Theta *
                     M_phitheta[j][m] * M_phitheta[k][m] * M_quadruleTheta->weight( m );
        }
    }

    return coeff;
}


}
