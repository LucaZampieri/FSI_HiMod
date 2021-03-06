#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// #include <lifev/himod/modalbasis/NSModalSpaceCircular.hpp>
#include "include/NSModalSpaceCircular.hpp"
#include <iomanip>
namespace LifeV
{

void NSModalSpaceCircular::
addSliceBC( const std::string& BCx, const Real& mux, const Real& Chix,
		    const std::string& BCr, const Real& mur, const Real& Chir,
		    const std::string& BCtheta, const Real& mutheta, const Real& Chitheta )
{
   M_xGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( BCx );
   M_xGenbasisRhoTheta -> setRho( M_Rho );
   M_xGenbasisRhoTheta -> setTheta( M_Theta );
   M_xGenbasisRhoTheta -> setMu( mux );
   M_xGenbasisRhoTheta -> setChi( Chix );
   M_xGenbasisRhoTheta -> setNumberModes( M_mx );

   M_rGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( BCr );
   M_rGenbasisRhoTheta -> setRho( M_Rho );
   M_rGenbasisRhoTheta -> setTheta( M_Theta );
   M_rGenbasisRhoTheta -> setMu( mur );
   M_rGenbasisRhoTheta -> setChi( Chir );
   M_rGenbasisRhoTheta -> setNumberModes( M_mr );

   M_thetaGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( BCtheta );
   M_thetaGenbasisRhoTheta -> setRho( M_Rho );
   M_thetaGenbasisRhoTheta -> setTheta( M_Theta );
   M_thetaGenbasisRhoTheta -> setMu( mutheta );
   M_thetaGenbasisRhoTheta -> setChi( Chitheta );
   M_thetaGenbasisRhoTheta -> setNumberModes( M_mtheta );

   // We choose Neumann BC for the pressure because we don't want to enforce any value on the boundary
   // ( which we do if we enforce Dirichlet BCs ).
   M_pGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( "neu" );
   M_pGenbasisRhoTheta -> setRho( M_Rho );
   M_pGenbasisRhoTheta -> setTheta( M_Theta );
   M_pGenbasisRhoTheta -> setMu( 1. );
   M_pGenbasisRhoTheta -> setChi( 0. );
   M_pGenbasisRhoTheta -> setNumberModes( M_mp );

   eigensProvider();
}

void NSModalSpaceCircular::
evaluateBasis()
{
    M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta, M_xEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta, M_rEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta, M_thetaEigenvalues,
    											M_quadruleRho, M_quadruleTheta );
    M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta, M_pEigenvalues, M_quadruleRho, M_quadruleTheta );
}

void NSModalSpaceCircular::
evaluateBasis( const std::string& xPoly,
				const std::string& rPoly,
    			const std::string& thetaPoly,
    			const std::string& pPoly,
    			const int trigonometricBasis )
{
	M_xGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( xPoly );
	M_xGenbasisRhoTheta -> setRho( M_Rho );
	M_xGenbasisRhoTheta -> setTheta( M_Theta );
	M_xGenbasisRhoTheta -> setNumberModes( M_mx );
	M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta, M_quadruleRho, M_quadruleTheta );

	M_rGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( rPoly );
	M_rGenbasisRhoTheta -> setRho( M_Rho );
	M_rGenbasisRhoTheta -> setTheta( M_Theta );
	M_rGenbasisRhoTheta -> setNumberModes( M_mr );
	M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta, M_quadruleRho, M_quadruleTheta );

	M_thetaGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( thetaPoly );
	M_thetaGenbasisRhoTheta -> setRho( M_Rho );
	M_thetaGenbasisRhoTheta -> setTheta( M_Theta );
	M_thetaGenbasisRhoTheta -> setNumberModes( M_mtheta );
	M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta, M_quadruleRho, M_quadruleTheta );

	M_pGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( pPoly );
	M_pGenbasisRhoTheta -> setRho( M_Rho );
	M_pGenbasisRhoTheta -> setTheta( M_Theta );
	M_pGenbasisRhoTheta -> setNumberModes( M_mp );
	M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta, M_quadruleRho, M_quadruleTheta );

}


void NSModalSpaceCircular::
showMe() const
{
    std::cout << "---- MODAL SPACE CIRCULAR SHOWME ---" << std::endl;
    std::cout << "Rho = "	<<	M_Rho		<< std::endl;
    std::cout << "Theta = "	<<	M_Theta		<< std::endl;
    std::cout << "M_mx = "	<<	M_mx		<< std::endl;
    std::cout << "M_mr = "	<<	M_mr		<< std::endl;
    std::cout << "M_mtheta = "	<<	M_mtheta	<< std::endl;
    std::cout << "M_mp = "	<<	M_mp		<< std::endl;
    std::cout << "------------------------------------" << std::endl;
}

void NSModalSpaceCircular::
eigensProvider()
{

    M_xGenbasisRhoTheta -> computeEigenvalues();
    for(UInt i = 0; i != ( M_xGenbasisRhoTheta->eigenValues() ).size(); ++i)
    {
    	M_xEigenvalues.push_back( ( M_xGenbasisRhoTheta->eigenValues() )[i] );
    }

    M_rGenbasisRhoTheta -> computeEigenvalues();
    for(UInt i = 0; i != ( M_rGenbasisRhoTheta->eigenValues() ).size(); ++i)
    {
    	M_rEigenvalues.push_back( ( M_rGenbasisRhoTheta->eigenValues() )[i] );
    }

    M_thetaGenbasisRhoTheta -> computeEigenvalues();
    for(UInt i = 0; i != ( M_thetaGenbasisRhoTheta->eigenValues() ).size(); ++i)
    {
    	M_thetaEigenvalues.push_back( ( M_thetaGenbasisRhoTheta->eigenValues() )[i] );
    }

	M_pGenbasisRhoTheta -> computeEigenvalues();
    for(UInt i = 0; i != ( M_pGenbasisRhoTheta->eigenValues() ).size(); ++i)
    {
    	M_pEigenvalues.push_back( ( M_pGenbasisRhoTheta->eigenValues() )[i] );
    }

    return;
}


//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpaceCircular::
xFourierCoefficients( const function_Type& g, const Real& t ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mx, 0.0 );

    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >		evaluate_g;

    // first resize
    evaluate_g.resize( M_quadruleRho->nbQuadPt() );

    // second resize
    for ( UInt k = 0; k != evaluate_g.size(); ++k )
    {
        evaluate_g[k].resize( M_quadruleTheta->nbQuadPt() );
    }
    // evaluation of g
    for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, 0, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t , 0 , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mx; ++k )
    {
        // extraction of the sub-indices
//        UInt p_k = ( M_eigenvalues[k] - 1 );

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt(2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
                FourCoeff[k] += evaluate_g[n][m] *
                                M_xphirho[k][n] * normrho *
                                M_xphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights

    }
    return FourCoeff;
}

//Fourier Coefficients in the case g and R are x-dependent
std::vector<Real> NSModalSpaceCircular::
xFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mx, 0.0 );

    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >		evaluate_g;

    // first resize
    evaluate_g.resize( M_quadruleRho->nbQuadPt() );

    // second resize
    for ( UInt k = 0; k != evaluate_g.size(); ++k )
    {
        evaluate_g[k].resize( M_quadruleTheta->nbQuadPt() );
    }
    // evaluation of g
    for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t , x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mx; ++k )
    {
        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt(2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        {
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
            	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;

                FourCoeff[k] += evaluate_g[n][m] *
                                M_xphirho[k][n] * normrho *
                                M_xphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * /*M_map->fJacobian()( t , x , inverseRhat , inverseThetahat , 0 ) **/
                                M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights
            }
        }
    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpaceCircular::
rFourierCoefficients( const function_Type& g, const Real& t ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mr, 0.0 );

    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >		evaluate_g;

    // first resize
    evaluate_g.resize( M_quadruleRho->nbQuadPt() );

    // second resize
    for ( UInt k = 0; k != evaluate_g.size(); ++k )
    {
        evaluate_g[k].resize( M_quadruleTheta->nbQuadPt() );
    }
    // evaluation of g
    for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
        {
            Real inverseRhat = M_map->inverseRhat()( t, 0, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
       	    Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, 0 , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mr; ++k )
    {
        // extraction of the sub-indeces
//        UInt p_k = ( M_eigenvalues[k] - 1 );

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
                FourCoeff[k] += evaluate_g[n][m] *
                                M_rphirho[k][n] * normrho *
                                M_rphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights

    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-dependent
std::vector<Real> NSModalSpaceCircular::
rFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mr, 0.0 );

    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >		evaluate_g;

    // first resize
    evaluate_g.resize( M_quadruleRho->nbQuadPt() );

    // second resize
    for ( UInt k = 0; k != evaluate_g.size(); ++k )
    {
        evaluate_g[k].resize( M_quadruleTheta->nbQuadPt() );
    }
    // evaluation of g
    for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
        {
            Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
       	    Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mr; ++k )
    {
        // extraction of the sub-indeces
//        UInt p_k = ( M_eigenvalues[k] - 1 );

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
            	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;

                FourCoeff[k] += evaluate_g[n][m] *
                                M_rphirho[k][n] * normrho *
                                M_rphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * /*M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                                M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights
            }

    }
    return FourCoeff;
}


//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpaceCircular::
thetaFourierCoefficients( const function_Type& g, const Real& t ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mtheta, 0.0 );

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
            Real inverseRhat = M_map->inverseRhat()( t, 0, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
            Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, 0 , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mtheta; ++k )
    {
        // extraction of the sub-indeces
//        UInt p_k = ( M_eigenvalues[k] - 1 );

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
                FourCoeff[k] += evaluate_g[n][m] *
                                M_thetaphirho[k][n] * normrho *
                                M_thetaphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights

    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpaceCircular::
thetaFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mtheta, 0.0 );

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
            Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
            Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mtheta; ++k )
    {
        // extraction of the sub-indeces
//        UInt p_k = ( M_eigenvalues[k] - 1 );

        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
            	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;

                FourCoeff[k] += evaluate_g[n][m] *
                                M_thetaphirho[k][n] * normrho *
                                M_thetaphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * /*M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                                M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights
            }
    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpaceCircular::
pFourierCoefficients( const function_Type& g, const Real& t ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mp, 0.0 );

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
            Real inverseRhat = M_map->inverseRhat()( t, 0, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
            Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, 0 , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mp; ++k )
    {
        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
                FourCoeff[k] += evaluate_g[n][m] *
                                M_pphirho[k][n] * normrho *
                                M_pphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights

    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-dependent
std::vector<Real> NSModalSpaceCircular::
pFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>					FourCoeff( M_mp, 0.0 );

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
            Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor( m, 0 ), 0);
            Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;
            evaluate_g[n][m] = g( t, x, inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mp; ++k )
    {
        Real normrho = 1.0 / M_Rho;
        Real normtheta = 1.0 / sqrt( 2. * M_PI );

        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0 );
            	Real inverseThetahat = M_quadruleTheta -> quadPointCoor( m, 0 ) * M_Theta;

                FourCoeff[k] += evaluate_g[n][m] *
                                M_pphirho[k][n] * normrho *
                                M_pphitheta[k][m] * normtheta *
                                M_quadruleRho->quadPointCoor( n, 0 ) *
                                M_Theta * /*M_map->fJacobian()( t, x, inverseRhat , inverseThetahat , 0 ) **/
                                M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );    // weights
            }
    }
    return FourCoeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpaceCircular::
xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_xphirho[k][n] * normrho * M_xphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
        }

    return coeff;
}

Real NSModalSpaceCircular::
xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_xphirho[k][n] * normrho * M_xphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * /* M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                     M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
        }

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpaceCircular::
rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt(2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_rphirho[k][n] * normrho * M_rphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * M_map->Jacobian()[n][m] *
                     M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

Real NSModalSpaceCircular::
rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt(2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_rphirho[k][n] * normrho * M_rphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * /*M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                     M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpaceCircular::
thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
	        Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

Real NSModalSpaceCircular::
thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
	        Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * /*M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                     M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpaceCircular::
pFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_pphirho[k][n] * normrho * M_pphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * M_map->Jacobian()[n][m] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

Real NSModalSpaceCircular::
pFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k, const bool& b ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
        	Real inverseRhat = M_map->inverseRhat()( t, x, M_quadruleRho -> quadPointCoor( n, 0 ), M_quadruleTheta -> quadPointCoor (m, 0), 0);
        	Real inverseThetahat = M_quadruleTheta -> quadPointCoor (m, 0) * M_Theta;
            coeff += f( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_pphirho[k][n] * normrho * M_pphitheta[k][m] * normtheta *
                     M_quadruleRho->quadPointCoor( n , 0 ) *
                     M_Theta * /*M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) **/
                     M_quadruleRho->weight( n ) * M_quadruleTheta->weight( m );
		}

    return coeff;
}

// ------------------		Compute Methods		---------------------------
// Note: normrho and normtheta are intended to be normalization constant.
//		The Jacobian links the physical domain with the reference one.
//		The further 2 * M_PI is a scale factor due to the numerical domain of integration:
//				int_0^(2*pi)	d_theta = int_0^1	2*pi*d_thetah

Real
NSModalSpaceCircular::
compute_r000xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    	{
            coeff += (  2*nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
                        M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
                        M_quadruleRho->quadPointCoor( n, 0 ) *
                        M_quadruleRho->quadPointCoor( n, 0 ) +

                        2*nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
                        M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
                        M_quadruleRho->quadPointCoor( n, 0 ) +

                        nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
                        M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
                        M_quadruleRho->quadPointCoor( n, 0 ) *
                        M_quadruleRho->quadPointCoor( n, 0 ) +

                        nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
                        M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
                        M_quadruleRho->quadPointCoor( n, 0 ) +

                        nu * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
                        M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
                        M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta +

                        alpha *
                        M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
                        M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
                        M_quadruleRho->quadPointCoor( n, 0 ) *
                        M_quadruleRho->quadPointCoor( n, 0 )
                    ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000xxx( const UInt& k, const UInt& j, const UInt& s ) const
{
    Real coeff( 0 );
    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
    					M_map->Dr()[n][h] *
		   				M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    			M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000xxr( const UInt& k, const UInt& j, const UInt& s ) const
{
    Real coeff( 0 );
    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_rphirho[s][n] * normrho * M_rphitheta[s][h] * normtheta *
		    			M_map->Jr()[n][h] *
		   				M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    			M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000xxt( const UInt& k, const UInt& j, const UInt& s ) const
{
    Real coeff( 0 );
    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_thetaphirho[s][n] * normrho * M_thetaphitheta[s][h] * normtheta *
		    			M_map->Jtheta()[n][h] *
		   				M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
		    			M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100xxx( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
		   				M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    			M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}
	return coeff;
}

Real
NSModalSpaceCircular::
compute_r001xx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_map->Dr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010xx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_map->Dr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 )
        					* M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100xx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
			{
				coeff += ( 2 * nu * M_map->Dr()[n][h] * ( M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
									M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) +

									M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
									M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta )
									) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
			}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r101xx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r110xx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( 2 * nu * M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 )
		    					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    	}
    }

    return coeff;
}

// ************************		rr		*********************************

Real
NSModalSpaceCircular::
compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
	    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
	{
	    coeff += ( nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

	    				nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

	    				nu * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
	    					M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta +

	    				2 * nu * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
	    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta +

	    				2 * nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

	   					2 * nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

						alpha *
							M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 )
	    					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
	}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000rrx( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
    					M_map->Dr()[n][h] *
		   				M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    			M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000rrr( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_rphirho[s][n] * normrho * M_rphitheta[s][h] * normtheta *
		    			M_map->Jr()[n][h] *
		   				M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    			M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000rrt( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_thetaphirho[s][n] * normrho * M_thetaphitheta[s][h] * normtheta *
		    			M_map->Jtheta()[n][h] *
		   				M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
		    			M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100rrx( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
		   				M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    			M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}
	return coeff;
}

Real
NSModalSpaceCircular::
compute_r001rr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 )
        				* M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010rr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100rr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );


		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( nu * M_map->Dr()[n][h] * (
		    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 ) +

		    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta )
		    					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r101rr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r110rr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// *********************	tt		************************************

Real
NSModalSpaceCircular::
compute_r000tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
	{
	    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
	{
	    coeff += ( nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					nu * M_map->Dr()[n][h] * M_map->Dr()[n][h] *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					nu * M_map->Jr()[n][h] * M_map->Jr()[n][h] *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 ) -

   					nu * M_map->Jr()[n][h] * M_map->Jtheta()[n][h]*
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					2 * nu * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
	   					M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta +

   					nu * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta -

   					nu * M_map->Jr()[n][h] * M_map->Jtheta()[n][h] *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) -

  					nu * M_map->Jr()[n][h] * M_map->Jtheta()[n][h] *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta +

					alpha *
						M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 )
	   					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000ttx( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
    					M_map->Dr()[n][h] *
		   				M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		    			M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000ttr( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_rphirho[s][n] * normrho * M_rphitheta[s][h] * normtheta *
		    			M_map->Jr()[n][h] *
		   				M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		    			M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000ttt( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_thetaphirho[s][n] * normrho * M_thetaphitheta[s][h] * normtheta *
		    			M_map->Jtheta()[n][h] *
		   				M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
		    			M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	}
    }
    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100ttx( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_xphirho[s][n] * normrho * M_xphitheta[s][h] * normtheta *
		   				M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		    			M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}
	return coeff;
}

Real
NSModalSpaceCircular::
compute_r001tt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] *
       					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010tt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100tt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( nu * M_map->Dr()[n][h] * (
		   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		   					M_quadruleRho->quadPointCoor( n, 0 ) +

		   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta )
		   					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
		   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r101tt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r110tt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
        				M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// ***********************		xr 		**********************************

Real
NSModalSpaceCircular::
compute_r000xr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jr()[n][h] * M_map->Dr()[n][h] * ( M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta  )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r001xr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010xr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// ***************************		xt		*********************************************

Real
NSModalSpaceCircular::
compute_r000xt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_map->Jtheta()[n][h] * ( M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r001xt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jtheta()[n][h] * M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010xt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jtheta()[n][h] * M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
        				M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
        				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// ****************************		rx		********************************************

Real
NSModalSpaceCircular::
compute_r000rx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_map->Jr()[n][h] * ( M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
       					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100rx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( nu * M_map->Jr()[n][h] * ( M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		   					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 ) +

		   					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		   					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta )
		   					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
		   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    return coeff;
}
// **************************		rt		*************************************************

Real
NSModalSpaceCircular::
compute_r000rt( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * ( M_map->Jr()[n][h] * M_map->Jtheta()[n][h] * M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
       					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_map->Jtheta()[n][h] * M_map->Jr()[n][h] *
       					M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta -

       					M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
       					M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta +

   					2 * M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
       					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000rtt( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( M_thetaphirho[s][n] * normrho * M_thetaphitheta[s][h] * normtheta *
    					M_map->Jtheta()[n][h] *
		   				M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    			M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}
	return coeff;
}

// ****************************		tx		***************************************************

Real
NSModalSpaceCircular::
compute_r000tx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Dr()[n][h] * M_map->Jtheta()[n][h] * M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100tx( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_map->Jtheta()[n][h] * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// ******************************		tr		*************************************************

Real
NSModalSpaceCircular::
compute_r000tr( const UInt& k, const UInt& j, const Real& nu ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * ( M_map->Jr()[n][h] * M_map->Jtheta()[n][h] * M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) -

       					M_map->Jtheta()[n][h] * M_map->Jtheta()[n][h] *
       					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta +

   					2 * M_map->Jtheta()[n][h] *
       					M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta )
       					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000trt( const UInt& k, const UInt& j, const UInt& s ) const
{
	Real coeff( 0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	{
    		coeff += ( - M_thetaphirho[s][n] * normrho * M_thetaphitheta[s][h] * normtheta *
    					M_map->Jtheta()[n][h] *
		   				M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		    			M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    			M_quadruleRho->quadPointCoor( n, 0 )
		    		) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}
	return coeff;
}

// ********************		xp		***************************************

Real
NSModalSpaceCircular::
compute_r000xp( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_map->Dr()[n][h] * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r100xp( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        			M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
        			) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
        			M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000rp( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_map->Jtheta()[n][h] *
        			M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
   					M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta +

					M_map->Jr()[n][h] *
					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
   					M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
   					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000tp( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( M_map->Jtheta()[n][h] * M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
		    			M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta
		    			) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000px( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_map->Dr()[n][h] * ( M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )

  					- M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta )
  					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
  					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r001px( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
   				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
   				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000pp( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta *
   		   M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta *
                   M_Rho * M_quadruleRho->quadPointCoor( n, 0 ) // this is the factor to avoid singularity... do we need it?
   		 ) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
   		 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r010px( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    	{
    	    coeff += ( - M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
    					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
    				) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
    				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    	}
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000pr( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_map->Jr()[n][h] * M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )

   					- M_map->Jr()[n][h] * M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta

  					- M_map->Jtheta()[n][h] * M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
   					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

Real
NSModalSpaceCircular::
compute_r000pt( const UInt& k, const UInt& j ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_map->Jtheta()[n][h] * M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
  					) * M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
  					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    return coeff;
}

// --------------------------------- OVERLOAD FOR X-DEPENDENCE ------------------------------------------------------
void
NSModalSpaceCircular::
compute_r000xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000xx ) const
{
    Real coeff1( 0.0 );
    Real coeff2( 0.0 );
    Real coeff3( 0.0 );
    Real coeff4( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    	{
		    coeff1 += M_map->Dr()[n][h] * M_map->Dr()[n][h] *
		              (  2 * nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->quadPointCoor( n, 0 ) +

						2 * nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 )
		    		  ) * M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

		   	coeff2 += (	nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->quadPointCoor( n, 0 ) +

		   					nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 )
		    		  ) * M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

		   	coeff3 +=	nu *
		    					M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta *
		    					M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

		   	coeff4 +=	1* //alpha *
		   						M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    	} // n for-loop
    } // h for-loop

    R000xx *= 0.0;
    R000xx += M_map->xDr();
    R000xx *= M_map->xDr();
    R000xx *= M_map->xJacobian();
    R000xx *= coeff1;

    vector_Type tmp( M_map->xJr() );
    tmp *= M_map->xJr();
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000xx += tmp;

    tmp = M_map->xJtheta();
    tmp *= M_map->xJtheta();
    tmp *= M_map->xJacobian();
    tmp *= coeff3;

    R000xx += tmp;

    tmp = M_map->xJacobian();
    tmp *= coeff4;

    R000xx += tmp;

    R000xx *= M_map->xR();
/*    R000xx = M_map->xR() * (
                            M_map->xDr() * M_map->xDr() * M_map->xJacobian() * coeff1 +
                            M_map->xJr() * M_map->xJr() * M_map->xJacobian() * coeff2 +
                            M_map->xJtheta() * M_map->xJtheta() * M_map->xJacobian() * coeff3 +
                            M_map->xJacobian() * coeff4
                            );*/
    return;
}

void
NSModalSpaceCircular::
compute_r001xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_map->Dr()[n][h] *
        					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R001xx = M_map->xDr();
    R001xx *= M_map->xJacobian();
    R001xx *= coeff;

    R001xx *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_map->Dr()[n][h]
        					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 )
        					* M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R010xx = M_map->xDr();
    R010xx *= M_map->xJacobian();
    R010xx *= coeff;

    R010xx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100xx ) const
{
    Real coeff( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
			{
				coeff +=  2 * nu * ( M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
									M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) +

									M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
									M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta
									) * M_map->Dr()[n][h] *
									M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
			}

	R100xx = M_map->xDr();
	R100xx *= M_map->xJacobian();
	R100xx *= coeff;

	R100xx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r101xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101xx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( 2 * nu * M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R101xx = M_map->xJacobian();
    R101xx *= coeff;

    R101xx *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r110xx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110xx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( 2 * nu * M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
		    					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 )
		    					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    	}
    }

    R110xx = M_map->xJacobian();
    R110xx *= coeff;

    R110xx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000rr ) const
{
    Real coeff1( 0.0 );
    Real coeff2( 0.0 );
    Real coeff3( 0.0 );
    Real coeff4( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
	    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
	{
	    coeff1 += ( nu *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

	    			nu *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 )
	    					) * M_map->Dr()[n][h] * M_map->Dr()[n][h]
	    					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

	    coeff2 += ( nu *
	    					M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta +

	    		    2 * nu *
	    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta
	    					)
	    					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

	    coeff3 += ( 2 * nu *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 ) +

	   				2 * nu *
	    					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 )
	    					)
	    					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

		coeff4 += 	alpha *
							M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
	    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_quadruleRho->quadPointCoor( n, 0 ) *
	    					M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
	}

    R000rr = M_map->xDr();
    R000rr *= M_map->xDr();
    R000rr *= M_map->xJacobian();
    R000rr *= coeff1;

    vector_Type tmp( M_map->xJtheta() );
    tmp *= M_map->xJtheta();
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000rr += tmp;

    tmp = M_map->xJr();
    tmp *= M_map->xJr();
    tmp *= M_map->xJacobian();
    tmp *= coeff3;

    R000rr += tmp;

    tmp *= M_map->xJacobian();
    tmp *= coeff4;

    R000rr += tmp;

    R000rr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r001rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001rr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_map->Dr()[n][h]
        				* M_Theta * M_quadruleRho->quadPointCoor( n, 0 )
        				* M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    R001rr = M_map->xDr();
    R001rr *= M_map->xJacobian();
    R001rr *= coeff;

    R001rr *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010rr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        					M_quadruleRho->quadPointCoor( n, 0 )
        					) * M_map->Dr()[n][h]
        					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R010rr = M_map->xDr();
    R010rr *= M_map->xJacobian();
    R010rr *= coeff;

    R010rr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100rr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );


		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += nu * (
		    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
		    					M_quadruleRho->quadPointCoor( n, 0 ) +

		    					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		    					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta
		    					) * M_map->Dr()[n][h]
		    					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
		    					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    R100rr = M_map->xDr();
    R100rr *= M_map->xJacobian();
    R100rr *= coeff;

    R100rr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r101rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101rr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R101rr = M_map->xJacobian();
    R101rr *= coeff;

    R101rr *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r110rr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110rr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
        				M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R110rr = M_map->xJacobian();
    R110rr *= coeff;

    R110rr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
										vector_Type& R000tt ) const
{
    Real coeff1( 0.0 );
    Real coeff2( 0.0 );
    Real coeff3( 0.0 );
    Real coeff4( 0.0 );
    Real coeff5( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
	for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
	{
	    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
	{
	    coeff1 += ( nu *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					nu *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 )
	   				) * M_map->Dr()[n][h] * M_map->Dr()[n][h]
	   				* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

   		coeff2 += (	nu *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) +

   					nu *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 )
	   				)
	   					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

   		coeff3 += (	- nu *
	   					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 )

	   				- nu *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 )

  					- nu *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
	   				)
	   					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

   		coeff4 += (	2 * nu *
	   					M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta +

   					nu *
	   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
	   				)
	   					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

		coeff5 += alpha *
						M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
	   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
	   					M_quadruleRho->quadPointCoor( n, 0 ) *
	   					M_quadruleRho->quadPointCoor( n, 0 )
	   					 * M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
		}
	}

	R000tt = M_map->xDr();
	R000tt *= M_map->xDr();
	R000tt *= M_map->xJacobian();
	R000tt *= coeff1;

	vector_Type tmp( M_map->xJr() );
	tmp *= M_map->xJr();
	tmp *= M_map->xJacobian();
	tmp *= coeff2;

	R000tt += tmp;

	tmp = M_map->xJr();
	tmp *= M_map->xJtheta();
	tmp *= M_map->xJacobian();
	tmp *= coeff3;

	R000tt += tmp;

	tmp = M_map->xJtheta();
	tmp *= M_map->xJtheta();
	tmp *= M_map->xJacobian();
	tmp *= coeff4;

	R000tt += tmp;

	tmp = M_map->xJacobian();
	tmp *= coeff5;

	R000tt += tmp;

	R000tt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r001tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001tt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu *
       					M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_map->Dr()[n][h] *
       					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R001tt = M_map->xDr();
    R001tt *= M_map->xJacobian();
    R001tt *= coeff;

    R001tt *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010tt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu *
                        M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_map->Dr()[n][h] *
       					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R010tt = M_map->xDr();
    R010tt *= M_map->xJacobian();
    R010tt *= coeff;

    R010tt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100tt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += nu * (
		   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		   					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
		   					M_quadruleRho->quadPointCoor( n, 0 ) +

		   					M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
		   					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
		   					) * M_map->Dr()[n][h] *
		   					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
		   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    R100tt = M_map->xDr();
    R100tt *= M_map->xJacobian();
    R100tt *= coeff;

    R100tt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r101tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R101tt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R101tt = M_map->xJacobian();
    R101tt *= coeff;

    R101tt *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r110tt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R110tt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
        				M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
        				M_quadruleRho->quadPointCoor( n, 0 )
        				) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R110tt = M_map->xJacobian();
    R110tt *= coeff;

    R110tt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000xr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += nu * (
                        M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta
       					) * M_map->Dr()[n][h]
       					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000xr = M_map->xDr();
    R000xr *= M_map->xJr();
    R000xr *= M_map->xJacobian();
    R000xr *= coeff;

    R000xr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r001xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu *
                        M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					)
       					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R001xr = M_map->xJr();
    R001xr *= M_map->xJacobian();
    R001xr *= coeff;

    R001xr *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010xr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xr ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu *
                        M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					)
       					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R010xr = M_map->xJr();
    R010xr *= M_map->xJacobian();
    R010xr *= coeff;

    R010xr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000xt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

	// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += nu * (
                        M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
       					) * M_map->Dr()[n][h]
       					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    R000xt = M_map->xDr();
    R000xt *= M_map->xJtheta();
    R000xt *= M_map->xJacobian();
    R000xt *= coeff;

    R000xt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r001xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R001xt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
       					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
       					)
       					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R001xt = M_map->xJtheta();
    R001xt *= M_map->xJacobian();
    R001xt *= coeff;

    R001xt *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010xt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R010xt ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_xphirho[k][n] * normrho * M_xdphitheta[k][h] * normtheta *
        				M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
        				)
        				* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R010xt = M_map->xJtheta();
    R010xt *= M_map->xJacobian();
    R010xt *= coeff;

    R010xt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000rx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000rx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += nu * ( M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
       					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 ) +

       					M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta
       					) * M_map->Dr()[n][h]
       					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000rx = M_map->xDr();
    R000rx *= M_map->xJr();
    R000rx *= M_map->xJacobian();
    R000rx *= coeff;

    R000rx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100rx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100rx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
		    for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += nu * ( M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		   					M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 ) +

		   					M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
		   					M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta
		   					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
		   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

	R100rx = M_map->xJr();
	R100rx *= M_map->xJacobian();
	R100rx *= coeff;

	R100rx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000rt( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000rt ) const
{
    Real coeff1( 0.0 );
    Real coeff2( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
        {
            coeff1 += nu * (
                            M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
           					M_thetadphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
           					M_quadruleRho->quadPointCoor( n, 0 ) +

           					M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
           					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta
           					)
           					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

           	coeff2 += nu * (
           	                - M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
           					M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta +

       					2 * M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
           					M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta
           					)
           					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
        }

    R000rt = M_map->xJr();
    R000rt *= M_map->xJtheta();
    R000rt *= M_map->xJacobian();
    R000rt *= coeff1;

    vector_Type tmp( M_map->xJtheta() );
    tmp *= M_map->xJtheta();
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000rt += tmp;

    R000rt *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000tx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000tx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta
       					) * M_map->Dr()[n][h] *
       					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000tx = M_map->xDr();
    R000tx *= M_map->xJtheta();
    R000tx *= M_map->xJacobian();
    R000tx *= coeff;

    R000tx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100tx( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R100tx ) const
{
    Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( nu * M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_xphirho[j][n] * normrho * M_xdphitheta[j][h] * normtheta
       					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
       					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R100tx = M_map->xJtheta();
    R100tx *= M_map->xJacobian();
    R100tx *= coeff;

    R100tx *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000tr( const UInt& k, const UInt& j, const Real& nu,
										vector_Type& R000tr ) const
{
    Real coeff1( 0.0 );
    Real coeff2( 0.0 );
    Real coeff3( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    // Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff1 +=  nu *
                        M_thetadphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta *
       					M_quadruleRho->quadPointCoor( n, 0 )
       					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

       	coeff2 += nu * (
       	                - M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta
       					)
       					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

   		coeff3 += 2 * nu *
       					M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
       					M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta
       					* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
    }

    R000tr = M_map->xJr();
    R000tr *= M_map->xJtheta();
    R000tr *= M_map->xJacobian();
    R000tr *= coeff1;

    vector_Type tmp( M_map->xJtheta() );
    tmp *= M_map->xJtheta();
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000tr += tmp;

    tmp = M_map->xJtheta();
    tmp *= M_map->xJacobian();
    tmp *= coeff3;

    R000tr += tmp;

    R000tr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000xp( const UInt& k, const UInt& j,
										vector_Type& R000xp ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_xdphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        		    M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta *
        		    M_quadruleRho->quadPointCoor( n, 0 )
				    ) * M_map->Dr()[n][h]
				    * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        			M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000xp = M_map->xDr();
    R000xp *= M_map->xJacobian();
    R000xp *= coeff;

    R000xp *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r100xp( const UInt& k, const UInt& j,
										vector_Type& R100xp ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
        			M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
        			) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
        			M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R100xp = M_map->xJacobian();
    R100xp *= coeff;

    R100xp *= M_map->xR();

    return ;
}

 void
NSModalSpaceCircular::
compute_r000rp( const UInt& k, const UInt& j,
										vector_Type& R000rp ) const
{
	Real coeff1( 0.0 );
	Real coeff2( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff1 +=  M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
   					M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta
   					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;

		coeff2 +=  M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
   					M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta *
   					M_quadruleRho->quadPointCoor( n, 0 )
   					 * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000rp = M_map->xJtheta();
    R000rp *= M_map->xJacobian();
    R000rp *= coeff1;

    vector_Type tmp( M_map->xJr() );
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000rp += tmp;

    R000rp *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000tp( const UInt& k, const UInt& j,
										vector_Type& R000tp ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
		    coeff += ( M_thetaphirho[k][n] * normrho * M_thetadphitheta[k][h] * normtheta *
		    			M_pphirho[j][n] * normrho * M_pphitheta[j][h] * normtheta
		    			) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
		    			M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
		}

    R000tp = M_map->xJtheta();
    R000tp *= M_map->xJacobian();
    R000tp *= coeff;

    R000tp *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000px( const UInt& k, const UInt& j,
										vector_Type& R000px ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += - ( M_xdphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )

  					+ M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
  					) * M_map->Dr()[n][h] *
  					M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
  					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000px = M_map->xDr();
    R000px *= M_map->xJacobian();
    R000px *= coeff;

    R000px *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r001px( const UInt& k, const UInt& j,
										vector_Type& R001px ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
   				) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
   				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R001px = M_map->xJacobian();
    R001px *= coeff;

    R001px *= M_map->xdR();

    return ;
}

void
NSModalSpaceCircular::
compute_r010px( const UInt& k, const UInt& j,
										vector_Type& R010px ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    {
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    	{
    	    coeff += ( - M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
    					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )
    				) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
    				M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    	}
    }

    R010px = M_map->xJacobian();
    R010px *= coeff;

    R010px *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000pr( const UInt& k, const UInt& j,
										vector_Type& R000pr ) const
{
	Real coeff1( 0.0 );
	Real coeff2( 0.0 );

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff1 += ( - M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta * M_quadruleRho->quadPointCoor( n, 0 )

   					- M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
   					)
   					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;

  		coeff2 += - M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
   					* M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
   					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000pr = M_map->xJr();
    R000pr *= M_map->xJacobian();
    R000pr *= coeff1;

    vector_Type tmp( M_map->xJtheta() );
    tmp *= M_map->xJacobian();
    tmp *= coeff2;

    R000pr += tmp;

    R000pr *= M_map->xR();

    return ;
}

void
NSModalSpaceCircular::
compute_r000pt( const UInt& k, const UInt& j,
										vector_Type& R000pt ) const
{
	Real coeff = 0.0;

    // If you have a function which is normalized in (0,1) with respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
    {
        coeff += ( - M_thetaphirho[j][n] * normrho * M_thetadphitheta[j][h] * normtheta *
   					M_pphirho[k][n] * normrho * M_pphitheta[k][h] * normtheta
  					) * M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
  					M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h ) ;
    }

    R000pt = M_map->xJtheta();
    R000pt *= M_map->xJacobian();
    R000pt *= coeff;

    R000pt *= M_map->xR();

    return ;
}

// ---------------------------------------------------------------------------------------------------------
//PHI
Real NSModalSpaceCircular::
compute_Phix( const UInt& k ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
		    	// note: r is taken into account 2 times: 1) multiplying of the beginning equation; 2) inner product.
				coeff += M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
					     M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

    return coeff;
}

void NSModalSpaceCircular::
compute_Phix( const UInt& k,
    			vector_Type& Phix ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
		    	// note: r is taken into account 2 times: 1) multiplying of the beginning equation; 2) inner product.
				coeff += M_xphirho[k][n] * normrho * M_xphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * /* M_Rho below */
					     M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

	Phix = M_map->xJacobian();
	Phix *= coeff;

	Phix *= M_map->xR();

    return ;
}

Real NSModalSpaceCircular::
compute_Phir( const UInt& k ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
				coeff += M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
					     M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

    return coeff;
}

void NSModalSpaceCircular::
compute_Phir( const UInt& k,
				vector_Type& Phir ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
				coeff += M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * /* M_Rho below */
					     M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

    Phir = M_map->xJacobian();
    Phir *= coeff;

    Phir *= M_map->xR();

    return ;
}

Real NSModalSpaceCircular::
compute_Phitheta( const UInt& k ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
				coeff += M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * M_Rho *
					     M_Theta * M_map->Jacobian()[n][h] * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

    return coeff;
}

void NSModalSpaceCircular::
compute_Phitheta( const UInt& k,
					vector_Type& Phitheta ) const
{
    Real coeff = 0;

    // If you have a function which is normalized in (0,1) in respect to L2 norm.
    // You have to consider other constant to obtain normalization in (0,L)
    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt( 2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
    	for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		    {
				coeff += M_thetaphirho[k][n] * normrho * M_thetaphitheta[k][h] * normtheta *
					     M_quadruleRho->quadPointCoor( n, 0 ) * /* M_Rho below */
					     M_Theta * M_quadruleRho->quadPointCoor( n, 0 ) *
						 M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
			}

    Phitheta = M_map->xJacobian();
    Phitheta *= coeff;

    Phitheta *= M_map->xR();

    return ;
}

// ------------------------------ New FSI --------------------------------------
void NSModalSpaceCircular::
addSliceBCFSI( const std::string& BCx, const Real& mux, const Real& Chix,
		    const std::string& BCtheta, const Real& mutheta, const Real& Chitheta )
{
   M_xGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( BCx );
   M_xGenbasisRhoTheta -> setRho( M_Rho );
   M_xGenbasisRhoTheta -> setTheta( M_Theta );
   M_xGenbasisRhoTheta -> setMu( mux );
   M_xGenbasisRhoTheta -> setChi( Chix );
   M_xGenbasisRhoTheta -> setNumberModes( M_mx );

	 // We choose Neumann BC for the pressure because we don't want to enforce any value on the boundary
   // ( which we do if we enforce Dirichlet BCs ).
   M_rGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( "neu" );
   M_rGenbasisRhoTheta -> setRho( M_Rho );
   M_rGenbasisRhoTheta -> setTheta( M_Theta );
   M_rGenbasisRhoTheta -> setMu( 1. );
   M_rGenbasisRhoTheta -> setChi( 0. );
   M_rGenbasisRhoTheta -> setNumberModes( M_mr );

   M_thetaGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( BCtheta );
   M_thetaGenbasisRhoTheta -> setRho( M_Rho );
   M_thetaGenbasisRhoTheta -> setTheta( M_Theta );
   M_thetaGenbasisRhoTheta -> setMu( mutheta );
   M_thetaGenbasisRhoTheta -> setChi( Chitheta );
   M_thetaGenbasisRhoTheta -> setNumberModes( M_mtheta );

   // We choose Neumann BC for the pressure because we don't want to enforce any value on the boundary
   // ( which we do if we enforce Dirichlet BCs ).
   M_pGenbasisRhoTheta = EducatedBasisFactory::instance().createObject( "neu" );
   M_pGenbasisRhoTheta -> setRho( M_Rho );
   M_pGenbasisRhoTheta -> setTheta( M_Theta );
   M_pGenbasisRhoTheta -> setMu( 1. );
   M_pGenbasisRhoTheta -> setChi( 0. );
   M_pGenbasisRhoTheta -> setNumberModes( M_mp );

   eigensProvider();
}

void NSModalSpaceCircular::
evaluateBasisFSI()
{
    M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta, M_xEigenvalues, M_quadruleRho, M_quadruleTheta );
		//M_rGenbasisRhoTheta -> evaluateBasis( M_rphirhoWall, M_rdphirhoWall, M_rphitheta, M_rdphitheta, M_rEigenvalues, M_quadruleRhoWall, M_quadruleTheta );
    M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta, M_rEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta, M_thetaEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta, M_pEigenvalues, M_quadruleRho, M_quadruleTheta );
}

void NSModalSpaceCircular::
evaluateBasisFSI( const std::string& xPoly, const std::string& rPoly, const std::string& thetaPoly, const std::string& pPoly, const int trigonometricBasis )
{
	M_xGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( xPoly );
	M_xGenbasisRhoTheta -> setRho( M_Rho );
	M_xGenbasisRhoTheta -> setTheta( M_Theta );
	M_xGenbasisRhoTheta -> setNumberModes( M_mx );
	M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta, M_quadruleRho, M_quadruleTheta );

	M_rGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( rPoly );
	M_rGenbasisRhoTheta -> setRho( M_Rho );
	M_rGenbasisRhoTheta -> setTheta( M_Theta );
	M_rGenbasisRhoTheta -> setNumberModes( M_mr );
	M_rGenbasisRhoTheta -> evaluateBasis( M_rphirhoWall, M_rdphirhoWall, M_rphitheta, M_rdphitheta, M_quadruleRhoWall, M_quadruleTheta );
	M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta, M_quadruleRho, M_quadruleTheta );

	M_thetaGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( thetaPoly );
	M_thetaGenbasisRhoTheta -> setRho( M_Rho );
	M_thetaGenbasisRhoTheta -> setTheta( M_Theta );
	M_thetaGenbasisRhoTheta -> setNumberModes( M_mtheta );
	M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta, M_quadruleRho, M_quadruleTheta );

	M_pGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( pPoly );
	M_pGenbasisRhoTheta -> setRho( M_Rho );
	M_pGenbasisRhoTheta -> setTheta( M_Theta );
	M_pGenbasisRhoTheta -> setNumberModes( M_mp );
	M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta, M_quadruleRho, M_quadruleTheta );
}

void NSModalSpaceCircular::compute_r000rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha, const Real& rho_s, const Real& h_s, const Real& e, vector_Type& R000rr ) const
{

		Real coeff1( 0.0 );
		Real coeff2( 0.0 );
		Real coeff3( 0.0 );
		Real coeff4( 0.0 );
		Real coeff5( 0.0 );

		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
		for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
				for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
				{
						coeff1 += ( nu *
							M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_quadruleRho->quadPointCoor( n, 0 ) +
							nu *
							M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 )
							) * M_map->Dr()[n][h] * M_map->Dr()[n][h]
							* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

						coeff2 += ( nu *
							M_rphirho[k][n] * normrho * M_rdphitheta[k][h] * normtheta *
							M_rphirho[j][n] * normrho * M_rdphitheta[j][h] * normtheta +
							2 * nu *
							M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta
							)
							* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

						coeff3 += ( 2 * nu *
							M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rdphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_quadruleRho->quadPointCoor( n, 0 ) +
							2 * nu *
							M_rdphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 )
							)
							* M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

						coeff4 += alpha *
							M_rphirho[k][n] * normrho * M_rphitheta[k][h] * normtheta *
							M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

				}
				/*
				coeff5 += (rho_s*h_s*alpha + e/alpha) *
					M_rphirhoWall[k][M_quadruleRhoWall->nbQuadPt()-1] * normrho * M_rphitheta[k][h] * normtheta *
					M_rphirhoWall[j][M_quadruleRhoWall->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
					M_Theta * M_quadruleTheta->weight( h );
				*/
				coeff5 += (rho_s*h_s*alpha + e/alpha) *
					M_rphirho[k][M_quadruleRho->nbQuadPt()-1] * normrho * M_rphitheta[k][h] * normtheta *
					M_rphirho[j][M_quadruleRho->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
					M_Theta * M_quadruleTheta->weight( h );
		}

		R000rr = M_map->xDr();
		R000rr *= M_map->xDr();
		R000rr *= M_map->xJacobian();
		R000rr *= coeff1;

		vector_Type tmp( M_map->xJtheta() );
		tmp *= M_map->xJtheta();
		tmp *= M_map->xJacobian();
		tmp *= coeff2;

		R000rr += tmp;

		tmp = M_map->xJr();
		tmp *= M_map->xJr();
		tmp *= M_map->xJacobian();
		tmp *= coeff3;

		R000rr += tmp;

		tmp = M_map->xJacobian();
		tmp *= coeff4;

		R000rr += tmp;

		tmp = M_map->xJacobianWall();
		tmp *= coeff5;

		R000rr += tmp;

		R000rr *= M_map->xR();

}

void NSModalSpaceCircular::compute_r00x( const UInt& j, const vector_Type& f, const vector_Type& u_old, const Real& alpha, vector_Type& R00x ) const
{

		Real coeff1;
		Real coeff2;

		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for (UInt m = 0; m < R00x.size(); ++m)
		{
				coeff1 = 0;
				coeff2 = 0;

				// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
				for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
				{
						for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
						{
								coeff1 += f[xcoord2index(m,n,h)] *
									M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

								coeff2 += alpha * u_old[xcoord2index(m,n,h)] *
									M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
						}
				}
				R00x[m] = (coeff1 + coeff2)*M_map->xJacobian()[m];
		}

		R00x *= M_map->xR();

}

void NSModalSpaceCircular::compute_r00r( const UInt& j, const vector_Type& f, const vector_Type& u_old, const vector_Type& urWall_old, const vector_Type& etar_old, const Real& alpha, const Real& rho_s, const Real& h_s, const Real& e, vector_Type& R00r ) const
{

		Real coeff1;
		Real coeff2;
		Real coeff3;
		Real coeff4;

		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho   = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for (UInt m = 0; m < R00r.size(); ++m)
		{
				coeff1 = 0;
				coeff2 = 0;
				coeff3 = 0;
				coeff4 = 0;

				// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
				for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
				{
						for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
						{
								coeff1 += f[rcoord2index(m,n,h)] *
									M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

								coeff2 += alpha * u_old[rcoord2index(m,n,h)] *
									M_rphirho[j][n] * normrho * M_rphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
						}
						/* using urwall
						coeff3 += rho_s*h_s*alpha * urWall_old[coord2indexWall(m,h)] *
							M_rphirhoWall[j][M_quadruleRhoWall->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
							M_Theta * M_quadruleTheta->weight( h );

						coeff4 += e * etar_old[m] *
							M_rphirhoWall[j][M_quadruleRhoWall->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
							M_Theta * M_quadruleTheta->weight( h );
							*/
							coeff3 += rho_s*h_s*alpha * u_old[rcoord2index(m,M_quadruleRho->nbQuadPt()-1,h)] *
								M_rphirho[j][M_quadruleRho->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
								M_Theta * M_quadruleTheta->weight( h );

							coeff4 += e * etar_old[m] *
								M_rphirho[j][M_quadruleRho->nbQuadPt()-1] * normrho * M_rphitheta[j][h] * normtheta *
								M_Theta * M_quadruleTheta->weight( h );
				}
				R00r[m] = ((coeff1 + coeff2)*M_map->xJacobian()[m] + (coeff3 - coeff4)*M_map->xJacobianWall()[m]);
		}

		R00r *= M_map->xR();

}

void NSModalSpaceCircular::compute_r00t( const UInt& j, const vector_Type& f, const vector_Type& u_old, const Real& alpha, vector_Type& R00t ) const
{

		Real coeff1;
		Real coeff2;

		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		for (UInt m = 0; m < R00t.size(); ++m)
		{
				coeff1 = 0;
				coeff2 = 0;

				// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
				for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
				{
						for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
						{
								coeff1 += f[thetacoord2index(m,n,h)] *
									M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );

								coeff2 += alpha * u_old[thetacoord2index(m,n,h)] *
									M_thetaphirho[j][n] * normrho * M_thetaphitheta[j][h] * normtheta *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_quadruleRho->quadPointCoor( n, 0 ) *
									M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
						}
				}
				R00t[m] = (coeff1 + coeff2)*M_map->xJacobian()[m];
		}

		R00t *= M_map->xR();

}

void NSModalSpaceCircular::compute_b0( const UInt& j, const Real& p1, Real& B0 ) const
{

		Real coeff1( 0.0 );


		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
		for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
				for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
				{
						coeff1 += p1 *
							M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_quadruleRho->quadPointCoor( n, 0 ) *
							M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
				}
		}

		B0 = coeff1*M_map->xJacobian()[0]*M_map->xR()[0];

}


void NSModalSpaceCircular::compute_bL( const UInt& j, const Real& p2, Real& BL ) const
{

		Real coeff1( 0.0 );

		// If you have a function which is normalized in (0,1) with respect to L2 norm.
		// You have to consider other constant to obtain normalization in (0,L)
		Real normrho = 1.0 / M_Rho;
		Real normtheta = 1.0 / sqrt( 2. * M_PI );

		// Note: quadPointCoord taken into account INSIDE round bracket in order to cancel the denominator and avoid singularity
		for ( UInt h = 0; h != M_quadruleTheta->nbQuadPt(); ++h )
		{
				for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
				{
						coeff1 += p2 *
							M_xphirho[j][n] * normrho * M_xphitheta[j][h] * normtheta * // phi_hat x,j
							M_quadruleRho->quadPointCoor( n, 0 ) * // r_hat
							M_quadruleRho->quadPointCoor( n, 0 ) * // r_hat
							M_Theta * M_quadruleRho->weight( n ) * M_quadruleTheta->weight( h );
				}
		}

		BL = coeff1*M_map->xJacobian()[M_map->xJacobian().size()-1]*M_map->xR()[M_map->xR().size()-1]; // verify if .back() exists (Andrea)

}

// -----------------------------------------------------------------------------

}
