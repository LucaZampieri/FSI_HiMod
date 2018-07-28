#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/himod/modalbasis/NSModalSpacePipe.hpp>
#include <iomanip>
namespace LifeV
{

void NSModalSpacePipe::
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

void NSModalSpacePipe::
evaluateBasis()
{
    M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta,
                                          M_xEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta,
                                          M_rEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta,
                                              M_thetaEigenvalues, M_quadruleRho, M_quadruleTheta );
    M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta,
                                          M_pEigenvalues, M_quadruleRho, M_quadruleTheta );
}

void NSModalSpacePipe::
evaluateBasis( const std::string& xPoly,
                const std::string& rPoly,
                const std::string& thetaPoly,
                const std::string& pPoly )
{
    M_xGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( xPoly );
    M_xGenbasisRhoTheta -> setRho( M_Rho );
    M_xGenbasisRhoTheta -> setTheta( M_Theta );
    M_xGenbasisRhoTheta -> setNumberModes( M_mx );
    M_xGenbasisRhoTheta -> evaluateBasis( M_xphirho, M_xdphirho, M_xphitheta, M_xdphitheta,
                                          M_quadruleRho, M_quadruleTheta );

    M_rGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( rPoly );
    M_rGenbasisRhoTheta -> setRho( M_Rho );
    M_rGenbasisRhoTheta -> setTheta( M_Theta );
    M_rGenbasisRhoTheta -> setNumberModes( M_mr );
    M_rGenbasisRhoTheta -> evaluateBasis( M_rphirho, M_rdphirho, M_rphitheta, M_rdphitheta,
                                          M_quadruleRho, M_quadruleTheta );
    
    M_thetaGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( thetaPoly );
    M_thetaGenbasisRhoTheta -> setRho( M_Rho );
    M_thetaGenbasisRhoTheta -> setTheta( M_Theta );
    M_thetaGenbasisRhoTheta -> setNumberModes( M_mtheta );
    M_thetaGenbasisRhoTheta -> evaluateBasis( M_thetaphirho, M_thetadphirho, M_thetaphitheta, M_thetadphitheta,
                                              M_quadruleRho, M_quadruleTheta );

    M_pGenbasisRhoTheta = UneducatedBasisFactory::instance().createObject( pPoly );
    M_pGenbasisRhoTheta -> setRho( M_Rho );
    M_pGenbasisRhoTheta -> setTheta( M_Theta );
    M_pGenbasisRhoTheta -> setNumberModes( M_mp );
    M_pGenbasisRhoTheta -> evaluateBasis( M_pphirho, M_pdphirho, M_pphitheta, M_pdphitheta,
                                          M_quadruleRho, M_quadruleTheta );

/*std::cout << "x-basis" << std::endl;
   for( UInt i=0; i<M_mx; ++i ) 
   {
       for( UInt k=0; k<32; ++k )
           std::cout << M_xphirho[i][k] << " ";
       std::cout << std::endl;
   }
std::cout << "p-basis" << std::endl;
    for( UInt i=0; i<M_mp; ++i ) 
   {
       for( UInt k=0; k<32; ++k )
           std::cout << M_pphirho[i][k] << " ";
       std::cout << std::endl;
   }*/
}


void NSModalSpacePipe::
showMe() const
{
    std::cout << "---- MODAL SPACE CIRCULAR SHOWME ---" << std::endl;
    std::cout << "M_mx = "    <<    M_mx        << std::endl;
    std::cout << "M_mr = "    <<    M_mr        << std::endl;
    std::cout << "M_mtheta = "<<    M_mtheta    << std::endl;
    std::cout << "M_mp = "    <<    M_mp        << std::endl;
    std::cout << "------------------------------------" << std::endl;
}

void NSModalSpacePipe::
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

std::vector<Real> NSModalSpacePipe::
xFourierCoefficients( const function_Type& g, const Real& t ) const
{
    return xFourierCoefficients( g,t,0 );
}
std::vector<Real> NSModalSpacePipe::
rFourierCoefficients( const function_Type& g, const Real& t ) const
{
    return rFourierCoefficients( g,t,0 );
}
std::vector<Real> NSModalSpacePipe::
thetaFourierCoefficients( const function_Type& g, const Real& t ) const
{
    return thetaFourierCoefficients( g,t,0 );
}
std::vector<Real> NSModalSpacePipe::
pFourierCoefficients( const function_Type& g, const Real& t ) const
{
    return pFourierCoefficients( g,t,0 );
}
//Fourier Coefficients in the case g and R are x-dependent
std::vector<Real> NSModalSpacePipe::
xFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>                    FourCoeff( M_mx, 0.0 );
    
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >        evaluate_g;

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
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            evaluate_g[n][m] = g( t , x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mx; ++k )
    {
        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
        {
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
                Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

                Real inverseThetahat = th * M_Theta;
                Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );

                FourCoeff[k] += evaluate_g[n][m] *
                                M_xphirho[k][n] * M_xphitheta[k][m] *
                                //M_map->fJacobian()( t , x , inverseRhat , inverseThetahat , 0 ) *
                                rh * M_quadruleRho->weight( n ) * //r dr 
                                M_Theta * M_quadruleTheta->weight( m );    // dtheta
            }
        }
    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-dependent
std::vector<Real> NSModalSpacePipe::
rFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>                    FourCoeff( M_mr, 0.0 );
    
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >        evaluate_g;

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
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            evaluate_g[n][m] = g( t, x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mr; ++k )
    {
        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
                Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

                Real inverseThetahat = th * M_Theta;
                Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
                
                FourCoeff[k] += evaluate_g[n][m] *
                                M_rphirho[k][n] * M_rphitheta[k][m] *
                                //M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) *
                                rh * M_quadruleRho->weight( n ) * //r dr
                                M_Theta * M_quadruleTheta->weight( m );    // d theta
            }

    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-independent
std::vector<Real> NSModalSpacePipe::
thetaFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>                    FourCoeff( M_mtheta, 0.0 );
    
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >        evaluate_g;

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
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            evaluate_g[n][m] = g( t, x , inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mtheta; ++k )
    {
        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
                Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

                Real inverseThetahat = th * M_Theta;
                Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
                
                FourCoeff[k] += evaluate_g[n][m] *
                                M_thetaphirho[k][n] * M_thetaphitheta[k][m] *
                                //M_map->fJacobian()( t, x , inverseRhat , inverseThetahat , 0 ) *
                                rh * M_quadruleRho->weight( n ) * //r dr
                                M_Theta * M_quadruleTheta->weight( m );    // dtheta
            }
    }
    return FourCoeff;
}

//Fourier Coefficients in the case g is x-dependent
std::vector<Real> NSModalSpacePipe::
pFourierCoefficients( const function_Type& g, const Real& t, const Real& x ) const
{
    // Initialization of FourCoeff as 0, length M_mtot
    std::vector<Real>                    FourCoeff( M_mp, 0.0 );
    
    // Loop to evaluate g on the quadrature nodes
    // a matrix 32x32 (NnodesXNnodes)
    std::vector<std::vector<Real> >        evaluate_g;

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
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            evaluate_g[n][m] = g( t, x, inverseRhat , inverseThetahat , 0 );
        }
    // Loop on the modes, for each mode compute the associated fourier coefficients
    // \f[ \int_{[0,Rho]\times[0,Theta]}  g(y,z) \phi_k drho dtheta \f]
    for ( UInt k = 0; k != M_mp; ++k )
    {
        //loop over quadrature nodes
        for ( UInt n = 0; n != M_quadruleRho -> nbQuadPt(); ++n )
            for ( UInt m = 0; m != M_quadruleTheta -> nbQuadPt(); ++m )
            {
                Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
                Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

                Real inverseThetahat = th * M_Theta;
                Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
                
                FourCoeff[k] += evaluate_g[n][m] *
                                M_pphirho[k][n] * M_pphitheta[k][m] *
                                //M_map->fJacobian()( t, x, inverseRhat , inverseThetahat , 0 ) *
                                rh * M_quadruleRho->weight( n ) * //r dr 
                                M_Theta * M_quadruleTheta->weight( m );    // dtheta
            }
    }
    return FourCoeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpacePipe::
xFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_xphirho[k][n] * M_xphitheta[k][m] *
                     rh * M_quadruleRho->weight( n ) * 
                     M_Theta * M_quadruleTheta->weight( m );
        }

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpacePipe::
rFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    Real normrho = 1.0 / M_Rho;
    Real normtheta = 1.0 / sqrt(2. * M_PI );

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_rphirho[k][n] * M_rphitheta[k][m] *
                     rh * M_quadruleRho->weight( n ) * 
                     M_Theta * M_quadruleTheta->weight( m );
        }

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpacePipe::
thetaFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            coeff += f ( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_thetaphirho[k][n] * M_thetaphitheta[k][m] *
                     rh * M_quadruleRho->weight( n ) * 
                     M_Theta * M_quadruleTheta->weight( m );
        }

    return coeff;
}

//Single Fourier coefficient in the case f is x-dependent
Real NSModalSpacePipe::
pFourierCoeffPointWise( const Real& t, const Real& x, const function_Type& f, const UInt& k ) const
{
    Real coeff = 0.0;

    for ( UInt n = 0; n != M_quadruleRho->nbQuadPt(); ++n )
        for ( UInt m = 0; m != M_quadruleTheta->nbQuadPt(); ++m )
        {    
            Real rh( M_quadruleRho -> quadPointCoor( n, 0 ) );   
            Real th( M_quadruleTheta -> quadPointCoor( m, 0 ) );   

            Real inverseThetahat = th * M_Theta;
            Real inverseRhat = M_map->inverseRhat()( t, x, rh, inverseThetahat, m );
            coeff += f( t, x , inverseRhat , inverseThetahat , 0 ) *
                     M_pphirho[k][n] * M_pphitheta[k][m] *
                     rh * M_quadruleRho->weight( n ) * 
                     M_Theta * M_quadruleTheta->weight( m );
        }

    return coeff;
}

// ------------------        Compute Methods        ---------------------------
//        The Jacobian links the physical domain with the reference one.
//        The further 2 * M_PI is a scale factor due to the numerical domain of integration:
//        int_0^(2*pi) d_theta = int_0^1 2*pi*d_thetah
// --------------------------------- X-THETA-DEPENDENCE ------------------------------------------------------
void
NSModalSpacePipe::
compute_r00xx( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
               vector_Type& R00xx ) const
{
    UInt dof( R00xx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    
    Real Idr110(0);
    Real Idr110rr(0);
    Real Idr101(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real Idr000(0);
    Real dr(0);
    Real r(0);

    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r =  qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110 += xdphirho[k][n] *
                  xdphirho[j][n] * 
                  r * dr;

        Idr110rr += xdphirho[k][n] *
                    xdphirho[j][n] * 
                    r * r * // This is for Dr^2
                    r * dr;
                  
        Idr101 += 1 / r *
                  xdphirho[k][n] *
                  xphirho[j][n] * 
                  r * dr;

        Idr101r += 1 / r *
                   xdphirho[k][n] *
                   xphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
                  
        Idr011 += 1 / r *
                  xphirho[k][n] *
                  xdphirho[j][n] * 
                  r * dr;

        Idr011r += 1 / r *
                  xphirho[k][n] *
                  xdphirho[j][n] * 
                  r * // This is for Dr
                  r * dr;
                  
        Idr002 += 1 / ( r * r ) *
                  xphirho[k][n] *
                  xphirho[j][n] * 
                  r * dr;
                  
        Idr000 += xphirho[k][n] *
                  xphirho[j][n] * 
                  r * dr;
    }
    
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R00xx[i] += ( 
                              Idr110rr *
                              2*nu * Dr[i][h] * Dr[i][h] *
                              xphitheta[k][h] *
                              xphitheta[j][h] +

                              Idr110 *
                              (
                                  nu * Jr[i][h] * Jr[i][h] +
                                  nu * Dthetar[i][h] * Dthetar[i][h]
                              ) *
                              xphitheta[k][h] *
                              xphitheta[j][h] +
                                   
                              Idr101r *
                              2*nu * Dr[i][h] * Dtheta[i][h] *
                              xphitheta[k][h] *
                              xdphitheta[j][h] +

                              Idr101 *
                              (
                                  nu * Jr[i][h] * Drtheta[i][h] +
                                  nu * Dthetar[i][h] * Jtheta[i][h]
                              ) *
                              xphitheta[k][h] *
                              xdphitheta[j][h] +
                               
                              Idr011r *    
                              2*nu * Dr[i][h] * Dtheta[i][h] *
                              xdphitheta[k][h] *
                              xphitheta[j][h] +

                              Idr011 *    
                              (
                                  nu * Jr[i][h] * Drtheta[i][h] +
                                  nu * Dthetar[i][h] * Jtheta[i][h]
                              ) *
                              xdphitheta[k][h] *
                              xphitheta[j][h] +
                                   
                              Idr002 *    
                              ( 2*nu * Dtheta[i][h] * Dtheta[i][h] +
                                  nu * Drtheta[i][h] * Drtheta[i][h] +
                                  nu * Jtheta[i][h] * Jtheta[i][h]
                              ) *
                              xdphitheta[k][h] *
                              xdphitheta[j][h] +
                                
                              Idr000 *
                              alpha * ( 
                                        xphitheta[k][h] *
                                        xphitheta[j][h]
                                      ) 
                            ) * Jacobian[i][h] *
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r01xx( const UInt& k, const UInt& j, const Real& nu,
                                        vector_Type& R01xx ) const
{
    UInt dof( R01xx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    
    Real Idr100r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);

    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100r += xdphirho[k][n] *
                   xphirho[k][n] *
                   r * // This is for Dr
                   r * dr;
                  
        Idr001 += 1 / qrRho.quadPointCoor( n, 0 ) *
                  xphirho[k][n] *
                  xphirho[j][n] *
                  r * dr;
    }    
            
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R01xx[i] += (  2*nu*( 
                                      Idr100r *
                                      Dr[i][h] *
                                      xphitheta[k][h] *
                                      xphitheta[j][h] +
                                      
                                      Idr001 *
                                      Dtheta[i][h] *
                                      xdphitheta[k][h] *
                                      xphitheta[j][h]
                                    )
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r10xx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10xx ) const
{
    UInt dof( R10xx.blockSize(0) );
    Real Theta(M_Theta);

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    
    Real Idr010r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);

    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );

        Idr010r += xphirho[k][n] * 
                   xdphirho[j][n] *
                   r * // This is for Dr
                   r * dr ;
                  
        Idr001 += 1 / qrRho.quadPointCoor( n, 0 ) *
                  xphirho[k][n] * 
                  xphirho[j][n] *
                  r * dr;
    }

    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R10xx[i] += (  2*nu*( 
                                      Idr010r *
                                      Dr[i][h] *
                                      xphitheta[k][h] *
                                      xphitheta[j][h] +
                                      
                                      Idr001 *
                                      Dtheta[i][h] *
                                      xphitheta[k][h] *
                                      xdphitheta[j][h]
                                    )
                            ) * Jacobian[i][h] * 
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r11xx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R11xx ) const
{
    UInt dof( R11xx.blockSize(0) );
    Real Theta( M_Theta);

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    
    Real Idr000(0);
    Real r(0);
    Real dr(0);
 
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr000 += xphirho[k][n] *
                  xphirho[j][n] *
                  r * dr;

    }
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R11xx[i] += (  2*nu*( Idr000 *
                                      xphitheta[k][h] *
                                      xphitheta[j][h]
                                    )
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r00xr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00xr ) const
{
    UInt dof( R00xr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr110r(0);
    Real Idr101(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);

    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )    
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110r += xdphirho[k][n] * 
                   rdphirho[j][n] * 
                   r * //This is for Dr
                   r * dr;
                  
        Idr101 += 1. / r *
                  xdphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
                  
        Idr011r += 1. / r *
                   xphirho[k][n] * 
                   rdphirho[j][n] * 
                   r * //This is for Dr
                   r * dr;
                  
        Idr002 += 1. / ( r * r ) *
                  xphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
    
    
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R00xr[i] += (  nu*( Idr110r *
                                    Jr[i][h] * Dr[i][h] *
                                    xphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    Idr101 * (
                                               Dtheta[i][h] * Jr[i][h] *
                                               xphitheta[k][h] *
                                               rdphitheta[j][h] +
                                               
                                               Dthetar[i][h] * Dtheta[i][h] *
                                               xphitheta[k][h] *
                                               rphitheta[j][h]
                                             ) +
                                    
                                    Idr011r *
                                    Drtheta[i][h] * Dr[i][h] *
                                    xdphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    Idr002 * (
                                               Drtheta[i][h] * Dtheta[i][h] *
                                               xdphitheta[k][h] *
                                               rdphitheta[j][h] +
                                               
                                               Dtheta[i][h] * Jtheta[i][h] *
                                               xdphitheta[k][h] *
                                               rphitheta[j][h]
                                             )
                                  )
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r01xr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01xr ) const
{
    UInt dof( R01xr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr100(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);

    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )    
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100 += xdphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
        
        Idr001 += 1. / r *
                  xphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
                  
    }
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R01xr[i] += (  nu*( Idr100 *
                                    Jr[i][h] *
                                    xphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    Idr001 *
                                    Drtheta[i][h] *
                                    xdphitheta[k][h] *
                                    rphitheta[j][h]
                                  )
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r00xt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00xt ) const
{
    UInt dof( R00xt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr110r(0);
    Real Idr101(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);
    
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110r += xdphirho[k][n] * 
                   thetadphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr101 += 1. / r *
                  xdphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
        
        Idr011r += 1. / r *
                   xphirho[k][n] * 
                   thetadphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  xphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R00xt[i] += (  nu*( Idr110r *
                                    Dthetar[i][h] * Dr[i][h] *
                                    xphitheta[k][h] *
                                    thetaphitheta[j][h] +
                                    
                                    Idr101 * (
                                               Dthetar[i][h] * Dtheta[i][h] *
                                               xphitheta[k][h] *
                                               thetadphitheta[j][h]
                                               
                                               - Dtheta[i][h] * Jr[i][h] *
                                               xphitheta[k][h] *
                                               thetaphitheta[j][h]
                                              ) +
                                    
                                    Idr011r *
                                    Dr[i][h] * Jtheta[i][h] *
                                    xdphitheta[k][h] *
                                    thetaphitheta[j][h] +
                                    
                                    Idr002 * (
                                               Dtheta[i][h] * Jtheta[i][h] *
                                               xdphitheta[k][h] *
                                               thetadphitheta[j][h]
                                               
                                               - Drtheta[i][h] * Dtheta[i][h] *
                                               xdphitheta[k][h] *
                                               thetaphitheta[j][h]
                                             )
                                  )
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r01xt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01xt ) const
{
    UInt dof( R01xt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr100(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r =  qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100 += xdphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
                  
        Idr001 += 1. / r *
                  xphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R01xt[i] += (  nu*( Idr100 *
                                    Dthetar[i][h] *
                                    xphitheta[k][h] *
                                    thetaphitheta[j][h] +
                                    
                                    Idr001 *
                                    Jtheta[i][h] *
                                    xdphitheta[k][h] *
                                    thetaphitheta[j][h]
                                  )
                            ) * Jacobian[i][h] *
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop
    
    return;
}

void
NSModalSpacePipe::
compute_r00xp( const UInt& k, const UInt& j,
               vector_Type& R00xp ) const
{

    UInt dof( R00xp.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr100r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100r += xdphirho[k][n] * 
                   pphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  xphirho[k][n] * 
                  pphirho[j][n] * 
                  r * dr;
    }
        
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R00xp[i] += -( Idr100r *
                              Dr[i][h] *
                              xphitheta[k][h] *
                              pphitheta[j][h] +
                               
                              Idr001 *
                              Dtheta[i][h] *
                              xdphitheta[k][h] *
                              pphitheta[j][h]
                            ) * Jacobian[i][h] *  
                                Theta * qrTheta.weight( h );
            } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10xp( const UInt& k, const UInt& j,
               vector_Type& R10xp ) const
{
    UInt dof( R10xp.blockSize(0) );
    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < M_quadruleRho->nbQuadPt(); ++n )
    {
        r = M_quadruleRho->quadPointCoor( n, 0 );
        dr = M_quadruleRho->weight( n );
        
        Idr000 += M_xphirho[k][n] * 
                  M_pphirho[j][n] * 
                  r * dr;
    }
    
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
            for ( h = 0; h < M_quadruleTheta->nbQuadPt(); ++h )
            {
                R10xp[i] += -( Idr000 * 
                               M_xphitheta[k][h] *
                               M_pphitheta[j][h]
                             ) * M_map->Jacobian()[i][h] *  
                                M_Theta * M_quadruleTheta->weight( h );
            } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00rr( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
               vector_Type& R00rr ) const
{

    UInt dof( R00rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr110rr(0);
    Real Idr110(0);
    Real Idr101(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110rr += rdphirho[k][n] * 
                    rdphirho[j][n] * 
                    r * r * // This is for Dr^2
                    r * dr;

        Idr110 += rdphirho[k][n] * 
                  rdphirho[j][n] * 
                  r * dr;
        
        Idr101 += 1. / r *
                  rdphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;

        Idr101r += 1. / r *
                   rdphirho[k][n] * 
                   rphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr011 += 1. / r *
                  rphirho[k][n] * 
                  rdphirho[j][n] * 
                  r * dr;

        Idr011r += 1. / r *
                   rphirho[k][n] * 
                   rdphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  rphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
        
        Idr000 += rphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00rr[i] += ( 
                                  Idr110rr *
                                  nu * Dr[i][h] * Dr[i][h] *
                                  rphitheta[k][h] *
                                  rphitheta[j][h] +

                                  Idr110 *
                                  (   
                                      nu * Dthetar[i][h] * Dthetar[i][h] + 
                                    2*nu * Jr[i][h] * Jr[i][h]
                                  ) *
                                  rphitheta[k][h] *
                                  rphitheta[j][h] +
                                         
                                  Idr101r *
                                  nu * Dr[i][h] * Dtheta[i][h] *
                                  rphitheta[k][h] *
                                  rdphitheta[j][h] +

                                  Idr101 *
                                  (
                                    (   
                                        nu * Dthetar[i][h] * Jtheta[i][h] +
                                      2*nu * Jr[i][h] * Drtheta[i][h]
                                    ) *
                                    rphitheta[k][h] *
                                    rdphitheta[j][h] +
                                    
                                    nu * Drtheta[i][h] * Dthetar[i][h] *
                                    rphitheta[k][h] *
                                    rphitheta[j][h]
                                  ) +
                                         
                                  Idr011r *
                                  nu * Dr[i][h] * Dtheta[i][h] *
                                  rdphitheta[k][h] *
                                  rphitheta[j][h] +

                                  Idr011 *
                                  (
                                    (   
                                        nu * Dthetar[i][h] * Jtheta[i][h] +
                                      2*nu * Jr[i][h] * Drtheta[i][h]
                                    ) *
                                    rdphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    nu * Dthetar[i][h] * Drtheta[i][h] *
                                    rphitheta[k][h] *
                                    rphitheta[j][h]
                                  ) +
                                  
                                  Idr002 *
                                  (
                                    (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                        nu * Jtheta[i][h] * Jtheta[i][h] +
                                      2*nu * Drtheta[i][h] * Drtheta[i][h]
                                    ) *
                                    rdphitheta[k][h] *
                                    rdphitheta[j][h] +
                                    
                                    (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                      2*nu * Jtheta[i][h] * Jtheta[i][h] +
                                        nu * Drtheta[i][h] * Drtheta[i][h]
                                    ) *
                                    rphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    nu * Jtheta[i][h] * Drtheta[i][h] *
                                    rdphitheta[k][h] *
                                    rphitheta[j][h] +
                                    
                                    nu * Jtheta[i][h] * Drtheta[i][h] *
                                    rphitheta[k][h] *
                                    rdphitheta[j][h]
                                  ) +
                                         
                                  Idr000 *
                                  alpha * (
                                            rphitheta[k][h] *
                                            rphitheta[j][h]
                                          )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r01rr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01rr ) const
{

    UInt dof( R01rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr100r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100r += rdphirho[k][n] * 
                   rphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
    
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R01rr[i] += nu * ( Idr100r *
                                       Dr[i][h] * 
                                       rphitheta[k][h] *
                                       rphitheta[j][h] +
                                       
                                       Idr001 *
                                       Dtheta[i][h] *
                                       rdphitheta[k][h] *
                                       rphitheta[j][h]
                                     ) * Jacobian[i][h] *  
                                         Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10rr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10rr ) const
{

    UInt dof( R10rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );

    Real Idr010r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010r += rphirho[k][n] * 
                   rdphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
        
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
            for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
            {
                R10rr[i] += nu * ( Idr010r *
                                   Dr[i][h] * 
                                   rphitheta[k][h] *
                                   rphitheta[j][h] +
                                   
                                   Idr001 *
                                   Dtheta[i][h] *
                                   rphitheta[k][h] *
                                   rdphitheta[j][h]
                                 ) * Jacobian[i][h] *
                                     Theta * qrTheta.weight( h );
            } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r11rr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R11rr ) const
{
    UInt dof( R11rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );

    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr000 += rphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
    
    // normrho and normtheta are included in the basis
    for ( i = 0; i < dof; ++i )
    {
        for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
        {
            R11rr[i] += nu * ( Idr000 * 
                               rphitheta[k][h] *
                               rphitheta[j][h]
                              ) * Jacobian[i][h] *  
                                  Theta * qrTheta.weight( h );
        } // n for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10rx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10rx ) const
{

    UInt dof( R10rx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr010(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010 += rphirho[k][n] * 
                  xdphirho[j][n] * 
                  r * dr;
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R10rx[i] += nu * ( Idr010 *
                                       Jr[i][h] * 
                                       rphitheta[k][h] *
                                       xphitheta[j][h] +
                                       
                                       Idr001 *
                                       Drtheta[i][h] *
                                       rphitheta[k][h] *
                                       xdphitheta[j][h]
                                     ) * Jacobian[i][h] *    
                                         Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00rx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00rx ) const
{
    UInt dof( R00rx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real Idr110r(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);
        
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110r += rdphirho[k][n] * 
                   xdphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr101r += 1. / r *
                   rdphirho[k][n] * 
                   xphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr011 += 1. / r *
                  rphirho[k][n] * 
                  xdphirho[j][n] * 
                  r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  rphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00rx[i] += nu * ( Idr110r *
                                       Jr[i][h] * Dr[i][h] *
                                       rphitheta[k][h] *
                                       xphitheta[j][h] +
                                       
                                       Idr101r *
                                       Dr[i][h] * Drtheta[i][h] *
                                       rphitheta[k][h] *
                                       xdphitheta[j][h] +
                                       
                                       Idr011 *
                                       ( 
                                         Jr[i][h] * Dtheta[i][h] *
                                         rdphitheta[k][h] *
                                         xphitheta[j][h] +
                                         
                                         Dthetar[i][h] * Dtheta[i][h] *
                                         rphitheta[k][h] *
                                         xphitheta[j][h]
                                       ) +
                                       
                                       Idr002 *
                                       (
                                         Dtheta[i][h] * Drtheta[i][h] *
                                         rdphitheta[k][h] *
                                         xdphitheta[j][h] +
                                         
                                         Dtheta[i][h] * Jtheta[i][h] *
                                         rphitheta[k][h] *
                                         xdphitheta[j][h]
                                       )
                                     ) * Jacobian[i][h] *  
                                         Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00rt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00rt ) const
{

    UInt dof( R00rt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr110(0);
    Real Idr101(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110 += rdphirho[k][n] * 
                  thetadphirho[j][n] * 
                  r * dr;
        
        Idr101 += 1. / r *
                  rdphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;

        Idr101r += 1. / r *
                   rdphirho[k][n] * 
                   thetaphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr011 += 1. / r *
                  rphirho[k][n] * 
                  thetadphirho[j][n] * 
                  r * dr;

        Idr011r += 1. / r *
                   rphirho[k][n] * 
                   thetadphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  rphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00rt[i] += ( nu * Idr110 *
                                       Jr[i][h] * Dthetar[i][h] *
                                       rphitheta[k][h] *
                                       thetaphitheta[j][h] +
                                       
                                       - Idr101r * 
                                       nu * Dtheta[i][h] * Dr[i][h] *
                                       rphitheta[k][h] *
                                       thetaphitheta[j][h] +

                                       Idr101 *
                                       (
                                         nu * Dthetar[i][h] * Drtheta[i][h] *
                                              rphitheta[k][h] *
                                              thetadphitheta[j][h] 
                                         
                                         - (
                                               nu * Dthetar[i][h] * Jtheta[i][h] +
                                             2*nu * Drtheta[i][h] * Jr[i][h]
                                           ) * rphitheta[k][h] *
                                               thetaphitheta[j][h] 
                                       ) +
                                       
                                       Idr011r *
                                       nu * Dtheta[i][h] * Dr[i][h] *
                                       rphitheta[k][h] *
                                       thetaphitheta[j][h] +

                                       Idr011 *
                                       (
                                         nu * Jr[i][h] * Jtheta[i][h] *
                                              rdphitheta[k][h] *
                                              thetaphitheta[j][h] +
                                              
                                         (
                                             nu * Drtheta[i][h] * Jr[i][h] +
                                           2*nu * Dthetar[i][h] * Jtheta[i][h]
                                         ) * rphitheta[k][h] *
                                             thetaphitheta[j][h] 
                                       ) +
                                       
                                       Idr002 *
                                       (
                                         nu * Drtheta[i][h] * Jtheta[i][h] *
                                              rdphitheta[k][h] *
                                              thetadphitheta[j][h]
                                         
                                         - (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                               nu * Jtheta[i][h] * Jtheta[i][h] +
                                             2*nu * Drtheta[i][h] * Drtheta[i][h]
                                           ) * rdphitheta[k][h] *
                                               thetaphitheta[j][h] +
                                               
                                         (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                           2*nu * Jtheta[i][h] * Jtheta[i][h] +
                                             nu * Drtheta[i][h] * Drtheta[i][h]  
                                           ) * rphitheta[k][h] *
                                               thetadphitheta[j][h]
                                               
                                         - nu * Drtheta[i][h] * Jtheta[i][h] *
                                           rphitheta[k][h] *
                                           thetaphitheta[j][h]
                                       )
                                     ) * Jacobian[i][h] *  
                                         Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10rt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10rt ) const
{
    UInt dof( R10rt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R10rt[i] += -nu * ( Idr001 *
                                        Dtheta[i][h] *
                                        rphitheta[k][h] *
                                        thetaphitheta[j][h]
                                      ) * Jacobian[i][h] *    
                                          Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r01rt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01rt ) const
{

    UInt dof( R01rt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R01rt[i] +=  nu * ( Idr001 *
                                        Dtheta[i][h] *
                                        rphitheta[k][h] *
                                        thetaphitheta[j][h]
                                      ) * Jacobian[i][h] *    
                                          Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00rp( const UInt& k, const UInt& j,
               vector_Type& R00rp ) const
{

    UInt dof( R00rp.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr100(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100 += rdphirho[k][n] * 
                  pphirho[j][n] * 
                  r * dr;
        
        Idr001 += 1. / r *
                  rphirho[k][n] * 
                  pphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00rp[i] += -( Idr100 *
                                  Jr[i][h] *
                                  rphitheta[k][h] *
                                  pphitheta[j][h] +
                                   
                                  Idr001 *
                                  (
                                    Drtheta[i][h] *
                                    rdphitheta[k][h] *
                                    pphitheta[j][h] +
                                      
                                    Jtheta[i][h] *
                                    rphitheta[k][h] *
                                    pphitheta[j][h]
                                  )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00tt( const UInt& k, const UInt& j, const Real& nu, const Real& alpha,
               vector_Type& R00tt ) const
{

    UInt dof( R00tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr110rr(0);
    Real Idr110(0);
    Real Idr101(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110rr += thetadphirho[k][n] * 
                    thetadphirho[j][n] * 
                    r * r * // This is for Dr^2
                    r * dr;

        Idr110 += thetadphirho[k][n] * 
                  thetadphirho[j][n] * 
                  r * dr;
        
        Idr101 += 1. / r *
                  thetadphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;

        Idr101r += 1. / r *
                  thetadphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * // This is for Dr
                  r * dr;
        
        Idr011 += 1. / r *
                  thetaphirho[k][n] * 
                  thetadphirho[j][n] * 
                  r * dr;

        Idr011r += 1. / r *
                   thetaphirho[k][n] * 
                   thetadphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  thetaphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
                  
        Idr000 += thetaphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00tt[i] += ( 
                                  Idr110rr *
                                  nu * Dr[i][h] * Dr[i][h] *
                                  thetaphitheta[k][h] *
                                  thetaphitheta[j][h] +

                                  Idr110 *
                                  (
                                      nu * Jr[i][h] * Jr[i][h] +
                                    2*nu * Dthetar[i][h] * Dthetar[i][h]
                                  ) * 
                                  thetaphitheta[k][h] *
                                  thetaphitheta[j][h] +
                                       
                                  Idr101r *
                                  nu * Dtheta[i][h] * Dr[i][h] *
                                  thetaphitheta[k][h] *
                                  thetadphitheta[j][h] +

                                  Idr101 *
                                  (
                                    (
                                        nu * Drtheta[i][h] * Jr[i][h] +
                                      2*nu * Dthetar[i][h] * Jtheta[i][h]
                                    ) * 
                                    thetaphitheta[k][h] *
                                    thetadphitheta[j][h]
                                    
                                    - nu * Jr[i][h] * Jtheta[i][h] *
                                           thetaphitheta[k][h] *
                                           thetaphitheta[j][h] 
                                  ) +
                                  
                                  Idr011r *
                                  nu * Dtheta[i][h] * Dr[i][h] *
                                  thetadphitheta[k][h] *
                                  thetaphitheta[j][h] +

                                  Idr011 *
                                  (
                                    (
                                        nu * Drtheta[i][h] * Jr[i][h] +
                                      2*nu * Dthetar[i][h] * Jtheta[i][h]
                                    ) * 
                                    thetadphitheta[k][h] *
                                    thetaphitheta[j][h] 
                                  
                                    - nu * Jr[i][h] * Jtheta[i][h] *
                                           thetaphitheta[k][h] *
                                           thetaphitheta[j][h] 
                                  ) +
                                  
                                  Idr002 *
                                  (
                                    (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                        nu * Drtheta[i][h] * Drtheta[i][h] +
                                      2*nu * Jtheta[i][h] * Jtheta[i][h]
                                    ) * 
                                    M_thetadphitheta[k][h] *
                                    M_thetadphitheta[j][h]
                                    
                                    - nu * Drtheta[i][h] * Jtheta[i][h] *
                                           thetadphitheta[k][h] *
                                           thetaphitheta[j][h] 
                                           
                                    - nu * Drtheta[i][h] * Jtheta[i][h] *
                                           thetaphitheta[k][h] *
                                           thetadphitheta[j][h] +
                                           
                                    (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                      2*nu * Drtheta[i][h] * Drtheta[i][h] +
                                        nu * Jtheta[i][h] * Jtheta[i][h]
                                    ) * 
                                    thetaphitheta[k][h] *
                                    thetaphitheta[j][h]
                                  ) +
                                     
                                  Idr000 *    
                                  alpha * (
                                            thetaphitheta[k][h] *
                                            thetaphitheta[j][h]
                                          )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10tt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10tt ) const
{
    UInt dof( R10tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr010r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010r += thetaphirho[k][n] * 
                   thetadphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R10tt[i] += ( nu * ( Idr010r *
                                         Dr[i][h] *
                                         thetaphitheta[k][h] *
                                         thetaphitheta[j][h] +
                                       
                                         Idr001 *
                                         Dtheta[i][h] *
                                         thetaphitheta[k][h] *
                                         thetadphitheta[j][h]
                                       )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r01tt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01tt ) const
{
    UInt dof( R01tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr100r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100r += thetadphirho[k][n] * 
                   thetaphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
        
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R01tt[i] += ( nu * ( Idr100r *
                                         Dr[i][h] *
                                         thetaphitheta[k][h] *
                                         thetaphitheta[j][h] +
                                       
                                         Idr001 *
                                         Dtheta[i][h] *
                                         thetadphitheta[k][h] *
                                         thetaphitheta[j][h]
                                       )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // n for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r11tt( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R11tt ) const
{

    UInt dof( R11tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr000 += thetaphirho[k][n] * 
                  thetaphirho[j][n] *
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R11tt[i] += ( nu * ( Idr000 * 
                                         thetaphitheta[k][h] *
                                         thetaphitheta[j][h]
                                       )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10tx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10tx ) const
{

    UInt dof( R10tx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr010(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010 += thetaphirho[k][n] * 
                  xdphirho[j][n] * 
                  r * dr;
        
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
        
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R10tx[i] += ( nu * ( Idr010 * 
                                         Dthetar[i][h] *
                                         thetaphitheta[k][h] *
                                         xphitheta[j][h] +
                                       
                                         Idr001 *
                                         Jtheta[i][h] *
                                         thetaphitheta[k][h] *
                                         xdphitheta[j][h]
                                       )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00tx( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00tx ) const
{

    UInt dof( R00tx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr110r(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110r += thetadphirho[k][n] * 
                   xdphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr101r += 1. / r *
                   thetadphirho[k][n] * 
                   xphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr011 += 1. / r *
                  thetaphirho[k][n] * 
                  xdphirho[j][n] * 
                  r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  thetaphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
    }
    
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00tx[i] += ( nu * ( Idr110r *
                                         Dr[i][h] * Dthetar[i][h] *
                                         thetaphitheta[k][h] *
                                         xphitheta[j][h] +
                                       
                                         Idr101r *
                                         Dr[i][h] * Jtheta[i][h] *
                                         thetaphitheta[k][h] *
                                         xdphitheta[j][h] +
                                         
                                         Idr011 *
                                         (
                                           Dthetar[i][h] * Dtheta[i][h] *
                                           thetadphitheta[k][h] *
                                           xphitheta[j][h]
                                           
                                          -Dtheta[i][h] * Jr[i][h] *
                                           thetaphitheta[k][h] *
                                           xphitheta[j][h]
                                         ) +
                                         
                                         Idr002 *
                                         (
                                           Dtheta[i][h] * Jtheta[i][h] *
                                           thetadphitheta[k][h] *
                                           xdphitheta[j][h]
                                           
                                          -Dtheta[i][h] * Drtheta[i][h] *
                                           thetaphitheta[k][h] *
                                           xdphitheta[j][h]
                                         )
                                       ) 
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00tr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R00tr ) const
{
    UInt dof( R00tr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real Idr110(0);
    Real Idr101(0);
    Real Idr101r(0);
    Real Idr011(0);
    Real Idr011r(0);
    Real Idr002(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr110 += thetadphirho[k][n] * 
                  rdphirho[j][n] * 
                  r * dr;
        
        Idr101 += 1. / r *
                  thetadphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;

        Idr101r += 1. / r *
                   thetadphirho[k][n] * 
                   rphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr011 += 1. / r *
                  thetaphirho[k][n] * 
                  rdphirho[j][n] * 
                  r * dr;

        Idr011r += 1. / r *
                   thetaphirho[k][n] * 
                   rdphirho[j][n] * 
                   r * // This is for Dr
                   r * dr;
        
        Idr002 += 1. / ( r * r ) *
                  thetaphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
        
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00tr[i] += ( Idr110 *
                                  nu * Dthetar[i][h] * Jr[i][h] *
                                       thetaphitheta[k][h] *
                                       rphitheta[j][h] +
                                       
                                  Idr101r *
                                  nu * Dtheta[i][h] * Dr[i][h] *
                                  thetaphitheta[k][h] *
                                  rphitheta[j][h] +

                                  Idr101 *
                                  (
                                    nu * Jtheta[i][h] * Jr[i][h] *
                                         thetaphitheta[k][h] *
                                         rdphitheta[j][h] +
                                  
                                    (
                                      2*nu * Jtheta[i][h] * Dthetar[i][h] +
                                        nu * Jr[i][h] * Drtheta[i][h]
                                    ) *
                                    thetaphitheta[k][h] *
                                    rphitheta[j][h]
                                  ) 
                                       
                                  - Idr011r *
                                  nu * Dtheta[i][h] * Dr[i][h] *
                                  thetaphitheta[k][h] *
                                  rphitheta[j][h] +

                                  Idr011 *
                                  (
                                    nu * Dthetar[i][h] * Drtheta[i][h] *
                                         thetadphitheta[k][h] *
                                         rphitheta[j][h]
                                         
                                   -(
                                        nu * Jtheta[i][h] * Dthetar[i][h] +
                                      2*nu * Jr[i][h] * Drtheta[i][h]
                                    ) *
                                    thetaphitheta[k][h] *
                                    rphitheta[j][h]
                                  ) +
                                       
                                  Idr002 *
                                  (
                                    nu * Jtheta[i][h] * Drtheta[i][h] *
                                         thetadphitheta[k][h] *
                                         rdphitheta[j][h]
                                         
                                    -(   nu * Dtheta[i][h] * Dtheta[i][h] +
                                         nu * Jtheta[i][h] * Jtheta[i][h] +
                                       2*nu * Drtheta[i][h] * Drtheta[i][h]
                                     ) *
                                     M_thetaphitheta[k][h] *
                                     M_rdphitheta[j][h] +
                                     
                                     (   nu * Dtheta[i][h] * Dtheta[i][h] +
                                       2*nu * Jtheta[i][h] * Jtheta[i][h] +
                                         nu * Drtheta[i][h] * Drtheta[i][h]
                                     ) *
                                     thetadphitheta[k][h] *
                                     rphitheta[j][h]
                                     
                                   -nu * Jtheta[i][h] * Drtheta[i][h] *
                                         thetaphitheta[k][h] *
                                         rphitheta[j][h]
                                  )
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r01tr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R01tr ) const
{

    UInt dof( R01tr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R01tr[i] +=  -nu * ( Idr001 *
                                         Dtheta[i][h] *
                                         thetaphitheta[k][h] *
                                         rphitheta[j][h]
                                      ) * Jacobian[i][h] *    
                                          Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r10tr( const UInt& k, const UInt& j, const Real& nu,
               vector_Type& R10tr ) const
{
    UInt dof( R10tr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R10tr[i] += nu * ( Idr001 *
                                        Dtheta[i][h] *
                                        thetaphitheta[k][h] *
                                        rphitheta[j][h]
                                      ) * Jacobian[i][h] *    
                                          Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00tp( const UInt& k, const UInt& j,
               vector_Type& R00tp ) const
{
    UInt dof( R00tp.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr100(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr100 += thetadphirho[k][n] * 
                  pphirho[j][n] * 
                  r * dr;
                  
        Idr001 += 1. / r *
                  thetaphirho[k][n] * 
                  pphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00tp[i] += ( -Idr100 *
                                   Dthetar[i][h] *
                                   thetaphitheta[k][h] *
                                   pphitheta[j][h]+
                                  
                                   Idr001 *
                                  (
                                   -Jtheta[i][h] *
                                    thetadphitheta[k][h] *
                                    pphitheta[j][h] +
                                    
                                    Drtheta[i][h] *
                                    thetaphitheta[k][h] *
                                    pphitheta[j][h]
                                  )
                                ) * Jacobian[i][h] *
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00px( const UInt& k, const UInt& j,
               vector_Type& R00px ) const
{
    UInt dof( R00px.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr010r(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010r += pphirho[k][n] * 
                   xdphirho[j][n] *
                   r * // This is for Dr
                   r * dr;
        
        Idr001 += 1. / r *
                  pphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00px[i] += -( Idr010r *
                                   Dr[i][h] *
                                   pphitheta[k][h] *
                                   xphitheta[j][h] +
                                  
                                   Idr001 *
                                   Dtheta[i][h] *
                                   pphitheta[k][h] *
                                   xdphitheta[j][h]
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r01px( const UInt& k, const UInt& j,
               vector_Type& R01px ) const
{
    UInt dof( R01px.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr000(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr000 += pphirho[k][n] * 
                  xphirho[j][n] * 
                  r * dr;
    }
        
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R01px[i] += -( Idr000 *
                                   pphitheta[k][h] *
                                   xphitheta[j][h]
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00pr( const UInt& k, const UInt& j,
               vector_Type& R00pr ) const
{
    UInt dof( R00pr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr010(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010 += pphirho[k][n] * 
                  rdphirho[j][n] * 
                  r * dr;
                  
        Idr001 += 1. / r *
                  pphirho[k][n] * 
                  rphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00pr[i] += -( Idr010 *
                                   Jr[i][h] *
                                   pphitheta[k][h] *
                                   rphitheta[j][h] +
                                  
                                   Idr001 *
                                   Drtheta[i][h] *
                                   pphitheta[k][h] *
                                   rdphitheta[j][h] +
                                   
                                   Idr001 *
                                   Jtheta[i][h] *
                                   pphitheta[k][h] *
                                   rphitheta[j][h]
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_r00pt( const UInt& k, const UInt& j,
               vector_Type& R00pt ) const
{
    UInt dof( R00pt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    MBMatrix_type pphirho( M_pphirho );
    MBMatrix_type pphitheta( M_pphitheta );
    
    Real Idr010(0);
    Real Idr001(0);
    Real r(0);
    Real dr(0);
    
    UInt n, i, h;

    for ( n = 0; n < qrRho.nbQuadPt(); ++n )
    {
        r = qrRho.quadPointCoor( n, 0 );
        dr = qrRho.weight( n );
        
        Idr010 += pphirho[k][n] * 
                  thetadphirho[j][n] * 
                  r * dr;
         
        Idr001 += 1. / r *
                  pphirho[k][n] * 
                  thetaphirho[j][n] * 
                  r * dr;
    }
    
        // normrho and normtheta are included in the basis
        for ( i = 0; i < dof; ++i )
        {
                for ( h = 0; h < qrTheta.nbQuadPt(); ++h )
                {
                    R00pt[i] += (-Idr010 *
                                  Dthetar[i][h] *
                                  pphitheta[k][h] *
                                  thetaphitheta[j][h]
                                  
                                 -Idr001 *
                                  Jtheta[i][h] *
                                  pphitheta[k][h] *
                                  thetadphitheta[j][h] +
                                  
                                  Idr001 *
                                  Drtheta[i][h] *
                                  pphitheta[k][h] *
                                  thetaphitheta[j][h]
                                ) * Jacobian[i][h] *  
                                    Theta * qrTheta.weight( h );
                } // h for-loop
        } // i for-loop

    return ;
}

// --------------- Non linear terms for NS -------------------
void
NSModalSpacePipe::
compute_b10rr( const UInt& k, const UInt& j,
               vector_Type& b10rr, const vector_Type& adv ) const
{

    UInt dof( b10rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );

    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
            
                b10rr[i] +=   adv[n + h * nquadRho + i * nquadTheta * nquadRho]
                            * rphirho[k][n] 
                            * rphirho[j][n]
                            * rphitheta[k][h]
                            * rphitheta[j][h]
                            * r * dr
                            * Jacobian[i][h] 
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b10xx( const UInt& k, const UInt& j,
               vector_Type& b10xx, const vector_Type& adv ) const
{

    UInt dof( b10xx.blockSize(0) );
    Real Theta(M_Theta);

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
    
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
                
                b10xx[i] += adv[n + h * nquadRho + i * nquadTheta * nquadRho]
                            * xphirho[k][n] 
                            * xphirho[j][n]
                            * xphitheta[k][h]
                            * xphitheta[j][h]
                            * r * dr
                            * Jacobian[i][h] 
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b10tt( const UInt& k, const UInt& j,
               vector_Type& b10tt, const vector_Type& adv ) const
{

    UInt dof( b10tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
        // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
            
                b10tt[i] += adv[n + h * nquadRho + i * nquadTheta * nquadRho]
                            * thetaphirho[k][n] 
                            * thetaphirho[j][n]
                            * thetaphitheta[k][h]
                            * thetaphitheta[j][h]
                            * r * dr
                            * Jacobian[i][h] 
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b00rr( const UInt& k, const UInt& j,
               vector_Type& b00rr, const vector_Type& adv ) const
{

    UInt dof( b00rr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rdphirho( M_rdphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type rdphitheta( M_rdphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
        // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );        
                
                UInt ix( n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt ir( dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt it( 2*dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                
                b00rr[i] += (
                              ( adv[ix] * Dtheta[i][h]
                              + adv[ir] * Drtheta[i][h]
                              + adv[it] * Jtheta[i][h]
                              ) / r
                              * rphirho[k][n] 
                              * rphirho[j][n]
                              * rdphitheta[k][h]
                              * rphitheta[j][h] +
                              
                              ( adv[ix] * Dr[i][h]*r
                              + adv[ir] * Jr[i][h]
                              + adv[it] * Dthetar[i][h]
                              )
                              * rdphirho[k][n] 
                              * rphirho[j][n]
                              * rphitheta[k][h]
                              * rphitheta[j][h]
                            )
                            * r * dr
                            * Jacobian[i][h]
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b00tr( const UInt& k, const UInt& j,
               vector_Type& b00tr, const vector_Type& adv ) const
{

    UInt dof( b00tr.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
            
                UInt ix( n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt ir( dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt it( 2*dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                
                b00tr[i] += - ( adv[ix] * Dtheta[i][h]
                              + adv[ir] * Drtheta[i][h]
                              + adv[it] * Jtheta[i][h]
                              ) / r
                              * thetaphirho[k][n] 
                              * rphirho[j][n]
                              * thetaphitheta[k][h]
                              * rphitheta[j][h]
                              * r * dr
                              * Jacobian[i][h] 
                              * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b00rt( const UInt& k, const UInt& j,
               vector_Type& b00rt, const vector_Type& adv ) const
{

    UInt dof( b00rt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
            
                UInt ix( n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt ir( dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt it( 2*dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                
                b00rt[i] +=  ( adv[ix] * Dtheta[i][h]
                              + adv[ir] * Drtheta[i][h]
                              + adv[it] * Jtheta[i][h]
                              ) / r
                              * rphirho[k][n] 
                              * thetaphirho[j][n]
                              * rphitheta[k][h]
                              * thetaphitheta[j][h]
                              * r * dr
                              * Jacobian[i][h] 
                              * M_Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b00tt( const UInt& k, const UInt& j,
               vector_Type& b00tt, const vector_Type& adv ) const
{

    UInt dof( b00tt.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetadphirho( M_thetadphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
    MBMatrix_type thetadphitheta( M_thetadphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
           
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );        
                
                UInt ix( n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt ir( dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt it( 2*dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                
                b00tt[i] += (
                              ( adv[ix] * Dtheta[i][h]
                              + adv[ir] * Drtheta[i][h]
                              + adv[it] * Jtheta[i][h]
                              ) / r
                              * thetaphirho[k][n] 
                              * thetaphirho[j][n]
                              * thetadphitheta[k][h]
                              * thetaphitheta[j][h] +
                              
                              ( adv[ix] * Dr[i][h]*r
                              + adv[ir] * Jr[i][h]
                              + adv[it] * Dthetar[i][h]
                              )
                              * thetadphirho[k][n] 
                              * thetaphirho[j][n]
                              * thetaphitheta[k][h]
                              * thetaphitheta[j][h]
                            )
                            * r * dr
                            * Jacobian[i][h]
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

void
NSModalSpacePipe::
compute_b00xx( const UInt& k, const UInt& j,
               vector_Type& b00xx, const vector_Type& adv ) const
{

    UInt dof( b00xx.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type Dr( M_map->Dr() );
    MBMatrix_type Jr( M_map->Jr() );
    MBMatrix_type Dthetar( M_map->Dthetar() );
    MBMatrix_type Drtheta( M_map->Drtheta() );
    MBMatrix_type Dtheta( M_map->Dtheta() );
    MBMatrix_type Jtheta( M_map->Jtheta() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xdphirho( M_xdphirho );
    MBMatrix_type xphitheta( M_xphitheta );
    MBMatrix_type xdphitheta( M_xdphitheta );
    
    Real r(0);
    Real dr(0);
    UInt nquadRho( qrRho.nbQuadPt() );
    UInt nquadTheta( qrTheta.nbQuadPt() );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {
        for ( UInt h = 0; h != nquadTheta; ++h )
        {
            for ( UInt n = 0; n != nquadRho; ++n )
            {
                r = qrRho.quadPointCoor( n, 0 );
                dr = qrRho.weight( n );
            
                UInt ix( n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt ir( dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                UInt it( 2*dof*nquadTheta*nquadRho + n + h*nquadRho + i*nquadTheta*nquadRho );
                
                b00xx[i] += (
                              ( adv[ix] * Dtheta[i][h]
                              + adv[ir] * Drtheta[i][h]
                              + adv[it] * Jtheta[i][h]
                              ) / r
                              * xphirho[k][n] 
                              * xphirho[j][n]
                              * xdphitheta[k][h]
                              * xphitheta[j][h] +
                              
                              ( adv[ix] * Dr[i][h]*r
                              + adv[ir] * Jr[i][h]
                              + adv[it] * Dthetar[i][h]
                              )
                              * xdphirho[k][n] 
                              * xphirho[j][n]
                              * xphitheta[k][h]
                              * xphitheta[j][h]
                            )
                            * r * dr
                            * Jacobian[i][h]
                            * Theta * qrTheta.weight( h );
            } // n for-loop
        } // h for-loop
    } // i for-loop

    return ;
}

// ---------------------  PHI -----------------------------
void NSModalSpacePipe::
compute_Phix( const UInt& k, vector_Type& Phix ) const
{
    UInt dof( Phix.blockSize(0) );
    Real Theta(M_Theta);

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type xphirho( M_xphirho );
    MBMatrix_type xphitheta( M_xphitheta );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {    
        for ( UInt n = 0; n != qrRho.nbQuadPt(); ++n )
            for ( UInt h = 0; h != qrTheta.nbQuadPt(); ++h )
            {
                // note: r is taken into account 2 times: 1) multiplying of the beginning equation; 2) inner product.
                Phix[i] += xphirho[k][n] * xphitheta[k][h] *
                           Jacobian[i][h] *
                           qrRho.quadPointCoor( n, 0 ) * qrRho.weight( n ) *
                           Theta * qrTheta.weight( h );
            }
    }

    return ;
}

void NSModalSpacePipe::
compute_Phir( const UInt& k, vector_Type& Phir ) const
{
    UInt dof( Phir.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type rphirho( M_rphirho );
    MBMatrix_type rphitheta( M_rphitheta );

    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {    
        for ( UInt n = 0; n != qrRho.nbQuadPt(); ++n )
            for ( UInt h = 0; h != qrTheta.nbQuadPt(); ++h )
            {
                Phir[i] += rphirho[k][n] * rphitheta[k][h] *
                           Jacobian[i][h] *
                           qrRho.quadPointCoor( n, 0 ) * qrRho.weight( n ) *
                           Theta * qrTheta.weight( h );
            }
    }

    return ;
}

void NSModalSpacePipe::
compute_Phitheta( const UInt& k, vector_Type& Phitheta ) const
{
    UInt dof( Phitheta.blockSize(0) );
    Real Theta( M_Theta );

    QuadratureRule qrRho( *M_quadruleRho );
    QuadratureRule qrTheta( *M_quadruleTheta );
    MBMatrix_type Jacobian( M_map->Jacobian() );
    MBMatrix_type thetaphirho( M_thetaphirho );
    MBMatrix_type thetaphitheta( M_thetaphitheta );
            
    // normrho and normtheta are included in the basis
    for ( UInt i = 0; i != dof; ++i )
    {    
        for ( UInt n = 0; n != qrRho.nbQuadPt(); ++n )
            for ( UInt h = 0; h != qrTheta.nbQuadPt(); ++h )
            {
                Phitheta[i] += thetaphirho[k][n] * thetaphitheta[k][h] *
                           Jacobian[i][h] *
                           qrRho.quadPointCoor( n, 0 ) * qrRho.weight( n ) *
                           Theta * qrTheta.weight( h );
            }
    }

    return ;
}

} // namespace
