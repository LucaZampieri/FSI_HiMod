#ifndef __HMAUTILITIESPIPE_HPP__
#define __HMAUTILITIESPIPE_HPP__

// ----------------------  UTILITY FUNCTIONS -------------------------------------

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
evaluateBase3DGrid( const vector_type& fun )
{
    UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();
    UInt ndofpFE = DataPressFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

    MapEpetra Map( nquadRho * nquadTheta * ( ndofuFE * 3 + ndofpFE ), Comm );
    vector_type fcoeff( Map, Unique );
    fcoeff *= 0.0;
 
    UInt lowerHalf = ( ndofuFE - 1 ) / 2;
    UInt index = 0;

    for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
    {
        if( remainder( s, 2 ) == 0 )
        {
            index = s / 2;
        }
        else
        {
            index = lowerHalf + ( s + 1 ) / 2;
        }
        for ( UInt j( 0 ); j != M_modalbasis->mx(); ++j ) // Ciclying over x-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                        fcoeff[nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +=
                                    fun[j * ndofuFE + index] *
                                    M_modalbasis->xphirho( j, nrho ) * M_modalbasis->xphitheta( j, ntheta );
                }
            }
        }
    }

    for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
    {
        if( remainder( s, 2 ) == 0 )
        {
            index = s / 2;
        }
        else
        {
            index = lowerHalf + ( s + 1 ) / 2;
        }
        for ( UInt j( 0 ); j != M_modalbasis->mr(); ++j ) // Ciclying over r-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                    fcoeff[ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +=
                                    fun[ndofuFE * M_modalbasis->mx() + j * ndofuFE + index] *
                                    M_modalbasis->rphirho( j, nrho ) * M_modalbasis->rphitheta( j, ntheta );
                }
            }
        }
    }
    
    for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
    {
        if( remainder( s, 2 ) == 0 )
        {
            index = s / 2;
        }
        else
        {
            index = lowerHalf + ( s + 1 ) / 2;
        }
        for ( UInt j( 0 ); j != M_modalbasis->mtheta(); ++j ) // Ciclying over theta-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                    fcoeff[2 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +=
                                    fun[ndofuFE  * ( M_modalbasis->mx() + M_modalbasis->mr() ) + j * ndofuFE + index] *
                                    M_modalbasis->thetaphirho( j, nrho ) * M_modalbasis->thetaphitheta( j, ntheta );
                }
            }
        }
    }
    
    for ( UInt s( 0 ); s != ndofpFE; ++s ) // on FE nodes
    {
        for ( UInt j( 0 ); j != M_modalbasis->mp(); ++j ) // Ciclying over p-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                    fcoeff[3*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho] +=
                                fun[ndofuFE*( M_modalbasis->mx()+M_modalbasis->mr()+M_modalbasis->mtheta() )+ j*ndofpFE + s] *
                                M_modalbasis->pphirho( j, nrho ) * M_modalbasis->pphitheta( j, ntheta );
                                                                      
                }
            }
        }
    }

    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
evaluateBase3DGrid( const function_Type& xVel, const function_Type& rVel, const function_Type& thetaVel, const function_Type& press, const Real& t )
{
    UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

    DOF DatauFESpace( M_velocityFespace->dof() );
    UInt ndofuFE = DatauFESpace.numTotalDof();
    DOF DatapFESpace( M_pressureFespace->dof() );
    UInt ndofpFE = DatapFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

    MapEpetra Map( nquadRho * nquadTheta * ( 3 * ndofuFE + ndofpFE ), Comm );
    vector_type fcoeff( Map, Unique );
    fcoeff *= 0.0;

    Real x, rho, theta, hStep;
    UInt sVel = 0;

    for ( UInt s( 0 ); s != ( ndofpFE - 1 ); ++s ) //on FE elements
    {
        // if you use P2 you have to modify this part of the code because vertices != nodes
        QuadratureRule interpQuad;
        interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ),
                                      M_velocityFespace->refFEPtr()->shape() );
        interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(),
                              std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
        CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
        interpCFE.update( M_velocityFespace->mesh()->element( s ), UPDATE_QUAD_NODES );

        x = interpCFE.quadNode( 0, 0 ); 
        hStep = ( interpCFE.quadNode( 1, 0 ) - x ) / 2;

        for( UInt w = 0; w != 2; ++w )
        {
            x += hStep * w;
            sVel = 2 * s + w; // Index of the velocity node
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) //on quadrature node nquadTheta
            {

                theta = M_modalbasis->Theta() * M_modalbasis->qrTheta()->quadPointCoor( ntheta, 0 );

                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) //on quadrature node nqadRho
                {
                    rho = M_modalbasis->fRho()(t,x,rho,theta,ntheta) * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 );
                    fcoeff[nrho + ntheta*nquadTheta + sVel*nquadTheta*nquadRho] = 
                                   xVel( t, x, rho, theta, 0 );
                    fcoeff[ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadTheta + sVel*nquadTheta*nquadRho] =
                                   rVel( t, x, rho, theta, 0 );
                    fcoeff[2*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadTheta + sVel*nquadTheta*nquadRho] = 
                                   thetaVel( t, x, rho, theta, 0 );

                    if ( s == (ndofpFE - 2) && w == 1 ) // one before the last node
                    {    
                        x = interpCFE.quadNode( 1, 0 );
                        fcoeff[nrho + ntheta*nquadRho + (sVel + 1)*nquadTheta*nquadRho] += 
                                      xVel( t, x, rho, theta, 0 );
                        fcoeff[ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + (sVel + 1)*nquadTheta*nquadRho] +=
                                      rVel( t, x, rho, theta, 0 );
                        fcoeff[2*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + (sVel + 1)*nquadTheta*nquadRho] +=
                                      thetaVel( t, x, rho, theta, 0 );
                        x -= hStep; 
                    }
                }
            }
        }
    }

    for ( UInt s( 0 ); s != ( ndofpFE - 1 ); ++s ) //on FE nodes
    {
        // if you use P2 you have to modify this part of the code
        QuadratureRule interpQuad;
        interpQuad.setDimensionShape( shapeDimension( M_pressureFespace->refFEPtr()->shape() ), M_pressureFespace->refFEPtr()->shape() );
        interpQuad.setPoints( M_pressureFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_pressureFespace->refFEPtr()->nbDof(), 0 ) );
        CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), getGeometricMap( *( M_pressureFespace->mesh() ) ), interpQuad );
        interpCFE.update( M_pressureFespace->mesh()->element( s ), UPDATE_QUAD_NODES );

        x = interpCFE.quadNode( 0, 0 ); 

        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) //on quadrature node nquadTheta
        {

            theta = M_modalbasis->Theta() * M_modalbasis->qrTheta()->quadPointCoor( ntheta, 0 );

            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) //on quadrature node nqadRho
            {
                rho = M_modalbasis->fRho()(t,x,rho,theta,ntheta) * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 );

                fcoeff[3*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho] = press( t,x,rho,theta,0 );

                if ( s == (ndofpFE - 2) )
                {
                    fcoeff[3*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + (s + 1)*nquadTheta*nquadRho] +=
                          press( t, interpCFE.quadNode( 1, 0 ), rho, theta, 0 ); 
                }
            }
        }
    }
    
    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
Real NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
normL2( const vector_type& fun, const std::string& solution )
{
    Real norm = 0.0;
    UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

    DOF DatauFESpace( M_velocityFespace->dof() );
    UInt ndofuFE = DatauFESpace.numTotalDof();
    DOF DatapFESpace( M_pressureFespace->dof() );
    UInt ndofpFE = DatapFESpace.numTotalDof();

    if( solution == "velocity" )
    {
        QuadratureRule interpQuad;
        interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
        interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
        CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
        interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );
    
        Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;
        Real w_x = 0.5 * h;
        Real w_rho;
        Real w_theta;
        UInt i(0);
        
        for( UInt c( 0 ); c != 3; ++c ) // velocity component
        {
            for( UInt s ( 0 ); s != ( ndofuFE - 1 ); ++s ) //ndofFE - velocity
            {
                i = (s%2==0)*s/2 + (s%2 != 0)*(s+ndofuFE)/2;
                for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) //on quadrature node nquadTheta
                {
                    for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) //on quadrature node nquadRho
                    {
                        w_rho = M_modalbasis->qrRho()->weight( nrho );
                        w_theta = M_modalbasis->qrTheta()->weight( ntheta );
                        norm +=     ( fun[c * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] *
                                      fun[c * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +
                                      fun[c * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] *
                                      fun[c * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] ) *
                                    M_modalbasis->map()->Jacobian()[i][ntheta] * // jacobian 
                                    M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * w_rho * // r dr 
                                    M_modalbasis->Theta() * w_theta * // dtheta
                                    w_x; // dx
                    }
                }
            }
        }
    }
    
    else if( solution == "pressure" )
    {
        QuadratureRule interpQuad;
        interpQuad.setDimensionShape( shapeDimension( M_pressureFespace->refFEPtr()->shape() ), M_pressureFespace->refFEPtr()->shape() );
        interpQuad.setPoints( M_pressureFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_pressureFespace->refFEPtr()->nbDof(), 0 ) );
        CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), getGeometricMap( *( M_pressureFespace->mesh() ) ), interpQuad );
        interpCFE.update( M_pressureFespace->mesh()->element (0), UPDATE_QUAD_NODES );
    
        Real h   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
        Real w_x = 0.5 * h;
        Real w_rho;
        Real w_theta;
        
        for( UInt s ( 0 ); s != ( ndofpFE - 1 ); ++s ) //ndofFE - pressure
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) //on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) //on quadrature node nquadRho
                {
                    w_rho = M_modalbasis->qrRho()->weight( nrho );
                    w_theta = M_modalbasis->qrTheta()->weight( ntheta );
                    norm += ( fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] *
                              fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +
                              fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] *
                              fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] ) *
                            w_rho * w_theta * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                            M_modalbasis->map()->Jacobian()[s][ntheta] * M_modalbasis->Theta() *
                            w_x;
                }
            }
        }
    }

    else
    {
        std::cout<<"Error: did you mean 'velocity' or 'pressure'?"<<std::endl;
    }
    return sqrt( norm );

}

template< typename mesh_type, typename matrix_type, typename vector_type>
Real NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
evaluateHiModFunc( const vector_ptrType& HMcoef, 
                   const Real& x, const Real& r, const Real& theta, const UInt& component )
{
    // FEM interpolation setup
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ),
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> (M_velocityFespace->refFEPtr()->nbDof(), 0) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );
    DOF DataFESpace( M_velocityFespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();
    UInt pndofFE = (ndofFE-1)/2 + 1;

    Real hp   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
    Real h   = hp/2;

    UInt i( std::floor( x / h ) );
    UInt ip( std::floor( x / hp ) );
    UInt ind_i;
    UInt ind_ip1;
    Real xp_i( ip*hp );
    Real xp_ip1( xp_i+hp );
    Real x_i;
    Real x_ip1;

    if( i%2==0 )
    {
        ind_i = i/2;
        ind_i = ( i==ndofFE-1? ndofFE-1 : ind_i );
        ind_ip1 = ( i==ndofFE-1? pndofFE-1 : ind_i+pndofFE );
        x_i = ( i==ndofFE-1? (i-1)*h : i*h );
        x_ip1 = x_i+h;
    }
    else
    {
        ind_i = (UInt)( (i+ndofFE)/2 );
        ind_ip1 = ind_i-pndofFE+1;
        x_i = i*h;
        x_ip1 = x_i+h;
    }

    if( ip == pndofFE-1 ) 
    {
        xp_ip1 = ip*h;
        xp_i = x_ip1-hp;
        ip = ip - 1; // right extremum
    }

    Real fem_i   = ( x_ip1 - x ) / h;
    Real fem_ip1 = ( x - x_i ) / h;
    Real pfem_i   = ( xp_ip1 - x ) / hp;
    Real pfem_ip1 = ( x - xp_i ) / hp;
    // ---------------------------------------------

    Real rh = r / M_modalbasis->fRho()( 0., x, 0., theta, 0 );
    rh = ( rh>1? 1.0:rh ); // rh>1 may happen with the FEM grid
    UInt mx( M_modalbasis->mx() );
    UInt mr( M_modalbasis->mr() );
    UInt mtheta( M_modalbasis->mtheta() );
    UInt mp( M_modalbasis->mp() );
    Real fxyz( 0. );

    switch( component )
    {
        case 0:
                {
                    Real ui( 0 );
                    Real uip1( 0 );

                    for ( UInt m( 0 ); m < mx; ++m )
                    {
                        ui = ( *HMcoef ) [ ind_i + m * ndofFE];
                        uip1 = ( *HMcoef ) [ ind_ip1 + m * ndofFE];
                        Real f = M_modalbasis->xGbRhoTheta()->basisFunction( m, rh, theta, 
                                                                               M_modalbasis->qrRho(), 
                                                                               M_modalbasis->xphirho() );
                        fxyz +=  ( ui * fem_i +  uip1 * fem_ip1 ) *
                                   M_modalbasis->xGbRhoTheta()->basisFunction( m, rh, theta, 
                                                                               M_modalbasis->qrRho(), 
                                                                               M_modalbasis->xphirho() );
                    }
                }break;

        case 1:
                {
                    Real ui( 0 );
                    Real uip1( 0 );

                    for ( UInt m( 0 ); m < mr; ++m )
                    {
                        ui = ( *HMcoef ) [ ind_i + (m+mx) * ndofFE];
                        uip1 = ( *HMcoef ) [ ind_ip1 + (m+mx) * ndofFE];

                        fxyz +=  ( ui * fem_i +  uip1 * fem_ip1 ) *
                                   M_modalbasis->rGbRhoTheta()->basisFunction( m, rh, theta, 
                                                                               M_modalbasis->qrRho(), 
                                                                               M_modalbasis->rphirho()
                                                                                ); 
                    }
                }break;

        case 2:
                {
                    Real ui( 0 );
                    Real uip1( 0 );

                    for ( UInt m( 0 ); m < mtheta; ++m )
                    {
                        ui = ( *HMcoef ) [ ind_i + (m+mx+mr) * ndofFE];
                        uip1 = ( *HMcoef ) [ ind_ip1 + (m+mx+mr) * ndofFE];

                        fxyz +=  ( ui * fem_i +  uip1 * fem_ip1 ) *
                                   M_modalbasis->thetaGbRhoTheta()->basisFunction( m, rh, theta, 
                                                                               M_modalbasis->qrRho(), 
                                                                               M_modalbasis->thetaphirho()
                                                                                ); 
                    }
                }break;

        case 3:
                {
                   Real pi( 0 );
                   Real pip1( 0 );

                   for ( UInt m( 0 ); m < mp; ++m )
                   {
                       pi = ( *HMcoef ) [ ip + m*pndofFE + (mx+mr+mtheta) * ndofFE];
                       pip1 = ( *HMcoef ) [ ip+1 + m*pndofFE + (mx+mr+mtheta) * ndofFE];

                       fxyz +=  ( pi * pfem_i +  pip1 * pfem_ip1 ) *
                                  M_modalbasis->pGbRhoTheta()->basisFunction( m, rh, theta, 
                                                                              M_modalbasis->qrRho(), 
                                                                              M_modalbasis->pphirho() 
                                                                               ); 
                   }
                }break;

        default:
                break;
    }

    return fxyz;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
std::vector<Real> NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
evaluateHiModFunc( const vector_ptrType& HMcoef,
                   const std::vector<Real>& x, const std::vector<Real>& y, const std::vector<Real>& z )
{

    // Cartesian to Polar coordinates and evaluation
    UInt n( y.size() );
    std::vector<Real> evals( 4*n, 0. );
    Real r( 0. );
    Real theta( 0. );

    for( UInt i( 0 ); i != n; ++i )
    {
        theta = std::atan( z[i]/y[i] ); // theta \in [-pi/2,pi/2]
        if( y[i] < 0 && z[i] >= 0 )
        {
            theta = theta + M_PI;
        }
        else if( y[i] <= 0 && z[i] < 0 )
        {
            theta = theta + M_PI;
        }
        else if( y[i] > 0 && z[i] < 0 )
        {
            theta = theta + 2*M_PI;
        }
        r = std::sqrt( y[i]*y[i] + z[i]*z[i] );
        evals[    i] = evaluateHiModFunc( HMcoef, x[i], r, theta, 0 );
        evals[  n+i] = evaluateHiModFunc( HMcoef, x[i], r, theta, 1 );
        evals[2*n+i] = evaluateHiModFunc( HMcoef, x[i], r, theta, 2 );
        evals[3*n+i] = evaluateHiModFunc( HMcoef, x[i], r, theta, 3 );
    }

    return evals;
    
}


#endif
