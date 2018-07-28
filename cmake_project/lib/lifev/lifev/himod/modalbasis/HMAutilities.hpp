#ifndef __HMAUTILITIES_HPP__
#define __HMAUTILITIES_HPP__

// ----------------------  UTILITY FUNCTIONS -------------------------------------

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
 
    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2 * M_PI );
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
                                    M_modalbasis->xphirho( j, nrho ) * normrho *
                                    M_modalbasis->xphitheta( j, ntheta ) * normtheta;
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
                                    M_modalbasis->rphirho( j, nrho ) * normrho *
                                    M_modalbasis->rphitheta( j, ntheta ) * normtheta;
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
                                    M_modalbasis->thetaphirho( j, nrho ) * normrho *
                                    M_modalbasis->thetaphitheta( j, ntheta ) * normtheta;
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
                                M_modalbasis->pphirho( j, nrho ) * normrho *
                                M_modalbasis->pphitheta( j, ntheta ) * normtheta;
                                                                      
                }
            }
        }
    }

    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
                    rho = M_modalbasis->Rho() * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 );
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
        interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
        interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
        CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
        interpCFE.update( M_velocityFespace->mesh()->element( s ), UPDATE_QUAD_NODES );

        x = interpCFE.quadNode( 0, 0 ); 

        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) //on quadrature node nquadTheta
        {

            theta = M_modalbasis->Theta() * M_modalbasis->qrTheta()->quadPointCoor( ntheta, 0 );

            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) //on quadrature node nqadRho
            {
                rho = M_modalbasis->Rho() * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 );

                fcoeff[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] = press( t, x, rho, theta, 0 );

                if ( s == (ndofpFE - 2) )
                {
                    fcoeff[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + (s + 1) *nquadTheta * nquadRho] +=
                                                                                                press( t, interpCFE.quadNode( 1, 0 ), rho, theta, 0 ); 
                }
            }
        }
    }
    
    return fcoeff;
}


template< typename mesh_type, typename matrix_type, typename vector_type>
Real NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
    
    for( UInt c( 0 ); c != 3; ++c ) // velocity component
    {
        for( UInt s ( 0 ); s != ( ndofuFE - 1 ); ++s ) //ndofFE - velocity
        {
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
                                w_rho * w_theta * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                M_modalbasis->map()->Jacobian()[nrho][ntheta] * M_modalbasis->Theta() *
                                w_x;
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
                    norm +=     ( fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] *
                                  fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho] +
                                  fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] *
                                  fun[3 * ndofuFE * nquadRho * nquadTheta + nrho + ntheta * nquadRho + ( s + 1 ) * nquadTheta * nquadRho] ) *
                                w_rho * w_theta * M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                M_modalbasis->map()->Jacobian()[nrho][ntheta] * M_modalbasis->Theta() *
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

/*
template< typename mesh_type, typename matrix_type, typename vector_type>
Real NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
evaluateHiModFunc( const vector_ptrType& fun, const Real& x, const Real& rho, const Real& theta )
{
    Real rhoh = rho / M_modalbasis->Rho();
    Real thetah = theta / M_modalbasis->Theta();
    std::vector<Real> ph;
    ph.push_back(rhoh);
    ph.push_back(thetah);

    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_fespace->refFEPtr()->refCoor(), std::vector<Real> (M_velocityFespace->refFEPtr()->nbDof(), 0) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_fespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );
    DOF DataFESpace( M_velocityFespace->dof() );
    UInt ndofFE = DataFESpace.numTotalDof();

    Real h   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );

    UInt i = std::floor( x / h );

    Real fem_i   = ( h * ( i + 1 ) - x ) / h;
    Real fem_ip1 = ( x - i * h ) / h;

    for ( UInt m( 0 ); m < M_modalbasis->mtot(); ++m )
    {
        return  ( ( *fun ) [ i + m * ndofFE] * fem_i +  ( *fun ) [ i + 1 + m * ndofFE] * fem_ip1 ) *
                M_modalbasis->gbRhoTheta()->evalSinglePoint( M_modalbasis->eigenvalues(m).lambda, m, ph ) / ( M_modalbasis->M_Rho ) * 1. / sqrt( M_modalbasis->M_Theta ); 
    }
}
*/

#endif
