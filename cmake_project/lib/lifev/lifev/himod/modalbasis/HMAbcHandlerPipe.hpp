#ifndef __HMABCHANDLERPIPE_HPP__
#define __HMABCHANDLERPIPE_HPP__

#include <lifev/core/LifeV.hpp>

#include <Epetra_SerialComm.h>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addBC( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const BCdata& bc, const Real& t )
{
    switch( bc.xInflow().BCType() )
    {
        case 0:
            addDirichletBC_xIn( systemMatrix, rhs, bc.xInflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_xIn( rhs, bc.xInflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Inflow BC available only." << std::endl;
    }
                
    switch( bc.rInflow().BCType() )
    {
        case 0:
            addDirichletBC_rIn( systemMatrix, rhs, bc.rInflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_rIn( rhs, bc.rInflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Inflow BC available only." << std::endl;
    }
    
    switch( bc.thetaInflow().BCType() )
    {
        case 0:
            addDirichletBC_thetaIn( systemMatrix, rhs, bc.thetaInflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_thetaIn( rhs, bc.thetaInflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Inflow BC available only." << std::endl;
    }
    
    switch( bc.xOutflow().BCType() )
    {
        case 0:
            addDirichletBC_xOut( systemMatrix, rhs, bc.xOutflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_xOut( rhs, bc.xOutflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
    }
    
    switch( bc.rOutflow().BCType() )
    {
        case 0:
            addDirichletBC_rOut( systemMatrix, rhs, bc.rOutflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_rOut( rhs, bc.rOutflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
    }
    
    switch( bc.thetaOutflow().BCType() )
    {
        case 0:
            addDirichletBC_thetaOut( systemMatrix, rhs, bc.thetaOutflow().BCData(), t );
            break;
        case 1:
            addNeumannBC_thetaOut( rhs, bc.thetaOutflow().BCData(), t );
            break;
        default:
            std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
    }
    
    diagonalizeBlocks( systemMatrix, rhs );
    
    return;
    
}

template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_xIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mx(),0.0);
    FCoefficients_g = M_modalbasis->xFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    
    // get h to scale the matrix 
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );
    Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;
        
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        dirNodes[j] = j * dof;
/*        systemMatrix->setCoefficient( j * dof, j * dof,
                                      1./(h*h) );
        rhs->setCoefficient( j * dof,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}


template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_rIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mr(),0.0);
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    
    // get h to scale the matrix 
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );
    Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;
        
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        dirNodes[j] = M_modalbasis->mx() * dof + j * dof;
/*        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + j * dof,
                                      M_modalbasis->mx() * dof + j * dof,
                                      1./(h*h) );
        rhs->setCoefficient( M_modalbasis->mx() * dof + j * dof,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_thetaIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mtheta(),0.0);
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    
    // get h to scale the matrix 
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );
    Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;
        
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        dirNodes[j] = M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof;
/*        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
                                      M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof, 
                                      1./(h*h) );
        rhs->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_xOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    // get Lx to compute the modal coefficients and h to scale the matrix
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );

    Real x( interpCFE.quadNode( 0, 0 ) ); 
    Real h( ( interpCFE.quadNode( 1, 0 ) - x ) / 2 );
        
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mx(),0.0);
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    FCoefficients_g = M_modalbasis->xFourierCoefficients( g, t, h * (udof-1) );
    
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        dirNodes[j] = j * udof + pdof - 1;
/*        systemMatrix->setCoefficient(  j * udof + pdof - 1,
                                       j * udof + pdof - 1,
                                       1./(h*h) );
        rhs->setCoefficient( j * udof + pdof - 1,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}

// non constant radius
template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_rOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );

    Real x( interpCFE.quadNode( 0, 0 ) ); 
    Real h( ( interpCFE.quadNode( 1, 0 ) - x ) / 2 );
    
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mr(),0.0);
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t, h * (udof-1) );
    
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        dirNodes[j] = M_modalbasis->mx() * udof + j * udof + pdof - 1;
/*        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
                                      M_modalbasis->mx() * udof + j * udof + pdof - 1,
                                      1./(h*h) );
        rhs->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}

// non constant radius
template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addDirichletBC_thetaOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), 
                                                  M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), 
                          std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );

    Real x( interpCFE.quadNode( 0, 0 ) ); 
    Real h( ( interpCFE.quadNode( 1, 0 ) - x ) / 2 );
    
    std::vector<Real> FCoefficients_g;
    std::vector<UInt> dirNodes(M_modalbasis->mtheta(),0.0);
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t, h * (udof-1) );
    
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        dirNodes[j] = M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1;
/*        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
                                      M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
                                      1./(h*h) );
        rhs->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
                             1./(h*h) * FCoefficients_g[j] );
*/    }
    systemMatrix->diagonalize( dirNodes, 1./(h*h), *rhs, FCoefficients_g );
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_xIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;

    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ),
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                    
                data += g( t, 0, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->xphirho( k, n ) * M_modalbasis->xphitheta( k, h ) * 
                        M_modalbasis->map()->Jacobian()[0][h] * 
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( k * dof,
                             ( *rhs )( k * dof ) + data );
        data = 0;
    }
    
    return;
    
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_rIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;

    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                    
                data += g( t, 0, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->rphirho( k, n ) * M_modalbasis->rphitheta( k, h ) *
                        M_modalbasis->map()->Jacobian()[0][h] * 
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( M_modalbasis->mx() + k * dof,
                             ( *rhs )( M_modalbasis->mx() + k * dof ) + data );
        data = 0;
    }
    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_thetaIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;

    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ),
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                    
                data += g( t, 0, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->thetaphirho( k, n ) * M_modalbasis->thetaphitheta( k, h ) *
                        M_modalbasis->map()->Jacobian()[0][h] * 
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( M_modalbasis->mx() + M_modalbasis->mr() + k * dof,
                             ( *rhs )( M_modalbasis->mx() + M_modalbasis->mr() + k * dof ) + data );
        data = 0;
    }
    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_xOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;

    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                DOF DatapFESpace( M_pressureFespace->dof() );
                UInt ndofpFE = DatapFESpace.numTotalDof();
                QuadratureRule interpQuad;
                interpQuad.setDimensionShape( shapeDimension( M_pressureFespace->refFEPtr()->shape() ), 
                                                              M_pressureFespace->refFEPtr()->shape() );
                interpQuad.setPoints( M_pressureFespace->refFEPtr()->refCoor(), 
                                      std::vector<Real> ( M_pressureFespace->refFEPtr()->nbDof(), 0 ) );
                CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), getGeometricMap( *( M_pressureFespace->mesh() ) ), 
                                     interpQuad );
                interpCFE.update( M_pressureFespace->mesh()->element (0), UPDATE_QUAD_NODES );
                Real hx   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
                Real Lx = hx * ( ndofpFE - 1 );
                    
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, Lx, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, h );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                    
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->xphirho( k, n ) * M_modalbasis->xphitheta( k, h ) *
                        M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, h ) *
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( ( k + 1 ) * dof - 1,
                             ( *rhs )( ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_rOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;

    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                DOF DatapFESpace( M_pressureFespace->dof() );
                UInt ndofpFE = DatapFESpace.numTotalDof();
                QuadratureRule interpQuad;
                interpQuad.setDimensionShape( shapeDimension( M_pressureFespace->refFEPtr()->shape() ), 
                                                              M_pressureFespace->refFEPtr()->shape() );
                interpQuad.setPoints( M_pressureFespace->refFEPtr()->refCoor(), 
                                      std::vector<Real> ( M_pressureFespace->refFEPtr()->nbDof(), 0 ) );
                                      CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), 
                                                           getGeometricMap( *( M_pressureFespace->mesh() ) ), interpQuad );
                interpCFE.update( M_pressureFespace->mesh()->element (0), UPDATE_QUAD_NODES );
                Real hx   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
                Real Lx = hx * ( ndofpFE - 1 );
                    
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, Lx, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, h );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                    
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->rphirho( k, n ) * M_modalbasis->rphitheta( k, h ) *
                        M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, h ) *
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addNeumannBC_thetaOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{    
    Real data = 0;
    
    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                DOF DatapFESpace( M_pressureFespace->dof() );
                UInt ndofpFE = DatapFESpace.numTotalDof();
                QuadratureRule interpQuad;
                interpQuad.setDimensionShape( shapeDimension( M_pressureFespace->refFEPtr()->shape() ), 
                                                              M_pressureFespace->refFEPtr()->shape() );
                interpQuad.setPoints( M_pressureFespace->refFEPtr()->refCoor(), 
                                      std::vector<Real> ( M_pressureFespace->refFEPtr()->nbDof(), 0 ) );
                CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), getGeometricMap( *( M_pressureFespace->mesh() ) ), 
                                     interpQuad );
                interpCFE.update( M_pressureFespace->mesh()->element (0), UPDATE_QUAD_NODES );
                Real hx   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
                Real Lx = hx * ( ndofpFE - 1 );
                    
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, Lx, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, h );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                    
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->thetaphirho( k, n ) * M_modalbasis->thetaphitheta( k, h ) *
                        M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, h ) *
                        rn * M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->Theta() * M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
diagonalizeBCmatrix( const matrix_ptrType& matrix, const vector_ptrType& rhs )
{
    UInt mx = M_modalbasis->mx();
    UInt mx_r = M_modalbasis->xGbRhoTheta()->mr1D();
    
    for ( UInt j = 0; j != mx_r; ++j )
    {
        matrix->setCoefficient( j, j,
                                    1. );
        rhs->setCoefficient( j, 0 );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
diagonalizeBlocks( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs )
{
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );

    Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;

    UInt dof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    
    UInt mx = M_modalbasis->mx();
    UInt mr = M_modalbasis->mr();
    UInt mtheta = M_modalbasis->mtheta();
    UInt mp = M_modalbasis->mp();

    for( UInt m(0); m!=mx; ++m )
    {
        if( (M_modalbasis->xGbRhoTheta()->eigenValues())[m].order == -1 )
        {
            UInt offset( m*dof );
            for ( UInt j = 0; j != dof; ++j )
            {
                systemMatrix->setCoefficient( offset+j, offset+j,
                                              1. / ( h * h ) );
                rhs->setCoefficient( offset+j, 0 );
            }
        }
    }
    
    for( UInt m(0); m!=mr; ++m )
    {
        if( (M_modalbasis->rGbRhoTheta()->eigenValues())[m].order == -1 )
        {
            UInt offset( (mx+m)*dof );
            for ( UInt j = 0; j != dof; ++j )
            {
                systemMatrix->setCoefficient( offset + j,
                                              offset + j,
                                              1. / ( h * h ) );
                rhs->setCoefficient( offset + j,
                                     0 );
            }
        }
    }
    
    for( UInt m(0); m!=mtheta; ++m )
    {
        if( (M_modalbasis->thetaGbRhoTheta()->eigenValues())[m].order == -1 )
        {
            UInt offset( (mx+mr+m)*dof );
            for ( UInt j = 0; j != dof; ++j )
            {
                systemMatrix->setCoefficient( offset + j,
                                              offset + j,
                                              1. / ( h * h ) );
                rhs->setCoefficient( offset + j,
                                     0 );
            }
        }
    }
    
    for( UInt m(0); m!=mp; ++m )
    {
        if( (M_modalbasis->pGbRhoTheta()->eigenValues())[m].order == -1 )
        {
            UInt offset( (mx+mr+mtheta)*dof+m*pdof );
            for ( UInt j = 0; j != pdof; ++j )
            {
                systemMatrix->setCoefficient( offset + j,
                                              offset + j,
                                              1. / ( h * h ) );
                rhs->setCoefficient( offset + j,
                                     0 );
            }
        }
    }
}

#endif
