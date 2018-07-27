#ifndef __HMABCHANDLER_HPP__
#define __HMABCHANDLER_HPP__

#include <lifev/core/LifeV.hpp>

#include <Epetra_SerialComm.h>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
	
	return;
	
}

// Uneducated version
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addBC( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const BCdata& bc, const Real& t, const int trigonometricBasis )
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
	
	switch( trigonometricBasis )
	{
		case 0:
				break;
		case 1:
				diagonalizeBlocks( systemMatrix, rhs );
				break;
	}
	
	return;
	
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addBC( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const BCdata& bc, const Real& t, const int trigonometricBasis, const bool& b )
{
	switch( bc.xInflow().BCType() )
	{
		case 0:
			addDirichletBC_xIn( systemMatrix, rhs, bc.xInflow().BCData(), t, b );
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
			addDirichletBC_rIn( systemMatrix, rhs, bc.rInflow().BCData(), t, b );
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
			addDirichletBC_thetaIn( systemMatrix, rhs, bc.thetaInflow().BCData(), t, b );
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
			addNeumannBC_xOut( rhs, bc.xOutflow().BCData(), t, b );
			break;
		default:
			std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
	}
	
	switch( bc.rOutflow().BCType() )
	{
		case 0:
			addDirichletBC_rOut( systemMatrix, rhs, bc.rOutflow().BCData(), t, b );
			break;
		case 1:
			addNeumannBC_rOut( rhs, bc.rOutflow().BCData(), t, b );
			break;
		default:
			std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
	}
	
	switch( bc.thetaOutflow().BCType() )
	{
		case 0:
			addDirichletBC_thetaOut( systemMatrix, rhs, bc.thetaOutflow().BCData(), t, b );
			break;
		case 1:
			addNeumannBC_thetaOut( rhs, bc.thetaOutflow().BCData(), t, b );
			break;
		default:
			std::cout << "Error: dir / neu Outflow BC available only." << std::endl;
	}
	
	switch( trigonometricBasis )
	{
		case 0:
				break;
		case 1:
				diagonalizeBlocks( systemMatrix, rhs );
				break;
	}
	
	return;
	
}

template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_xIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->xFourierCoefficients( g, t );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        systemMatrix->setCoefficient( j * dof, j * dof,
        							1e+30 );
        rhs->setCoefficient( j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_xIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->xFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        systemMatrix->setCoefficient( j * dof, j * dof,
        							1e+30 );
        rhs->setCoefficient( j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}


template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_rIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + j * dof,
        							M_modalbasis->mx() * dof + j * dof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * dof + j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_rIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + j * dof,
        							M_modalbasis->mx() * dof + j * dof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * dof + j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}


template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_thetaIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
        							M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_thetaIn( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t, 0 );
    UInt dof = M_etufespace->dof().numTotalDof();
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
        							M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * dof + M_modalbasis->mr() * dof + j * dof,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_xOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->xFourierCoefficients( g, t );
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        systemMatrix->setCoefficient(  j * udof + pdof - 1,
        							   j * udof + pdof - 1, 1e+30 );
        rhs->setCoefficient( j * udof + pdof - 1,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_rOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t );
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
        								M_modalbasis->mx() * udof + j * udof + pdof - 1, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
        					1e+30 * FCoefficients_g[j] );
    }
}

// non constant radius
template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_rOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );

    Real x( interpCFE.quadNode( 0, 0 ) ); 
	Real hStep( ( interpCFE.quadNode( 1, 0 ) - x ) / 2 );
    
    std::vector<Real> FCoefficients_g;
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    FCoefficients_g = M_modalbasis->rFourierCoefficients( g, t, hStep * (udof-1) );
    
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
                                      M_modalbasis->mx() * udof + j * udof + pdof - 1,
        							1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_thetaOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t )
{
    std::vector<Real> FCoefficients_g;
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t );
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
        								M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
        					1e+30 * FCoefficients_g[j] );
    }
}

// non constant radius
template< typename mesh_type, typename matrix_type, typename vector_type >
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addDirichletBC_thetaOut( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element( 0 ), UPDATE_QUAD_NODES );

    Real x( interpCFE.quadNode( 0, 0 ) ); 
	Real hStep( ( interpCFE.quadNode( 1, 0 ) - x ) / 2 );
    
    std::vector<Real> FCoefficients_g;
    UInt udof = M_etufespace->dof().numTotalDof();
    UInt pdof = M_etpfespace->dof().numTotalDof();
    FCoefficients_g = M_modalbasis->thetaFourierCoefficients( g, t, hStep * (udof-1) );
    
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
                                      M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
        							1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1,
        					1e+30 * FCoefficients_g[j] );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_xIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
    UInt dof = M_etufespace->dof().numTotalDof();
    
    for( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        for ( UInt n = 0; n != M_modalbasis->qrRho()->nbQuadPt(); ++n )
            for ( UInt h = 0; h != M_modalbasis->qrTheta()->nbQuadPt(); ++h )
            {
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
                //Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ),
                //                                                       thetah, 0 );
                Real inverseRhat = M_modalbasis->Rho()*rn; 
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
					
                data += g( t, 0, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->xphirho( k, n ) * normrho * M_modalbasis->xphitheta( k, h ) * normtheta * inverseRhat * 
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( k * dof,
                             ( *rhs )( k * dof ) + data );
        data = 0;
    }
	
    return;
	
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_rIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                        M_modalbasis->rphirho( k, n ) * normrho * M_modalbasis->rphitheta( k, h ) * normtheta * inverseRhat *
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( M_modalbasis->mx() + k * dof,
                             ( *rhs )( M_modalbasis->mx() + k * dof ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_thetaIn( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                        M_modalbasis->thetaphirho( k, n )*normrho * M_modalbasis->thetaphitheta( k, h )*normtheta * inverseRhat *
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( M_modalbasis->mx() + M_modalbasis->mr() + k * dof,
                             ( *rhs )( M_modalbasis->mx() + M_modalbasis->mr() + k * dof ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_xOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->xphirho( k, n ) * normrho * M_modalbasis->xphitheta( k, h ) * normtheta * inverseRhat *
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( ( k + 1 ) * dof - 1,
                                 ( *rhs )( ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_xOut( const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->xphirho( k, n ) * normrho * M_modalbasis->xphitheta( k, h ) * normtheta * inverseRhat *
                        rn * M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, 0 ) * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
            }
        rhs->setCoefficient( ( k + 1 ) * dof - 1,
                             ( *rhs )( ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_rOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                CurrentFE interpCFE( *( M_pressureFespace->refFEPtr() ), getGeometricMap( *( M_pressureFespace->mesh() ) ), 
                                     interpQuad );
                interpCFE.update( M_pressureFespace->mesh()->element (0), UPDATE_QUAD_NODES );
                Real hx   = interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 );
                Real Lx = hx * ( ndofpFE - 1 );
				    
                Real thetah = M_modalbasis->qrTheta()->quadPointCoor( h, 0 );
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->rphirho( k, n ) * normrho * M_modalbasis->rphitheta( k, h ) * normtheta * inverseRhat *
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_rOut( const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->rphirho( k, n ) * normrho * M_modalbasis->rphitheta( k, h ) * normtheta * inverseRhat *
                        rn * M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, 0 ) * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_thetaOut( const vector_ptrType& rhs, const function_Type& g, const Real& t )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                Real inverseRhat = M_modalbasis->map()->inverseRhat()( t, 0, M_modalbasis->qrRho()->quadPointCoor( n, 0 ), 
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->thetaphirho( k, n )*normrho * M_modalbasis->thetaphitheta( k, h )*normtheta * inverseRhat *
                        rn * M_modalbasis->map()->Jacobian()[n][h] * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNeumannBC_thetaOut( const vector_ptrType& rhs, const function_Type& g, const Real& t, const bool& b )
{	
    // HYPOTHESIS: g is x and t-independent
    Real data = 0;

    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2. * M_PI );
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
                                                                       thetah, 0 );
                Real inverseThetahat = M_modalbasis->Theta() * thetah;
                Real rn = M_modalbasis->qrRho()->quadPointCoor( n, 0 );
					
                data += g( t, Lx, inverseRhat, inverseThetahat, 0 ) *
                        M_modalbasis->thetaphirho( k, n )*normrho * M_modalbasis->thetaphitheta( k, h )*normtheta * inverseRhat *
                        rn * M_modalbasis->map()->fJacobian()( t, Lx, inverseRhat, inverseThetahat, 0 ) * M_modalbasis->Theta() *
                        M_modalbasis->qrRho()->weight( n ) *
                        M_modalbasis->qrTheta()->weight( h );
        }
        rhs->setCoefficient( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1,
                             ( *rhs )( M_modalbasis->mx() + M_modalbasis->mr() + ( k + 1 ) * dof - 1 ) + data );
        data = 0;
    }
	
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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
	UInt mx_r = M_modalbasis->xGbRhoTheta()->mr1D();
	
	UInt mr = M_modalbasis->mr();
	UInt mr_r = M_modalbasis->rGbRhoTheta()->mr1D();
	
	UInt mtheta = M_modalbasis->mtheta();
	UInt mtheta_r = M_modalbasis->thetaGbRhoTheta()->mr1D();
	
	UInt mp = M_modalbasis->mp();
	UInt mp_r = M_modalbasis->pGbRhoTheta()->mr1D();

    for ( UInt j = 0; j != dof*mx_r; ++j )
    {
        systemMatrix->setCoefficient( j, j,
        							1. / ( h * h ) );
		rhs->setCoefficient( j, 0 );
    }
    
    for ( UInt j = 0; j != dof*mr_r; ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * dof + j,
        								M_modalbasis->mx() * dof + j,
    	    							1. / ( h * h ) );
    	rhs->setCoefficient( M_modalbasis->mx() * dof + j,
    	    				0 );
    }
    
    for ( UInt j = 0; j != dof*mtheta_r; ++j )
    {
        systemMatrix->setCoefficient( ( M_modalbasis->mx() + M_modalbasis->mr() ) * dof + j,
        								( M_modalbasis->mx() + M_modalbasis->mr() ) * dof + j,
	        							1. / ( h * h ) );
		rhs->setCoefficient( ( M_modalbasis->mx() + M_modalbasis->mr() ) * dof + j,
	        				0 );
    }
    
    for ( UInt j = 0; j != pdof*mp_r; ++j )
    {
        systemMatrix->setCoefficient( ( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * dof + j,
        								( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * dof + j,
    	    							1. / ( h * h ) );
		rhs->setCoefficient( ( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * dof + j,
	        				0 );
    }
}

#endif
