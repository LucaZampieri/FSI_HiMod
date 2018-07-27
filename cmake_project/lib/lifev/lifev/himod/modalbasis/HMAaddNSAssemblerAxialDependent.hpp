#ifndef __HMAADDNSASSEMBLERAXIALDEPENDENT_HPP__
#define __HMAADDNSASSEMBLERAXIALDEPENDENT_HPP__

/*
    Assemble the matrix of the problem, for the moment with x- dependent radius
*/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addStokesProblem( const matrix_ptrType& systemMatrix, const Real& nu, ReferenceMap& refMap, const Real& t, const Real& alpha, const bool& b )
{

	DOF DatauFESpace( M_velocityFespace->dof() );
    UInt ndofuFE = DatauFESpace.numTotalDof();
    DOF DatapFESpace( M_pressureFespace->dof() );
    UInt ndofpFE = DatapFESpace.numTotalDof();
    
    std::vector<Real> Reval( ndofuFE );
    std::vector<Real> dReval( ndofuFE );
    
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );
    CurrentFE interpCFE( *( M_velocityFespace->refFEPtr() ), getGeometricMap( *( M_velocityFespace->mesh() ) ), interpQuad );
    interpCFE.update( M_velocityFespace->mesh()->element (0), UPDATE_QUAD_NODES );

    Real h   = ( interpCFE.quadNode( 1, 0 ) - interpCFE.quadNode( 0, 0 ) ) / 2;
    
    // Evaluate the reference map and the radius along the mainstream
    // M_modalbasis->map()->evaluateAxialMap( h, ndofuFE, M_modalbasis->fRho(), M_modalbasis->fdRho(), M_velocityFespace->map() );
    
    // xx-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice

            vector_ptrType R000xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R101xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R110xx( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000xx( k, j, nu, alpha, *R000xx );     
            M_modalbasis->compute_r001xx( k, j, nu, *R001xx ); 
            M_modalbasis->compute_r010xx( k, j, nu, *R010xx );
            M_modalbasis->compute_r100xx( k, j, nu, *R100xx );  
            M_modalbasis->compute_r101xx( k, j, nu, *R101xx );
            M_modalbasis->compute_r110xx( k, j, nu, *R110xx );
            
            //Assemble the (j,k)1D FEMproblem
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000xx ) * phi_i * phi_j +
                            value( M_etufespace, *R001xx ) * phi_i * phi_j
                            + value( M_etufespace, *R010xx ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R100xx ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R101xx ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R110xx ) * dot( grad( phi_i ), grad( phi_j ) ) 
                            
                          )                         
	                     >> ( systemMatrix->block( j, k ) );
            }
        }
    }
    
	// rr-block
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice
            vector_ptrType R000rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R101rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R110rr( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000rr( k, j, nu, alpha, *R000rr );     
            M_modalbasis->compute_r001rr( k, j, nu, *R001rr ); 
            M_modalbasis->compute_r010rr( k, j, nu, *R010rr );
            M_modalbasis->compute_r100rr( k, j, nu, *R100rr );  
            M_modalbasis->compute_r101rr( k, j, nu, *R101rr );
            M_modalbasis->compute_r110rr( k, j, nu, *R110rr );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000rr ) * phi_i * phi_j
                            + value( M_etufespace, *R001rr ) * phi_i * phi_j
                            + value( M_etufespace, *R010rr ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R100rr ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R101rr ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R110rr ) * dot( grad( phi_i ), grad( phi_j ) ) 
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + j, M_modalbasis->mx() + k ) );
            }
        }
    }
    
    // tt-block
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            vector_ptrType R000tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R101tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R110tt( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000tt( k, j, nu, alpha, *R000tt );     
            M_modalbasis->compute_r001tt( k, j, nu, *R001tt ); 
            M_modalbasis->compute_r010tt( k, j, nu, *R010tt );
            M_modalbasis->compute_r100tt( k, j, nu, *R100tt );  
            M_modalbasis->compute_r101tt( k, j, nu, *R101tt );
            M_modalbasis->compute_r110tt( k, j, nu, *R110tt );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000tt ) * phi_i * phi_j 
                            + value( M_etufespace, *R001tt ) * phi_i * phi_j
                            + value( M_etufespace, *R010tt ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R100tt ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R101tt ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R110tt ) * dot( grad( phi_i ), grad( phi_j ) ) 
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
            }
        }
    }
    
    // xr-block
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            vector_ptrType R000xr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001xr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010xr( new vector_type( M_velocityFespace->map(), Repeated ) );


            M_modalbasis->compute_r000xr( k, j, nu, *R000xr );     
            M_modalbasis->compute_r001xr( k, j, nu, *R001xr ); 
            M_modalbasis->compute_r010xr( k, j, nu, *R010xr );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000xr ) * phi_i * phi_j
                            + value( M_etufespace, *R001xr ) * phi_i * phi_j
                            + value( M_etufespace, *R010xr ) * dot( grad( phi_i ), value( oneVector ) ) * phi_j 
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + j, k ) );
            }
        }
    }
    
    // xt-block
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            vector_ptrType R000xt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001xt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010xt( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000xt( k, j, nu, *R000xt );     
            M_modalbasis->compute_r001xt( k, j, nu, *R001xt ); 
            M_modalbasis->compute_r010xt( k, j, nu, *R010xt );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000xt ) * phi_i * phi_j
                            + value( M_etufespace, *R001xt ) * phi_i * phi_j
                            + value( M_etufespace, *R010xt ) * dot( grad( phi_i ), value( oneVector ) * phi_j )
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, k ) );
            }
        }
    }
    
    // rx-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000rx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100rx( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000rx( k, j, nu, *R000rx );     
            M_modalbasis->compute_r100rx( k, j, nu, *R100rx );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000rx ) * phi_i * phi_j
                            + value( M_etufespace, *R100rx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
	                     >> ( systemMatrix->block( j, M_modalbasis->mx() + k ) );
            }
        }
    }
    
    // rt-block
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000rt( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000rt( k, j, nu, *R000rt );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R000rt ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, M_modalbasis->mx() + k ) );
            }
        }
    }
    
    // tx-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000tx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100tx( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000tx( k, j, nu, *R000tx );     
            M_modalbasis->compute_r100tx( k, j, nu, *R100tx );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000tx ) * phi_i * phi_j
                            + value( M_etufespace, *R100tx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
	                     >> ( systemMatrix->block( j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
            }
        }
    }
    
    // tr-block
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000tr( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000tr( k, j, nu, *R000tr );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R000tr ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
            }
        }
    }

    // px-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
        {
            //j refers to the test function, k refers to the solution
            vector_ptrType R000px( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R001px( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R010px( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000px( k, j, *R000px );
            M_modalbasis->compute_r001px( k, j, *R001px ); 
            M_modalbasis->compute_r010px( k, j, *R010px );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            value( M_etufespace, *R000px ) * phi_i * phi_j
                            + value( M_etufespace, *R001px ) * phi_i * phi_j
                            + value( M_etufespace, *R010px ) * dot( grad( phi_i ), value( oneVector ) ) * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( j, M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
            }
        }
    }

    // pr-block
    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000pr( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000pr( k, j, *R000pr );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                           value( M_etufespace, *R000pr ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + j, M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
            }
        }
    }
    
    // ptheta-block
    for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000pt( new vector_type( M_velocityFespace->map(), Repeated ) );
            
            M_modalbasis->compute_r000pt( k, j, *R000pt );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                           value( M_etufespace, *R000pt ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j,
	                     							M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
            }
        }
    }
    
    // xp-block
    for ( UInt j = 0; j != M_modalbasis->mp(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000xp( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_ptrType R100xp( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000xp( k, j, *R000xp );     
            M_modalbasis->compute_r100xp( k, j, *R100xp );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R000xp ) * phi_i * phi_j
                            + value( M_etufespace, *R100xp ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + j, k ) );
            }
        }
    }
    
    // rp-block
    for ( UInt j = 0; j != M_modalbasis->mp(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000rp( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000rp( k, j, *R000rp );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R000rp ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + j,
	                     							M_modalbasis->mx() + k ) );
            }
        }
    }
    
    // thetap-block
    for ( UInt j = 0; j != M_modalbasis->mp(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
        {
            //j refers to the test function, k refers to the solution

            vector_ptrType R000tp( new vector_type( M_velocityFespace->map(), Repeated ) );

            M_modalbasis->compute_r000tp( k, j, *R000tp );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R000tp ) * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + j,
	                     							M_modalbasis->mx() + M_modalbasis->mr() + k ) );
            }
        }
    }
    
    // pp-block
    for ( UInt j = 0; j != M_modalbasis->mp(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
        {
            //j refers to the test function, k refers to the solution

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etpfespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etpfespace,
                           M_etpfespace,
                           
                            0 * phi_i * phi_j
                            
                          )                         
	                     >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + j,
	                     							M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
            }
        }
    }
    
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addrhs( const vector_ptrType& rhs, const Real& fx, const Real& fr, const Real& ftheta, const bool& b )
{
    //Cycling on x-block
    for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        vector_ptrType Phix( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->compute_Phix( k, *Phix );
        
        *Phix *= fx;

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phix ) * phi_i
                      )
                    >> ( rhs->block( k ) );
        }

    }
    
    //Cycling on r-block
    for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        vector_ptrType Phir( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->compute_Phir( k, *Phir );
        
        *Phir *= fr;

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phir ) * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + k ) );
        }

    }
    
    //Cycling on theta-block
    for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        vector_ptrType Phitheta( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->compute_Phitheta( k, *Phitheta );
        
        *Phitheta *= ftheta;

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phitheta ) * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + k ) );
        }

    }
    
    //Cycling on p-block : null force
    for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etpfespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etpfespace,
                        0 * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
        }

    }
    
    return;
}


#endif
