#ifndef __HMAADDNSASSEMBLERPIPE_HPP__
#define __HMAADDNSASSEMBLERPIPE_HPP__

/*********************************************************************************
*                               STOKES PROBLEM                                   *
**********************************************************************************/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addStokesProblem( const matrix_ptrType& systemMatrix, const Real& nu, ReferenceMap& refMap, const Real& t, const Real& alpha )
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

    // Enable OpenMP nested parallel regions
    omp_set_nested(1);
    omp_lock_t( lockMS );
    omp_init_lock( &lockMS );

    int j, k;
    int mx( M_modalbasis->mx() );
    int mr( M_modalbasis->mr() );
    int mtheta( M_modalbasis->mtheta() );
    int mp( M_modalbasis->mp() );
    //NSModalSpacePipe ms( *M_modalbasis );

    #pragma omp parallel sections
    {
    // xx-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) default(shared) private( j,k )
    for ( j = 0; j < mx; ++j )
    {
        for ( k = 0; k < mx; ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice
            omp_set_lock( &lockMS );
            vector_ptrType R00xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00xx *= 0;
            M_modalbasis->compute_r00xx( k, j, nu, alpha, *R00xx );     

            omp_set_lock( &lockMS );
            vector_ptrType R01xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01xx *= 0;
            M_modalbasis->compute_r01xx( k, j, nu, *R01xx ); 

            omp_set_lock( &lockMS );
            vector_ptrType R10xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10xx *= 0;
            M_modalbasis->compute_r10xx( k, j, nu, *R10xx );

            omp_set_lock( &lockMS );
            vector_ptrType R11xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R11xx *= 0;
            M_modalbasis->compute_r11xx( k, j, nu, *R11xx );
            
            //Assemble the (j,k)1D FEMproblem
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00xx ) * phi_i * phi_j
                            + value( M_etufespace, *R01xx ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R10xx ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R11xx ) * dot( grad( phi_i ), grad( phi_j ) ) 
                          )                         
                         >> ( systemMatrix->block( j, k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "xx block closed." << std::endl; 
    } //section

    // rr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice
            omp_set_lock( &lockMS );
            vector_ptrType R00rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00rr *= 0;
            M_modalbasis->compute_r00rr( k, j, nu, alpha, *R00rr );     

            omp_set_lock( &lockMS );
            vector_ptrType R01rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01rr *= 0;
            M_modalbasis->compute_r01rr( k, j, nu, *R01rr ); 

            omp_set_lock( &lockMS );
            vector_ptrType R10rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10rr *= 0;
            M_modalbasis->compute_r10rr( k, j, nu, *R10rr );

            omp_set_lock( &lockMS );
            vector_ptrType R11rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R11rr *= 0;
            M_modalbasis->compute_r11rr( k, j, nu, *R11rr );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00rr ) * phi_i * phi_j
                            + value( M_etufespace, *R01rr ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R10rr ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R11rr ) * dot( grad( phi_i ), grad( phi_j ) ) 
                            
                          )                         
                         >> ( systemMatrix->block( mx + j, mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "rr block closed." << std::endl; 
    } //section
    
    // tt-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice
            omp_set_lock( &lockMS );
            vector_ptrType R00tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00tt *= 0;
            M_modalbasis->compute_r00tt( k, j, nu, alpha, *R00tt );     

            omp_set_lock( &lockMS );
            vector_ptrType R01tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01tt *= 0;
            M_modalbasis->compute_r01tt( k, j, nu, *R01tt ); 

            omp_set_lock( &lockMS );
            vector_ptrType R10tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10tt *= 0;
            M_modalbasis->compute_r10tt( k, j, nu, *R10tt );

            omp_set_lock( &lockMS );
            vector_ptrType R11tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R11tt *= 0;
            M_modalbasis->compute_r11tt( k, j, nu, *R11tt );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00tt ) * phi_i * phi_j 
                            + value( M_etufespace, *R01tt ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *R10tt ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + value( M_etufespace, *R11tt ) * dot( grad( phi_i ), grad( phi_j ) ) 
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j,
                                                       mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "tt block closed." << std::endl; 
    } //section

    // xr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mx; ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            omp_set_lock( &lockMS );
            vector_ptrType R00xr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00xr *= 0;
            M_modalbasis->compute_r00xr( k, j, nu, *R00xr );      

            omp_set_lock( &lockMS );
            vector_ptrType R01xr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01xr *= 0;
            M_modalbasis->compute_r01xr( k, j, nu, *R01xr );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00xr ) * phi_i * phi_j
                            + value( M_etufespace, *R01xr ) * dot( grad( phi_i ), value( oneVector ) ) * phi_j 
                            
                          )                         
                         >> ( systemMatrix->block( mx + j, k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "xr block closed." << std::endl; 
    }
    
    // xt-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mx; ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            omp_set_lock( &lockMS );
            vector_ptrType R00xt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00xt *= 0;
            M_modalbasis->compute_r00xt( k, j, nu, *R00xt );      

            omp_set_lock( &lockMS );
            vector_ptrType R01xt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01xt *= 0;
            M_modalbasis->compute_r01xt( k, j, nu, *R01xt );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00xt ) * phi_i * phi_j
                            + value( M_etufespace, *R01xt ) * dot( grad( phi_i ), value( oneVector ) * phi_j )
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j, k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "xt block closed." << std::endl; 
    } // section
    
    // rx-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mx; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00rx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00rx *= 0;
            M_modalbasis->compute_r00rx( k, j, nu, *R00rx );     

            omp_set_lock( &lockMS );
            vector_ptrType R10rx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10rx *= 0;
            M_modalbasis->compute_r10rx( k, j, nu, *R10rx );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00rx ) * phi_i * phi_j
                            + value( M_etufespace, *R10rx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
                         >> ( systemMatrix->block( j, mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "rx block closed." << std::endl;
    } // section
    
    // rt-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00rt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00rt *= 0;
            M_modalbasis->compute_r00rt( k, j, nu, *R00rt );

            omp_set_lock( &lockMS );
            vector_ptrType R01rt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01rt *= 0;
            M_modalbasis->compute_r01rt( k, j, nu, *R01rt );

            omp_set_lock( &lockMS );
            vector_ptrType R10rt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10rt *= 0;
            M_modalbasis->compute_r10rt( k, j, nu, *R10rt );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R00rt ) * phi_i * phi_j +
                           value( M_etufespace, *R10rt ) * phi_i * dot( grad( phi_j ), value( oneVector ) ) +
                           value( M_etufespace, *R01rt ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j, mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "rt block closed." << std::endl; 
    } // section
    
    // tx-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mx; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00tx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00tx *= 0;
            M_modalbasis->compute_r00tx( k, j, nu, *R00tx );     

            omp_set_lock( &lockMS );
            vector_ptrType R10tx( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10tx *= 0;
            M_modalbasis->compute_r10tx( k, j, nu, *R10tx );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00tx ) * phi_i * phi_j
                            + value( M_etufespace, *R10tx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
                         >> ( systemMatrix->block( j, mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "tx block closed." << std::endl; 
    } //section
    
    // tr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00tr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00tr *= 0;
            M_modalbasis->compute_r00tr( k, j, nu, *R00tr );

            omp_set_lock( &lockMS );
            vector_ptrType R01tr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01tr *= 0;
            M_modalbasis->compute_r01tr( k, j, nu, *R01tr );

            omp_set_lock( &lockMS );
            vector_ptrType R10tr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10tr *= 0;
            M_modalbasis->compute_r10tr( k, j, nu, *R10tr );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R00tr ) * phi_i * phi_j +
                           value( M_etufespace, *R10tr ) * phi_i * dot( grad( phi_j ), value( oneVector ) ) +
                           value( M_etufespace, *R01tr ) * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            
                          )                         
                         >> ( systemMatrix->block( mx + j, mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    std::cout << "tr block closed." << std::endl; 
    }

    // px-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mx; ++j )
    {
        for ( k = 0; k < mp; ++k )
        {
            //j refers to the test function, k refers to the solution
            omp_set_lock( &lockMS );
            vector_ptrType R00px( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00px *= 0;
            M_modalbasis->compute_r00px( k, j, *R00px ); 

            omp_set_lock( &lockMS );
            vector_ptrType R01px( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R01px *= 0;
            M_modalbasis->compute_r01px( k, j, *R01px );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            value( M_etufespace, *R00px ) * phi_i * phi_j
                            + value( M_etufespace, *R01px ) * dot( grad( phi_i ), value( oneVector ) ) * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( j, mx + mr + mtheta + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "px block closed." << std::endl; 
    } //section

    // pr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mp; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00pr( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00pr *= 0;

            M_modalbasis->compute_r00pr( k, j, *R00pr );
            
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                           value( M_etufespace, *R00pr ) * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + j,
                                                       mx + mr + mtheta + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "pr block closed." << std::endl; 
    } //section
    
    // ptheta-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mp; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00pt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00pt *= 0;
            
            M_modalbasis->compute_r00pt( k, j, *R00pt );
            
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                           value( M_etufespace, *R00pt ) * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j,
                                                       mx + mr + mtheta + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "pt block closed." << std::endl; 
    } //section
    
    // xp-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mp; ++j )
    {
        for ( k = 0; k < mx; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00xp( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00xp *= 0;
            M_modalbasis->compute_r00xp( k, j, *R00xp );     

            omp_set_lock( &lockMS );
            vector_ptrType R10xp( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R10xp *= 0;
            M_modalbasis->compute_r10xp( k, j, *R10xp );  

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            value( M_etufespace, *R00xp ) * phi_i * phi_j
                            + value( M_etufespace, *R10xp ) * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + mtheta + j, k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "xp block closed." << std::endl; 
    } //section
    
    // rp-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mp; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00rp( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00rp *= 0;
            M_modalbasis->compute_r00rp( k, j, *R00rp );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R00rp ) * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + mtheta + j,
                                                       mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "rp block closed." << std::endl; 
    } //section
    
    // thetap-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mp; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType R00tp( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *R00tp *= 0;
            M_modalbasis->compute_r00tp( k, j, *R00tp );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                           value( M_etufespace, *R00tp ) * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + mtheta + j,
                                                       mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "tp block closed." << std::endl; 
    } //section
    
    // pp-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mp; ++j )
    {
        for ( k = 0; k < mp; ++k )
        {
            //j refers to the test function, k refers to the solution

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etpfespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etpfespace,
                           M_etpfespace,
                           
                            0 * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + mtheta + j,
                                                       mx + mr + mtheta + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
   std::cout << "pp block closed." << std::endl; 
    } //section
    } //sections
   omp_destroy_lock( &lockMS );
    return;
}

/*********************************************************************************
*                             NON-LINEAR TERM (NS)                               *
**********************************************************************************/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addAdvection( const matrix_ptrType& systemMatrix, const vector_type& adv )
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
    
    // Enable OpenMP nested parallel regions
    omp_set_nested(1);
    omp_lock_t( lockMS );
    omp_init_lock( &lockMS );

    int j, k;
    int mx( M_modalbasis->mx() );
    int mr( M_modalbasis->mr() );
    int mtheta( M_modalbasis->mtheta() );
    int mp( M_modalbasis->mp() );

    #pragma omp parallel sections
    {
    // xx-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) default(shared) private( j,k )
    for ( j = 0; j < mx; ++j )
    {
        for ( k = 0; k < mx; ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice

            omp_set_lock( &lockMS );
            vector_ptrType b00xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my00Adv( adv );
            omp_unset_lock( &lockMS );
            *b00xx *= 0;
            M_modalbasis->compute_b00xx( k, j, *b00xx, my00Adv );

            omp_set_lock( &lockMS );
            vector_ptrType b10xx( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my10Adv( adv );
            omp_unset_lock( &lockMS );
            *b10xx *= 0;
            M_modalbasis->compute_b10xx( k, j, *b10xx, my10Adv );
            
            //Assemble the (j,k)1D FEMproblem
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *b00xx ) * phi_i * phi_j
                            + value( M_etufespace, *b10xx ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                          )                         
                         >> ( systemMatrix->block( j, k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    } //section
    
    // rr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice
            omp_set_lock( &lockMS );
            vector_ptrType b00rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my00Adv( adv );
            omp_unset_lock( &lockMS );
            *b00rr *= 0;
            M_modalbasis->compute_b00rr( k, j, *b00rr, my00Adv );

            omp_set_lock( &lockMS );
            vector_ptrType b10rr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my10Adv( adv );
            omp_unset_lock( &lockMS );
            *b10rr *= 0;
            M_modalbasis->compute_b10rr( k, j, *b10rr, my10Adv );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *b00rr ) * phi_i * phi_j
                            + value( M_etufespace, *b10rr ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            
                          )                         
                         >> ( systemMatrix->block( mx + j, mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    } //section
    
    // tt-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution

            //For each blocks compute the integral coefficient on the slice

            omp_set_lock( &lockMS );
            vector_ptrType b00tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my00Adv( adv );
            omp_unset_lock( &lockMS );
            *b00tt *= 0;
            M_modalbasis->compute_b00tt( k, j, *b00tt, my00Adv );

            omp_set_lock( &lockMS );
            vector_ptrType b10tt( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my10Adv( adv );
            omp_unset_lock( &lockMS );
            *b10tt *= 0;
            M_modalbasis->compute_b10tt( k, j, *b10tt, my10Adv );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *b00tt ) * phi_i * phi_j 
                            + value( M_etufespace, *b10tt ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j,
                                                       mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    } //section
    
    // xr-block
    // xt-block
    // rx-block
    // rt-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mtheta; ++j )
    {
        for ( k = 0; k < mr; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_type my00Adv( adv );
            vector_ptrType b00rt( new vector_type( M_velocityFespace->map(), Repeated ) );
            omp_unset_lock( &lockMS );
            *b00rt *= 0;
            M_modalbasis->compute_b00rt( k, j, *b00rt, my00Adv );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *b00rt ) * phi_i * phi_j 
                            
                          )                         
                         >> ( systemMatrix->block( mx + mr + j, mx + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    } //section
    
    // tx-block
    // tr-block
    #pragma omp section
    {
    #pragma omp parallel for collapse(2) schedule( static ) private( j,k )
    for ( j = 0; j < mr; ++j )
    {
        for ( k = 0; k < mtheta; ++k )
        {
            //j refers to the test function, k refers to the solution

            omp_set_lock( &lockMS );
            vector_ptrType b00tr( new vector_type( M_velocityFespace->map(), Repeated ) );
            vector_type my00Adv( adv );
            omp_unset_lock( &lockMS );
            *b00tr *= 0;
            M_modalbasis->compute_b00tr( k, j, *b00tr, my00Adv );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                omp_set_lock( &lockMS );
                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                           value( M_etufespace, *b00tr ) * phi_i * phi_j
                            
                          )                         
                         >> ( systemMatrix->block( mx + j, mx + mr + k ) );
                omp_unset_lock( &lockMS );
            }
        }
    }
    } //section
    } //sections

    // px-block
    // pr-block
    // ptheta-block
    // xp-block
    // rp-block
    // thetap-block
    // pp-block
   omp_destroy_lock( &lockMS );
    
    return;
}

/*********************************************************************************
*                                     RHS                                        *
**********************************************************************************/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addrhs( const vector_ptrType& rhs, const Real& fx, const Real& fr, const Real& ftheta )
{
    // Enable OpenMP nested parallel regions
    omp_set_nested(1);
    omp_lock_t( lockMS );
    omp_init_lock( &lockMS );

    int k;
    int mx( M_modalbasis->mx() );
    int mr( M_modalbasis->mr() );
    int mtheta( M_modalbasis->mtheta() );
    int mp( M_modalbasis->mp() );

    #pragma omp parallel sections
    {
    #pragma omp section
    {
    //Cycling on x-block
    #pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mx; ++k )
    {
        omp_set_lock( &lockMS );
        vector_ptrType Phix( new vector_type( M_velocityFespace->map(), Repeated ) );
        omp_unset_lock( &lockMS );
        M_modalbasis->compute_Phix( k, *Phix );
        
        *Phix *= fx;

        {
            using namespace ExpressionAssembly;
            
            omp_set_lock( &lockMS );
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phix ) * phi_i
                      )
                    >> ( rhs->block( k ) );
            omp_unset_lock( &lockMS );
        }

    }
    } //section
    
    //Cycling on r-block
    #pragma omp section
    {
    #pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mr; ++k )
    {
        omp_set_lock( &lockMS );
        vector_ptrType Phir( new vector_type( M_velocityFespace->map(), Repeated ) );
        omp_unset_lock( &lockMS );
        M_modalbasis->compute_Phir( k, *Phir );
        
        *Phir *= fr;

        {
            using namespace ExpressionAssembly;
            
            omp_set_lock( &lockMS );
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phir ) * phi_i
                      )
                    >> ( rhs->block( mx + k ) );
            omp_unset_lock( &lockMS );
        }

    }
    } //section
    
    #pragma omp section
    {
    //Cycling on theta-block
    #pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mtheta; ++k )
    {
        omp_set_lock( &lockMS );
        vector_ptrType Phitheta( new vector_type( M_velocityFespace->map(), Repeated ) );
        omp_unset_lock( &lockMS );
        M_modalbasis->compute_Phitheta( k, *Phitheta );
        
        *Phitheta *= ftheta;

        {
            using namespace ExpressionAssembly;
            
            omp_set_lock( &lockMS );
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *Phitheta ) * phi_i
                      )
                    >> ( rhs->block( mx + mr + k ) );
            omp_unset_lock( &lockMS );
        }

    }
    } //section
    
    //Cycling on p-block : null force
    #pragma omp section
    {
    #pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mp; ++k )
    {
        {
            using namespace ExpressionAssembly;
            
            omp_set_lock( &lockMS );
            integrate ( elements( M_etpfespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etpfespace,
                        0 * phi_i
                      )
                    >> ( rhs->block( mx + mr + mtheta + k ) );
            omp_unset_lock( &lockMS );
        }

    }
    } //section
    } //sections
    
   omp_destroy_lock( &lockMS );
    return;
}


/*********************************************************************************
*                                     RHS                                        *
**********************************************************************************/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
addrhs( const vector_ptrType& rhs, const vector_ptrType& f_interpolated, const vector_ptrType& uold_interpolated )
{

    // Enable OpenMP nested parallel regions
    /*omp_set_nested(1);
    omp_lock_t( lockMS );
    omp_init_lock( &lockMS );
*/
    int k, s;
    int mx( M_modalbasis->mx() );
    int mr( M_modalbasis->mr() );
    int mtheta( M_modalbasis->mtheta() );
    int mp( M_modalbasis->mp() );
    int numTotalDof( M_velocityFespace->dofPtr()->numTotalDof() );
    int pnumTotalDof( M_pressureFespace->dofPtr()->numTotalDof() );

    //#pragma omp parallel sections private( k,s )
    {
    //#pragma omp section
    {
    //Cycling on x-blocks
    //#pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mx; ++k )
    {
        //omp_setlock( &lockMS );
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        UInt fInd( f_interpolated->block( k )->firstIndex() );
        UInt uInd( uold_interpolated->block( k )->firstIndex() );
        vector_type f( *f_interpolated );
        vector_type u( *uold_interpolated );
        //omp_unsetlock( &lockMS );

        //#pragma omp parallel for schedule( static ) default(shared) private( s )
        for ( s = 0; s < numTotalDof; ++s )
        {
            rhs_k[s] = f[s + fInd ] + 
                       u[s + uInd];
        }

        {
            using namespace ExpressionAssembly;
            //omp_setlock( &lockMS );
            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( k ) );
            //omp_unsetlock( &lockMS );
        }

    }
    } //section
    
    //#pragma omp section
    {
    //Cycling on r-blocks
    //#pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mr; ++k )
    {
        //omp_setlock( &lockMS );
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        UInt fInd( f_interpolated->block( mx+k )->firstIndex() );
        UInt uInd( uold_interpolated->block( mx+k )->firstIndex() );
        vector_type f( *f_interpolated );
        vector_type u( *uold_interpolated );
        //omp_unsetlock( &lockMS );

        //#pragma omp parallel for schedule( static ) default(shared) private( s )
        for ( s = 0; s < numTotalDof; ++s )
        {
            rhs_k[s] = f[s + fInd ] +
                       u[s + uInd ];
        }


        {
            using namespace ExpressionAssembly;

            //omp_setlock( &lockMS );
            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( mx + k ) );
            //omp_unsetlock( &lockMS );
        }

    }
    } //section
    
    //#pragma omp section
    {
    //Cycling on theta-blocks
    //#pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mtheta; ++k )
    {
        //omp_setlock( &lockMS );
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        UInt fInd( f_interpolated->block( mx+mr+k )->firstIndex() );
        UInt uInd( uold_interpolated->block( mx+mr+k )->firstIndex() );
        vector_type f( *f_interpolated );
        vector_type u( *uold_interpolated );
        //omp_unsetlock( &lockMS );

        //#pragma omp parallel for schedule( static ) default(shared) private( s )
        for ( s = 0; s < numTotalDof; ++s )
        {
            rhs_k[s] = f[s + fInd ] +
                       u[s + uInd ];
        }


        {
            using namespace ExpressionAssembly;

            //omp_setlock( &lockMS );
            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( mx + mr + k ) );
            //omp_unsetlock( &lockMS );
        }

    }
    } //section
    
    //#pragma omp section
    {
    //Cycling on p-blocks (null force)
    //#pragma omp parallel for schedule( static ) default(shared) private( k )
    for ( k = 0; k < mp; ++k )
    {
        //omp_setlock( &lockMS );
        VectorEpetra rhs_k( M_pressureFespace->map(), Repeated );
        //omp_unsetlock( &lockMS );

        //#pragma omp parallel for schedule( static ) default(shared) private( s )
        for ( s = 0; s < pnumTotalDof; ++s )
        {
            rhs_k[s] = 0;
        }

       {
            using namespace ExpressionAssembly;

            //omp_setlock( &lockMS );
            integrate  ( elements( M_etpfespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etpfespace,
                         value( M_etpfespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( mx + mr + mtheta + k ) );
            //omp_unsetlock( &lockMS );
        }
    }
    } //section
    } //sections
   //omp_destroy_lock( &lockMS );

    return;
}
#endif
