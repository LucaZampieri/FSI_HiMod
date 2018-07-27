#ifndef __HMAADDNSASSEMBLER_HPP__
#define __HMAADDNSASSEMBLER_HPP__

//AddNSAssembler
/*
    Assemble the matrix of the problem, for the moment with costant coefficients (mu,beta,gamma)
*/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addStokesProblem( const matrix_ptrType& systemMatrix, const Real& nu, ReferenceMap& refMap, const Real& t, const Real& alpha )
{
    // xx-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice

            VectorSmall<6> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xx( k, j, nu, alpha );     
            Coeff[1] = M_modalbasis->compute_r001xx( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010xx( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r100xx( k, j, nu );  
            Coeff[4] = M_modalbasis->compute_r101xx( k, j, nu );
            Coeff[5] = M_modalbasis->compute_r110xx( k, j, nu );
            
            //Assemble the (j,k)1D FEMproblem
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + ( Coeff[3] * M_modalbasis->Rho() + Coeff[4] * M_modalbasis->dRho() ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[5] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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
            VectorSmall<6> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rr( k, j, nu, alpha );     
            Coeff[1] = M_modalbasis->compute_r001rr( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010rr( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r100rr( k, j, nu );  
            Coeff[4] = M_modalbasis->compute_r101rr( k, j, nu );
            Coeff[5] = M_modalbasis->compute_r110rr( k, j, nu );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + ( Coeff[3] * M_modalbasis->Rho() + Coeff[4] * M_modalbasis->dRho() ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[5] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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

            VectorSmall<6> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tt( k, j, nu, alpha );     
            Coeff[1] = M_modalbasis->compute_r001tt( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010tt( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r100tt( k, j, nu );  
            Coeff[4] = M_modalbasis->compute_r101tt( k, j, nu );
            Coeff[5] = M_modalbasis->compute_r110tt( k, j, nu );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + ( Coeff[3] * M_modalbasis->Rho() + Coeff[4] * M_modalbasis->dRho() ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[5] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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

            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xr( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r001xr( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010xr( k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) ) * phi_j 
                            
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

            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xt( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r001xt( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010xt( k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) * phi_j )
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rx( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r100rx( k, j, nu ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rt( k, j, nu );     

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tx( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r100tx( k, j, nu ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tr( k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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
            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000px( k, j );
            Coeff[1] = M_modalbasis->compute_r001px( k, j );
            Coeff[2] = M_modalbasis->compute_r010px( k, j );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) ) * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000pr( k, j );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000pt( k, j );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xp( k, j );     
            Coeff[1] = M_modalbasis->compute_r100xp( k, j ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rp( k, j );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tp( k, j );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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
pressureMassMatrix( const matrix_ptrType& massMatrix, ReferenceMap& refMap )
{
    // pp-block
    for ( UInt j = 0; j != M_modalbasis->mp(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
        {
            //j refers to the test function, k refers to the solution

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> coeff;
                coeff[0] = M_modalbasis->compute_r000pp( k, j );

                integrate( elements( M_etpfespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etpfespace,
                           M_etpfespace,
                           
                           coeff[0] * phi_i * phi_j
                            
                          )                         
	                  >> ( massMatrix->block( j, k ) );
            } // end namespace
        } // end k-loop
    } // end j-loop
   
    massMatrix->globalAssemble();
 
    return;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addNavierStokesProblem( const matrix_ptrType& systemMatrix, const Real& nu, ReferenceMap& refMap, const Real& t, const Real& alpha,
						const vector_type& advection )
{
    UInt Nelements( M_velocityFespace->mesh()->numElements() );
    UInt udof( 2 * Nelements + 1 );
    
    MapEpetra rowMap( M_velocityFespace->map() );
    matrix_ptrType xMatrixAdvection( new matrix_type( rowMap ) );
    matrix_ptrType rMatrixAdvection( new matrix_type( rowMap ) );
    matrix_ptrType tMatrixAdvection( new matrix_type( rowMap ) );
    
    *xMatrixAdvection *= 0.0;
    *rMatrixAdvection *= 0.0;
    *tMatrixAdvection *= 0.0;
    
    for( UInt i( 0 ); i != udof; ++i )
    {
        for( UInt k( 0 ); k != M_modalbasis->mx(); ++k )
            xMatrixAdvection->addToCoefficient( i, k, advection[i+k*udof] );
            
        for( UInt k( 0 ); k != M_modalbasis->mr(); ++k )
            rMatrixAdvection->addToCoefficient( i, k, advection[i+(k+M_modalbasis->mx())*udof] );

        for( UInt k( 0 ); k != M_modalbasis->mtheta(); ++k )
            tMatrixAdvection->addToCoefficient( i, k, advection[i+(k+M_modalbasis->mx()+M_modalbasis->mr())*udof] );
    }

    xMatrixAdvection->globalAssemble();
    rMatrixAdvection->globalAssemble();
    tMatrixAdvection->globalAssemble();
    xMatrixAdvection->spy("xAdv");
    rMatrixAdvection->spy("rAdv");
    tMatrixAdvection->spy("tAdv");

    // xx-block
    for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
    {
        for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
        {
            //j refers to the test function, k refers to the solution
            //For each blocks compute the integral coefficient on the slice
            VectorSmall<4> Coeff;

			boost::shared_ptr<VectorEpetra> r100xx_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			boost::shared_ptr<VectorEpetra> r000xx_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			
			interpolateR100xx( r100xx_interpolated, xMatrixAdvection, k, j, nu );
			interpolateR000xx( r000xx_interpolated, xMatrixAdvection, rMatrixAdvection, tMatrixAdvection, k, j, nu, alpha );
     
            Coeff[0] = M_modalbasis->compute_r001xx( k, j, nu ); 
            Coeff[1] = M_modalbasis->compute_r010xx( k, j, nu );
            Coeff[2] = M_modalbasis->compute_r101xx( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r110xx( k, j, nu );

            //Assemble the (j,k)1D FEMproblem
            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *r000xx_interpolated ) * M_modalbasis->Rho() * phi_i * phi_j
                            + Coeff[0] * M_modalbasis->dRho() * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *r100xx_interpolated ) * M_modalbasis->Rho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[2] * M_modalbasis->dRho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[3] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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
            VectorSmall<4> Coeff;
            
            boost::shared_ptr<VectorEpetra> r100rr_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			boost::shared_ptr<VectorEpetra> r000rr_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			
			interpolateR100rr( r100rr_interpolated, xMatrixAdvection, k, j, nu );
			interpolateR000rr( r000rr_interpolated, xMatrixAdvection, rMatrixAdvection, tMatrixAdvection, k, j, nu, alpha );

            Coeff[0] = M_modalbasis->compute_r001rr( k, j, nu ); 
            Coeff[1] = M_modalbasis->compute_r010rr( k, j, nu );
            Coeff[2] = M_modalbasis->compute_r101rr( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r110rr( k, j, nu );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *r000rr_interpolated ) * M_modalbasis->Rho() * phi_i * phi_j
                            + Coeff[0] * M_modalbasis->dRho() * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *r100rr_interpolated ) * M_modalbasis->Rho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[2] * M_modalbasis->dRho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[3] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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

            VectorSmall<4> Coeff;
            
            boost::shared_ptr<VectorEpetra> r100tt_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			boost::shared_ptr<VectorEpetra> r000tt_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );
			
			interpolateR100tt( r100tt_interpolated, xMatrixAdvection, k, j, nu );
			interpolateR000tt( r000tt_interpolated, xMatrixAdvection, rMatrixAdvection, tMatrixAdvection, k, j, nu, alpha );

            Coeff[0] = M_modalbasis->compute_r001tt( k, j, nu ); 
            Coeff[1] = M_modalbasis->compute_r010tt( k, j, nu );
            Coeff[2] = M_modalbasis->compute_r101tt( k, j, nu );
            Coeff[3] = M_modalbasis->compute_r110tt( k, j, nu );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            value( M_etufespace, *r000tt_interpolated ) * M_modalbasis->Rho() * phi_i * phi_j
                            + Coeff[0] * M_modalbasis->dRho() * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_j * dot( grad( phi_i ), value( oneVector ) )
                            + value( M_etufespace, *r100tt_interpolated ) * M_modalbasis->Rho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[2] * M_modalbasis->dRho() * dot( grad( phi_j ), value( oneVector ) ) * phi_i
                            + Coeff[3] * M_modalbasis->Rho() * dot( grad( phi_i ), grad( phi_j ) ) 
                            
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

            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xr( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r001xr( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010xr( k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) ) * phi_j 
                            
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

            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xt( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r001xt( k, j, nu ); 
            Coeff[2] = M_modalbasis->compute_r010xt( k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) * phi_j )
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rx( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r100rx( k, j, nu ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

			boost::shared_ptr<VectorEpetra> r000rt_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );            

            interpolateR000rt( r000rt_interpolated, tMatrixAdvection, k, j, nu );  

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( value( M_etufespace, *r000rt_interpolated ) * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tx( k, j, nu );     
            Coeff[1] = M_modalbasis->compute_r100tx( k, j, nu ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

            boost::shared_ptr<VectorEpetra> r000tr_interpolated( new VectorEpetra( M_velocityFespace->map(), Repeated ) );            

            interpolateR000tr( r000tr_interpolated, tMatrixAdvection, k, j, nu );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etufespace,
                           M_etufespace,
                           
                            ( value( M_etufespace, *r000tr_interpolated ) * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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
            VectorSmall<3> Coeff;

            Coeff[0] = M_modalbasis->compute_r000px( k, j );
            Coeff[1] = M_modalbasis->compute_r001px( k, j );
            Coeff[2] = M_modalbasis->compute_r010px( k, j );

            {
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() + Coeff[1] * M_modalbasis->dRho() ) * phi_i * phi_j
                            + Coeff[2] * M_modalbasis->Rho() * dot( grad( phi_i ), value( oneVector ) ) * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000pr( k, j );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000pt( k, j );
            
			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_pressureFespace->qr(),
                           M_etufespace,
                           M_etpfespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<2> Coeff;

            Coeff[0] = M_modalbasis->compute_r000xp( k, j );     
            Coeff[1] = M_modalbasis->compute_r100xp( k, j ); 

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            + Coeff[1] * M_modalbasis->Rho() * phi_i * dot( grad( phi_j ), value( oneVector ) )
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000rp( k, j );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

            VectorSmall<1> Coeff;

            Coeff[0] = M_modalbasis->compute_r000tp( k, j );

			{
                using namespace ExpressionAssembly;

                VectorSmall<1> oneVector;
                oneVector[0] = 1.0;

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
                           M_etpfespace,
                           M_etufespace,
                           
                            ( Coeff[0] * M_modalbasis->Rho() ) * phi_i * phi_j
                            
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

                integrate( elements( M_etufespace->mesh() ),
                           M_velocityFespace->qr(),
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
addrhs( const vector_ptrType& rhs, const Real& fx, const Real& fr, const Real& ftheta )
{
    //Cycling on x-block
    for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        Real Coeff = M_modalbasis->compute_Phix( k );

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        fx * Coeff * phi_i
                      )
                    >> ( rhs->block( k ) );
        }

    }
    
    //Cycling on r-block
    for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        Real Coeff = M_modalbasis->compute_Phir( k );

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        fr * Coeff * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + k ) );
        }

    }
    
    //Cycling on theta-block
    for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        Real Coeff = M_modalbasis->compute_Phitheta( k );

        {
            using namespace ExpressionAssembly;
            
            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        ftheta * Coeff * phi_i
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

//
//1) Verify efficency of assignment on the little rhs_k
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addrhs( const vector_ptrType& rhs, const vector_ptrType& f_interpolated, const vector_ptrType& uold_interpolated )
{

    //Cycling on x-blocks
    for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        for ( UInt s( 0 ); s != M_velocityFespace->dofPtr()->numTotalDof(); ++s )
        {
            rhs_k[s] = ( *f_interpolated ) [s + f_interpolated->block( k )->firstIndex()] + 
						( *uold_interpolated ) [s + uold_interpolated->block( k )->firstIndex()];
        }

        {
            using namespace ExpressionAssembly;

            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( k ) );
        }

    }
    
    //Cycling on r-blocks
    for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        for ( UInt s( 0 ); s != M_velocityFespace->dofPtr()->numTotalDof(); ++s )
        {
            rhs_k[s] = ( *f_interpolated ) [s + f_interpolated->block( M_modalbasis->mx() + k )->firstIndex()] +
                       ( *uold_interpolated ) [s + uold_interpolated->block( M_modalbasis->mx() + k )->firstIndex()];
        }


        {
            using namespace ExpressionAssembly;

            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( M_modalbasis->mx() + k ) );
        }

    }
    
    //Cycling on theta-blocks
    for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        VectorEpetra rhs_k( M_velocityFespace->map(), Repeated );
        for ( UInt s( 0 ); s != M_velocityFespace->dofPtr()->numTotalDof(); ++s )
        {
            rhs_k[s] = ( *f_interpolated )[s + f_interpolated->block( M_modalbasis->mx()+M_modalbasis->mr()+k )->firstIndex()] +
                       ( *uold_interpolated ) [s + uold_interpolated->block( M_modalbasis->mx()+M_modalbasis->mr()+k )->firstIndex()];
        }


        {
            using namespace ExpressionAssembly;

            integrate  ( elements( M_etufespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etufespace,
                         value( M_etufespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + k ) );
        }

    }
    
    //Cycling on p-blocks (null force)
    for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
    {
	VectorEpetra rhs_k( M_pressureFespace->map(), Repeated );
        for ( UInt s( 0 ); s != M_pressureFespace->dofPtr()->numTotalDof(); ++s )
        {
            rhs_k[s] = 0;
        }

       {
            using namespace ExpressionAssembly;

            integrate  ( elements( M_etpfespace->mesh() ),
                         M_velocityFespace->qr(),
                         M_etpfespace,
                         value( M_etpfespace, rhs_k ) * phi_i
                       )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
        }

    }

    return;
}

#endif
