#ifndef __HMAADDNSASSEMBLERFSI_HPP__
#define __HMAADDNSASSEMBLERFSI_HPP__

/*
    Assemble the matrix of the problem, for the moment with x- dependent radius
*/
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addStokesProblemFSI( const matrix_ptrType& systemMatrix, const Real& mu, const Real& rho_f, const Real& dt, ReferenceMap& refMap, const Real& t )
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

            M_modalbasis->computeFSI_r000xx( k, j, mu, rho_f, dt, *R000xx );
            M_modalbasis->computeFSI_r001xx( k, j, mu, *R001xx );
            M_modalbasis->computeFSI_r010xx( k, j, mu, *R010xx );
            M_modalbasis->computeFSI_r100xx( k, j, mu, *R100xx );
            M_modalbasis->computeFSI_r101xx( k, j, mu, *R101xx );
            M_modalbasis->computeFSI_r110xx( k, j, mu, *R110xx );

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

            M_modalbasis->computeFSI_r000rr( k, j, mu, rho_f, dt, *R000rr );
            M_modalbasis->computeFSI_r001rr( k, j, mu, *R001rr );
            M_modalbasis->computeFSI_r010rr( k, j, mu, *R010rr );
            M_modalbasis->computeFSI_r100rr( k, j, mu, *R100rr );
            M_modalbasis->computeFSI_r101rr( k, j, mu, *R101rr );
            M_modalbasis->computeFSI_r110rr( k, j, mu, *R110rr );

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

            M_modalbasis->computeFSI_r000tt( k, j, mu, rho_f, dt, *R000tt );
            M_modalbasis->computeFSI_r001tt( k, j, mu, *R001tt );
            M_modalbasis->computeFSI_r010tt( k, j, mu, *R010tt );
            M_modalbasis->computeFSI_r100tt( k, j, mu, *R100tt );
            M_modalbasis->computeFSI_r101tt( k, j, mu, *R101tt );
            M_modalbasis->computeFSI_r110tt( k, j, mu, *R110tt );

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


            M_modalbasis->computeFSI_r000xr( k, j, mu, *R000xr );
            M_modalbasis->computeFSI_r001xr( k, j, mu, *R001xr );
            M_modalbasis->computeFSI_r010xr( k, j, mu, *R010xr );

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

            M_modalbasis->computeFSI_r000xt( k, j, mu, *R000xt );
            M_modalbasis->computeFSI_r001xt( k, j, mu, *R001xt );
            M_modalbasis->computeFSI_r010xt( k, j, mu, *R010xt );

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

            M_modalbasis->computeFSI_r000rx( k, j, mu, *R000rx );
            M_modalbasis->computeFSI_r100rx( k, j, mu, *R100rx );

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

            M_modalbasis->computeFSI_r000rt( k, j, mu, *R000rt );

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

            M_modalbasis->computeFSI_r000tx( k, j, mu, *R000tx );
            M_modalbasis->computeFSI_r100tx( k, j, mu, *R100tx );

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

            M_modalbasis->computeFSI_r000tr( k, j, mu, *R000tr );

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

            M_modalbasis->computeFSI_r000px( k, j, *R000px );
            M_modalbasis->computeFSI_r001px( k, j, *R001px );
            M_modalbasis->computeFSI_r010px( k, j, *R010px );

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

            M_modalbasis->computeFSI_r000pr( k, j, *R000pr );

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

            M_modalbasis->computeFSI_r000pt( k, j, *R000pt );

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

            M_modalbasis->computeFSI_r000xp( k, j, *R000xp );
            M_modalbasis->computeFSI_r100xp( k, j, *R100xp );

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

            M_modalbasis->computeFSI_r000rp( k, j, *R000rp );

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

            M_modalbasis->computeFSI_r000tp( k, j, *R000tp );

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

//------------------- Parte relativa alla velocità ALE -------------------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addALEProblemFSI( const matrix_ptrType& systemMatrix, const Real& rho_f, const vector_type& wr, ReferenceMap& refMap, const Real& t )
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
/*
		// xx-block
		for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
				{
						//j refers to the test function, k refers to the solution
						//For each blocks compute the integral coefficient on the slice

						vector_ptrType W000xx( new vector_type( M_velocityFespace->map(), Repeated ) );
						vector_ptrType W100xx( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000xx( k, j, wx, rho_f, *W000xx );
						M_modalbasis->computeFSI_w100xx( k, j, wx, rho_f, *W100xx );

						//Assemble the (j,k)1D FEMproblem
						{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000xx ) * phi_i * phi_j +
														+ value( M_etufespace, *W100xx ) * dot( grad( phi_j ), value( oneVector ) ) * phi_i

													)
											 >> ( systemMatrix->block( j, k ) );
						}
				}
		}
*/
	// rr-block
		for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
				{
						//j refers to the test function, k refers to the solution

						//For each blocks compute the integral coefficient on the slice
						vector_ptrType W000rr( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000rr( k, j, wr, rho_f, *W000rr );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000rr ) * phi_i * phi_j

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

						vector_ptrType W000tt( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000tt( k, j, wr/*, wtheta*/, rho_f, *W000tt );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000tt ) * phi_i * phi_j

													)
											 >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
						}
				}
		}
/*
		// xr-block
		for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
				{
						//j refers to the test function, k refers to the solution

						//For each blocks compute the integral coefficient on the slice

						vector_ptrType W000xr( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000xr( k, j, wx, rho_f, *W000xr );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000xr ) * phi_i * phi_j

													)
											 >> ( systemMatrix->block( M_modalbasis->mx() + j, k ) );
						}
				}
		}
*/
/*
		// xt-block
		for ( UInt j = 0; j != M_modalbasis->mtheta(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
				{
						//j refers to the test function, k refers to the solution

						//For each blocks compute the integral coefficient on the slice

						vector_ptrType W000xt( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000xt( k, j, wx, rho_f, *W000xt );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000xt ) * phi_i * phi_j

													)
											 >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, k ) );
						}
				}
		}
*/
		// rx-block
		for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
				{
						//j refers to the test function, k refers to the solution

						vector_ptrType W000rx( new vector_type( M_velocityFespace->map(), Repeated ) );
						vector_ptrType W100rx( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000rx( k, j, wr, rho_f, *W000rx );
						M_modalbasis->computeFSI_w100rx( k, j, wr, rho_f, *W100rx );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000rx ) * phi_i * phi_j
														+ value( M_etufespace, *W100rx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )

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

						vector_ptrType W000rt( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000rt( k, j, wr/*, wtheta*/, rho_f, *W000rt );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

													 value( M_etufespace, *W000rt ) * phi_i * phi_j

													)
											 >> ( systemMatrix->block( M_modalbasis->mx() + M_modalbasis->mr() + j, M_modalbasis->mx() + k ) );
						}
				}
		}
/*
		// tx-block
		for ( UInt j = 0; j != M_modalbasis->mx(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
				{
						//j refers to the test function, k refers to the solution

						vector_ptrType W000tx( new vector_type( M_velocityFespace->map(), Repeated ) );
						vector_ptrType W100tx( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000tx( k, j, wtheta, rho_f, *W000tx );
						M_modalbasis->computeFSI_w100tx( k, j, wtheta, rho_f, *W100tx );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

														value( M_etufespace, *W000tx ) * phi_i * phi_j
														+ value( M_etufespace, *W100tx ) * phi_i * dot( grad( phi_j ), value( oneVector ) )

													)
											 >> ( systemMatrix->block( j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
						}
				}
		}
*/
/*
		// tr-block
		for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
		{
				for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
				{
						//j refers to the test function, k refers to the solution

						vector_ptrType W000tr( new vector_type( M_velocityFespace->map(), Repeated ) );

						M_modalbasis->computeFSI_w000tr( k, j, wtheta, rho_f, *W000tr );

			{
								using namespace ExpressionAssembly;

								VectorSmall<1> oneVector;
								oneVector[0] = 1.0;

								integrate( elements( M_etufespace->mesh() ),
													 M_velocityFespace->qr(),
													 M_etufespace,
													 M_etufespace,

													 value( M_etufespace, *W000tr ) * phi_i * phi_j

													)
											 >> ( systemMatrix->block( M_modalbasis->mx() + j, M_modalbasis->mx() + M_modalbasis->mr() + k ) );
						}
				}
		}
*/
		return;
}

//------------------------------------------------------------------------------

//---------------- Parte relativa all'integrale di bordo laterale --------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addWallProblemFSI( const matrix_ptrType& systemMatrix, const Real& rho_s, const Real& h_s, const Real& dt, const Real& E, const Real& csi, const Real& R0, ReferenceMap& refMap, const Real& t )
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

		// rr-block
	    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
	    {
	        for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
	        {
	            //j refers to the test function, k refers to the solution

	            //For each blocks compute the integral coefficient on the slice
	            vector_ptrType B000rr( new vector_type( M_velocityFespace->map(), Repeated ) );

	            M_modalbasis->computeFSI_b000rr( k, j, rho_s, h_s, dt, E, csi, R0, *B000rr );

				{
	                using namespace ExpressionAssembly;

	                VectorSmall<1> oneVector;
	                oneVector[0] = 1.0;

	                integrate( elements( M_etufespace->mesh() ),
	                           M_velocityFespace->qr(),
	                           M_etufespace,
	                           M_etufespace,

	                            value( M_etufespace, *B000rr ) * phi_i * phi_j

	                          )
		                     >> ( systemMatrix->block( M_modalbasis->mx() + j, M_modalbasis->mx() + k ) );
	            }
	        }
	    }

			return;
}

//------------------------------------------------------------------------------

//---------------- Parte relativa all'integrale di volume del rhs --------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addrhsFSI( const vector_ptrType& rhs, const vector_type& f, const vector_type& u_old, const Real& rho_f, const Real& dt )
{
    //Cycling on x-block
    for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        vector_ptrType F00x( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->computeFSI_f00x( k, f, u_old, rho_f, dt, *F00x );

        {
            using namespace ExpressionAssembly;

            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *F00x ) * phi_i
                      )
                    >> ( rhs->block( k ) );
        }

    }

    //Cycling on r-block
    for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
    {
        vector_ptrType F00r( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->computeFSI_f00r( k, f, u_old, rho_f, dt, *F00r );

        {
            using namespace ExpressionAssembly;

            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *F00r ) * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + k ) );
        }

    }

    //Cycling on theta-block
    for ( UInt k = 0; k != M_modalbasis->mtheta(); ++k )
    {
        vector_ptrType F00t( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->computeFSI_f00t( k, f, u_old, rho_f, dt, *F00t );

        {
            using namespace ExpressionAssembly;

            integrate ( elements( M_etufespace->mesh() ),
                        M_velocityFespace->qr(),
                        M_etufespace,
                        value( M_etufespace, *F00t ) * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + k ) );
        }

    }

    //Cycling on p-block : null force
    for ( UInt k = 0; k != M_modalbasis->mp(); ++k )
    {
        {
            using namespace ExpressionAssembly;

            integrate ( elements( M_etpfespace->mesh() ),
                        M_pressureFespace->qr(),
                        M_etpfespace,
                        0 * phi_i
                      )
                    >> ( rhs->block( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() + k ) );
        }

    }

    return;
}

//------------------------------------------------------------------------------

//---------------- Parte relativa all'integrale di bordo del rhs ---------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addWallrhsFSI( const vector_ptrType& rhs, const vector_type& urWall_old, const vector_type& etar_old, const Real& rho_s, const Real& h_s, const Real& dt, const Real& E, const Real& csi, const Real& R0 )
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

	//Cycling on r-block
	for ( UInt k = 0; k != M_modalbasis->mr(); ++k )
	{
			vector_ptrType B00r( new vector_type( M_velocityFespace->map(), Repeated ) );
			M_modalbasis->computeFSI_b00r( k, urWall_old, etar_old, rho_s, h_s, dt, E, csi, R0, *B00r );

			{
					using namespace ExpressionAssembly;

					integrate ( elements( M_etufespace->mesh() ),
											M_velocityFespace->qr(),
											M_etufespace,
											value( M_etufespace, *B00r ) * phi_i
										)
									>> ( rhs->block( M_modalbasis->mx() + k ) );
			}

	}

	return;
}

//------------------------------------------------------------------------------

//------------ Parte relativa agli integrali in inflow e outflow ---------------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addBCxInOutFSI( const vector_ptrType& rhs, const Real& p1, const Real& p2 )
{
		UInt dof = M_etufespace->dof().numTotalDof();

		for( UInt k = 0; k != M_modalbasis->mx(); ++k )
		{
				Real P00x1;
				M_modalbasis->computeFSI_p00x1( k, p1, P00x1 );
				Real P00x2;
				M_modalbasis->computeFSI_p00x2( k, p2, P00x2 );

				rhs->setCoefficient( k * dof,
														 ( *rhs )( k * dof ) + P00x1 );
				rhs->setCoefficient( ( k + 1 ) * dof - 1,
										         ( *rhs )( ( k + 1 ) * dof - 1 ) + P00x2 );
		}
}

//------------------------------------------------------------------------------

//------ Parte relativa alle condizioni al bordo di Dirichlet omogenee ---------

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addBCrInOutFSI( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs )
{
    UInt udof = M_etufespace->dof().numTotalDof();
		UInt pdof = M_etpfespace->dof().numTotalDof();

    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + j * udof, M_modalbasis->mx() * udof + j * udof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + j * udof, 0 );

				systemMatrix->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1, M_modalbasis->mx() * udof + j * udof + pdof - 1, 1e+30 );
				rhs->setCoefficient( M_modalbasis->mx() * udof + j * udof + pdof - 1, 0 );
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addBCthetaInOutFSI( const matrix_ptrType& systemMatrix, const vector_ptrType& rhs )
{
    UInt udof = M_etufespace->dof().numTotalDof();
		UInt pdof = M_etpfespace->dof().numTotalDof();

    for ( UInt j = 0; j != M_modalbasis->mr(); ++j )
    {
        systemMatrix->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof, M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof, 1e+30 );
        rhs->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof, 0 );

				systemMatrix->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1, M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1, 1e+30 );
				rhs->setCoefficient( M_modalbasis->mx() * udof + M_modalbasis->mr() * udof + j * udof + pdof - 1, 0 );
    }
}

//------------------------------------------------------------------------------

/*
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
addrhs( const vector_ptrType& rhs, const Real& fx, const Real& fr, const Real& ftheta, const bool& b )
{
    //Cycling on x-block
    for ( UInt k = 0; k != M_modalbasis->mx(); ++k )
    {
        vector_ptrType Phix( new vector_type( M_velocityFespace->map(), Repeated ) );
        M_modalbasis->computeFSI_Phix( k, *Phix );

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
        M_modalbasis->computeFSI_Phir( k, *Phir );

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
        M_modalbasis->computeFSI_Phitheta( k, *Phitheta );

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
*/

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
evaluateForce3DGridFSI( const function_Type& fx, const function_Type& fr, const function_Type& ftheta, const Real& t, const grid_type& grid )
{
    UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

    MapEpetra Map( nquadRho * nquadTheta * ndofuFE * 3, Comm );
    vector_type fcoeff( Map, Unique );
    fcoeff *= 0.0;

    for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
    {
      for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
      {
        for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
        {
          fcoeff[nrho + ntheta*nquadRho + s*nquadTheta*nquadRho] = fx( t, std::get<0>(grid[s][nrho][ntheta]),
																																			    std::get<1>(grid[s][nrho][ntheta]),
																																			    std::get<2>(grid[s][nrho][ntheta]), 1 );
				 	fcoeff[ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho] = fr( t, std::get<0>(grid[s][nrho][ntheta]),
																																						 												    std::get<1>(grid[s][nrho][ntheta]),
																																						 												    std::get<2>(grid[s][nrho][ntheta]), 1 );
					fcoeff[2*ndofuFE*nquadRho*nquadTheta + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho] = ftheta( t, std::get<0>(grid[s][nrho][ntheta]),
																																							 														    std::get<1>(grid[s][nrho][ntheta]),
																																							 													 	    std::get<2>(grid[s][nrho][ntheta]), 1 );


        }
      }
    }
		return fcoeff;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
evaluateBase3DGridFSI( const vector_type& fun )
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
evaluateBaseWallGridFSI( const vector_type& fun )
{
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();

    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

    MapEpetra Map( nquadTheta * ndofuFE, Comm );
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
        for ( UInt j( 0 ); j != M_modalbasis->mr(); ++j ) // Ciclying over r-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                  fcoeff[ntheta + s * nquadTheta] +=
                  fun[ndofuFE * M_modalbasis->mx() + j * ndofuFE + index] *
                                  M_modalbasis->rphirhoWall( j ) * normrho *
                                  M_modalbasis->rphitheta( j, ntheta ) * normtheta;
            }
        }
    }

    return fcoeff;
}

template< typename mesh_type, typename matrix_type, typename vector_type>
vector_type NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
evaluateInitialVelocityWallFSI( const function_Type& ur0, const grid_type& grid, const Real& R )
{
	UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();

	DOF DataVelFESpace( M_velocityFespace->dof() );
	UInt ndofuFE = DataVelFESpace.numTotalDof();

	boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

	MapEpetra Map( nquadTheta * ndofuFE, Comm );
	vector_type fcoeff( Map, Unique );
	fcoeff *= 0.0;

	for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
	{
		for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
		{
			fcoeff[ntheta + s*nquadTheta] = ur0( 0, std::get<0>(grid[s][0][ntheta]), R, std::get<2>(grid[s][0][ntheta]), 1 );
		}
	}
	return fcoeff;
}

#endif
