#ifndef __HMAINTERPOLATIONPIPE_HPP__
#define __HMAINTERPOLATIONPIPE_HPP__

//Note:
//1) Added get methods in FESpace.hpp to get res FE pointer
//2) There's not a function to add a coefficient to a VectorEpetraStructured, it should be into VectorStructuredView.hpp like
//   in MatrixStructuredView.hpp
//3) Problem with vect.setCoefficient if you use parallel mode?
//4) MeshPartitioner.hpp contiene dei float!!!!
// ----------------- R( x ) -------------------------------
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
interpolate( const function_Type& fx, const function_Type& fr, const function_Type& ftheta,
            const vector_ptrType& f_interpolated, const Real& t )
{

    // First, we build a "quadrature" that consists in the nodes (0 weight)
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape( shapeDimension( M_velocityFespace->refFEPtr()->shape() ), M_velocityFespace->refFEPtr()->shape() );
    interpQuad.setPoints( M_velocityFespace->refFEPtr()->refCoor(), std::vector<Real> ( M_velocityFespace->refFEPtr()->nbDof(), 0 ) );

    // Then, we define a currentFE with nodes on the reference nodes
    CurrentFE interpCFE( M_velocityFespace->refFE(), getGeometricMap( *M_velocityFespace->mesh() ), interpQuad );

    // Some constants
    UInt totalNumberElements( M_velocityFespace->mesh()->numElements() );
    UInt numberLocalDof( M_velocityFespace->dofPtr()->numLocalDof() );//Local to the element
    UInt numberTotalDof( M_velocityFespace->dofPtr()->numTotalDof() );//Total
    UInt numberTotalPDof( M_pressureFespace->dofPtr()->numTotalDof() );//Total

    // Storage for the values
    std::vector<Real> nodalValues( numberLocalDof, 0 );
    std::vector<Real> FEValues( numberLocalDof, 0 );

    //Do the loop over the x-frequencies
    for( UInt k( 0 ); k != M_modalbasis->mx() ; ++k )
    {
        // Do the loop over the cells
        for( UInt iterElement( 0 ); iterElement != totalNumberElements; ++iterElement )
        {
            // We update the CurrentFE so that we get the coordinates of the nodes
            interpCFE.update( M_velocityFespace->mesh()->element( iterElement ), UPDATE_QUAD_NODES );

            // Loop over the degrees of freedom (= quadrature nodes)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Store the nodal value of the fourier coefficient
                nodalValues[iterDof] =  M_modalbasis->xFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), fx, k );
            }

            // Transform the nodal values in FE values
            FEValues = M_velocityFespace->refFEPtr()->nodalToFEValues( nodalValues ); //In P1 case it's useless

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Find the ID of the considered DOF
                ID globalDofID( M_velocityFespace->dofPtr()->localToGlobalMap( iterElement, iterDof ) );
                // Compute the value of the function and set it
                //#pragma GCC diagnostic ignored "-Wconversion"
                f_interpolated->setCoefficient( globalDofID + k * numberTotalDof, FEValues[iterDof] );
                //#pragma GCC diagnostic warning "-Wconversion"
            }
        }
    }
    
    //Do the loop over the r-frequencies
    for( UInt k( 0 ); k != M_modalbasis->mr() ; ++k )
    {
        // Do the loop over the cells
        for( UInt iterElement( 0 ); iterElement != totalNumberElements; ++iterElement )
        {
            // We update the CurrentFE so that we get the coordinates of the nodes
            interpCFE.update( M_velocityFespace->mesh()->element( iterElement ), UPDATE_QUAD_NODES );

            // Loop over the degrees of freedom (= quadrature nodes)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Store the nodal value of the fourier coefficient
                nodalValues[iterDof] =  M_modalbasis->rFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), fr, k );
            }

            // Transform the nodal values in FE values
            FEValues = M_velocityFespace->refFEPtr()->nodalToFEValues( nodalValues ); //In P1 case it's uselesss

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Find the ID of the considered DOF
                ID globalDofID( M_velocityFespace->dofPtr()->localToGlobalMap( iterElement, iterDof ) );
                // Compute the value of the function and set it
                //#pragma GCC diagnostic ignored "-Wconversion"
                f_interpolated->setCoefficient( globalDofID + ( M_modalbasis->mx() + k ) * numberTotalDof, FEValues[iterDof] );
                //#pragma GCC diagnostic warning "-Wconversion"
            }
        }
    }
    
    //Do the loop over the theta-frequencies
    for( UInt k( 0 ); k != M_modalbasis->mtheta() ; ++k )
    {
        // Do the loop over the cells
        for( UInt iterElement( 0 ); iterElement != totalNumberElements; ++iterElement )
        {
            // We update the CurrentFE so that we get the coordinates of the nodes
            interpCFE.update( M_velocityFespace->mesh()->element( iterElement ), UPDATE_QUAD_NODES );

            // Loop over the degrees of freedom (= quadrature nodes)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Store the nodal value of the fourier coefficient
                nodalValues[iterDof] =  M_modalbasis->thetaFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), ftheta, k );
            }

            // Transform the nodal values in FE values
            FEValues = M_velocityFespace->refFEPtr()->nodalToFEValues( nodalValues ); //In P1 case it's uselesss

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Find the ID of the considered DOF
                ID globalDofID( M_velocityFespace->dofPtr()->localToGlobalMap( iterElement, iterDof ) );
                // Compute the value of the function and set it
                //#pragma GCC diagnostic ignored "-Wconversion"
                f_interpolated->setCoefficient( globalDofID + ( M_modalbasis->mx() + M_modalbasis->mr() + k ) * numberTotalDof, FEValues[iterDof] );
                //#pragma GCC diagnostic warning "-Wconversion"
            }
        }
    }
    
    //Do the loop over the p-frequencies
    for( UInt k( 0 ); k != M_modalbasis->mp() ; ++k )
    {
        // Do the loop over the cells
        for( UInt iterElement( 0 ); iterElement != totalNumberElements; ++iterElement )
        {
            // We update the CurrentFE so that we get the coordinates of the nodes
            interpCFE.update( M_pressureFespace->mesh()->element( iterElement ), UPDATE_QUAD_NODES );

            // Then on the dimension of the FESpace (scalar field vs vectorial field)
            for( UInt iterDof( 0 ); iterDof != numberLocalDof; ++iterDof )
            {
                // Find the ID of the considered DOF
                ID globalDofID( M_pressureFespace->dofPtr()->localToGlobalMap( iterElement, iterDof ) );
                // Null force for pressure component
                //#pragma GCC diagnostic ignored "-Wconversion"
                f_interpolated->setCoefficient( globalDofID + k * numberTotalPDof +
                                                 ( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * numberTotalDof, 0 );
                //#pragma GCC diagnostic warning "-Wconversion"
            }
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 2>::
interpolate( const vector_type& fxrtheta, const vector_ptrType& f_interpolated )
{
    UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();
    
    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();
    UInt ndofpFE = DataPressFESpace.numTotalDof();
    
    Real Theta( M_modalbasis->Theta() );
    int mx( M_modalbasis->mx() );
    int mr( M_modalbasis->mr() );
    int mtheta( M_modalbasis->mtheta() );
    int mp( M_modalbasis->mp() );

    const QuadratureRule* qrRho( M_modalbasis->qrRho() );
    const QuadratureRule* qrTheta( M_modalbasis->qrTheta() );
    MBMatrix_type Jacobian( M_modalbasis->map()->Jacobian() );
    MBMatrix_type xphirho( M_modalbasis->xphirho() );
    MBMatrix_type xphitheta( M_modalbasis->xphitheta() );
    MBMatrix_type rphirho( M_modalbasis->rphirho() );
    MBMatrix_type rphitheta( M_modalbasis->rphitheta() );
    MBMatrix_type thetaphirho( M_modalbasis->thetaphirho() );
    MBMatrix_type thetaphitheta( M_modalbasis->thetaphitheta() );
    MBMatrix_type pphirho( M_modalbasis->pphirho() );
    MBMatrix_type pphitheta( M_modalbasis->pphitheta() );

    UInt index( 0 );
    
    for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes
    {
        if( s < ndofpFE )
        {
            index = 2 * s;
        }
        else
        {
            index = 2 * s - ndofuFE;
        }
        for ( UInt j( 0 ); j != mx; ++j ) // Ciclying over x-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                        ( *f_interpolated )( s + j * ndofuFE )  += fxrtheta( nrho + ntheta*nquadRho + index*nquadTheta*nquadRho ) *
                                                                   xphirho[j][nrho] * xphitheta[j][ntheta] * // phi_k
                                                                   Jacobian[s][ntheta] * // Jacobian
                                                                   qrRho->quadPointCoor( nrho, 0 ) * qrRho->weight( nrho ) * // r dr
                                                                   Theta * qrTheta->weight( ntheta ); // dtheta
                }
            }
        } // x-modes

        for ( UInt j( 0 ); j != mr; ++j ) // Ciclying over r-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                        ( *f_interpolated )( s + ( j+mx ) * ndofuFE )  +=
                        fxrtheta( ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho + index*nquadTheta*nquadRho ) * // u_old
                        rphirho[j][nrho] * rphitheta[j][ntheta] * // phi_k
                        Jacobian[s][ntheta] * // Jacobian
                        qrRho->quadPointCoor( nrho, 0 ) * qrRho->weight( nrho ) * // r dr
                        Theta *  qrTheta->weight( ntheta ); // dtheta
                }
            }
        } // r-modes

        for ( UInt j( 0 ); j != mtheta; ++j ) // Ciclying over theta-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                        ( *f_interpolated )( s + ( j+mx+mr ) * ndofuFE )  +=
                        fxrtheta( 2*ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho + index*nquadTheta*nquadRho ) *
                        thetaphirho[j][nrho] * thetaphitheta[j][ntheta] * // phi_k
                        Jacobian[s][ntheta] * // jacobian
                        qrRho->quadPointCoor( nrho, 0 ) * qrRho->weight( nrho ) * // r dr
                        Theta * qrTheta->weight( ntheta ); // dtheta
                }
            }
        } // theta-modes
    } // x-nodes
    
    for ( UInt s( 0 ); s != ndofpFE; ++s ) // on FE p-nodes
    {
        for ( UInt j( 0 ); j != mp; ++j ) // Ciclying over p-modes contribute
        {
            for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
            {
                for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
                {
                        ( *f_interpolated )( s + j*ndofpFE + ( mx+mr+mtheta ) * ndofuFE )  +=
                        fxrtheta( 3*ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho ) *
                        pphirho[j][nrho] * pphitheta[j][ntheta] * // phi_k
                        Jacobian[s][ntheta] * // jacobian
                        qrRho->quadPointCoor( nrho, 0 ) * qrRho->weight( nrho ) * // r dr
                        Theta * qrTheta->weight( ntheta ); // dtheta
                }
            }
        } // p-modes
    }// x-nodes
}

#endif
