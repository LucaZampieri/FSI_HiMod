#ifndef __HMAINTERPOLATION_HPP__
#define __HMAINTERPOLATION_HPP__

	template< typename mesh_type, typename matrix_type, typename vector_type>  
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR100xx( boost::shared_ptr<VectorEpetra>& R100xx,
						const matrix_ptrType& advection,
						const UInt& k,
						const UInt& j,
						const Real& nu )
	{
	    VectorEpetra R100xxx( R100xx->map() );
	    R100xxx *= 0.0;
	    
		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R100xxx[s] = M_modalbasis->compute_r100xxx( k, j, s );
		}
		*R100xx = ( (*advection) * R100xxx );
		*R100xx += M_modalbasis->compute_r100xx( k, j, nu );
		R100xx->globalAssemble();
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR000xx( boost::shared_ptr<VectorEpetra>& R000xx,
						const matrix_ptrType& xAdvection,
						const matrix_ptrType& rAdvection,
						const matrix_ptrType& tAdvection,
						const UInt& k,
						const UInt& j,
						const Real& nu,
						const Real& alpha )
	{
	    VectorEpetra R000xxx( R000xx->map() );
	    VectorEpetra R000xxr( R000xx->map() );
	    VectorEpetra R000xxt( R000xx->map() );
	    R000xxx *= 0.0;
	    R000xxr *= 0.0;
	    R000xxt *= 0.0;
	
	
		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R000xxx[s] = M_modalbasis->compute_r000xxx( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mr(); ++s )
		{
			R000xxr[s] = M_modalbasis->compute_r000xxr( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mtheta(); ++s )
		{
			R000xxt[s] = M_modalbasis->compute_r000xxt( k, j, s );
		}
		
		*R000xx = ( (*xAdvection) * R000xxx );
		*R000xx += ( (*rAdvection) * R000xxr );
		*R000xx += ( (*tAdvection) * R000xxt );
		*R000xx += M_modalbasis->compute_r000xx( k, j, nu, alpha );
        
		R000xx->globalAssemble();
		
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR100rr( boost::shared_ptr<VectorEpetra>& R100rr,
						const matrix_ptrType& advection,
						const UInt& k,
						const UInt& j,
						const Real& nu )
	{
		VectorEpetra R100rrx( R100rr->map() );
        R100rrx *= 0.0;

		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R100rrx[s] = M_modalbasis->compute_r100rrx( k, j, s );
		}
		*R100rr = ( (*advection) * R100rrx );
		*R100rr += M_modalbasis->compute_r100rr( k, j, nu );
		R100rr->globalAssemble();
		
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR000rr( boost::shared_ptr<VectorEpetra>& R000rr,
						const matrix_ptrType& xAdvection,
						const matrix_ptrType& rAdvection,
						const matrix_ptrType& tAdvection,
						const UInt& k,
						const UInt& j,
						const Real& nu,
						const Real& alpha )
	{
		VectorEpetra R000rrx( R000rr->map() );
	    VectorEpetra R000rrr( R000rr->map() );
	    VectorEpetra R000rrt( R000rr->map() );
	    R000rrx *= 0.0;
	    R000rrr *= 0.0;
	    R000rrt *= 0.0;
	
		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R000rrx[s] = M_modalbasis->compute_r000rrx( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mr(); ++s )
		{
			R000rrr[s] = M_modalbasis->compute_r000rrr( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mtheta(); ++s )
		{
			R000rrt[s] = M_modalbasis->compute_r000rrt( k, j, s );
		}
		
		*R000rr = ( (*xAdvection) * R000rrx );
		*R000rr += ( (*rAdvection) * R000rrr );
		*R000rr += ( (*tAdvection) * R000rrt );
		*R000rr += M_modalbasis->compute_r000rr( k, j, nu, alpha );

		R000rr->globalAssemble();
		
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR100tt( boost::shared_ptr<VectorEpetra>& R100tt,
						const matrix_ptrType& advection,
						const UInt& k,
						const UInt& j,
						const Real& nu )
	{
		VectorEpetra R100ttx( R100tt->map() );
		R100ttx *= 0.0;
	    
		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R100ttx[s] = M_modalbasis->compute_r100ttx( k, j, s );
		}
		*R100tt = ( (*advection) * R100ttx );
		*R100tt += M_modalbasis->compute_r100tt( k, j, nu );
		R100tt->globalAssemble();
		
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR000tt( boost::shared_ptr<VectorEpetra>& R000tt,
						const matrix_ptrType& xAdvection,
						const matrix_ptrType& rAdvection,
						const matrix_ptrType& tAdvection,
						const UInt& k,
						const UInt& j,
						const Real& nu,
						const Real& alpha )
	{
		VectorEpetra R000ttx( R000tt->map() );
	    VectorEpetra R000ttr( R000tt->map() );
	    VectorEpetra R000ttt( R000tt->map() );
	    R000ttx *= 0.0;
	    R000ttr *= 0.0;
	    R000ttt *= 0.0;
	
		for( UInt s = 0; s != M_modalbasis->mx(); ++s )
		{
			R000ttx[s] = M_modalbasis->compute_r000ttx( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mr(); ++s )
		{
			R000ttr[s] = M_modalbasis->compute_r000ttr( k, j, s );
		}
		for( UInt s = 0; s != M_modalbasis->mtheta(); ++s )
		{
			R000ttt[s] = M_modalbasis->compute_r000ttt( k, j, s );
		}
		
		*R000tt = ( (*xAdvection) * R000ttx );
		*R000tt += ( (*rAdvection) * R000ttr );
		*R000tt += ( (*tAdvection) * R000ttt );
		*R000tt += M_modalbasis->compute_r000tt( k, j, nu, alpha );

		R000tt->globalAssemble();
		
		return;
	}

	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR000rt( boost::shared_ptr<VectorEpetra>& R000rt,
						const matrix_ptrType& advection,
						const UInt& k,
						const UInt& j,
						const Real& nu )
	{
		VectorEpetra R000rtt( R000rt->map() );
		R000rtt *= 0.0;
	    
		for( UInt s = 0; s != M_modalbasis->mtheta(); ++s )
		{
			R000rtt[s] = M_modalbasis->compute_r000rtt( k, j, s );
		}
		*R000rt = ( (*advection) * R000rtt );
		*R000rt += M_modalbasis->compute_r000rt( k, j, nu );

		R000rt->globalAssemble();
		
		return;
	}
    						  
	template< typename mesh_type, typename matrix_type, typename vector_type>
	void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
	interpolateR000tr( boost::shared_ptr<VectorEpetra>& R000tr,
						const matrix_ptrType& advection,
						const UInt& k,
						const UInt& j,
						const Real& nu )
	{
		VectorEpetra R000trt( R000tr->map() );
		R000trt *= 0.0;
	    
		for( UInt s = 0; s != M_modalbasis->mtheta(); ++s )
		{
			R000trt[s] = M_modalbasis->compute_r000trt( k, j, s );
		}
		*R000tr = ( (*advection) * R000trt );
		*R000tr += M_modalbasis->compute_r000tr( k, j, nu );

		R000tr->globalAssemble();
		
		return;
	}


//Note:
//1) Added get methods in FESpace.hpp to get res FE pointer
//2) There's not a function to add a coefficient to a VectorEpetraStructured, it should be into VectorStructuredView.hpp like
//   in MatrixStructuredView.hpp
//3) Problem with vect.setCoefficient if you use parallel mode?
//4) MeshPartitioner.hpp contiene dei float!!!!
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
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

// ----------------- R( x ) -------------------------------
template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
interpolate( const function_Type& fx, const function_Type& fr, const function_Type& ftheta,
			const vector_ptrType& f_interpolated, const Real& t, const bool& b )
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
                nodalValues[iterDof] =  M_modalbasis->xFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), fx, k, b );
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
                nodalValues[iterDof] =  M_modalbasis->rFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), fr, k, b );
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
                nodalValues[iterDof] =  M_modalbasis->thetaFourierCoeffPointWise( t, interpCFE.quadNode( iterDof, 0 ), ftheta, k, b );
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
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
interpolate( const vector_type& fxrtheta, const vector_ptrType& f_interpolated )
{
	UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();
    
    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();
    UInt ndofpFE = DataPressFESpace.numTotalDof();
    
    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2 * M_PI );
    
    UInt index( 0 );
    
	for ( UInt s( 0 ); s != ndofuFE; ++s ) // on FE nodes (P1 first, then P2)
	{
            if( s < ndofpFE )
            {
                index = 2 * s; // Progressive index (alternating P1-P2 FE nodes)
            }
            else
            {
                index = 2 * s - ndofuFE;
            }
	    for ( UInt j( 0 ); j != M_modalbasis->mx(); ++j ) // Ciclying over x-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
	                (*f_interpolated)( s + j*ndofuFE )  += fxrtheta( nrho + ntheta*nquadRho + index * nquadTheta * nquadRho ) *
                                                               M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                                               M_modalbasis->Rho() * // u_old * R * rhat
	                                                       M_modalbasis->xphirho( j, nrho ) * normrho *
	                                                       M_modalbasis->xphitheta( j, ntheta ) * normtheta * // phi_k
                                                               M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                                               M_modalbasis->Theta() *
                                                               M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
	                                                       M_modalbasis->qrRho()->weight( nrho ) *
                                                               M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mr(); ++j ) // Ciclying over r-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
                        ( *f_interpolated )( s + ( j + M_modalbasis->mx() ) * ndofuFE )  +=
                                        fxrtheta( ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho +
                                                  index*nquadTheta*nquadRho ) *
                                        M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->Rho() * // u_old * R * rhat
                                        M_modalbasis->rphirho( j, nrho ) * normrho *
                                        M_modalbasis->rphitheta( j, ntheta ) * normtheta * // phi_k
                                        M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                        M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
                                        M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mtheta(); ++j ) // Ciclying over theta-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
                        ( *f_interpolated )( s + ( j + M_modalbasis->mx() + M_modalbasis->mr() ) * ndofuFE )  +=
                                        fxrtheta( 2*ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho + 
                                                     index*nquadTheta*nquadRho ) *
                                        M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->Rho() * // u_old * R * rhat
	                                M_modalbasis->thetaphirho( j, nrho ) * normrho *
	                                M_modalbasis->thetaphitheta( j, ntheta ) * normtheta * // phi_k
		             		M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
                                        M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
	for ( UInt s( 0 ); s != ndofpFE; ++s ) // on FE p-nodes
	{
	    for ( UInt j( 0 ); j != M_modalbasis->mp(); ++j ) // Ciclying over p-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + j*ndofpFE +
                                                     ( M_modalbasis->mx()+M_modalbasis->mr()+M_modalbasis->mtheta() )*ndofuFE )  +=
                                       fxrtheta( 3*ndofuFE*nquadTheta*nquadRho + nrho + ntheta*nquadRho + s*nquadTheta*nquadRho ) *
                                       M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->Rho() * // u_old * R * rhat
                                       M_modalbasis->pphirho( j, nrho ) * normrho *
                                       M_modalbasis->pphitheta( j, ntheta ) * normtheta * // phi_k
                                       M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
                                       M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
                                       M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
plainInterpolate( const vector_type& fxrtheta, const vector_ptrType& f_interpolated )
{
	UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();
    
    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();
    UInt ndofpFE = DataPressFESpace.numTotalDof();
    
    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2 * M_PI );
    
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
	    for ( UInt j( 0 ); j != M_modalbasis->mx(); ++j ) // Ciclying over x-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + j * ndofuFE )  += fxrtheta( nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
	                                                    M_modalbasis->xphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->xphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mr(); ++j ) // Ciclying over r-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + ( j + M_modalbasis->mx() ) * ndofuFE )  +=
		                								fxrtheta( ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
	                                                    M_modalbasis->rphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->rphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mtheta(); ++j ) // Ciclying over theta-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + ( j + M_modalbasis->mx() + M_modalbasis->mr() ) * ndofuFE )  +=
		                								fxrtheta( 2 * ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
	                                                    M_modalbasis->thetaphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->thetaphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
	for ( UInt s( 0 ); s != ndofpFE; ++s ) // on FE p-nodes
	{
	    for ( UInt j( 0 ); j != M_modalbasis->mp(); ++j ) // Ciclying over p-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + ( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * ndofuFE + j * ndofpFE )  +=
		                								fxrtheta( 3 * ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho ) *
	                                                    M_modalbasis->pphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->pphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->Jacobian()[nrho][ntheta] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void NSHiModAssembler<mesh_type, matrix_type, vector_type, 1>::
interpolate( const vector_type& fxrtheta, const vector_ptrType& f_interpolated, const bool& b )
{
	UInt nquadRho = M_modalbasis->qrRho()->nbQuadPt();
    UInt nquadTheta = M_modalbasis->qrTheta()->nbQuadPt();
    
    DOF DataVelFESpace( M_velocityFespace->dof() );
    DOF DataPressFESpace( M_pressureFespace->dof() );
    UInt ndofuFE = DataVelFESpace.numTotalDof();
    UInt ndofpFE = DataPressFESpace.numTotalDof();
    
    Real normrho = 1.0 / M_modalbasis->Rho();
    Real normtheta = 1.0 / sqrt( 2 * M_PI );
    
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
	    for ( UInt j( 0 ); j != M_modalbasis->mx(); ++j ) // Ciclying over x-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + j * ndofuFE )  += fxrtheta( nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->map()->xR()[s] * // u_old * R * rhat
	                                                    M_modalbasis->xphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->xphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->xJacobian()[s] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mr(); ++j ) // Ciclying over r-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + ( j + M_modalbasis->mx() ) * ndofuFE )  +=
		                								fxrtheta( ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->map()->xR()[s] * // u_old * R * rhat
	                                                    M_modalbasis->rphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->rphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->xJacobian()[s] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
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
	    for ( UInt j( 0 ); j != M_modalbasis->mtheta(); ++j ) // Ciclying over theta-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + ( j + M_modalbasis->mx() + M_modalbasis->mr() ) * ndofuFE )  +=
		                								fxrtheta( 2 * ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + index * nquadTheta * nquadRho ) *
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->map()->xR()[s] * // u_old * R * rhat
	                                                    M_modalbasis->thetaphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->thetaphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->xJacobian()[s] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
	
	for ( UInt s( 0 ); s != ndofpFE; ++s ) // on FE p-nodes
	{
	    for ( UInt j( 0 ); j != M_modalbasis->mp(); ++j ) // Ciclying over p-modes contribute
	    {
	        for ( UInt ntheta( 0 ); ntheta != nquadTheta; ++ntheta ) // on quadrature node nquadTheta
	        {
	            for ( UInt nrho( 0 ); nrho != nquadRho; ++nrho ) // on quadrature node nqadRho
	            {
		                ( *f_interpolated )( s + j * ndofpFE + ( M_modalbasis->mx() + M_modalbasis->mr() + M_modalbasis->mtheta() ) * ndofuFE )  +=
		                								fxrtheta( 3 * ndofuFE * nquadTheta * nquadRho + nrho + ntheta * nquadRho + s * nquadTheta * nquadRho ) *
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) * M_modalbasis->map()->xR()[s] * // u_old * R * rhat
	                                                    M_modalbasis->pphirho( j, nrho ) * normrho *
	                                                    M_modalbasis->pphitheta( j, ntheta ) * normtheta * // phi_k
		                								M_modalbasis->qrRho()->quadPointCoor( nrho, 0 ) *
	                                                    M_modalbasis->Theta() * M_modalbasis->map()->xJacobian()[s] * // rhat * jacobian
	                                                    M_modalbasis->qrRho()->weight( nrho ) * M_modalbasis->qrTheta()->weight( ntheta );
	            }
	        }
	    }
	}
}

#endif
