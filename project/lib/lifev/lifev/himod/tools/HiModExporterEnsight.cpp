#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/himod/tools/HiModExporterEnsight.hpp>

namespace LifeV
{

std::vector< std::vector<Real> >				HiModExporterEnsight::M_points;
std::vector< std::vector<UInt> >				HiModExporterEnsight::M_tetrahedra;


HiModExporterEnsight::
HiModExporterEnsight( std::string fileName,
                      const UInt& mx, const UInt& mr, const UInt& mtheta, const UInt& mp,
                      const NSModalSpaceCircular& MB,
                      const Real& h, const UInt& Nelements,
                      const UInt& Ntimesteps, const Real& dt ):
                      M_uMeshSize( h ), M_Nelements( Nelements ), 
                      M_mx( mx ), M_mr( mr ), M_mtheta( mtheta ), M_mp( mp ),
                      M_quadrulerho( MB.qrRho() ), M_quadruletheta( MB.qrTheta() ),
                      M_xBasis( MB.xGbRhoTheta() ), M_rBasis( MB.rGbRhoTheta() ), 
                      M_tBasis( MB.thetaGbRhoTheta() ), M_pBasis( MB.pGbRhoTheta() )
{
    readGeometry();
    writeGeoFile();
    writeCaseFile( fileName, Ntimesteps, dt );

    M_xRadialBasis.resize( mx, vector_Type( MB.qrRho()->nbQuadPt(), 0 ) );
    M_rRadialBasis.resize( mr, vector_Type( MB.qrRho()->nbQuadPt(), 0 ) );
    M_tRadialBasis.resize( mtheta, vector_Type( MB.qrRho()->nbQuadPt(), 0 ) );
    M_pRadialBasis.resize( mp, vector_Type( MB.qrRho()->nbQuadPt(), 0 ) );
							    
    for( UInt i( 0 ); i != mx; ++i )
        for( UInt j( 0 ); j != MB.qrRho()->nbQuadPt(); ++j )
            M_xRadialBasis[i][j] = MB.xphirho( i, j );
							            
    for( UInt i( 0 ); i != mr; ++i )
        for( UInt j( 0 ); j != MB.qrRho()->nbQuadPt(); ++j )
            M_rRadialBasis[i][j] = MB.rphirho( i, j );
							            
    for( UInt i( 0 ); i != mtheta; ++i )
        for( UInt j( 0 ); j != MB.qrRho()->nbQuadPt(); ++j )
            M_tRadialBasis[i][j] = MB.thetaphirho( i, j );
							            
    for( UInt i( 0 ); i != mp; ++i )
        for( UInt j( 0 ); j != MB.qrRho()->nbQuadPt(); ++j )
            M_pRadialBasis[i][j] = MB.pphirho( i, j );
}

void
HiModExporterEnsight::
readGeometry()
{
	std::ifstream coord( "PointsCoord.txt" );
	std::ifstream tetra( "TetraLabels.txt" );
	
	vector_Type 		allPoints;
	vector_UInt_Type	allTetra;
	
	std::copy( std::istream_iterator<Real>( coord ),
       		    std::istream_iterator<Real>(),
       		    std::back_inserter<vector_Type >( allPoints ) );
       		    
	std::copy( std::istream_iterator<Real>( tetra ),
       		    std::istream_iterator<Real>(),
       		    std::back_inserter<vector_UInt_Type >( allTetra ) );

	// Each point has three coordinates
	UInt Npoints( allPoints.size() / 3 );
	// Netgen file has 12 columns
	UInt Ntetra( allTetra.size() / 12 );

	// M_points is zeros( Npoints, 3 );
	M_points.resize( Npoints, vector_Type( 3, 0. ) );
	// Each tetrahedron is identified by 10 points (vertices+middle points)
	M_tetrahedra.resize( Ntetra, vector_UInt_Type( 10, 0 ) );

	for( UInt i( 0 ); i != Npoints; ++i )
	{
		M_points[i][0] = allPoints[i*3+0];
		M_points[i][1] = allPoints[i*3+1];
		M_points[i][2] = allPoints[i*3+2];
	}
	
	for( UInt i( 0 ); i != Ntetra; ++i )
	{
		for( UInt node( 0 ); node != 10; ++node )
		{
			// skip the first two columns
			M_tetrahedra[i][node] = allTetra[i*12+2+node];
		}
	}
}

void
HiModExporterEnsight::
writeGeoFile()
{
	std::ofstream geo("geometry.geo");

	std::stringstream Npoints;
	Npoints << M_points.size();

	UInt digitNpoints( static_cast<UInt>( log10( (double) M_points.size() ) ) + 1 );
    
    std::stringstream space;
    space.str("");
        		
    for( UInt t( 0 ); t != 8 - digitNpoints; ++t )
    {
        space << " ";
    }
    
	geo << "Geometry file" << std::endl;
	geo << "Generated by LifeV" << std::endl;
	geo << "node id given" << std::endl;
	geo << "element id given" << std::endl;
	geo << "coordinates" << std::endl;
	geo << space.str() << Npoints.str() << std::endl;

	for( UInt s( 1 ); s != digitNpoints; ++s )
	{
		// set tabs
		std::stringstream tabs;
		tabs << "";
		
		for( UInt i( 0 ); i != digitNpoints - s; ++i )
		{
			tabs << " ";
		}
		
		std::stringstream x;
		std::stringstream y;
		std::stringstream z;
	
		for( UInt i( pow( 10, s - 1 ) - 1 ); i != pow( 10, s ) - 1; ++i )
		{
			UInt index( i );
		    geo << space.str() << tabs.str() << i + 1;
		
			if( M_points[index][0] >= 0 )
			{
				x << " ";
			}
			x.setf( std::ios_base::scientific, std::ios_base::floatfield );
			x.precision( 5 );
			x << M_points[index][0];
		
			if( M_points[index][1] >= 0 )
			{
				y << " ";
			}
			y.setf( std::ios_base::scientific, std::ios_base::floatfield );
			y.precision( 5 );
			y << M_points[index][1];
		
			if( M_points[index][2] >= 0 )
			{
				z << " ";
			}
			z.setf( std::ios_base::scientific, std::ios_base::floatfield );
			z.precision( 5 );
			z << M_points[index][2];
		
			geo << x.str() << y.str() << z.str() << std::endl;
		
			x.str("");
			y.str("");
			z.str("");
		}
		
		tabs.str("");
	}

	for( UInt i( pow( 10, digitNpoints - 1 ) - 1 ); i != M_points.size(); ++i )
	{
		std::stringstream x;
		std::stringstream y;
		std::stringstream z;
				
		    geo << space.str() << i + 1;

		    UInt index( i );
			if( M_points[index][0] >= 0 )
			{
				x << " ";
			}
			x.setf( std::ios_base::scientific, std::ios_base::floatfield );
			x.precision( 5 );
			x << M_points[index][0];

			if( M_points[index][1] >= 0 )
			{
				y << " ";
			}
			y.setf( std::ios_base::scientific, std::ios_base::floatfield );
			y.precision( 5 );
			y << M_points[index][1];

			if( M_points[index][2] >= 0 )
			{
				z << " ";
			}
			z.setf( std::ios_base::scientific, std::ios_base::floatfield );
			z.precision( 5 );
			z << M_points[index][2];

			geo << x.str() << y.str() << z.str() << std::endl;
		
			x.str("");
			y.str("");
			z.str("");
	}
	
	// Write tetra -----------------------------------------------------------------------
	std::stringstream Ntetra;
	Ntetra << M_tetrahedra.size();

	UInt digitNtetra( static_cast<UInt>( log10( (double) M_tetrahedra.size() ) ) + 1 );
	
	space.str("");
	for( UInt t( 0 ); t != 8 - digitNtetra; ++t )
    {
        space << " ";
    }
    
	geo << "part\t\t1" << std::endl;
	geo << "full geometry" << std::endl;
	geo << "tetra10" << std::endl;
	geo << space.str() << Ntetra.str() << std::endl;
	
	for( UInt s( 1 ); s != digitNtetra; ++s )
	{
		// set tabs
		std::stringstream tabs;
		tabs << "";
		
		for( UInt i( 0 ); i != digitNpoints - s; ++i )
		{
			tabs << " ";
		}
		
		for( UInt i( pow( 10, s - 1 ) - 1 ); i != pow( 10, s ) - 1; ++i )
		{
			UInt index( i );
		    geo << space.str() << tabs.str() << i + 1;
		    
		    // Print vertices of i-th tetra
		    for( UInt v( 0 ); v != 10; ++v )
		    {
		        std::stringstream coltabs;
        		coltabs.str("");
        		
		        for( UInt t( 0 ); t != 8 - (UInt)( log10( (double) M_tetrahedra[i][v] ) ) + 1; ++t )
		        {
		            coltabs << " ";
		        }
		        
		        geo << coltabs.str() << M_tetrahedra[i][v];
		    }
		    
		    geo << std::endl;
		}
	}
	
	for( UInt i( pow( 10, digitNtetra - 1 ) - 1 ); i != M_tetrahedra.size(); ++i )
	{
	    std::stringstream x;
        std::stringstream y;
		std::stringstream z;
				
		geo << space.str() << i + 1;

        // Print vertices of i-th tetra
		    for( UInt v( 0 ); v != 10; ++v )
		    {
		        std::stringstream coltabs;
        		coltabs.str("");
        		
		        for( UInt t( 0 ); t != 8 - (UInt)( log10( (double) M_tetrahedra[i][v] ) ) + 1; ++t )
		        {
		            coltabs << " ";
		        }
		        
		        geo << coltabs.str() << M_tetrahedra[i][v];
		    }
		    
		    geo << std::endl;

	}	
	
	geo.close();
}

void
HiModExporterEnsight::
writeCaseFile( std::string& fileName, const UInt& Ntimesteps, const Real& dt )
{ // IMPORTANT: INITIAL SOLUTION NOT INCLUDED
    fileName += ".case";
    std::ofstream fileCase( fileName.c_str() );

	std::stringstream Tsteps;
	Tsteps << Ntimesteps + 1;
	
	fileCase << "FORMAT" << std::endl;
	fileCase << "type: ensight" << std::endl;

	fileCase << "GEOMETRY" << std::endl;
	fileCase << "model: 1 geometry.geo" << std::endl;

	fileCase << "VARIABLE" << std::endl;
	fileCase << "vector per node: 1 velocity velocity.***.vct" << std::endl;
//	fileCase << "scalar per node: 1 pressure pressure.***.scl" << std::endl;
	
	fileCase << "TIME" << std::endl;
	fileCase << "time set: 1" << std::endl;
	fileCase << "number of steps: " << Tsteps.str() << std::endl;
	fileCase << "filename start number: 0" << std::endl;
	fileCase << "filename increment: 1" << std::endl;
	fileCase << "time values: " << std::endl;
	fileCase << "0 ";
	
	Tsteps.str("");
	for( UInt tIter( 0 ); tIter != Ntimesteps; ++tIter )
	{
	    Tsteps << ( tIter + 1 ) * dt;
	    fileCase << Tsteps.str() << " ";
	    Tsteps.str("");
	}
	
	fileCase.close();
	
	return;
}

void
HiModExporterEnsight::
writeSolution( const VectorEpetraStructured& sinCoeff, const VectorEpetraStructured& cosCoeff,
               const UInt& timeIter, const UInt& Nelements ) const
{
	const UInt udof( 2 * Nelements + 1 );
	const UInt pdof( Nelements + 1 );

	std::stringstream fileIter;
    std::string filenameVel( "velocity." );
//    std::string filenamePress( "pressure." );

    UInt nDigits( static_cast<UInt>( log10( (double) timeIter ) ) + 1 );

    if( nDigits > 3 )
    {
    	std::cout << "Error: too many time steps." <<std::endl;
    	return;
    }

    for( UInt s( 0 ); s != 3 - nDigits; ++s )
    {
        fileIter << "0";
    }
    fileIter << timeIter;

    filenameVel += fileIter.str();
    filenameVel += ".vct";
    // filenamePress += fileIter.str(); filenamePress += ".scl";

	vector_Type interpSin( ( M_mx + M_mr + M_mtheta ) * M_points.size(), 0.0 );
	vector_Type interpCos( ( M_mx + M_mr + M_mtheta ) * M_points.size(), 0.0 );
	vector_Type gridSinSolution( 3 * M_points.size(), 0 );
	vector_Type gridCosSolution( 3 * M_points.size(), 0 );
	vector_Type gridSolution( 3 * M_points.size(), 0 );

	matrix_pointType xbasis( M_points.size(), vector_Type( M_mx, 0.0 ) );
	matrix_pointType rbasis( M_points.size(), vector_Type( M_mr, 0.0 ) );
	matrix_pointType tbasis( M_points.size(), vector_Type( M_mtheta, 0.0 ) );
    matrix_pointType pbasis( M_points.size(), vector_Type( M_mp, 0.0 ) );

	interpolate( sinCoeff, interpSin );
	interpolate( cosCoeff, interpCos );

	evaluateModalBasis( xbasis, rbasis, tbasis, pbasis, sine );
	evaluateSolution( interpSin, xbasis, rbasis, tbasis, pbasis, gridSinSolution );
	
    evaluateModalBasis( xbasis, rbasis, tbasis, pbasis, cosine );
	evaluateSolution( interpCos, xbasis, rbasis, tbasis, pbasis, gridCosSolution );

	std::transform( gridSinSolution.begin(), gridSinSolution.end(), gridCosSolution.begin(), gridSolution.begin(), std::plus<Real>() );
	
	// Swap to cartesian coordinates
	vector_Type gridCartesianSolution( gridSolution );
	for(  UInt ix = 0; ix!= M_points.size(); ++ix )
	{
	    Real theta( atan( M_points[ix][3] / M_points[ix][2] ) );    
	    gridCartesianSolution[M_points.size()+ix] = gridSolution[M_points.size()+ix] * cos( theta ) - gridSolution[M_points.size()*2+ix] * sin( theta );
	    gridCartesianSolution[2*M_points.size()+ix] = gridSolution[M_points.size()+ix] * sin( theta ) + gridSolution[M_points.size()*2+ix] * cos( theta );
	}
	
	// write the evaluation of the solution on the mesh
	std::ofstream velocity( filenameVel.c_str(), std::ofstream::out );
//    std::ofstream pressure( filenamePress );
		
		
    std::stringstream vx;
    std::stringstream vy;
    std::stringstream vz;
    velocity << "Vector per node" << std::endl;
    for( UInt ix = 0; ix!= M_points.size(); ++ix )
    {
        if( gridCartesianSolution[ix] >= 0 )
        {
            vx << " ";
        }
        else
        {
            vx << "";
        }
        vx.setf( std::ios_base::scientific, std::ios_base::floatfield );
        vx.precision( 5 );
        if( gridCartesianSolution[ix] < pow( 10, -99 ) && gridCartesianSolution[ix] >= 0 )
        {
            vx << 0.;
        }
        else if( gridCartesianSolution[ix] > - pow( 10, -99 ) && gridCartesianSolution[ix] <= 0 )
        {
            vx << -0.;
        }
        else
        {
            vx << gridCartesianSolution[ix];
        }
       			
        if( gridCartesianSolution[ix+M_points.size()] >= 0 )
        {
            vy << " ";
        }
        else
        {
            vy << "";
        }
        vy.setf( std::ios_base::scientific, std::ios_base::floatfield );
        vy.precision( 5 );
        if( gridCartesianSolution[ix + M_points.size()] < pow( 10, -99 ) && 
            gridCartesianSolution[ix + M_points.size()] >= 0 )
        {
            vy << 0.;
        }
        else if( gridCartesianSolution[ix + M_points.size()] > - pow( 10, -99 ) && 
                 gridCartesianSolution[ix + M_points.size()] <= 0 )
        {
            vy << -0.;
        }
        else
        {
            vy << gridCartesianSolution[ix + M_points.size()];
        }
	            
        if( gridCartesianSolution[ix+M_points.size()*2] >= 0 )
        {
            vz << " ";
        }
        else
        {
            vz << "";
        }
        vz.setf( std::ios_base::scientific, std::ios_base::floatfield );
        vz.precision( 5 );
        if( gridCartesianSolution[ix + M_points.size()*2] < pow( 10, -99 ) &&
            gridCartesianSolution[ix + M_points.size()*2] >= 0 )
        {
            vz << 0.;
        }
        else if( gridCartesianSolution[ix + M_points.size()*2] > - pow( 10, -99 ) && 
                 gridCartesianSolution[ix + M_points.size()*2] <= 0 )
        {
            vy << -0.;
        }
        else
        {
            vz << gridCartesianSolution[ix + M_points.size()*2];
        }
	            
        velocity << vx.str() << vy.str() << vz.str();
				
        if( remainder( ix+1, 2 ) == 0 )
        {
            velocity << std::endl;
        }
				
        vx.str("");
        vy.str("");
        vz.str("");
    }
	
/*	pressure << "Scalar per node" << std::endl;
	for( UInt ix = 3 * M_points.size(); ix!= 4 * M_points.size(); ++ix)
	{
				pressure <<
				gridSolution[ix] << std::endl;
	}
*/			
	
    velocity.close();	
//	pressure.close();
}


void
HiModExporterEnsight::
interpolate( const VectorEpetraStructured& coeff, vector_Type& interpolation ) const
{
    UInt index;
    UInt udof( 2 * M_Nelements + 1 );

        for( UInt i( 0 ); i != M_points.size(); ++i )
        {
            // Identify left node
            UInt s( static_cast<UInt>( M_points[i][0] / M_uMeshSize ) );
            if( remainder( s, 2 ) == 0 )
            {
                index = s / 2;
            }
            else
            {
                index = ( s + 1 ) / 2 + M_Nelements;
            }
            
            if( s == 2 * M_Nelements ) // right boundary
            {
                // x-velocity
                for( UInt k( 0 ); k != M_mx; ++k )
                {
                    interpolation[i+k*M_points.size()] = coeff[index+k*udof];
                }
                
                // r-velocity
                for( UInt k( 0 ); k != M_mr; ++k )
                {
                    interpolation[i+(M_mx+k)*M_points.size()] = coeff[index+(M_mx+k)*udof];
                }
                
                // theta-velocity
                for( UInt k( 0 ); k != M_mtheta; ++k )
                {
                    interpolation[i+(M_mx+M_mr+k)*M_points.size()] = coeff[index+(M_mx+M_mr+k)*udof];
                }
            }
            else if( s == 2 * M_Nelements - 1 ) // middle point of the last interval
            {
                // x-velocity
                for( UInt k( 0 ); k != M_mx; ++k )
                {
                    interpolation[i+k*M_points.size()] = ( coeff[M_Nelements+k*udof] - coeff[index+k*udof] ) / ( M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+k*udof];
                }
                
                // r-velocity
                for( UInt k( 0 ); k != M_mr; ++k )
                {
                    interpolation[i+(M_mx+k)*M_points.size()] = ( coeff[M_Nelements+(M_mx+k)*udof] - coeff[index+(M_mx+k)*udof] ) / ( M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+(M_mx+k)*udof];
                }
                
                // theta-velocity
                for( UInt k( 0 ); k != M_mtheta; ++k )
                {
                    interpolation[i+(M_mx+M_mr+k)*M_points.size()] = ( coeff[M_Nelements+(M_mx+M_mr+k)*udof] - coeff[index+(M_mx+M_mr+k)*udof] ) / ( M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+(M_mx+M_mr+k)*udof];
                }
            }
            else
            {
                // x-velocity
                for( UInt k( 0 ); k != M_mx; ++k )
                {
                    interpolation[i+k*M_points.size()] = ( coeff[index+1+k*udof] - coeff[index+k*udof] ) / ( 2 * M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+k*udof];
                }
                
                // r-velocity
                for( UInt k( 0 ); k != M_mr; ++k )
                {
                    interpolation[i+(M_mx+k)*M_points.size()] = ( coeff[index+1+(M_mx+k)*udof] - coeff[index+(M_mx+k)*udof] ) / ( 2 * M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+(M_mx+k)*udof];
                }
                
                // theta-velocity
                for( UInt k( 0 ); k != M_mtheta; ++k )
                {
                    interpolation[i+(M_mx+M_mr+k)*M_points.size()] = ( coeff[index+1+(M_mx+M_mr+k)*udof] - coeff[index+(M_mx+M_mr+k)*udof] ) / ( 2 * M_uMeshSize ) *
                                                ( M_points[i][0] - s * M_uMeshSize ) +
                                                coeff[index+(M_mx+M_mr+k)*udof];
                }
            }
            
        }

    return; 
}

void
HiModExporterEnsight::
evaluateModalBasis( matrix_pointType& xbasis,
                         matrix_pointType& rbasis,
                         matrix_pointType& tbasis,
                         matrix_pointType& pbasis,
                         const int& trigonometricBasis ) const
{	
    vector_Type x( 2, 0. );
    for( UInt i( 0 ); i != M_points.size(); ++i )
    {
        x[0] = M_points[i][1]; // y 
        x[1] = M_points[i][2]; // z
        
        for( UInt j( 0 ); j != M_mx; ++j )
        {
            xbasis[i][j] = M_xBasis->evalSinglePoint( j, x, M_quadrulerho, M_xRadialBasis, trigonometricBasis );
        }
        
        for( UInt j( 0 ); j != M_mr; ++j )
        {
            rbasis[i][j] = M_rBasis->evalSinglePoint( j, x, M_quadrulerho, M_rRadialBasis, trigonometricBasis );
        }
        
        for( UInt j( 0 ); j != M_mtheta; ++j )
        {
            tbasis[i][j] = M_tBasis->evalSinglePoint( j, x, M_quadrulerho, M_tRadialBasis, trigonometricBasis );
        }
        
        for( UInt j( 0 ); j != M_mp; ++j )
        {
            pbasis[i][j] = M_pBasis->evalSinglePoint( j, x, M_quadrulerho, M_pRadialBasis, trigonometricBasis );
        }
    }
    return;
}

void
HiModExporterEnsight::
evaluateSolution( const vector_Type& interp,
	                        const matrix_pointType& xbasis,
	                        const matrix_pointType& rbasis,
	                        const matrix_pointType& tbasis,
	                        const matrix_pointType& pbasis,
	                        vector_Type& gridSolution ) const
{

    UInt N( M_points.size() );
    
    for( UInt i( 0 ); i != N; ++i )
    {
        for( UInt k( 0 ); k != M_mx; ++k )
        {
            gridSolution[i] += interp[i+k*N] * xbasis[i][k];
        }
        for( UInt k( 0 ); k != M_mr; ++k )
        {
            gridSolution[i+N] += interp[i+(M_mx+k)*N] * rbasis[i][k];
        }
        for( UInt k( 0 ); k != M_mtheta; ++k )
        {
            gridSolution[i+2*N] += interp[i+(M_mx+M_mr+k)*N] * tbasis[i][k];
        }
    }
    
    return;
}

} // end namespace
