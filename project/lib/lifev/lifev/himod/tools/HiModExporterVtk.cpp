#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/himod/tools/HiModExporterVtk.hpp>

namespace LifeV
{

void
HiModExporterVtk::
writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D, const UInt& timeIter ) const
{
    const UInt udof( 2 * M_Nelements + 1 );
    const UInt pdof( M_Nelements + 1 );
    
    const UInt nquadRho( M_quadrulerho->nbQuadPt() );
    const UInt nquadTheta( M_quadruletheta->nbQuadPt() );
    
    Real Theta( 2 * M_PI );

    std::stringstream fileIter;
    std::string filenameVel( fileName );
    filenameVel += "_velocity.";
    std::string filenamePress( fileName );
    filenamePress += "_pressure.";

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
    filenameVel += ".vtk";

    filenamePress += fileIter.str();
    filenamePress += ".vtk";
    
    // Swap to cartesian coordinates
    VectorEpetraStructured gridCartesianSolution( solution_3D );
    for(  UInt ix = 0; ix!= udof; ++ix )
    {
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                UInt ny( ( ix + udof ) * nquadTheta * nquadRho + itheta * nquadRho + irho );
                UInt nz( ( ix + 2 * udof ) * nquadTheta * nquadRho + itheta * nquadRho + irho );
                Real theta( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta );    
                
                gridCartesianSolution[ny] = solution_3D[ny] * cos( theta ) - solution_3D[nz] * sin( theta );
                gridCartesianSolution[nz] = solution_3D[ny] * sin( theta ) + solution_3D[nz] * cos( theta );
            }
        }
    }
    
    
    // ______________________    velocity    _____________________________
    std::ofstream velocity( filenameVel.c_str(), std::ofstream::out );

    velocity << "# vtk DataFile Version 3.0" << std::endl;
    velocity << "LifeV output" << std::endl;
    velocity << "ASCII" << std::endl;
    velocity << "DATASET UNSTRUCTURED_GRID" << std::endl;
    velocity << "POINTS " << udof * nquadRho * nquadTheta << " float" << std::endl;

    for( UInt ix( 0 ); ix != udof; ++ix )
    {     
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                velocity <<
                ix * M_uMeshSize <<" "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * cos( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << " "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * sin( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << std::endl;
                                    
            }
        }
    }
    
    velocity << "CELLS " << ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << " " 
                            << ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) * 9 * ( udof - 1 ) << std::endl;

    UInt shift( nquadRho * nquadTheta );                            
    for( UInt s( 0 ); s != udof - 1; ++s )
    {
        // radial polygons
        for( UInt h( 0 ); h != nquadTheta - 1; ++h )
        {
            for( UInt n( 0 ); n != nquadRho - 1; ++n )
            {
                UInt i( s * shift + h * nquadRho + n );

                velocity << 8 << " " << i << " " << i + 1 << " " << i + 1 + nquadRho << " " << i + nquadRho << " "
                                     << i + shift << " " << i + 1 + shift << " " << i + 1 + nquadRho + shift << " " << i + nquadRho + shift << std::endl;
            }
        }
            // restoration of periodicity
        for( UInt n( 0 ); n != nquadRho - 1; ++n )
        {
            UInt i( s * shift + ( nquadTheta - 1 ) * nquadRho + n );

            velocity << 8 << " " << i << " " << i + 1 << " " << n + 1 + s * shift << " " << n + s * shift << " "
                                 << i + shift << " " << i + 1 + shift << " " << n + 1 + shift + s * shift << " " << n + shift + s * shift << std::endl;
        }
    }
    
    // inner cylindrical core
    for( UInt s( 0 ); s != udof - 1; ++s )
    {
        for( UInt v( 0 ); v != ( nquadTheta - 4 ) / 2 + 1; ++v )
        {
            UInt i( s * shift );
            velocity << 8 << " " << i << " " << i + ( 2 * v + 1 ) * nquadRho << " " << i + ( 2 * v + 2 ) * nquadRho << " " << i + ( 2 * v + 3 ) * nquadRho << " "
                                 << i + shift << " " << i + ( 2 * v + 1 ) * nquadRho + shift << " "
                                 << i + ( 2 * v + 2 ) * nquadRho + shift << " " << i + ( 2 * v + 3 ) * nquadRho + shift << std::endl;
        }
    }
      
    velocity << "CELL_TYPES " << ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << std::endl;
    for( UInt h( 0 ); h != ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ); ++h )
    {
        velocity << 12 << std::endl;
    }
    
    velocity << "POINT_DATA " << udof * nquadRho * nquadTheta << std::endl;
    velocity << "FIELD FieldData 1" << std::endl;
    velocity << "Velocity " << 3 << " " << udof * nquadRho * nquadTheta << " double" << std::endl;
    for( UInt ix( 0 ); ix != udof; ++ix )
    {     
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta ) 
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                velocity <<
                gridCartesianSolution[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << " " <<
                gridCartesianSolution[ (udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << " " <<
                gridCartesianSolution[ (2*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << std::endl;                    
            }
        }
    }
            
    velocity.close();

    // ______________________    pressure    _____________________________
    
    std::ofstream pressure( filenamePress.c_str(), std::ofstream::out );
    
    pressure << "# vtk DataFile Version 3.0" << std::endl;
    pressure << "LifeV output" << std::endl;
    pressure << "ASCII" << std::endl;
    pressure << "DATASET UNSTRUCTURED_GRID" << std::endl;
    pressure << "POINTS " << pdof * nquadRho * nquadTheta << " float" << std::endl;

    
    for( UInt ix( 0 ); ix != pdof; ++ix )
            {     
                for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
                {  
                    for( UInt irho( 0 ); irho != nquadRho; ++irho )
                    {
                        pressure <<
                        ix * 2 * M_uMeshSize <<" "<<
                        M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * cos( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << " "<<
                        M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * sin( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << std::endl;
                                            
                    }
                }
            }
    pressure << "CELLS " << ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << " " 
                            << ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) * 9 * ( pdof - 1 ) << std::endl;
                            
    for( UInt s( 0 ); s != pdof - 1; ++s )
    {
        // radial polygons
        for( UInt h( 0 ); h != nquadTheta - 1; ++h )
        {
            for( UInt n( 0 ); n != nquadRho - 1; ++n )
            {
                UInt i( s * shift + h * nquadRho + n );

                pressure << 8 << " " << i << " " << i + 1 << " " << i + 1 + nquadRho << " " << i + nquadRho << " "
                                     << i + shift << " " << i + 1 + shift << " " << i + 1 + nquadRho + shift << " " << i + nquadRho + shift << std::endl;
            }
        }
            // restoration of periodicity
        for( UInt n( 0 ); n != nquadRho - 1; ++n )
        {
            UInt i( s * shift + ( nquadTheta - 1 ) * nquadRho + n );

            pressure << 8 << " " << i << " " << i + 1 << " " << n + 1 + s * shift << " " << n + s * shift << " "
                                 << i + shift << " " << i + 1 + shift << " " << n + 1 + shift + s * shift << " " << n + shift + s * shift << std::endl;
        }
    }
    
    // inner cylindrical core
    for( UInt s( 0 ); s != pdof - 1; ++s )
    {
        for( UInt v( 0 ); v != ( nquadTheta - 4 ) / 2 + 1; ++v )
        {
            UInt i( s * shift );
            pressure << 8 << " " << i << " " << i + ( 2 * v + 1 ) * nquadRho << " " << i + ( 2 * v + 2 ) * nquadRho << " " << i + ( 2 * v + 3 ) * nquadRho << " "
                                 << i + shift << " " << i + ( 2 * v + 1 ) * nquadRho + shift << " "
                                 << i + ( 2 * v + 2 ) * nquadRho + shift << " " << i + ( 2 * v + 3 ) * nquadRho + shift << std::endl;
        }
    }
      
    pressure << "CELL_TYPES " << ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << std::endl;
    for( UInt h( 0 ); h != ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ); ++h )
    {
        pressure << 12 << std::endl;
    }
    
            
            pressure << "POINT_DATA " << pdof * nquadRho * nquadTheta << std::endl;
            pressure << "SCALARS pressure double 1" << std::endl;
            pressure << "LOOKUP_TABLE default" << std::endl;
            
            for( UInt ix( 0 ); ix != pdof; ++ix )
            {     
                for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta ) 
                {  
                    for( UInt irho( 0 ); irho != nquadRho; ++irho )
                    {
                        pressure << solution_3D[ (3*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << std::endl;
                    }
                }
            }            
    
    pressure.close();
}

void
HiModExporterVtk::
writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D, const UInt& timeIter, const bool xDependent ) const
{
    const UInt udof( 2 * M_Nelements + 1 );
    const UInt pdof( M_Nelements + 1 );
    
    const UInt nquadRho( M_quadrulerho->nbQuadPt() );
    const UInt nquadTheta( M_quadruletheta->nbQuadPt() );
    
    Real Theta( 2 * M_PI );

    std::stringstream fileIter;
    std::string filenameVel( fileName );
    filenameVel += "_velocity.";
    std::string filenamePress( fileName );
    filenamePress += "_pressure.";

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
    filenameVel += ".vtk";

    filenamePress += fileIter.str();
    filenamePress += ".vtk";
    
    // Swap to cartesian coordinates
    VectorEpetraStructured gridCartesianSolution( solution_3D );
    for(  UInt ix = 0; ix!= udof; ++ix )
    {
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                UInt ny( ( ix + udof ) * nquadTheta * nquadRho + itheta * nquadRho + irho );
                UInt nz( ( ix + 2 * udof ) * nquadTheta * nquadRho + itheta * nquadRho + irho );
                Real theta( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta );    
                
                gridCartesianSolution[ny] = solution_3D[ny] * cos( theta ) - solution_3D[nz] * sin( theta );
                gridCartesianSolution[nz] = solution_3D[ny] * sin( theta ) + solution_3D[nz] * cos( theta );
            }
        }
    }
    
    // ______________________    velocity    _____________________________
    std::ofstream velocity( filenameVel.c_str(), std::ofstream::out );

    velocity << "# vtk DataFile Version 3.0" << std::endl;
    velocity << "LifeV output" << std::endl;
    velocity << "ASCII" << std::endl;
    velocity << "DATASET UNSTRUCTURED_GRID" << std::endl;
    velocity << "POINTS " << udof * nquadRho * nquadTheta << " float" << std::endl;

    for( UInt ix( 0 ); ix != udof; ++ix )
    {     
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                Real x(ix * M_uMeshSize);
                Real z(M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta);
                velocity <<
                ix * M_uMeshSize <<" "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_fR( 0, x, 0, z, itheta ) * cos( z ) << " "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_fR( 0, x, 0, z, itheta ) * sin( z ) << std::endl;
                                    
            }
        }
    }
    
    velocity << "CELLS " << ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << " " 
                            << ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) * 9 * ( udof - 1 ) << std::endl;

    UInt shift( nquadRho * nquadTheta );                            
    for( UInt s( 0 ); s != udof - 1; ++s )
    {
        // radial polygons
        for( UInt h( 0 ); h != nquadTheta - 1; ++h )
        {
            for( UInt n( 0 ); n != nquadRho - 1; ++n )
            {
                UInt i( s * shift + h * nquadRho + n );

                velocity << 8 << " " << i << " " << i + 1 << " " << i + 1 + nquadRho << " " << i + nquadRho << " "
                                     << i + shift << " " << i + 1 + shift << " " << i + 1 + nquadRho + shift << " " << i + nquadRho + shift << std::endl;
            }
        }
            // restoration of periodicity
        for( UInt n( 0 ); n != nquadRho - 1; ++n )
        {
            UInt i( s * shift + ( nquadTheta - 1 ) * nquadRho + n );

            velocity << 8 << " " << i << " " << i + 1 << " " << n + 1 + s * shift << " " << n + s * shift << " "
                                 << i + shift << " " << i + 1 + shift << " " << n + 1 + shift + s * shift << " " << n + shift + s * shift << std::endl;
        }
    }
    
    // inner cylindrical core
    for( UInt s( 0 ); s != udof - 1; ++s )
    {
        for( UInt v( 0 ); v != ( nquadTheta - 4 ) / 2 + 1; ++v )
        {
            UInt i( s * shift );
            velocity << 8 << " " << i << " " << i + ( 2 * v + 1 ) * nquadRho << " " << i + ( 2 * v + 2 ) * nquadRho << " " << i + ( 2 * v + 3 ) * nquadRho << " "
                                 << i + shift << " " << i + ( 2 * v + 1 ) * nquadRho + shift << " "
                                 << i + ( 2 * v + 2 ) * nquadRho + shift << " " << i + ( 2 * v + 3 ) * nquadRho + shift << std::endl;
        }
    }
      
    velocity << "CELL_TYPES " << ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << std::endl;
    for( UInt h( 0 ); h != ( udof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ); ++h )
    {
        velocity << 12 << std::endl;
    }
    
    velocity << "POINT_DATA " << udof * nquadRho * nquadTheta << std::endl;
    velocity << "FIELD FieldData 1" << std::endl;
    velocity << "Velocity " << 3 << " " << udof * nquadRho * nquadTheta << " double" << std::endl;
    for( UInt ix( 0 ); ix != udof; ++ix )
    {     
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta ) 
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                velocity <<
                gridCartesianSolution[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << " " <<
                gridCartesianSolution[ (udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << " " <<
                gridCartesianSolution[ (2*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << std::endl;                    
            }
        }
    }
            
    velocity.close();

    // ______________________    pressure    _____________________________
    
    std::ofstream pressure( filenamePress.c_str(), std::ofstream::out );
    
    pressure << "# vtk DataFile Version 3.0" << std::endl;
    pressure << "LifeV output" << std::endl;
    pressure << "ASCII" << std::endl;
    pressure << "DATASET UNSTRUCTURED_GRID" << std::endl;
    pressure << "POINTS " << pdof * nquadRho * nquadTheta << " float" << std::endl;

    
    for( UInt ix( 0 ); ix != pdof; ++ix )
    {     
        for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
        {  
            for( UInt irho( 0 ); irho != nquadRho; ++irho )
            {
                Real x(ix * 2*M_uMeshSize);
                Real z(M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta);
                pressure <<
                x <<" "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_fR( 0, x, 0, z, itheta ) * cos( z ) << " "<<
                M_quadrulerho->quadPointCoor( irho, 0 ) * M_fR( 0, x, 0, z, itheta ) * sin( z ) << std::endl;
            }
        }
    }
    pressure << "CELLS " << ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << " " 
                            << ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) * 9 * ( pdof - 1 ) << std::endl;
                            
    for( UInt s( 0 ); s != pdof - 1; ++s )
    {
        // radial polygons
        for( UInt h( 0 ); h != nquadTheta - 1; ++h )
        {
            for( UInt n( 0 ); n != nquadRho - 1; ++n )
            {
                UInt i( s * shift + h * nquadRho + n );

                pressure << 8 << " " << i << " " << i + 1 << " " << i + 1 + nquadRho << " " << i + nquadRho << " "
                                     << i + shift << " " << i + 1 + shift << " " << i + 1 + nquadRho + shift << " " << i + nquadRho + shift << std::endl;
            }
        }
            // restoration of periodicity
        for( UInt n( 0 ); n != nquadRho - 1; ++n )
        {
            UInt i( s * shift + ( nquadTheta - 1 ) * nquadRho + n );

            pressure << 8 << " " << i << " " << i + 1 << " " << n + 1 + s * shift << " " << n + s * shift << " "
                                 << i + shift << " " << i + 1 + shift << " " << n + 1 + shift + s * shift << " " << n + shift + s * shift << std::endl;
        }
    }
    
    // inner cylindrical core
    for( UInt s( 0 ); s != pdof - 1; ++s )
    {
        for( UInt v( 0 ); v != ( nquadTheta - 4 ) / 2 + 1; ++v )
        {
            UInt i( s * shift );
            pressure << 8 << " " << i << " " << i + ( 2 * v + 1 ) * nquadRho << " " << i + ( 2 * v + 2 ) * nquadRho << " " << i + ( 2 * v + 3 ) * nquadRho << " "
                                 << i + shift << " " << i + ( 2 * v + 1 ) * nquadRho + shift << " "
                                 << i + ( 2 * v + 2 ) * nquadRho + shift << " " << i + ( 2 * v + 3 ) * nquadRho + shift << std::endl;
        }
    }
      
    pressure << "CELL_TYPES " << ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << std::endl;
    for( UInt h( 0 ); h != ( pdof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ); ++h )
    {
        pressure << 12 << std::endl;
    }
    
            
            pressure << "POINT_DATA " << pdof * nquadRho * nquadTheta << std::endl;
            pressure << "SCALARS pressure double 1" << std::endl;
            pressure << "LOOKUP_TABLE default" << std::endl;
            
            for( UInt ix( 0 ); ix != pdof; ++ix )
            {     
                for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta ) 
                {  
                    for( UInt irho( 0 ); irho != nquadRho; ++irho )
                    {
                        pressure << solution_3D[ (3*udof + ix) * nquadTheta * nquadRho + itheta * nquadRho + irho ] << std::endl;
                    }
                }
            }            
    
    pressure.close();
}

void
HiModExporterVtk::
writeSolution( std::string fileName, const VectorEpetraStructured& solution_3D ) const
{
    const UInt dof( M_Nelements + 1 );
    
    const UInt nquadRho( M_quadrulerho->nbQuadPt() );
    const UInt nquadTheta( M_quadruletheta->nbQuadPt() );
    
    Real Theta( 2 * M_PI );

    std::string filenameVtk( fileName );

    filenameVtk += ".vtk";
    
    std::ofstream solution( filenameVtk.c_str(), std::ofstream::out );
    
    solution << "# vtk DataFile Version 3.0" << std::endl;
    solution << "LifeV output" << std::endl;
    solution << "ASCII" << std::endl;
    solution << "DATASET UNSTRUCTURED_GRID" << std::endl;
    solution << "POINTS " << dof * nquadRho * nquadTheta << " float" << std::endl;

    
    for( UInt ix( 0 ); ix != dof; ++ix )
            {     
                for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta )
                {  
                    for( UInt irho( 0 ); irho != nquadRho; ++irho )
                    {
                        solution <<
                        ix * M_uMeshSize <<" "<<
                        M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * cos( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << " "<<
                        M_quadrulerho->quadPointCoor( irho, 0 ) * M_R * sin( M_quadruletheta->quadPointCoor( itheta, 0 ) * Theta ) << std::endl;
                                            
                    }
                }
            }
    solution << "CELLS " << ( dof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << " " 
                            << ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) * 9 * ( dof - 1 ) << std::endl;

    UInt shift( nquadRho * nquadTheta );        
    for( UInt s( 0 ); s != dof - 1; ++s )
    {
        // radial polygons
        for( UInt h( 0 ); h != nquadTheta - 1; ++h )
        {
            for( UInt n( 0 ); n != nquadRho - 1; ++n )
            {
                UInt i( s * shift + h * nquadRho + n );

                solution << 8 << " " << i << " " << i + 1 << " " << i + 1 + nquadRho << " " << i + nquadRho << " "
                                     << i + shift << " " << i + 1 + shift << " " << i + 1 + nquadRho + shift << " " << i + nquadRho + shift << std::endl;
            }
        }
            // restoration of periodicity
        for( UInt n( 0 ); n != nquadRho - 1; ++n )
        {
            UInt i( s * shift + ( nquadTheta - 1 ) * nquadRho + n );

            solution << 8 << " " << i << " " << i + 1 << " " << n + 1 + s * shift << " " << n + s * shift << " "
                                 << i + shift << " " << i + 1 + shift << " " << n + 1 + shift + s * shift << " " << n + shift + s * shift << std::endl;
        }
    }
    
    // inner cylindrical core
    for( UInt s( 0 ); s != dof - 1; ++s )
    {
        for( UInt v( 0 ); v != ( nquadTheta - 4 ) / 2 + 1; ++v )
        {
            UInt i( s * shift );
            solution << 8 << " " << i << " " << i + ( 2 * v + 1 ) * nquadRho << " " << i + ( 2 * v + 2 ) * nquadRho << " " << i + ( 2 * v + 3 ) * nquadRho << " "
                                 << i + shift << " " << i + ( 2 * v + 1 ) * nquadRho + shift << " "
                                 << i + ( 2 * v + 2 ) * nquadRho + shift << " " << i + ( 2 * v + 3 ) * nquadRho + shift << std::endl;
        }
    }
      
    solution << "CELL_TYPES " << ( dof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ) << std::endl;
    for( UInt h( 0 ); h != ( dof - 1 ) * ( ( nquadRho - 1 ) * nquadTheta + ( nquadTheta - 4 ) / 2 + 1 ); ++h )
    {
        solution << 12 << std::endl;
    }
    
            
            solution << "POINT_DATA " << dof * nquadRho * nquadTheta << std::endl;
            solution << "SCALARS solution double 1" << std::endl;
            solution << "LOOKUP_TABLE default" << std::endl;
            
            for( UInt ix( 0 ); ix != dof; ++ix )
            {     
                for( UInt itheta( 0 ); itheta != nquadTheta; ++itheta ) 
                {  
                    for( UInt irho( 0 ); irho != nquadRho; ++irho )
                    {
                        solution << solution_3D[ ix * nquadTheta * nquadRho + itheta * nquadRho + irho ] << std::endl;
                    }
                }
            }            
    
    solution.close();
}

} // end namespace
