#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/himod/tools/ReferenceMap.hpp>
#include <functional>
#include <boost/bind.hpp>

namespace LifeV
{

    // Constructor for patient-specific geometries
    ReferenceMap::
    ReferenceMap( const std::string& Rfile, const std::string& dxRfile, const std::string& dtRfile,
                  const MapEpetra& epetraMap, const UInt& dof, const Real& h,
                  const QuadratureRule* quadrho, const QuadratureRule* quadtheta ) :
    M_quadruleRho( quadrho ), M_quadruleTheta( quadtheta )
    {
        // We actually do not need these, but they have to be instantiated because they are const members
        M_xJr = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xJtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xDr = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xDtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xDrtheta = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xDthetar = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xJacobian = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xR = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        M_xdR = static_cast<vector_ptrType>( new vector_Type( epetraMap, Repeated ) );
        
        UInt dimTheta = M_quadruleTheta->nbQuadPt();
               
        // Read geometry from file
        readMatrix(   Rfile, dof, dimTheta, M_R );
        readMatrix( dxRfile, dof, dimTheta, M_dxR );
        readMatrix( dtRfile, dof, dimTheta, M_dtR );
        
        // Instantiate Map containers
        M_Jr.resize( dof );
        M_Jtheta.resize( dof );
        M_Dr.resize( dof );
        M_Dtheta.resize( dof );
        M_Drtheta.resize( dof );
        M_Dthetar.resize( dof );
        M_Jacobian.resize( dof );
        
        for( UInt i = 0; i != dof; ++i )
        {
            M_Jr[i].resize( dimTheta );
            M_Jtheta[i].resize( dimTheta );
            M_Dr[i].resize( dimTheta );
            M_Dtheta[i].resize( dimTheta );
            M_Drtheta[i].resize( dimTheta );
            M_Dthetar[i].resize( dimTheta );
            M_Jacobian[i].resize( dimTheta );
        }
        
        UInt index(0);
        // Initialize Map containers
        for( UInt i = 0; i != dof; ++i )
        {
            
            for( UInt j = 0; j != M_quadruleTheta->nbQuadPt(); ++j )
            {

                 M_Jr[i][j]       = 1. / M_R[i][j];
                 M_Jtheta[i][j]   = M_Jr[i][j]; //1. / M_R[index][j];
                 M_Dr[i][j]       = -M_dxR[i][j] * M_Jr[i][j]; // -M_dxR / M_R[index][j];
                 M_Dtheta[i][j]   = 0;
                 M_Drtheta[i][j]  = 0;
                 M_Dthetar[i][j]  = -M_dtR[i][j] * M_Jr[i][j] * M_Jr[i][j]; // -M_dtR / M_R^2;
                 M_Jacobian[i][j] = fabs(1. / ( M_Jr[i][j]*M_Jtheta[i][j] - M_Drtheta[i][j]*M_Dthetar[i][j] ));
            }
        }
        
        // Wrap functions
        M_Rfunctor.setData(M_R);
        M_Rfunctor.setH(h);
        M_Rfunctor.setScaling(0);
        M_Rfunctor.setThetaQuadRule( *M_quadruleTheta );
        M_functionR = boost::bind( &MapFunctor::f, &M_Rfunctor, _1, _2, _3, _4, _5 );
        
        M_RhatFunctor.setData(M_R);
        M_RhatFunctor.setH(h);
        M_RhatFunctor.setScaling(1);
        M_RhatFunctor.setThetaQuadRule( *M_quadruleTheta );
        M_inverseRhat = boost::bind( &MapFunctor::f, &M_RhatFunctor, _1, _2, _3, _4, _5 );
        
        M_Jfunctor.setData(M_Jacobian);
        M_Jfunctor.setH(h);
        M_Jfunctor.setScaling(0);
        M_Jfunctor.setThetaQuadRule( *M_quadruleTheta );
        M_functionJacobian = boost::bind( &MapFunctor::f, &M_Jfunctor, _1, _2, _3, _4, _5 );
        
        /* REMARK: To initialize these functions, more functors need to be added as private members
        MapFunctor dxRfunctor( M_dxR, h );
        M_functionDxR = boost::bind( &MapFunctor::f, &dxRfunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor dtRfunctor( M_dtR, h );
        M_functionDtR = boost::bind( &MapFunctor::f, &dtRfunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor JrFunctor( M_Jr, h );
        M_functionJr = boost::bind( &MapFunctor::f, &JrFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor JthetaFunctor( M_Jtheta, h );
        M_functionJtheta = boost::bind( &MapFunctor::f, &JthetaFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor DrFunctor( M_Dr, h );
        M_functionDr = boost::bind( &MapFunctor::f, &DrFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor DthetaFunctor( M_Jr, h );
        M_functionDtheta = boost::bind( &MapFunctor::f, &DthetaFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor DrthetaFunctor( M_Drtheta, h );
        M_functionDrtheta = boost::bind( &MapFunctor::f, &DrthetaFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor DthetarFunctor( M_Dthetar, h );
        M_functionDthetar = boost::bind( &MapFunctor::f, &DthetarFunctor, _1, _2, _3, _4, _5 );
        
        MapFunctor jFunctor( M_Jr, h );
        M_functionJacobian = boost::bind( &MapFunctor::f, &jFunctor, _1, _2, _3, _4, _5 );
        */
    }
    
    void ReferenceMap::
    readMatrix( const std::string fileName, const UInt& nR, const UInt& nC, MBMatrix_type& M )
    {
        // Define file and resize containers
        std::ifstream file( fileName.c_str() );
        std::vector<Real> Mvec;
        M.resize( nR );
        
        for( UInt i( 0 ); i != nR; ++ i )
        {
            M[i].resize( nC );
        }
        
        // Read matrix in vector form
        std::copy( std::istream_iterator<Real> (file),
                   std::istream_iterator<Real> (),
                   std::back_inserter< std::vector<Real> > (Mvec) );
                   
        for( UInt k(0); k != nR*nC; ++k )
        {
            UInt col( k/nR );
            UInt index( k%nR );
            UInt row;
            
            if( index%2 == 0 )
            {
                row = index/2;
            }
            else
            {
                row = (index+nR)/2;
            }
            
            M[row][col] = Mvec[k];
        }
        
        return;
    }
                       
    void ReferenceMap::evaluateMap( const Real& t )
    {
        // Set the dimensions
        UInt dimRho = M_quadruleRho->nbQuadPt();
        UInt dimTheta = M_quadruleTheta->nbQuadPt();
        
        // Resize loops
        M_Jr.resize( dimRho );
        M_Jtheta.resize( dimRho );
        M_Dr.resize( dimRho );
        M_Dtheta.resize( dimRho );
        M_Drtheta.resize( dimRho );
        M_Dthetar.resize( dimRho );
        M_Jacobian.resize( dimRho );
        
        for( UInt i = 0; i != dimRho; ++i )
        {
            M_Jr[i].resize( dimTheta );
            M_Jtheta[i].resize( dimTheta );
            M_Dr[i].resize( dimTheta );
            M_Dtheta[i].resize( dimTheta );
            M_Drtheta[i].resize( dimTheta );
            M_Dthetar[i].resize( dimTheta );
            M_Jacobian[i].resize( dimTheta );
        }

        // Evaluation loops
        for( UInt i = 0; i != dimRho; ++i )
        {
            Real rh( M_quadruleRho->quadPointCoor( i, 0 ) );
            for( UInt j = 0; j != dimTheta; ++j )
            {
                Real theta( 2*M_PI*M_quadruleTheta->quadPointCoor( j, 0 ) );

                M_Jr[i][j] = M_functionJr( t, 0, rh, theta, 0 );
                M_Jtheta[i][j] = M_functionJtheta( t, 0, rh, theta, 0 );
                M_Dr[i][j] = M_functionDr( t, 0, rh, theta, 0 );
                M_Dtheta[i][j] = M_functionDtheta( t, 0, rh, theta, 0 );
                M_Drtheta[i][j] = M_functionDrtheta( t, 0, rh, theta, 0 );
                M_Dthetar[i][j] = M_functionDthetar( t, 0, rh, theta, 0 );
                M_Jacobian[i][j] = M_functionJacobian( t, 0, rh, theta, 0 );

            }
        }
    } // end evaluateMap
    
    void ReferenceMap::evaluateFtensor( const Real& t, const Real& h, const UInt& dof )
    {
        // Set the dimensions
        UInt dimRho = M_quadruleRho->nbQuadPt();
        UInt dimTheta = M_quadruleTheta->nbQuadPt();
    
        // Resize loops
        M_Jr.resize( dof );
        M_Jtheta.resize( dof );
        M_Dr.resize( dof );
        M_Dtheta.resize( dof );
        M_Drtheta.resize( dof );
        M_Dthetar.resize( dof );
        M_Jacobian.resize( dof );
            
        for( UInt i( 0 ); i != dof; ++ i )
        {
            M_Jr[i].resize( dimTheta );
            M_Jtheta[i].resize( dimTheta );
            M_Dr[i].resize( dimTheta );
            M_Dtheta[i].resize( dimTheta );
            M_Drtheta[i].resize( dimTheta );
            M_Dthetar[i].resize( dimTheta );
            M_Jacobian[i].resize( dimTheta );
        }
            
        UInt index(0);
        Real vertex(0);
        Real theta(0);

        // Evaluation loops
        for( UInt i = 0; i != dof; ++i )
        {
            if( i < ( dof - 1 ) / 2 + 1 ) // Pressure DOF
            {
                index = 2 * i;
            }
            else // Velocity DOF
            {
                index = 2 * i - dof;
            }
            for( UInt j = 0; j != dimTheta; ++j )
            {
                 vertex = index * h;
                 theta = 2*M_PI * M_quadruleTheta->quadPointCoor( j, 0 );

                 M_Jr[i][j]       = M_functionJr( t, vertex, 0, theta, 0 );
                 M_Jtheta[i][j]   = M_functionJtheta( t, vertex, 0, theta, 0 );
                 M_Dr[i][j]       = M_functionDr( t, vertex, 0, theta, 0 );
                 M_Dtheta[i][j]   = M_functionDtheta( t, vertex, 0, theta, 0 );
                 M_Drtheta[i][j]  = M_functionDrtheta( t, vertex, 0, theta, 0 );
                 M_Dthetar[i][j]  = M_functionDthetar( t, vertex, 0, theta, 0 );
                 M_Jacobian[i][j] = M_functionJacobian( t, vertex, 0, theta, 0 );
            }
        }
    } // end evaluateFtensor
    
    void ReferenceMap::evaluateAxialMap( const Real& h, const UInt& dof, const function_Type& Radius, const function_Type& dRadius )
    {
        const Real t( 0 );
        
        UInt index( 0 );
        
        for( UInt i( 0 ); i != dof; ++ i )
        {
            if( i < ( dof - 1 ) / 2 + 1 )
            {
                index = 2 * i; // Pressure DOF
            }
            else
            {
                index = 2 * i - dof; // Velocity DOF
            }
            (*M_xJr)[i] = M_functionJr( t, index * h , 0, 0, 0 );
            (*M_xJtheta)[i] = M_functionJtheta( t, index * h, 0, 0, 0 );
            (*M_xDr)[i] = M_functionDr( t, index * h, 1, 0, 0 );
            (*M_xDtheta)[i] = M_functionDtheta( t, index * h, 0, 0, 0 );
            (*M_xDrtheta)[i] = M_functionDrtheta( t, index * h, 0, 0, 0 );
            (*M_xDthetar)[i] = M_functionDthetar( t, index * h, 0, 0, 0 );
            (*M_xJacobian)[i] = M_functionJacobian( t, index * h, 0, 0, 0 );
            (*M_xR)[i] = Radius( t, index * h, 0, 0, 0 );
            (*M_xdR)[i] = dRadius( t, index * h, 0, 0, 0 );
        } // end i-for
        // Evaluation loop for the radial component of Dr
        UInt dimRho = M_quadruleRho->nbQuadPt();
        UInt dimTheta = M_quadruleTheta->nbQuadPt();
        
        M_Dr.resize( dimRho, std::vector<Real>( dimTheta ) );
        for( UInt i = 0; i != dimRho; ++i )
        {
            for( UInt j = 0; j != dimTheta; ++j )
                {
                    M_Dr[i][j] = M_quadruleRho->quadPointCoor( i, 0 );
                }
        }
    }// end evaluateAxialMap
} // end namespace
