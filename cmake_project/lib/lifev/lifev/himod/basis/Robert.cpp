#include <lifev/himod/basis/Robert.hpp>
#include <cmath>

namespace LifeV
{
    void
    Robert::
    evaluateBasis(     MBMatrix_type& phirho,
                       MBMatrix_type& dphirho,
                       MBMatrix_type& phitheta,
                       MBMatrix_type& dphitheta,
                       const QuadratureRule* quadrulerho,     
                       const QuadratureRule* quadruletheta )
    {
        int trigonometricBasis(0);
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
        MapEpetra Map( M_mtot, Comm );

        //Using the map it is possible to define the system matrix
        M_massMatrix = static_cast< boost::shared_ptr<matrix_Type>  > ( new matrix_Type( Map ) );

        std::vector<UInt> block_row( 1, M_mtot ); // one block with mtot cells
        std::vector<UInt> block_col( 1, M_mtot );

        M_massMatrix->setBlockStructure( block_row, block_col );
        *M_massMatrix *= 0.0;
        
        M_eigenvalues.resize( M_mtot );
        phirho.resize( M_mtot );
        dphirho.resize( M_mtot );
        phitheta.resize( M_mtot );
        dphitheta.resize( M_mtot );

        for ( UInt i = 0; i != M_mtot; ++i )
        {
            phirho[i].resize( quadrulerho->nbQuadPt() );
            dphirho[i].resize( quadrulerho->nbQuadPt() );
            phitheta[i].resize( quadruletheta->nbQuadPt() );
            dphitheta[i].resize( quadruletheta->nbQuadPt() );

        }
        
        MBVector_type normR( M_mtot, 0 );
        
        UInt p( 0 );
        
        for ( UInt pt = 0; pt != mtheta1D; ++pt )
        {
            for ( UInt pr = 0; pr != mr1D; ++pr )
            {
                M_eigenvalues[pt*mr1D+pr].lambda = pt*mr1D+pr;
                M_eigenvalues[pt*mr1D+pr].order = pt;
                M_eigenvalues[pt*mr1D+pr].index = pr;
                    
                // We need to distinguish in order to fulfil the parity theorem
                if( remainder( pt, 2 ) == 0 )
                {
                    p = 2 * pr;
                }
                else
                {
                    p = 2 * pr + 1;
                }
                // compute the norm
                for( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                {
                    Real x = quadrulerho->quadPointCoor( n, 0 );
                    normR[mr1D*pt + pr] += ( pow( quadruletheta->quadPointCoor( n, 0 ), 2 * pt ) *
                                            cos( p * acos( x ) ) * cos( p * acos( x ) ) 
                                            ) *
                                            M_Rho * quadrulerho->quadPointCoor( n, 0 ) * 
                                            M_Rho * quadrulerho->weight( n );
                } // end norm-for
                for ( UInt n = 0; n != quadrulerho->nbQuadPt(); ++n )
                {
                        Real x = quadrulerho->quadPointCoor( n, 0 );
                    // NOTA BENE: LE FUNZIONI SONO MOLTIPLICATE PER UN M_Rho IN PIÙ PER COMPENSARE IL normrho CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                        phirho[mr1D*pt + pr][n] = M_Rho *
                                                    pow( x, pt ) *
                                                    cos( p * acos( x ) ) / sqrt( normR[mr1D*pt + pr] );
                                                    
                        dphirho[mr1D*pt + pr][n] = M_Rho * ( 
                                            pt * pow( x, pt - 1 ) *
                                            cos( p * acos( x ) ) +
                                
                                            pow( x, pt ) *
                                            p * sin( p * acos( x ) ) /
                                            sqrt( 1 - x * x )
                                            ) / sqrt( normR[mr1D*pt + pr] );
                } // end n-for
            } // end pr-for
        }//end pt-for
        
        gramSchmidtOrthogonalization( phirho, dphirho,
                                        quadrulerho,
                                        mr1D,
                                        mtheta1D );
                                        
        checkOrthogonality( phirho, quadrulerho, mr1D, mtheta1D );
        
        Real normrho( 1. / M_Rho );
        Real normtheta( 1. / sqrt( 2 * M_PI ) );    
        Real mass( 0 );    
        
        switch( trigonometricBasis )
        {
            case 0:

                    for ( UInt pt = 0; pt != mtheta1D; ++pt )
                    {
                        Real correction( 1 );
                        if( pt==0 )
                        {
                            correction = 1. / sqrt( 2 );
                        }
                            
                        for ( UInt pr = 0; pr != mr1D; ++pr )
                        {
                            for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                            {
                                // NOTA BENE: M_normtheta È sqrt(2) PER COMPENSARE IL normtheta CHE COMPARE NEGLI INTEGRALI DI ModalSpace e HMA.
                                // la norma di sin + cos, infatti, è sqrt(2*pi), mentre quella di sin o di cos è solo sqrt(pi).
                                phitheta[mr1D*pt + pr][n]  = M_normTheta * correction * ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );

                                dphitheta[mr1D*pt + pr][n] = M_normTheta * pt * correction *
                                                                  ( - sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                            } // end n-for
                        } // end pr-for
                    } // end pt-for
                    
                    mass = 0;
                    
                    // Fill in mass matrix
                    // Cycle on rows
                    for ( UInt row = 0; row != M_mtot; ++row )
                    {
                        // Cycle on columns
                        for ( UInt col = 0; col != M_mtot; ++col )
                        {
                            // Compute integral
                            for ( UInt nt = 0; nt != quadruletheta->nbQuadPt(); ++nt )
                            {
                                for ( UInt nr = 0; nr != quadrulerho->nbQuadPt(); ++nr )
                                {
                                    mass += phirho[row][nr] * normrho *
                                            phitheta[row][nt] * normtheta *
                                            phirho[col][nr] * normrho *
                                            phitheta[col][nt] * normtheta *
                                            M_Rho * quadrulerho->quadPointCoor( nr, 0 ) *
                                            M_Rho * quadrulerho->weight( nr ) *
                                            M_Theta * quadruletheta->weight( nt );
                                }
                            }
                            M_massMatrix->addToCoefficient( row, col, mass );
                            mass = 0;
                        }
                    }

                    break;                            
                            
            case 1:
                    for ( UInt pt = 0; pt != mtheta1D; ++pt )
                    {
                        for ( UInt pr = 0; pr != mr1D; ++pr )
                        {
                            for ( UInt n = 0; n != quadruletheta->nbQuadPt(); ++n )
                            {
                                phitheta[mr1D*pt + pr][n]  = M_normTheta * ( sin ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
        
                                dphitheta[mr1D*pt + pr][n] = M_normTheta * pt *
                                                      ( cos ( pt * 2 * M_PI * quadruletheta->quadPointCoor( n, 0 ) ) );
                            } // end n-for
                        } // end pr-for
                    } // end pt-for
                            
                    mass = 0;
                    
                    // Fill in mass matrix
                    // Cycle on rows
                    for ( UInt row = 0; row != M_mtot; ++row )
                    {
                        // Cycle on columns on double indices
                        for ( UInt col = 0; col != M_mtot; ++col )
                        {
                            // Compute integral
                            for ( UInt nt = 0; nt != quadruletheta->nbQuadPt(); ++nt )
                            {
                                for ( UInt nr = 0; nr != quadrulerho->nbQuadPt(); ++nr )
                                {
                                    mass += phirho[row][nr] * normrho *
                                            phitheta[row][nt] * normtheta *
                                            phirho[col][nr] * normrho *
                                            phitheta[col][nt] * normtheta *
                                            M_Rho * quadrulerho->quadPointCoor( nr, 0 ) *
                                            M_Rho * quadrulerho->weight( nr ) *
                                            M_Theta * quadruletheta->weight( nt );
                                } // end nr-for
                            } // end nt-for
                            M_massMatrix->addToCoefficient(row, col, mass );
                            mass = 0;
                        } // end col-for
                    } // end row-for
                            
                    break;

            default:
                    break;
        }
            
        M_massMatrix->globalAssemble();
        
        return;
    }
    
    /*Real
    Robert::
    evalSinglePoint( const UInt index, const MBVector_type& p,
                       const QuadratureRule* quadrulerho, const int trigonometricBasis ) const
    {
        UInt mr1D( 0 ), mtheta1D( 0 );
        rectangularTruncation( M_mtot, mr1D, mtheta1D );
        
        UInt rIndex = remainder( index, mr1D );
        UInt thetaIndex = static_cast<UInt> ( index / mr1D );
        
        Real rhoh = p[0];
        Real thetah = p[1];
        
        switch( trigonometricBasis )
        {
            case 0:
                    return pow( rhoh, thetaIndex ) * cos( rIndex * acos( rhoh ) ) * cos( thetaIndex * 2 * M_PI * thetah );
            case 1:
                    return pow( rhoh, thetaIndex ) * cos( rIndex * acos( rhoh ) ) * sin( thetaIndex * 2 * M_PI * thetah );
            default:
                    return 0;
        }
    }*/
} //End LifeV namespace

