#include<lifev/himod/util/ConvergenceErrorHandler.hpp>

#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>

#pragma GCC diagnostic ignored "-Wconversion"
#include <lifev/himod/modalbasis/HiModAssembler.hpp>
#pragma GCC diagnostic warning "-Wconversion"

namespace LifeV
{

ConvergenceErrorHandler::ConvergenceErrorHandler (const UInt& m, const UInt& n, const UInt& mstart, const UInt& step)
{

    M_Nit = n;
    M_error.resize (m);
    UInt M = mstart;
    for (UInt i (0); i < m ; ++i)
    {
        M_error[i].resize (n + 1);
        M_error[i][0] = M;
        M *= step;
    }
}


void
ConvergenceErrorHandler::convergeFile (std::string filename)
{
    std::fstream fileOut;
    fileOut.open (filename.c_str(), std::fstream::out);

    if (fileOut.is_open() == false)
    {
        std::cerr << "File not opened" << std::endl;
        exit (1);
    }

    fileOut << M_Nit << std::endl;

    for (UInt i (0); i < M_error.size(); ++i)
    {
        for (UInt j (0); j < M_error[i].size(); ++j)
        {
            fileOut << M_error[i][j] << '\t';
        }

        fileOut <<i+1<< std::endl;
    }

    fileOut.close();
}

/*
    Add error in the correct cell of the matrix ( mode , kind of spatial discretization )
*/
void
ConvergenceErrorHandler::addError (const Real& err, const UInt& mod, const UInt& n)
{
    M_error[mod][n + 1] = err;
}

/*!
    Compute normL2_ues with higher precision (Nel elements in the 1d FESpace)
*/
Real
ConvergenceErrorHandler::computeNormUes (const UInt& Nel, const UInt& M, const function_type& ues, const Real& Lx, const Real& Ly, const Real& Lz,
                                         const std::string& down,  const std::string& up,
                                         const std::string& left,  const std::string& right)
{
    boost::shared_ptr<Epetra_Comm> Comm( new Epetra_SerialComm );

    boost::shared_ptr< mesh_type > fullMeshPtr_aus( new mesh_type );

    regularMesh1D ( *fullMeshPtr_aus, 0, Nel , false, Lx, 0.0 );

    boost::shared_ptr<FESpace< mesh_type, MapEpetra > > uSpace_aus
    ( new FESpace< mesh_type, MapEpetra > ( fullMeshPtr_aus, "P1", 1, Comm ) );

    boost::shared_ptr<ModalSpaceRectangular> MB_aus( new ModalSpaceRectangular( Ly, Lz, M, &quadRuleLobSeg64pt, &quadRuleLobSeg64pt ) );
    MB_aus->addSliceBCY( left, right );
    MB_aus->addSliceBCZ( down, up );
    MB_aus->evaluateBasis();

    HiModAssembler<mesh_type , matrix_type  , vector_type, Rectangular > HM_aus( uSpace_aus, MB_aus, Comm );

    UInt nquadY = MB_aus->qrY().nbQuadPt();
    UInt nquadZ = MB_aus->qrZ().nbQuadPt();

    MapEpetra Map_3D_aus( nquadY * nquadZ * (Nel + 1), Comm );

    vector_type ues_3D_aus( Map_3D_aus, Unique );

    ues_3D_aus = HM_aus.evaluateBase3DGrid( ues );

    Real norm_ues;

    norm_ues = HM_aus.normL2 (ues_3D_aus);

    return norm_ues;
}

}//End lifev namespace
