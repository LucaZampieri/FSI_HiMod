#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>

#include <lifev/core/LifeV.hpp>

#include <lifev/himod/util/GeneralTest.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/himod/modalbasis/ModalSpaceRectangular.hpp>

#include <lifev/himod/modalbasis/HiModAssembler.hpp>
#include <lifev/himod/util/ConvergenceErrorHandler.hpp>
#include <lifev/himod/util/CaseTest.hpp>
#include <lifev/himod/util/CaseTestStruct.hpp>


namespace LifeV
{

GeneralTest::GeneralTest (const QuadratureRule* quady, const QuadratureRule* quadz) : M_quady (quady), M_quadz (quadz)
{
    M_CaseTestSwitch["Camini"] = Camini;
    M_CaseTestSwitch["DDDD"] = DDDD;
    M_CaseTestSwitch["DDDD_Adv"] = DDDD_Adv;
    M_CaseTestSwitch["DDDD_ADR"] = DDDD_ADR;
    M_CaseTestSwitch["RRRR"] = RRRR;
    M_CaseTestSwitch["DRDR"] = DRDR;
    M_CaseTestSwitch["RR2D"] = RR2D;
    M_CaseTestSwitch["BDRR"] = BDRR;
    M_CaseTestSwitch["DD2D"] = DD2D;
    M_CaseTestSwitch["BDDD"] = BDDD;
}

void
GeneralTest::setDatafile (GetPot const& DataFile)
{
    M_DataFile = DataFile;

    CaseTestStruct CurrentCase;
    std::string CaseName = M_DataFile ("himod/casetest", "");
    switch (M_CaseTestSwitch[CaseName])
    {
        case Camini:
            CurrentCase = CaseTest::Camini;
            break;
        case DDDD:
            CurrentCase = CaseTest::DDDD;
            break;
        case DDDD_Adv:
            CurrentCase = CaseTest::DDDD_Adv;
            break;
        case DDDD_ADR:
            CurrentCase = CaseTest::DDDD_ADR;
            break;
        case RRRR:
            CurrentCase = CaseTest::RRRR;
            break;
        case DRDR:
            CurrentCase = CaseTest::DRDR;
            break;
        case RR2D:
            CurrentCase = CaseTest::RR2D;
            break;
        case BDRR:
            CurrentCase = CaseTest::BDRR;
            break;
        case DD2D:
            CurrentCase = CaseTest::DD2D;
            break;
        case BDDD:
            CurrentCase = CaseTest::BDDD;
            break;
        default:
            std::cout << "This case test does not exist, please add it in CaseTest.hpp" << std::endl;
            break;
    }

    setExport (DataFile ("himod/export/nx", 30), DataFile ("himod/export/ny", 30), DataFile ("himod/export/nz", 30) );
    setData (DataFile ("himod/m", 30), DataFile ("himod/Nel", 30) );

    setCoefficients (CurrentCase.mu, CurrentCase.beta, CurrentCase.sigma,  CurrentCase.chiy, CurrentCase.chiz);
    setBC (CurrentCase.UP, CurrentCase.DOWN, CurrentCase.LEFT, CurrentCase.RIGHT);
    setFunctions (CurrentCase.f, CurrentCase.g);
    setDomain (CurrentCase.Lx, CurrentCase.Ly, CurrentCase.Lz);
    setCaseName (CurrentCase.CaseName);
}

void
GeneralTest::setExport (UInt const& nx, UInt const& ny, UInt const& nz)
{
    M_nx_grid =   nx;
    M_ny_grid =   ny;
    M_nz_grid =   nz;
}


void
GeneralTest::setCoefficients (Real const& mu, TreDvector_type const& beta, Real const& sigma,  Real const& chiy, Real const& chiz)
{
    M_mu        =   mu;
    M_beta      =   beta;
    M_sigma     =   sigma;
    M_chiy       =   chiy;
    M_chiz       =   chiz;
}

void
GeneralTest::setBC (std::string bcUp, std::string bcDown, std::string bcLeft, std::string bcRight)
{
    M_bcDown        =   bcDown;
    M_bcUp          =   bcUp;
    M_bcLeft        =   bcLeft;
    M_bcRight       =   bcRight;
}


void
GeneralTest::setFunctions (function_type const& ForceTerm, function_type const& DirichletInflow)
{
    M_ForceTerm         =   ForceTerm;
    M_DirichletInflow   =   DirichletInflow;
}

void
GeneralTest::setData (   UInt const& m,
                         UInt const& Nel)
{
    M_m =   m;
    M_Nel =   Nel;
}

void
GeneralTest::setDomain (Real const& lx, Real const& ly, Real const& lz)
{
    M_lx    =   lx;
    M_ly    =   ly;
    M_lz    =   lz;
}

void
GeneralTest::run()
{

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);


    //! creating the 1D-mesh
    boost::shared_ptr< mesh_type > MeshPtr (new mesh_type);

    regularMesh1D ( *MeshPtr, 0, M_Nel, false, M_lx , 0.0);

    //! creating the 1D-fespace
    boost::shared_ptr<FESpace< mesh_type, MapEpetra > > uSpace ( new FESpace< mesh_type, MapEpetra > (MeshPtr, "P1", 1, Comm) );

    //Checking the number of DoF in the fespace

    DOF DataFESpace (uSpace->dof() );
    UInt numdof = DataFESpace.numTotalDof();
    //std::cout << "Degrees of freedom FEM = " << numdof << std::endl;


    //! Creation the Modal basis
    boost::shared_ptr<ModalSpaceRectangular> MB (new ModalSpaceRectangular (M_ly, M_lz, M_m, M_quady, M_quadz) );
    MB->addSliceBCY (M_bcLeft, M_bcRight , M_mu , M_chiy);
    MB->addSliceBCZ (M_bcDown, M_bcUp    , M_mu , M_chiz);
    MB->evaluateBasis();

    //! Creation of the HiMod space
    HiModAssembler< mesh_type , matrix_type , vector_type, Rectangular > HM (uSpace, MB, Comm);

    //! We need a map in order to define the system matrix and the other vectors
    MapEpetra Map (M_m * numdof, Comm);

    //! Declaration of the matrix pointer
    boost::shared_ptr<matrix_type> systemMatrix (new matrix_type ( Map ) );

    //! Setting of block structure of the matrix
    std::vector<UInt> block_row (MB->mtot(), DataFESpace.numTotalDof() );
    std::vector<UInt> block_col (MB->mtot(), DataFESpace.numTotalDof() );
    systemMatrix->setBlockStructure (block_row, block_col);

    //! Initialization of the matrix with all zeros
    *systemMatrix *= 0.0;

    //!  Now we actually build the matrix
    HM.addADRProblem (systemMatrix, M_mu, M_beta, M_sigma);
    //! We project the force term on the modal basis and then interpolate the fourier coefficients on the FEspace.
    boost::shared_ptr<vector_type> f_interpolated (new vector_type ( Map, Repeated ) );
    *f_interpolated *= 0.0;
    f_interpolated->setBlockStructure (block_row);
    HM.interpolate (M_ForceTerm, f_interpolated);

    //! Assembling of the rhs with the f_interpolated
    boost::shared_ptr<vector_type> rhs (new vector_type ( Map, Repeated ) );
    *rhs *= 0.0;
    rhs->setBlockStructure (block_row);
    HM.addrhs (rhs, f_interpolated);

    //! Now we add boundary condition on inflow
    HM.addDirichletBC_In (systemMatrix, rhs, M_DirichletInflow );

    // We can finally close the rhs and the matrix
    rhs->globalAssemble();
    systemMatrix->globalAssemble();

    //! We start now the solution phase by setting the preconditioner.
    boost::shared_ptr<prec_Type> precPtr (new prec_Type);
    precPtr->setDataFromGetPot (M_DataFile, "prec");

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    LinearSolver linearSolver;
    linearSolver.setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList );
    linearSolver.setPreconditioner ( precPtr );

    //! the actual solution of the system
    boost::shared_ptr<vector_type> solution ( new vector_type ( Map, Unique ) );

    linearSolver.setOperator ( systemMatrix );
    linearSolver.setRightHandSide ( rhs );

    linearSolver.solve ( solution );

    HM.exporterStructuredVTK (  M_nx_grid, M_ny_grid, M_nz_grid, solution,
                                M_DataFile, M_caseName, M_DataFile ("himod/export/solution", "solution") );
    HM.exporterFunctionVTK (        M_nx_grid, M_ny_grid, M_nz_grid, M_ForceTerm,
                                    M_DataFile, M_caseName, M_DataFile ("himod/export/forceterm", "forceterm") );

    std::cout << " ###########################################" << std::endl;
    std::cout << " Dimension of the modal basis = " << M_m << std::endl;
    std::cout << " Number of element in x-direction = " << M_Nel << std::endl;
    std::cout << " ###########################################" << std::endl;

}



}//End of lifeV namespace

