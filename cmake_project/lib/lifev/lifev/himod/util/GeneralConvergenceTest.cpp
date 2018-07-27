#include <Epetra_ConfigDefs.h>
#include <Epetra_SerialComm.h>

#include <lifev/core/LifeV.hpp>

#include <lifev/himod/util/GeneralConvergenceTest.hpp>

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

GeneralConvergenceTest::GeneralConvergenceTest (const QuadratureRule* quady, const QuadratureRule* quadz) : M_quady (quady), M_quadz (quadz)
{
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
GeneralConvergenceTest::setDatafile (GetPot const& DataFile)
{
    M_DataFile = DataFile;
    CaseTestStruct CurrentCase;
    std::string CaseName = M_DataFile ("himod/casetest", "");

    switch (M_CaseTestSwitch[CaseName])
    {
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
            std::cout << "This case test does not exist, please add it" << std::endl;
            break;
    }

    setCoefficients (CurrentCase.mu, CurrentCase.beta, CurrentCase.sigma,  CurrentCase.chiy, CurrentCase.chiz);
    setBC (CurrentCase.UP, CurrentCase.DOWN, CurrentCase.LEFT, CurrentCase.RIGHT);
    setFunctions (CurrentCase.ues, CurrentCase.f, CurrentCase.g);
    setDomain (CurrentCase.Lx, CurrentCase.Ly, CurrentCase.Lz);
    setCaseName (CurrentCase.CaseName);

	setLoopData (	M_DataFile("himod/LoopData/m_start",2),
                    M_DataFile("himod/LoopData/Nel_start",10),
				M_DataFile("himod/LoopData/n_space_it",4),
				M_DataFile("himod/LoopData/n_modal_it",5),
				M_DataFile("himod/LoopData/m_step",2),
				M_DataFile("himod/LoopData/Nel_step",2));

	setExport(	M_DataFile("himod/export/makeoutput",false),
				M_DataFile("himod/export/makegraph",true));

	setChrono(	M_DataFile("himod/chrono",false));
}


void
GeneralConvergenceTest::setCoefficients (Real const& mu, TreDvector_type const& beta, Real const& sigma,  Real const& chiy, Real const& chiz)
{
    M_mu        =   mu;
    M_beta      =   beta;
    M_sigma     =   sigma;
    M_chiy       =   chiy;
    M_chiz       =   chiz;
}

void
GeneralConvergenceTest::setBC (std::string bcUp, std::string bcDown, std::string bcLeft, std::string bcRight)
{
    M_bcDown        =   bcDown;
    M_bcUp          =   bcUp;
    M_bcLeft        =   bcLeft;
    M_bcRight       =   bcRight;
}


void
GeneralConvergenceTest::setFunctions (function_type const& ExactSolution, function_type const& ForceTerm, function_type const& DirichletInflow)
{
    M_ExactSolution     =   ExactSolution;
    M_ForceTerm         =   ForceTerm;
    M_DirichletInflow   =   DirichletInflow;
}
void
GeneralConvergenceTest::setLoopData (   UInt const& m_start,
                                        UInt const& Nel_start,
                                        UInt const& n_space_it,
                                        UInt const& n_modal_it,
                                        UInt const& m_step,
                                        UInt const& Nel_step)
{
    M_m_start       =   m_start;
    M_m_step        =   m_step;
    M_Nel_start =   Nel_start;
    M_Nel_step  =   Nel_step;
    M_n_space_it    =   n_space_it;
    M_n_modal_it    =   n_modal_it;
}

void
GeneralConvergenceTest::setDomain (Real const& lx, Real const& ly, Real const& lz)
{
    M_lx    =   lx;
    M_ly    =   ly;
    M_lz    =   lz;
}

void GeneralConvergenceTest::setCaseName (const std::string& CN)
{
    M_caseName = CN;
}

void GeneralConvergenceTest::setExport (bool const& exp_VTK, bool const& exp_graph)
{
    M_export = exp_VTK;
    M_graph  = exp_graph;
}

void GeneralConvergenceTest::setChrono(bool const& chrono)
	{
		M_chrono = chrono;
	}



//----------------------------- RUN method ------------------------------------------------------------

void
GeneralConvergenceTest::run()
{

    UInt Nelem  = M_Nel_start;
    UInt m    = M_m_start;
    UInt Nstep  = M_Nel_step;
    UInt mstep  = M_m_step;

    Real norm_ues;
    norm_ues = ConvergenceErrorHandler::computeNormUes (100 , 50, M_ExactSolution , M_lx, M_ly, M_lz, M_bcDown, M_bcUp, M_bcLeft, M_bcRight );

    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);

    ConvergenceErrorHandler error (M_n_modal_it, M_n_space_it, m, mstep);

	LifeChrono HiModchrono;
	ofstream tempi;
	UInt iter_time(1);
	if(M_chrono) tempi.open("AssembleTime.dat");
    
	for (UInt it_space (0); it_space < M_n_space_it; ++it_space)
    {


        //! creating the 1D-mesh
        boost::shared_ptr< mesh_type > MeshPtr (new mesh_type);
        regularMesh1D ( *MeshPtr, 0, Nelem, false, M_lx , 0.0);

        //! creating the 1D-fespace
        boost::shared_ptr<FESpace< mesh_type, MapEpetra > > uSpace ( new FESpace< mesh_type, MapEpetra > (MeshPtr, "P1", 1, Comm) );

        //Checking the number of DoF in the fespace

        DOF DataFESpace (uSpace->dof() );
        UInt numdof = DataFESpace.numTotalDof();
        //std::cout << "Degrees of freedom FEM = " << numdof << std::endl;

        m = M_m_start;
        for (UInt it_modal (0); it_modal < M_n_modal_it; ++it_modal)
        {

		if(M_chrono) HiModchrono.start();
          
		  //! Creation the Modal basis
            boost::shared_ptr<ModalSpaceRectangular> MB (new ModalSpaceRectangular (M_ly, M_lz, m, M_quady, M_quadz) );
            MB->addSliceBCY (M_bcLeft, M_bcRight , M_mu , M_chiy);
            MB->addSliceBCZ (M_bcDown, M_bcUp    , M_mu , M_chiz);
            MB->evaluateBasis();
            //MB->showMe();

            //! Creation of the HiMod space
            HiModAssembler< mesh_type , matrix_type , vector_type, Rectangular > HM (uSpace, MB, Comm);

            if (M_export && ( (it_modal + it_space) == 0) )
            {
                HM.exporterFunctionVTK (    M_DataFile ("himod/export/nx", 10),
                                            M_DataFile ("himod/export/ny", 10),
                                            M_DataFile ("himod/export/nz", 10),
                                            M_ExactSolution, M_DataFile,
                                            M_caseName, M_DataFile ("himod/export/exactsolution", "UES") );
            }
            //! We need a map in order to define the system matrix and the other vectors
            MapEpetra Map (m * numdof, Comm);

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

		if(M_chrono) 
			{
				HiModchrono.stop();
				tempi<<iter_time<<'\t'<<m<<'\t'<<HiModchrono.diff()<<std::endl;
				iter_time++;
			}

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

            if (M_export && it_space == M_n_space_it - 1)
            {
                std::string name = "m=";
                std::ostringstream aus;
                aus << m;
                name += aus.str();
                HM.exporterStructuredVTK (  M_DataFile ("himod/export/nx", 10),
                                            M_DataFile ("himod/export/ny", 10),
                                            M_DataFile ("himod/export/nz", 10),
                                            solution, M_DataFile, M_caseName, name);
            }

            //Compute ERRORS and try a test.
            UInt nquadY = MB->qrY().nbQuadPt();
            UInt nquadZ = MB->qrZ().nbQuadPt();
            MapEpetra Map_3D (nquadY * nquadZ * numdof, Comm);

            vector_type solution_3D (Map_3D, Unique);
            solution_3D = HM.evaluateBase3DGrid (*solution);

            vector_type ues_3D (Map_3D, Unique);
            ues_3D = HM.evaluateBase3DGrid (M_ExactSolution);

            vector_type err_3D (Map_3D, Unique);
            err_3D = solution_3D;
            err_3D += ues_3D * (-1.);

            Real normL2 = HM.normL2 (err_3D);

            std::cout << " ###########################################" << std::endl;

            std::cout << " Dimension of the modal basis = " << m << std::endl;
            std::cout << " Number of element in x-direction = " << Nelem << std::endl;
            std::cout << " Error L2, normalized: " << normL2 / norm_ues << std::endl;
            std::cout << " ###########################################" << std::endl;

            error.addError (normL2 / norm_ues, it_modal, it_space);
            m = m * mstep;
        }
        Nelem *= Nstep;
    }

#pragma GCC diagnostic ignored "-Wunused-result"
    system ("rm -f Errors.dat");
    error.convergeFile ("Errors.dat");
    system("gnuplot --persist graph.gplot");
#pragma GCC diagnostic warning "-Wunused-result"

}



}//End of lifeV namespace

