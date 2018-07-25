//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief ML preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-2006
 */

#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerML::PreconditionerML ( boost::shared_ptr<Epetra_Comm> comm ) :
    super(),
    M_comm ( comm ),
    M_operator(),
    M_preconditioner(),
    M_analyze (false),
    M_visualizationDataAvailable (false)
{

}

PreconditionerML::~PreconditionerML()
{
    M_preconditioner.reset();
    M_operator.reset();
}


// ===================================================
// Methods
// ===================================================
Int
PreconditionerML::buildPreconditioner ( operator_type& matrix )
{

    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_operator.
    //to avoid the risk of dandling pointers always deallocate M_preconditioner first and then M_operator
    M_preconditioner.reset();
    M_operator = matrix;

    M_precType = M_list.get ( "prec type", "undefined??" );
    M_precType += "_ML";

    // <one-level-postsmoothing> / <two-level-additive>
    // <two-level-hybrid> / <two-level-hybrid2>

    M_preconditioner.reset ( new prec_raw_type ( * (M_operator->matrixPtr() ), this->parametersList(), true ) );

    if ( M_analyze )
    {
        ML_Epetra::MultiLevelPreconditioner* prec;
        prec = dynamic_cast<ML_Epetra::MultiLevelPreconditioner*> ( M_preconditioner.get() );
        Int NumPreCycles = 5;
        Int NumPostCycles = 1;
        Int NumMLCycles = 10;
        prec->AnalyzeHierarchy ( true, NumPreCycles, NumPostCycles, NumMLCycles );
    }

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
PreconditionerML::resetPreconditioner()
{
    //the Trilinos::MultiLevelPreconditioner unsafely access to the area of memory co-owned by M_operator.
    //to avoid the risk of dandling pointers always deallocate M_preconditioner first and then M_operator

    M_preconditioner.reset();
    M_operator.reset();

    this->M_preconditionerCreated = false;
}

void
PreconditionerML::createMLList ( list_Type& list,
                                 const GetPot&       dataFile,
                                 const std::string&  section,
                                 const std::string&  subSection,
                                 const bool&         verbose )
{
    list.setName ( "ML paramters list" );

    std::string defList = dataFile ( (section + "/" + subSection + "/default_parameter_list").data(), "SA" );
    if ( defList != "none" )
    {
        ML_Epetra::SetDefaults ( defList, list );
    }

    bool found;

    bool MLPrintParameterList = dataFile ( (section + "/displayList").data(), false, found );

    Int MLOutput             = dataFile ( (section + "/" + subSection + "/MLOuput").data(), 0, found);
    if ( found )
    {
        list.set ( "ML output", MLOutput );
    }

    Int printUnused          = dataFile ( (section + "/" + subSection + "/print_unused").data(), -2, found);
    if ( found )
    {
        list.set ( "print unused", printUnused );
    }

    Int PDEEquations         = dataFile ( (section + "/" + subSection + "/pde_equations").data(), 1, found);
    if ( found )
    {
        list.set ( "PDE equations", PDEEquations );
    }

    Int CycleApplications    = dataFile ( (section + "/" + subSection + "/cycle_applications").data(), 1, found);
    if ( found )
    {
        list.set ( "cycle applications", CycleApplications );
    }

    Int MaxLevels            = dataFile ( (section + "/" + subSection + "/max_levels").data(), 2, found);
    if ( found )
    {
        list.set ( "max levels", MaxLevels );
    }

    std::string IncOrDec     = dataFile ( (section + "/" + subSection + "/inc_or_dec").data(), "increasing", found);
    if ( found )
    {
        list.set ( "increasing or decreasing", IncOrDec );
    }

    std::string precType     = dataFile ( (section + "/" + subSection + "/prec_type").data(), "MGV", found);
    if ( found )
    {
        list.set ( "prec type", precType );
    }

    // Int NumProjectedModes    = dataFile( (section + "/" + subSection + "/number_of_prejected_modes").data(), 0 );
    // if ( found ) list.set( "ML print parameter list", MLPrintParameterList );

    std::string eigenAnalysisType = dataFile ( (section + "/" + subSection + "/eigne-analysis/type").data(), "cg", found );
    if ( found )
    {
        list.set ( "eigen-analysis: type", eigenAnalysisType );
    }

    Int eigenAnalysisIterations   = dataFile ( (section + "/" + subSection + "/eigne-analysis/iterations").data(), 10, found );
    if ( found )
    {
        list.set ( "eigen-analysis: iterations", eigenAnalysisIterations );
    }

    // Aggregation options

    std::string AggregationType                  = dataFile ( (section + "/" + subSection + "/aggregation/type").data(), "Uncoupled", found );
    if ( found )
    {
        list.set ( "aggregation: type", AggregationType );
    }

    bool AggregationBlockScaling          = dataFile ( (section + "/" + subSection + "/aggregation/block_scaling").data(), false, found );
    if ( found )
    {
        list.set ( "aggregation: block scaling", AggregationBlockScaling );
    }

    Real        AggregationThreshold             = dataFile ( (section + "/" + subSection + "/aggregation/threshold").data(), 0.0 , found );
    if ( found )
    {
        list.set ( "aggregation: threshold", AggregationThreshold );
    }

    Real        AggregationDampingFactor         = dataFile ( (section + "/" + subSection + "/aggregation/damping_factor").data(), 4. / 3. , found );
    if ( found )
    {
        list.set ( "aggregation: damping factor", AggregationDampingFactor );
    }

    Int AggregationSmoothingSweeps               = dataFile ( (section + "/" + subSection + "/aggregation/smoothing_sweeps").data(), 1, found );
    if ( found )
    {
        list.set ( "aggregation: smoothing sweeps", AggregationSmoothingSweeps );
    }

    Int AggregationGlobalAggregates              = dataFile ( (section + "/" + subSection + "/aggregation/global_aggregates").data(), 1, found );
    if ( found )
    {
        list.set ( "aggregation: global aggregates", AggregationGlobalAggregates );
    }

    Int AggregationLocalAggregates               = dataFile ( (section + "/" + subSection + "/aggregation/local_aggregates").data(), 1, found );
    if ( found )
    {
        list.set ( "aggregation: local aggregates", AggregationLocalAggregates );
    }

    Int AggregationNodesPerAggregate             = dataFile ( (section + "/" + subSection + "/aggregation/nodes_per_aggregate").data(), 1, found );
    if ( found )
    {
        list.set ( "aggregation: nodes per aggregate",  AggregationNodesPerAggregate );
    }

    Int AggregationNextLevelAggregatesPerProcess = dataFile ( (section + "/" + subSection + "/aggregation/next-level_aggregates_per_process").data(), 128, found );
    if ( found )
    {
        list.set ( "aggregation: next level aggregates per process", AggregationNextLevelAggregatesPerProcess );
    }

    bool AggregationUseTentativeRestriction      = dataFile ( (section + "/" + subSection + "/aggregation/tentative_restriction").data(), false, found );
    if ( found )
    {
        list.set ( "aggregation: use tentative restriction", AggregationUseTentativeRestriction );
    }

    bool AggregationSymmetrize                   = dataFile ( (section + "/" + subSection + "/aggregation/symmetrize").data(), false, found );
    if ( found )
    {
        list.set ( "aggregation: symmetrize", AggregationSymmetrize );
    }

    Int AggregationNumLevelTypes                 = dataFile ( (section + "/" + subSection + "/aggregation/level_type").data(), 0, found );
    if ( found )
    {
        for (Int i (0); i < AggregationNumLevelTypes; ++i)
        {
            std::string levelIndex             = dataFile ( (section + "/" + subSection + "/aggregation/level_type").data(), "0", 2 * i + 1 );
            std::string levelParamValue = dataFile ( (section + "/" + subSection + "/aggregation/level_type").data(), "METIS", 2 * i + 2 );
            list.set ( ("aggregation: type (level " + levelIndex + ")").data(), levelParamValue );
        }
    }

    bool   EnergyMinimizationEnable              = dataFile ( (section + "/" + subSection + "/energy_minimization/enable").data(), false, found );
    if ( found )
    {
        list.set ( "energy minimization: enable", EnergyMinimizationEnable );
    }

    Int    EnergyMinimizationType                = dataFile ( (section + "/" + subSection + "/energy_minimization/type").data(), 2, found );
    if ( found )
    {
        list.set ( "energy minimization: type", EnergyMinimizationType );
    }

    Real   EnergyMinimizationDropTol             = dataFile ( (section + "/" + subSection + "/energy_minimization/droptol").data(), 0., found );
    if ( found )
    {
        list.set ( "energy minimization: droptol", EnergyMinimizationDropTol );
    }

    bool   EnergyMinimizationCheap               = dataFile ( (section + "/" + subSection + "/energy_minimization/cheap").data(), false, found );
    if ( found )
    {
        list.set ( "energy minimization: cheap", EnergyMinimizationCheap );
    }

    // Smoothing parameters

    std::string SmootherType                = dataFile ( (section + "/" + subSection + "/smoother/type").data(), "IFPACK", found );
    if ( found )
    {
        list.set ( "smoother: type", SmootherType );
    }

    Int SmootherSweeps                      = dataFile ( (section + "/" + subSection + "/smoother/sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "smoother: sweeps", SmootherSweeps );
    }

    Real   SmootherDampingFactor            = dataFile ( (section + "/" + subSection + "/smoother/damping_factor").data(), 1.0, found );
    if ( found )
    {
        list.set ( "smoother: damping factor", SmootherDampingFactor );
    }

    std::string SmootherPreOrPost           = dataFile ( (section + "/" + subSection + "/smoother/pre_or_post").data(), "both", found );
    if ( found )
    {
        list.set ( "smoother: pre or post", SmootherPreOrPost );
    }

    Real   SmootherChebyshevAlpha           = dataFile ( (section + "/" + subSection + "/smoother/Chebyshev_alpha").data(), 20., found );
    if ( found )
    {
        list.set ( "smoother: Chebyshev alpha", SmootherChebyshevAlpha );
    }

    bool SmootherHiptmairEfficientSymmetric = dataFile ( (section + "/" + subSection + "/smoother/Hiptmair_efficient_symmetric").data(), true, found );
    if ( found )
    {
        list.set ( "smoother: Hiptmair efficient symmetric", SmootherHiptmairEfficientSymmetric );
    }

    std::string SmootherIfpackType           = dataFile ( (section + "/" + subSection + "/smoother/ifpack_type").data(), "ILU", found );
    if ( found )
    {
        list.set ( "smoother: ifpack type", SmootherIfpackType );
    }

    Int SmootherIfpackOverlap                = dataFile ( (section + "/" + subSection + "/smoother/ifpack_overlap").data(), 1, found );
    if ( found )
    {
        list.set ( "smoother: ifpack overlap", SmootherIfpackOverlap );
    }

    Int SmootherNumLevelTypes               = dataFile ( (section + "/" + subSection + "/smoother/level_type").data(), 0, found );
    if ( found )
    {
        for (Int i (0); i < SmootherNumLevelTypes; ++i)
        {
            std::string levelIndex          = dataFile ( (section + "/" + subSection + "/smoother/level_type").data(), "0", 2 * i + 1 );
            std::string levelParamValue     = dataFile ( (section + "/" + subSection + "/smoother/level_type").data(), "IFPACK", 2 * i + 2 );
            list.set ( ("smoother: type (level " + levelIndex + ")").data(), levelParamValue );
        }
    }

    Int SmootherNumLevelSweeps              = dataFile ( (section + "/" + subSection + "/smoother/level_sweeps").data(), 0, found );
    if ( found )
    {
        for (Int i (0); i < SmootherNumLevelSweeps; ++i)
        {
            std::string levelIndex          = dataFile ( (section + "/" + subSection + "/smoother/level_sweeps").data(), "0", 2 * i + 1 );
            Int levelParamValue             = dataFile ( (section + "/" + subSection + "/smoother/level_sweeps").data(), 1, 2 * i + 2 );
            list.set ( ("smoother: sweeps (level " + levelIndex + ")").data(), levelParamValue );
        }
    }

    // subsmoother parameter

    std::string SubSmootherType             = dataFile ( (section + "/" + subSection + "/subsmoother/type").data(), "Chebyshev", found );
    if ( found )
    {
        list.set ( "subsmoother: type", SubSmootherType );
    }

    Real   SubSmootherChebyshevAlpha        = dataFile ( (section + "/" + subSection + "/subsmoother/Chebyshev_alpha").data(), 20., found );
    if ( found )
    {
        list.set ( "subsmoother: Chebyshev alpha", SubSmootherChebyshevAlpha );
    }

    //    Real   SubSmootherSGSDampingFactor      = dataFile((section + "/" + subSection + "/subsmoothers/SGS_damping_factor").data(), 1., found );
    Int SubSmootherEdgeSweeps               = dataFile ( (section + "/" + subSection + "/subsmoother/edge_sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "subsmoother: edge sweeps", SubSmootherEdgeSweeps );
    }

    Int SubSmootherNodeSweeps               = dataFile ( (section + "/" + subSection + "/subsmoother/node_sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "subsmoother: node sweeps", SubSmootherNodeSweeps );
    }


    // Coarsest Grid Parameters

    Int CoarseMaxSize                 = dataFile ( (section + "/" + subSection + "/coarse/max_size").data(), 128, found );
    if ( found )
    {
        list.set ( "coarse: max size",  CoarseMaxSize );
    }

    std::string CoarseType            = dataFile ( (section + "/" + subSection + "/coarse/type").data(), "Chebyshev", found );
    if ( found )
    {
        list.set ( "coarse: type", CoarseType );
    }

    std::string CoarsePreOrPost       = dataFile ( (section + "/" + subSection + "/coarse/pre_or_post").data(), "post", found );
    if ( found )
    {
        list.set ( "coarse: pre or post", CoarsePreOrPost );
    }

    Int CoarseSweeps                  = dataFile ( (section + "/" + subSection + "/coarse/sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "coarse: sweeps", CoarseSweeps );
    }

    Real   CoarseDampingFactor        = dataFile ( (section + "/" + subSection + "/coarse/damping_factor").data(), 1.0, found );
    if ( found )
    {
        list.set ( "coarse: damping factor", CoarseDampingFactor );
    }

    std::string CoarseSubsmootherType = dataFile ( (section + "/" + subSection + "/coarse/subsmoother_type").data(), "Chebyshev", found );
    if ( found )
    {
        list.set ( "coarse: subsmoother type", CoarseSubsmootherType );
    }

    Int CoarseNodeSweeps              = dataFile ( (section + "/" + subSection + "/coarse/node_sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "coarse: node sweeps", CoarseNodeSweeps );
    }

    Int CoarseEdgeSweeps              = dataFile ( (section + "/" + subSection + "/coarse/edge_sweeps").data(), 2, found );
    if ( found )
    {
        list.set ( "coarse: edge sweeps", CoarseEdgeSweeps );
    }

    Real   CoarseChebyshevAlpha       = dataFile ( (section + "/" + subSection + "/coarse/Chebyshev_alpha").data(), 30., found );
    if ( found )
    {
        list.set ( "coarse: Chebyshev alpha", CoarseChebyshevAlpha );
    }

    Int CoarseMaxProcesses            = dataFile ( (section + "/" + subSection + "/coarse/max_processes").data(), -1, found );
    if ( found )
    {
        list.set ( "coarse: max processes", CoarseMaxProcesses );
    }


    // Load-balancing Options
    Int RepartitionEnable              = dataFile ( (section + "/" + subSection + "/repartition/enable").data(), 0, found );
    if ( found )
    {
        list.set ( "repartition: enable", RepartitionEnable );
    }

    std::string RepartitionPartitioner = dataFile ( (section + "/" + subSection + "/repartition/partitioner").data(), "ParMETIS", found );
    if ( found )
    {
        list.set ( "repartition: partitioner", RepartitionPartitioner );
    }

    Real   RepartitionMaxMinRatio      = dataFile ( (section + "/" + subSection + "/repartition/max_min_ratio").data(), 1.3, found );
    if ( found )
    {
        list.set ( "repartition: max min ratio", RepartitionMaxMinRatio );
    }

    Int RepartitionMinPerProc          = dataFile ( (section + "/" + subSection + "/repartition/min_per_proc").data(), 512, found );
    if ( found )
    {
        list.set ( "repartition: min per proc", RepartitionMinPerProc );
    }

    Real   RepartitionNodeMaxMinRatio  = dataFile ( (section + "/" + subSection + "/repartition/node_max_min_ratio").data(), 1.3, found );
    if ( found )
    {
        list.set ( "repartition: node max min ratio", RepartitionNodeMaxMinRatio );
    }

    Int RepartitionNodeMinPerProc      = dataFile ( (section + "/" + subSection + "/repartition/node_min_per_proc").data(), 170, found );
    if ( found )
    {
        list.set ( "repartition: node min per proc", RepartitionNodeMinPerProc );
    }

    Int RepartitionZoltanDimensions    = dataFile ( (section + "/" + subSection + "/repartition/Zoltan_dimensions").data(), 2, found );
    if ( found )
    {
        list.set ( "repartition: Zoltan dimensions", RepartitionZoltanDimensions );
    }

    if ( MLPrintParameterList && verbose )
    {
        std::cout << "ML parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

void PreconditionerML::showMe ( std::ostream& output ) const
{
    output << "showMe must be implemented for the PreconditionerML class" << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerML::setDataFromGetPot ( const GetPot&      dataFile,
                                      const std::string& section )
{
    bool verbose = M_comm->MyPID() == 0;

    M_analyze = dataFile ( (section + "/" + "ML" + "/analyze_smoother" ).data(), false); // To be moved in createMLList

    // ML List
    createMLList ( M_list, dataFile, section, "ML", verbose );

    // visualization
    bool found (false);
    bool enableViz    = dataFile ( (section + "/" + "ML" + "/visualization/enable").data(), false, found );
    if ( found )
    {
        /*
           If the visualization is desired and we have set the required data,
           we set the following variables.
           (see Trilinos::ML manual for more details)
         */
        if (M_visualizationDataAvailable)
        {
            M_list.set ( "viz: enable", enableViz);
            M_list.set ( "viz: output format", "vtk");
            M_list.set ("x-coordinates", & ( (*M_xCoord) [0]) );
            M_list.set ("y-coordinates", & ( (*M_yCoord) [0]) );
            M_list.set ("z-coordinates", & ( (*M_zCoord) [0]) );
        }
        else
        {
            std::cout << "Warning: Visualization options are not available if you have not use setVerticesCoordinates first!" << std::endl;
        }
    }

    // IfPack list
    bool MLPrintParameterList = dataFile ( (section + "/displayList").data(), false );
    if ( MLPrintParameterList && verbose )
    {
        std::cout << "Smoother: ";
    }
    list_Type& SmootherIFSubList = M_list.sublist ( "smoother: ifpack list" );
    PreconditionerIfpack::createIfpackList ( SmootherIFSubList, dataFile, section, "ML", verbose );
}

void
PreconditionerML::setVerticesCoordinates (boost::shared_ptr<std::vector<Real> > xCoord,
                                          boost::shared_ptr<std::vector<Real> > yCoord,
                                          boost::shared_ptr<std::vector<Real> > zCoord)
{
    M_xCoord = xCoord;
    M_yCoord = yCoord;
    M_zCoord = zCoord;
    M_visualizationDataAvailable = true;
}


// ===================================================
// Get Methods
// ===================================================
Real
PreconditionerML::condest()
{
    return 0.;
}

Preconditioner::prec_raw_type*
PreconditionerML::preconditioner()
{
    return M_preconditioner.get();
}

} // namespace LifeV
