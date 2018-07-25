//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
    @brief Test for PartitionIO class - cut and write

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>
    @date 10-05-2012

    Partition a mesh using a single (MPI) process and save mesh parts
    to an HDF5 file.
 */

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAS_HDF5


#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>


#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/filter/ExporterHDF5Mesh3D.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>

using namespace LifeV;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

int main (int argc, char** argv)
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    typedef RegionMesh<LinearTetra> mesh_Type;

    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    if (comm->NumProc() != 1)
    {
        std::cout << "This test needs to be run "
                  << "with a single process. Aborting."
                  << std::endl;
        return (EXIT_FAILURE);
    }

    GetPot commandLine (argc, argv);
    string dataFileName = commandLine.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    const UInt numElements (dataFile ("mesh/nelements", 10) );
    const UInt numParts (dataFile ("test/num_parts", 4) );
    const std::string partsFileName (dataFile ("test/hdf5_file_name", "cube.h5") );
    const std::string ioClass (dataFile ("test/io_class", "new") );

    std::cout << "Number of elements in mesh: " << numElements << std::endl;
    std::cout << "Number of parts: " << numParts << std::endl;
    std::cout << "Name of HDF5 container: " << partsFileName << std::endl;

    boost::shared_ptr<mesh_Type> fullMeshPtr (new mesh_Type ( comm ) );
    regularMesh3D (*fullMeshPtr, 1, numElements, numElements, numElements,
                   false, 2.0, 2.0, 2.0, -1.0, -1.0, -1.0);

    MeshPartitioner<mesh_Type> meshPart;
    meshPart.setup (numParts, comm);

    meshPart.attachUnpartitionedMesh (fullMeshPtr);
    meshPart.doPartitionGraph();
    meshPart.doPartitionMesh();

    // Release the original mesh from the MeshPartitioner object and
    // delete the RegionMesh object
    meshPart.releaseUnpartitionedMesh();
    fullMeshPtr.reset();

    // Write mesh parts to HDF5 container
    if (! ioClass.compare ("old") )
    {
        ExporterHDF5Mesh3D<mesh_Type> HDF5Output (dataFile,
                                                  meshPart.meshPartition(),
                                                  partsFileName,
                                                  comm->MyPID() );
        HDF5Output.addPartitionGraph (meshPart.elementDomains(), comm);
        HDF5Output.addMeshPartitionAll (meshPart.meshPartitions(), comm);
        HDF5Output.postProcess (0);
        HDF5Output.closeFile();
    }
    else
    {
        boost::shared_ptr<Epetra_MpiComm> mpiComm =
            boost::dynamic_pointer_cast<Epetra_MpiComm> (comm);
        PartitionIO<mesh_Type> partitionIO (partsFileName, mpiComm);
        partitionIO.write (meshPart.meshPartitions() );
    }

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* LIFEV_HAS_HDF5 */

    return (EXIT_SUCCESS);
}
