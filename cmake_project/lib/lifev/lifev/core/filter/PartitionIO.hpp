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
  @brief Class that handles I/O of mesh parts (for offline partitioning mode)

  @date 10-05-2012
  @author Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef PARTITION_IO_H_
#define PARTITION_IO_H_

#include <algorithm>


#include <Epetra_config.h>

#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

#include <Epetra_MpiComm.h>


#include<lifev/core/LifeV.hpp>
#include <lifev/core/filter/HDF5IO.hpp>

namespace LifeV
{

/*!
  @brief Class that handles I/O of mesh parts (for offline partitioning mode)
  @author Radu Popescu radu.popescu@epfl.ch

  This class is used to write the mesh parts produced by (offline) mesh
  partitioning into a single HDF5 container. This part is done with a single
  MPI process (usually on a workstation).

  The class is later used during an online simulation (multiple MPI processes,
  on a cluster, supercomputer etc.) to read the mesh parts.

  Usage:
      - call the constructor (a default one is also available)
      - if the default constructor was used, call setup method to initialize
        the object
      - call write method to store a set of mesh parts in the HDF5 file
      - call read method to load the mesh part associate with the current
        MPI process from the HDF5 file

  Creating, opening and closing the HDF5 file is done automatically by the
  object.

  Description of the storage format:
  N - number of mesh parts

  The HDF5 container contains 6 tables (num. of rows is the leading dimension):

  1. stats - N x 15 (unsigned integer) - each row belongs to a mesh part and
     contains the following values, in order:
         - N
         - graph part size (unused)
         - number of points
         - number of boundary points
         - number of vertices
         - number of boundary vertices
         - number of global vertices
         - number of edges
         - number of boundary edges
         - number of global edges
         - number of faces
         - number of boundary faces
         - number of global faces
         - number of volumes
         - number of global volumes
  2. point_ids - (3 * N) x (maximum number of points in mesh parts) (unsigned
     integer) - each mesh part uses three consecutive rows, with the following
     contents:
         - point markers
         - point global ids
         - point flags
  3. point_coords - (3 * N) x (maximum number of points in mesh parts) (double)
     each mesh part uses three consecutive rows, with the following contents:
         - point x coordinates
         - point y coordinates
         - point z coordinates
  4. edges - (5 * N) x (maximum number of edges in mesh parts) (unsigned
     integer) - each mesh part uses three consecutive rows, with the following
     contents:
         - local ids of first point of edges
         - local ids of second point of edges
         - edge markers
         - edge global ids
         - edge flags
  5. faces ((7 + num of face nodes) * N) x (maximum number of faces in mesh
     parts) (unsigned integer) - each mesh part uses (7 + num of face nodes)
     consecutive rows, with the following contents:
         - one row for the local ids of each face node of the faces in the part
         - face markers
         - face global ids
         - id of first neighbour element
         - id of second neighbour element
         - position of first neighbour element
         - position of second neighbour element
         - face flags
  6. elements ((3 + num of element nodes) * N) x (maximum number of elements
     in mesh parts) (unsigned integer) - each mesh part uses (2 + num of
     element nodes) consecutive rows, with the following contents:
         - one row for the local ids of each element node of the elements in
           the mesh part
         - element markers
         - element global ids
         - element flag
*/
template<typename MeshType>
class PartitionIO
{
public:
    //! @name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef std::vector<meshPtr_Type> meshParts_Type;
    typedef boost::shared_ptr<meshParts_Type> meshPartsPtr_Type;
    // TODO: enable the following two typedefs if I decide to also
    // write the graph to the HDF5 container
    // typedef std::vector<std::vector<Int> > table_Type;
    // typedef boost::shared_ptr<table_Type> tablePtr_Type;
    typedef boost::shared_ptr<Epetra_MpiComm> commPtr_Type;
    //@}

    //! \name Constructors & Destructors
    //@{
    //! Default empty constructor
    PartitionIO() {}

    //! Constructor
    /*!
     * Non-empty constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param transposeInFile (default false) write/read tables transposed
     *        in the HDF5. This option is available because the transposed
     *        format may be faster on certain machines. Leave as default
     *        unless certain that it works for a given machine.
     */
    PartitionIO (const std::string& fileName,
                 const commPtr_Type& comm,
                 const bool transposeInFile = false);

    //! Empty destructor
    virtual ~PartitionIO() {}
    //@}

    //! \name Public Methods
    //@{
    //! Initialization method
    /*!
     * Initialization method that is used with the default constructor
     * \param fileName the name of the HDF5 file to be used
     * \param comm pointer to Epetra_Comm
     * \param transposeInFile (default false) write/read tables transposed
     *        in the HDF5. This option is available because the transposed
     *        format may be faster on certain machines. Leave as default
     *        unless certain that it works for a given machine.
     */
    void setup (const std::string& fileName,
                const commPtr_Type& comm,
                const bool transposeInFile = false);
    //! Write method
    /*!
     * Call this method to write the mesh parts to disk
     * \param meshParts pointer to a vector containing pointers to the mesh
     *        parts (RegionMesh objects). This pointer is released after
     *        writing, so the mesh parts can be cleared from memory
     */
    void write (const meshPartsPtr_Type& meshParts);
    //! Read method
    /*!
     * Call this method to read from the HDF5 file the mesh part associated
     * with the rank of the current MPI process
     * \param meshPart pointer to a RegionMesh that will contain the mesh
     *        part. If the RegionMesh has been initialized and contains any
     *        state, it will be destroyed before reading
     */
    void read (meshPtr_Type& meshPart);
    //@}

private:
    // Copy constructor and assignment operator are disabled
    PartitionIO (const PartitionIO&);
    PartitionIO& operator= (const PartitionIO&);

    //! Private Methods
    //@{
    // Methods for writing
    void writeStats();
    void writePoints();
    void writeEdges();
    void writeFaces();
    void writeElements();
    // Methods for reading
    void readStats();
    void readPoints();
    void readEdges();
    void readFaces();
    void readElements();
    //@}

    //! Private Data Members
    //@{
    commPtr_Type M_comm;
    UInt M_myRank;
    std::string M_fileName;
    bool M_transposeInFile;
    meshPartsPtr_Type M_meshPartsOut;
    meshPtr_Type M_meshPartIn;
    // Mesh geometry
    UInt M_elementNodes;
    UInt M_faceNodes;
    UInt M_maxNumPoints;
    UInt M_maxNumEdges;
    UInt M_maxNumFaces;
    UInt M_maxNumElements;
    UInt M_numParts;
    // Counters when loading the mesh part
    UInt M_numPoints;
    UInt M_numEdges;
    UInt M_numFaces;
    UInt M_numElements;
    //HDF5 I/O filter
    HDF5IO M_HDF5IO;
    // Buffers for reading/writing
    std::vector<UInt> M_uintBuffer;
    std::vector<Real> M_realBuffer;
    //@}
}; // class PartitionIO

template<typename MeshType>
inline PartitionIO<MeshType>::PartitionIO (const std::string& fileName,
                                           const commPtr_Type& comm,
                                           const bool transposeInFile) :
    M_comm (comm),
    M_fileName (fileName),
    M_transposeInFile (transposeInFile),
    M_maxNumPoints (0),
    M_maxNumEdges (0),
    M_maxNumFaces (0),
    M_maxNumElements (0)
{
    M_elementNodes = MeshType::elementShape_Type::S_numPoints;
    M_faceNodes = MeshType::elementShape_Type::GeoBShape::S_numPoints;

    M_myRank = M_comm->MyPID();
}

template<typename MeshType>
inline void PartitionIO<MeshType>::setup (const std::string& fileName,
                                          const commPtr_Type& comm,
                                          const bool transposeInFile)
{
    M_comm = comm;
    M_fileName = fileName;
    M_transposeInFile = transposeInFile;
    M_maxNumPoints = 0;
    M_maxNumEdges = 0;
    M_maxNumFaces = 0;
    M_maxNumElements = 0;

    M_elementNodes = MeshType::elementShape_Type::S_numPoints;
    M_faceNodes = MeshType::elementShape_Type::GeoBShape::S_numPoints;
}

template<typename MeshType>
void PartitionIO<MeshType>::write (const meshPartsPtr_Type& meshParts)
{
    M_meshPartsOut = meshParts;
    M_numParts = M_meshPartsOut->size();

    M_HDF5IO.openFile (M_fileName, M_comm, false);
    writeStats();
    writePoints();
    writeEdges();
    writeFaces();
    writeElements();
    M_HDF5IO.closeFile();

    M_meshPartsOut.reset();
}

template<typename MeshType>
void PartitionIO<MeshType>::read (meshPtr_Type& meshPart)
{
    meshPart.reset();
    M_meshPartIn.reset (new mesh_Type);

    M_HDF5IO.openFile (M_fileName, M_comm, true);
    readStats();
    readPoints();
    readEdges();
    readFaces();
    readElements();
    M_HDF5IO.closeFile();

    meshPart = M_meshPartIn;
    M_meshPartIn.reset();
}

template<typename MeshType>
void PartitionIO<MeshType>::writeStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = M_numParts;
        currentSpaceDims[1] = 15;
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = 15;
        currentSpaceDims[1] = M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
    }

    // Create new table
    M_HDF5IO.createTable ("stats", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (15);
    for (UInt i = 0; i < M_numParts; ++i)
    {
        mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
        M_uintBuffer[0] = M_numParts;
        // Next one is unused. I'll put it to keep compatibility
        // in case I decide to write the graphs
        M_uintBuffer[1] = 0;
        M_uintBuffer[2] = currentPart.numPoints();
        if (M_uintBuffer[2] > M_maxNumPoints)
        {
            M_maxNumPoints = M_uintBuffer[2];
        }
        M_uintBuffer[3] = currentPart.numBPoints();
        M_uintBuffer[4] = currentPart.numVertices();
        M_uintBuffer[5] = currentPart.numBVertices();
        M_uintBuffer[6] = currentPart.numGlobalVertices();
        M_uintBuffer[7] = currentPart.numEdges();
        if (M_uintBuffer[7] > M_maxNumEdges)
        {
            M_maxNumEdges = M_uintBuffer[7];
        }
        M_uintBuffer[8] = currentPart.numBEdges();
        M_uintBuffer[9] = currentPart.numGlobalEdges();
        M_uintBuffer[10] = currentPart.numFaces();
        if (M_uintBuffer[10] > M_maxNumFaces)
        {
            M_maxNumFaces = M_uintBuffer[10];
        }
        M_uintBuffer[11] = currentPart.numBFaces();
        M_uintBuffer[12] = currentPart.numGlobalFaces();
        M_uintBuffer[13] = currentPart.numVolumes();
        if (M_uintBuffer[13] > M_maxNumElements)
        {
            M_maxNumElements = M_uintBuffer[13];
        }
        M_uintBuffer[14] = currentPart.numGlobalVolumes();

        hsize_t currentOffset[2];
        if (! M_transposeInFile)
        {
            currentOffset[0] = i;
            currentOffset[1] = 0;
        }
        else
        {
            currentOffset[0] = 0;
            currentOffset[1] = i;
        }
        M_HDF5IO.write ("stats", H5T_NATIVE_UINT, currentCount, currentOffset,
                        &M_uintBuffer[0]);
    }

    M_HDF5IO.closeTable ("stats");
}

template<typename MeshType>
void PartitionIO<MeshType>::writePoints()
{
    // Write points (N = number of parts)
    // Two tables: one is (3 * N) x max_num_points of int
    //             second is (3 * N) x max_num_points of real
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = 3 * M_numParts;
        currentSpaceDims[1] = M_maxNumPoints;
        currentCount[0] = 3;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumPoints;
        currentSpaceDims[1] = 3 * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 3;
    }

    M_HDF5IO.createTable ("point_ids", H5T_STD_U32BE, currentSpaceDims);
    M_HDF5IO.createTable ("point_coords", H5T_IEEE_F64BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_realBuffer.resize (currentCount[0] * currentCount[1], 0);
    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numPoints(); ++j)
            {
                M_uintBuffer[j] = currentPart.pointList[j].markerID();
                M_uintBuffer[stride + j] = currentPart.point (j).id();
                M_uintBuffer[2 * stride + j] =
                    static_cast<int> (currentPart.pointList[j].flag() );
                M_realBuffer[j] = currentPart.pointList[j].x();
                M_realBuffer[stride + j] = currentPart.pointList[j].y();
                M_realBuffer[2 * stride + j] = currentPart.pointList[j].z();
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                            currentOffset, &M_realBuffer[0]);
        }
    }
    else
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numPoints(); ++j)
            {
                M_uintBuffer[stride * j] = currentPart.pointList[j].markerID();
                M_uintBuffer[stride * j + 1] = currentPart.point (j).id();
                M_uintBuffer[stride * j + 2] =
                    static_cast<int> (currentPart.pointList[j].flag() );
                M_realBuffer[stride * j] = currentPart.pointList[j].x();
                M_realBuffer[stride * j + 1] = currentPart.pointList[j].y();
                M_realBuffer[stride * j + 2] = currentPart.pointList[j].z();
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("point_ids", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
            M_HDF5IO.write ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                            currentOffset, &M_realBuffer[0]);
        }
    }
    M_realBuffer.resize (0);

    M_HDF5IO.closeTable ("point_ids");
    M_HDF5IO.closeTable ("point_coords");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeEdges()
{
    // Write edges (N = number of parts)
    // Table is (5 * N) x max_num_edges of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = 5 * M_numParts;
        currentSpaceDims[1] = M_maxNumEdges;
        currentCount[0] = 5;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumEdges;
        currentSpaceDims[1] = 5 * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 5;
    }

    M_HDF5IO.createTable ("edges", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numEdges(); ++j)
            {
                M_uintBuffer[j] = currentPart.edgeList[j].point (0).localId();
                M_uintBuffer[stride + j] =
                    currentPart.edgeList[j].point (1).localId();
                M_uintBuffer[2 * stride + j] =
                    currentPart.edgeList[j].markerID();
                M_uintBuffer[3 * stride + j] = currentPart.edgeList[j].id();
                M_uintBuffer[4 * stride + j] =
                    static_cast<int> (currentPart.edgeList[j].flag() );
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numEdges(); ++j)
            {
                M_uintBuffer[stride * j] =
                    currentPart.edgeList[j].point (0).localId();
                M_uintBuffer[stride * j + 1] =
                    currentPart.edgeList[j].point (1).localId();
                M_uintBuffer[stride * j + 2] =
                    currentPart.edgeList[j].markerID();
                M_uintBuffer[stride * j + 3] = currentPart.edgeList[j].id();
                M_uintBuffer[stride * j + 4] =
                    static_cast<int> (currentPart.edgeList[j].flag() );
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("edges", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("edges");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeFaces()
{
    // Write faces (N = number of parts)
    // Table is ((7 + num_face_nodes * N) x max_num_faces of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (7 + M_faceNodes) * M_numParts;
        currentSpaceDims[1] = M_maxNumFaces;
        currentCount[0] = 7 + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumFaces;
        currentSpaceDims[1] = (7 + M_faceNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 7 + M_faceNodes;
    }

    M_HDF5IO.createTable ("faces", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numFaces(); ++j)
            {
                for (UInt k = 0; k < M_faceNodes; ++k)
                {
                    M_uintBuffer[k * stride + j] =
                        currentPart.faceList[j].point (k).localId();
                }
                M_uintBuffer[M_faceNodes * stride + j] =
                    currentPart.faceList[j].markerID();
                M_uintBuffer[ (M_faceNodes + 1) * stride + j] =
                    currentPart.faceList[j].id();
                M_uintBuffer[ (M_faceNodes + 2) * stride + j] =
                    currentPart.faceList[j].firstAdjacentElementIdentity();
                M_uintBuffer[ (M_faceNodes + 3) * stride + j] =
                    currentPart.faceList[j].secondAdjacentElementIdentity();
                M_uintBuffer[ (M_faceNodes + 4) * stride + j] =
                    currentPart.faceList[j].firstAdjacentElementPosition();
                M_uintBuffer[ (M_faceNodes + 5) * stride + j] =
                    currentPart.faceList[j].secondAdjacentElementPosition();
                M_uintBuffer[ (M_faceNodes + 6) * stride + j] =
                    static_cast<int> (currentPart.faceList[j].flag() );
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numFaces(); ++j)
            {
                for (UInt k = 0; k < M_faceNodes; ++k)
                {
                    M_uintBuffer[stride * j + k] =
                        currentPart.faceList[j].point (k).localId();
                }
                M_uintBuffer[stride * j + M_faceNodes] =
                    currentPart.faceList[j].markerID();
                M_uintBuffer[stride * j + M_faceNodes + 1] =
                    currentPart.faceList[j].id();
                M_uintBuffer[stride * j + M_faceNodes + 2] =
                    currentPart.faceList[j].firstAdjacentElementIdentity();
                M_uintBuffer[stride * j + M_faceNodes + 3] =
                    currentPart.faceList[j].secondAdjacentElementIdentity();
                M_uintBuffer[stride * j + M_faceNodes + 4] =
                    currentPart.faceList[j].firstAdjacentElementPosition();
                M_uintBuffer[stride * j + M_faceNodes + 5] =
                    currentPart.faceList[j].secondAdjacentElementPosition();
                M_uintBuffer[stride * j + M_faceNodes + 6] =
                    static_cast<int> (currentPart.faceList[j].flag() );
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("faces", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("faces");
}

template<typename MeshType>
void PartitionIO<MeshType>::writeElements()
{
    // Write elements (N = number of parts)
    // Table is ((3 + num_element_nodes * N) x max_num_elements of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    if (! M_transposeInFile)
    {
        currentSpaceDims[0] = (3 + M_elementNodes) * M_numParts;
        currentSpaceDims[1] = M_maxNumElements;
        currentCount[0] = 3 + M_elementNodes;
        currentCount[1] = currentSpaceDims[1];
    }
    else
    {
        currentSpaceDims[0] = M_maxNumElements;
        currentSpaceDims[1] = (3 + M_elementNodes) * M_numParts;
        currentCount[0] = currentSpaceDims[0];;
        currentCount[1] = 3 + M_elementNodes;
    }

    M_HDF5IO.createTable ("elements", H5T_STD_U32BE, currentSpaceDims);

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numVolumes(); ++j)
            {
                for (UInt k = 0; k < M_elementNodes; ++k)
                {
                    M_uintBuffer[k * stride + j] =
                        currentPart.volumeList[j].point (k).localId();
                }
                M_uintBuffer[M_elementNodes * stride + j] =
                    currentPart.volumeList[j].markerID();
                M_uintBuffer[ (M_elementNodes + 1) * stride + j] =
                    currentPart.volumeList[j].id();
                M_uintBuffer[ (M_elementNodes + 2) * stride + j] =
                    static_cast<int> (currentPart.volumeList[j].flag() );
            }

            hsize_t currentOffset[2] = {i* currentCount[0], 0};
            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }
    else
    {
        for (UInt i = 0; i < M_numParts; ++i)
        {
            mesh_Type& currentPart = (* (*M_meshPartsOut) [i]);
            for (UInt j = 0; j < currentPart.numVolumes(); ++j)
            {
                for (UInt k = 0; k < M_elementNodes; ++k)
                {
                    M_uintBuffer[stride * j + k] =
                        currentPart.volumeList[j].point (k).localId();
                }
                M_uintBuffer[stride * j + M_elementNodes] =
                    currentPart.volumeList[j].markerID();
                M_uintBuffer[stride * j + M_elementNodes + 1] =
                    currentPart.volumeList[j].id();
                M_uintBuffer[stride * j + M_elementNodes + 2] =
                    static_cast<int> (currentPart.volumeList[j].flag() );
            }

            hsize_t currentOffset[2] = {0, i* currentCount[1]};
            M_HDF5IO.write ("elements", H5T_NATIVE_UINT, currentCount,
                            currentOffset, &M_uintBuffer[0]);
        }
    }

    M_HDF5IO.closeTable ("elements");
}

template<typename MeshType>
void PartitionIO<MeshType>::readStats()
{
    // Write mesh partition stats (N = number of parts)
    // This is an N x 15 table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("stats", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 1;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 1;
        currentOffset[0] = 0;
        currentOffset[1] = M_myRank;
    }

    // Fill buffer
    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("stats", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("stats");

    // Insert stats into mesh partition object
    M_numPoints = M_uintBuffer[2];
    M_meshPartIn->setMaxNumPoints (M_uintBuffer[2], true);
    M_meshPartIn->setNumBPoints (M_uintBuffer[3]);

    // Vertices
    M_meshPartIn->setNumVertices (M_uintBuffer[4]);
    M_meshPartIn->setNumBVertices (M_uintBuffer[5]);
    M_meshPartIn->setNumGlobalVertices (M_uintBuffer[6]);

    // Edges
    M_numEdges = M_uintBuffer[7];
    M_meshPartIn->setNumEdges (M_uintBuffer[7]);
    M_meshPartIn->setMaxNumEdges (M_uintBuffer[7], true);
    M_meshPartIn->setNumBEdges (M_uintBuffer[8]);
    M_meshPartIn->setMaxNumGlobalEdges (M_uintBuffer[9]);

    // Faces
    M_numFaces = M_uintBuffer[10];
    M_meshPartIn->setNumFaces (M_uintBuffer[10]);
    M_meshPartIn->setMaxNumFaces (M_uintBuffer[10], true);
    M_meshPartIn->setNumBFaces (M_uintBuffer[11]);
    M_meshPartIn->setMaxNumGlobalFaces (M_uintBuffer[12]);

    // Volumes
    M_numElements = M_uintBuffer[13];
    M_meshPartIn->setMaxNumVolumes (M_uintBuffer[13], true);
    M_meshPartIn->setMaxNumGlobalVolumes (M_uintBuffer[14]);
}

template<typename MeshType>
void PartitionIO<MeshType>::readPoints()
{
    // Read mesh points (N = number of parts)
    // There are two tables: a (3 * N) x max_num_points table of int and
    // a (3 * N) x max_num_points table of real
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    //hid_t realDataset = H5Dopen(M_fileId, "point_coords", H5P_DEFAULT);

    M_HDF5IO.openTable ("point_ids", currentSpaceDims);
    M_HDF5IO.openTable ("point_coords", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 3;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 3;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_realBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("point_ids", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);
    M_HDF5IO.read ("point_coords", H5T_NATIVE_DOUBLE, currentCount,
                   currentOffset, &M_realBuffer[0]);

    M_HDF5IO.closeTable ("point_ids");
    M_HDF5IO.closeTable ("point_coords");

    // Insert points into the mesh partition object
    M_meshPartIn->pointList.reserve (M_numPoints);
    M_meshPartIn->_bPoints.reserve (M_meshPartIn->numBPoints() );

    typename MeshType::point_Type* pp = 0;

    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt j = 0; j < M_numPoints; ++j)
        {
            pp = & ( M_meshPartIn->addPoint ( false, false ) );
            pp->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[2 * stride + j]) );
            pp->setMarkerID (M_uintBuffer[j]);
            pp->x() = M_realBuffer[j];
            pp->y() = M_realBuffer[stride + j];
            pp->z() = M_realBuffer[2 * stride + j];
            pp->setId (M_uintBuffer[stride + j]);
        }
    }
    else
    {
        for (UInt j = 0; j < M_numPoints; ++j)
        {
            pp = & ( M_meshPartIn->addPoint ( false, false ) );
            pp->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[stride * j + 2]) );
            pp->setMarkerID (M_uintBuffer[stride * j]);
            pp->x() = M_realBuffer[stride * j];
            pp->y() = M_realBuffer[stride * j + 1];
            pp->z() = M_realBuffer[stride * j + 2];
            pp->setId (M_uintBuffer[stride * j + 1]);
        }
    }
    M_realBuffer.resize (0);
}

template<typename MeshType>
void PartitionIO<MeshType>::readEdges()
{
    // Read mesh edges (N = number of parts)
    // Read a (5 * N) x max_num_edges table of int and
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("edges", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 5;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 5;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("edges", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("edges");

    M_meshPartIn->edgeList.reserve (M_numEdges);

    typename MeshType::edge_Type* pe;

    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt j = 0; j < M_numEdges; ++j)
        {
            pe = & (M_meshPartIn->addEdge (false) );
            pe->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[4 * stride + j]) );
            pe->setId (M_uintBuffer[3 * stride + j]);
            pe->setPoint (0, M_meshPartIn->point (M_uintBuffer[j]) );
            pe->setPoint (1, M_meshPartIn->point (M_uintBuffer[stride + j]) );
            pe->setMarkerID (M_uintBuffer[2 * stride + j]);
        }
    }
    else
    {
        for (UInt j = 0; j < M_numEdges; ++j)
        {
            pe = & (M_meshPartIn->addEdge (false) );
            pe->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[stride * j + 4]) );
            pe->setId (M_uintBuffer[stride * j + 3]);
            pe->setPoint (0, M_meshPartIn->point (M_uintBuffer[stride * j]) );
            pe->setPoint (1, M_meshPartIn->point (M_uintBuffer[stride * j + 1]) );
            pe->setMarkerID (M_uintBuffer[stride * j + 2]);
        }
    }
}

template<typename MeshType>
void PartitionIO<MeshType>::readFaces()
{
    // read mesh faces (N = number of parts)
    // Read a ((7 + num_face_points) * N) x max_num_faces table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("faces", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 7 + M_faceNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 7 + M_faceNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("faces", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("faces");

    typename MeshType::face_Type* pf = 0;

    M_meshPartIn->faceList.reserve (M_numFaces);

    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt j = 0; j < M_numFaces; ++j)
        {
            pf = & (M_meshPartIn->addFace (false) );
            pf->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[ (6 + M_faceNodes)
                                                      * stride + j]) );
            pf->setId (M_uintBuffer[ (M_faceNodes + 1) * stride + j]);
            pf->firstAdjacentElementIdentity() =
                M_uintBuffer[ (M_faceNodes + 2) * stride + j];
            pf->secondAdjacentElementIdentity() =
                M_uintBuffer[ (M_faceNodes + 3) * stride + j];
            pf->firstAdjacentElementPosition() =
                M_uintBuffer[ (M_faceNodes + 4) * stride + j];
            pf->secondAdjacentElementPosition() =
                M_uintBuffer[ (M_faceNodes + 5) * stride + j];
            pf->setMarkerID (M_uintBuffer[M_faceNodes * stride + j]);
            for (UInt k = 0; k < M_faceNodes; ++k)
            {
                pf->setPoint (k, M_meshPartIn->point (
                                  M_uintBuffer[k * stride + j]) );
            }
        }
    }
    else
    {
        for (UInt j = 0; j < M_numFaces; ++j)
        {
            pf = & (M_meshPartIn->addFace (false) );
            pf->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[stride * j
                                                     + M_faceNodes + 6]) );
            pf->setId (M_uintBuffer[stride * j + M_faceNodes + 1]);
            pf->firstAdjacentElementIdentity() =
                M_uintBuffer[stride * j + M_faceNodes + 2];
            pf->secondAdjacentElementIdentity() =
                M_uintBuffer[stride * j + M_faceNodes + 3];
            pf->firstAdjacentElementPosition() =
                M_uintBuffer[stride * j + M_faceNodes + 4];
            pf->secondAdjacentElementPosition() =
                M_uintBuffer[stride * j + M_faceNodes + 5];
            pf->setMarkerID (M_uintBuffer[ (7 + M_faceNodes) * j + M_faceNodes]);
            for (UInt k = 0; k < M_faceNodes; ++k)
            {
                pf->setPoint (k, M_meshPartIn->point (
                                  M_uintBuffer[stride * j + k]) );
            }
        }
    }

    M_meshPartIn->setLinkSwitch ("HAS_ALL_FACETS");
    M_meshPartIn->setLinkSwitch ("FACETS_HAVE_ADIACENCY");
}

template<typename MeshType>
void PartitionIO<MeshType>::readElements()
{
    // Read mesh elements (N = number of parts)
    // Read a ((3 + num_element_points) * N) x max_num_elements table of int
    hsize_t currentSpaceDims[2];
    hsize_t currentCount[2];
    hsize_t currentOffset[2];

    M_HDF5IO.openTable ("elements", currentSpaceDims);

    if (! M_transposeInFile)
    {
        currentCount[0] = 3 + M_elementNodes;
        currentCount[1] = currentSpaceDims[1];
        currentOffset[0] = currentCount[0] * M_myRank;
        currentOffset[1] = 0;
    }
    else
    {
        currentCount[0] = currentSpaceDims[0];
        currentCount[1] = 3 + M_elementNodes;
        currentOffset[0] = 0;
        currentOffset[1] = currentCount[1] * M_myRank;
    }

    M_uintBuffer.resize (currentCount[0] * currentCount[1], 0);
    M_HDF5IO.read ("elements", H5T_NATIVE_UINT, currentCount, currentOffset,
                   &M_uintBuffer[0]);

    M_HDF5IO.closeTable ("elements");

    M_meshPartIn->volumeList.reserve (M_numElements);

    typename MeshType::volume_Type* pv = 0;

    UInt stride = currentCount[1];
    if (! M_transposeInFile)
    {
        for (UInt j = 0; j < M_numElements; ++j)
        {
            pv = & (M_meshPartIn->addVolume() );
            pv->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[ (M_elementNodes + 2) * stride + j]) );
            pv->setId (M_uintBuffer[ (M_elementNodes + 1) * stride + j]);
            pv->setLocalId (j);
            for (UInt k = 0; k < M_elementNodes; ++k)
            {
                pv->setPoint (k, M_meshPartIn->point (
                                  M_uintBuffer[k * stride + j]) );
            }
            pv->setMarkerID (M_uintBuffer[M_elementNodes * stride + j]);
        }
    }
    else
    {
        for (UInt j = 0; j < M_numElements; ++j)
        {
            pv = & (M_meshPartIn->addVolume() );
            pv->replaceFlag (
                static_cast<flag_Type> (M_uintBuffer[stride * j + M_elementNodes + 2]) );
            pv->setId (M_uintBuffer[stride * j + M_elementNodes + 1]);
            pv->setLocalId (j);
            for (UInt k = 0; k < M_elementNodes; ++k)
            {
                pv->setPoint (k, M_meshPartIn->point (
                                  M_uintBuffer[stride * j + k]) );
            }
            pv->setMarkerID (M_uintBuffer[stride * j + M_elementNodes]);
        }
    }

    M_meshPartIn->updateElementEdges (false, false);
    M_meshPartIn->updateElementFaces (false, false);
}

} /* namespace LifeV */

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */
#endif /* PARTITION_IO_H_ */
