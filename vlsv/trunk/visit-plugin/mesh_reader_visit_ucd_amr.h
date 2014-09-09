/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2014 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MESH_READER_VISIT_UCD_AMR_H
#define MESH_READER_VISIT_UCD_AMR_H

#include <unordered_map>
#include <mesh_reader.h>
#include <duplicate_node_elimination.h>

#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>

namespace vlsvplugin {

   class VisitUCDAMRReader: public MeshReader {
    public:
      VisitUCDAMRReader();
      virtual ~VisitUCDAMRReader();
      
      virtual bool readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output);
      virtual bool readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output);
      
    protected:
      bool xPeriodic;
      bool yPeriodic;
      bool zPeriodic;
      float* crds_node_x;
      float* crds_node_y;
      float* crds_node_z;
      uint64_t maxRefinementLevel;
      bool nodeCoordinateArraysRead;
      uint64_t N_nodes_x;
      uint64_t N_nodes_y;
      uint64_t N_nodes_z;

      virtual void cartesianNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					 uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
					 vtkUnstructuredGrid* ugrid);
      virtual void cartesianNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
				       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
				       vtkUnstructuredGrid* ugrid);
      virtual void cylindricalNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					   uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
					   vtkUnstructuredGrid* ugrid);
      virtual void cylindricalNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					   uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
					   vtkUnstructuredGrid* ugrid);
      virtual void sphericalNodeLookup(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
				       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
				       vtkUnstructuredGrid* ugrid);
      
      virtual void insertCartesianNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox);
      virtual void insertCartesianNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox);
      virtual void insertCylindricalNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					    uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox);
      virtual void insertCylindricalNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					    uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox);
      virtual void insertSphericalNodes(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
					uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox);
      virtual bool readNodeCoordinateArrays(vlsv::Reader* vlsvReader,const std::string& meshName);
   };

} // namespace plugin

#endif
