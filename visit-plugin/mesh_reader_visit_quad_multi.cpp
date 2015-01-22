/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013 Finnish Meteorological Institute
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

#include <mesh_reader_visit_quad_multi.h>
#include <mesh_metadata_visit_quad_multi.h>
#include <duplicate_node_elimination.h>

#include <typeinfo>
#include <unordered_map>

#include <DebugStream.h>
#include <avtGhostData.h>

#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>

using namespace std;

namespace vlsvplugin {
   VisitQuadMultiMeshReader::VisitQuadMultiMeshReader(): MeshReader() { }
   
   VisitQuadMultiMeshReader::~VisitQuadMultiMeshReader() { }
   
   bool VisitQuadMultiMeshReader::readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitQuadMultiMeshReader::readMesh called, domain: " << domain << endl;
      output = NULL;
      
      // Check that VLSVReader exists:
      if (vlsvReader == NULL) {
         debug2 << "VLSV\t\t ERROR: VLSVReader is NULL" << endl;
         return false;
      }
      
      // Check that metadata is not NULL:
      if (md == NULL) {
         debug2 << "VLSV\t\t ERROR: MeshMetadata object is NULL" << endl;
         return false;
      }
      
      // Check that given metadata is of correct type:
      VisitQuadMultiMeshMetadata* const metadata = dynamic_cast<VisitQuadMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitQuadMultiMeshMedata" << endl;
         return false;
      }
      
      debug4 << "VLSV\t\t arraysize:  " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize: " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:   " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:   " << metadata->getDatatype() << endl;

      // Get domain offset arrays:
      const uint64_t* cellOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;      
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfo(vlsvReader,domain,cellOffsets,ghostOffsets,variableOffsets) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
         return false;
      }
      const uint64_t N_cells = cellOffsets[domain+1]-cellOffsets[domain];
      debug4 << "VLSV\t\t N_cells: " << N_cells << endl;
      
      // Get mesh bounding box:
      const float* const bbox = metadata->getMeshBoundingBox();
      if (bbox == NULL) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain mesh bounding box" << endl;
         return false;
      }
      
      // Read domain's cells' (i,j,k) indices:
      int32_t* cells = NULL;
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",metadata->getName()));
      if (vlsvReader->read("MESH",attribs,cellOffsets[domain],cellOffsets[domain+1]-cellOffsets[domain],cells) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read cell indices" << endl;
         delete [] cells; cells = NULL;
         return false;
      }
      
      // For each cell, attempt to insert its eight nodes to unordered_map. 
      // The insertion will fail if the node already exists in the map, in which case 
      // counter is not increased. Counter, i.e. the value of unordered_map for given
      // NodeIndices, tells node's index.
      vtkIdType counter = 0;
      int32_t* ptr = cells;
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual> nodeIndices;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t i=0; i<N_cells; ++i) {
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+0,ptr[1]+1,ptr[2]+0),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+0,ptr[1]+0,ptr[2]+0),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+1,ptr[1]+0,ptr[2]+0),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+1,ptr[1]+1,ptr[2]+0),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+0,ptr[1]+1,ptr[2]+1),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+0,ptr[1]+0,ptr[2]+1),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+1,ptr[1]+0,ptr[2]+1),counter));
         if (result.second == true) ++counter;
         result = nodeIndices.insert(make_pair(NodeIndices(ptr[0]+1,ptr[1]+1,ptr[2]+1),counter));
         if (result.second == true) ++counter;
         ptr += 3;
      }
      
      // unordered_map now contains all unique nodes. Its size is equal to number
      // of nodes in this mesh piece:
      debug4 << "VLSV\t\t domain has " << nodeIndices.size() << " unique nodes" << endl;
      const size_t N_uniqueNodes = nodeIndices.size();

      // Create vtkPoints object and copy node coordinates to it:
      vtkPoints* coordinates = vtkPoints::New();
      coordinates->SetNumberOfPoints(N_uniqueNodes);
      float* pointer = reinterpret_cast<float*>(coordinates->GetVoidPointer(0));

      // Get the transformation matrix
      const double* transform = metadata->getTransform();
      
      for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator 
	   it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
         const vtkIdType position = it->second;
         if (3*position+2 >= N_uniqueNodes*3) {
            cerr << "position exceeds array size!" << endl;
            exit(1);
         }
         
         float crds[3];
         crds[0] = bbox[0] + it->first.i * bbox[3];
         crds[1] = bbox[1] + it->first.j * bbox[4];
         crds[2] = bbox[2] + it->first.k * bbox[5];
         
         pointer[3*position+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
         pointer[3*position+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
         pointer[3*position+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
      }

      // Create vtkUnstructuredGrid:
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(coordinates);      
      coordinates->Delete();
      ugrid->Allocate(N_cells);

      // Add all cells' connectivity information to vtkUnstructuredGrid:
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      ptr = cells;
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      for (uint64_t i=0; i<N_cells; ++i) {
         // Find indices of all eight nodes:
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+0,ptr[1]+1,ptr[2]+0));
         vertices[0] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+0,ptr[1]+0,ptr[2]+0));
         vertices[1] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+1,ptr[1]+0,ptr[2]+0));
         vertices[2] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+1,ptr[1]+1,ptr[2]+0));
         vertices[3] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+0,ptr[1]+1,ptr[2]+1));
         vertices[4] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+0,ptr[1]+0,ptr[2]+1));
         vertices[5] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+1,ptr[1]+0,ptr[2]+1));
         vertices[6] = nodeIt->second;
         nodeIt = nodeIndices.find(NodeIndices(ptr[0]+1,ptr[1]+1,ptr[2]+1));
         vertices[7] = nodeIt->second;
         ptr += 3;
         
         // Insert cell:
         ugrid->InsertNextCell(cellType,8,vertices);
      }
      delete [] cells; cells = NULL;

      // Determine correct values for real and ghost (internal to problem) zones:
      unsigned char cellIsReal = 0;
      unsigned char cellIsGhost = 0;
      avtGhostData::AddGhostZoneType(cellIsGhost,DUPLICATED_ZONE_INTERNAL_TO_PROBLEM);
      
      // Create an array that flags each zone either as real or internal ghost:
      vtkUnsignedCharArray* ghostZones = vtkUnsignedCharArray::New();
      ghostZones->SetName("avtGhostZones");
      ghostZones->Allocate(N_cells);      
      const uint64_t N_ghosts = ghostOffsets[domain+1]-ghostOffsets[domain];
      for (uint64_t i=0; i<N_cells-N_ghosts; ++i) ghostZones->InsertNextValue(cellIsReal);
      for (uint64_t i=N_cells-N_ghosts; i<N_cells; ++i) ghostZones->InsertNextValue(cellIsGhost);
      
      // Copy ghost cell information to vtkUnstructuredGrid:
      ugrid->GetCellData()->AddArray(ghostZones);
      #ifdef NEW_VTK_API
         vtkStreamingDemandDrivenPipeline::SetUpdateGhostLevel(ugrid->GetInformation(), 0);
      #else
         ugrid->SetUpdateGhostLevel(0);
      #endif
      ghostZones->Delete();
	
      output = ugrid;
      return true;
   }
   
   bool VisitQuadMultiMeshReader::readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output) {
      debug2 << "VLSV\t VisitQuadMultiMeshReader::readVariable called, domain: " << domain << endl;
      output = NULL;

      // Check that VLSVReader exists:
      if (vlsvReader == NULL) {
         debug2 << "VLSV\t\t ERROR: VLSVReader is NULL" << endl;
         return false;
      }

      // Check that metadata is not NULL:
      if (md == NULL) {
         debug2 << "VLSV\t\t ERROR: MeshMetadata object is NULL" << endl;
         return false;
      }
      
      // Check that given mesh metadata is of correct type:
      VisitQuadMultiMeshMetadata* const metadata = dynamic_cast<VisitQuadMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitQuadMultiMeshMedata" << endl;
         return false;
      }
    
      // Get domain offset arrays:
      const uint64_t* cellOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfo(vlsvReader,domain,cellOffsets,ghostOffsets,variableOffsets) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
         return false;
      }
      const uint64_t N_totalCells = cellOffsets[domain+1] - cellOffsets[domain];
      const uint64_t N_ghosts     = ghostOffsets[domain+1] - ghostOffsets[domain];
      const uint64_t N_cells      = N_totalCells - N_ghosts;
      const uint64_t components   = vmd.vectorSize;

      // Create a vtkFloatArray for variable data:
      bool success = true;
      vtkFloatArray* rv = vtkFloatArray::New();
      rv->SetNumberOfComponents(vmd.vectorSize);
      rv->SetNumberOfTuples(N_totalCells);
      float* variableData = rv->GetPointer(0);
      
      // Read variable values from domain's real cells:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",vmd.name));
      attribs.push_back(make_pair("mesh",md->getName()));
      if (vlsvReader->read("VARIABLE",attribs,variableOffsets[domain],N_cells,variableData,false) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's real cell variable data" << endl;
         success = false;
      }

      // Read array that tell which domains contain ghost cell data:
      list<pair<string,string> > meshAttribs;
      meshAttribs.push_back(make_pair("mesh",md->getName()));
      uint64_t* ghostDomains = NULL;
      if (vlsvReader->read("MESH_GHOST_DOMAINS",meshAttribs,ghostOffsets[domain],N_ghosts,ghostDomains,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_DOMAINS array" << endl;
         success = false;
      }
      
      // Read array that tells local IDs of ghost cell data in each domain:
      uint64_t* ghostLocalIDs = NULL;
      if (vlsvReader->read("MESH_GHOST_LOCALIDS",meshAttribs,ghostOffsets[domain],N_ghosts,ghostLocalIDs,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_LOCALIDS array" << endl;
         success = false;
      }
      
      // Read variable values for domain's ghost cells:
      if (success == true) {
         float* ptr = variableData + N_cells*components;
         for (uint64_t i=0; i<N_ghosts; ++i) {
            const uint64_t ghostDomainID    = ghostDomains[i];
            const uint64_t ghostValueOffset = variableOffsets[ghostDomainID] + ghostLocalIDs[i];	 
            if (vlsvReader->read("VARIABLE",attribs,ghostValueOffset,1,ptr,false) == false) {
               debug2 << "VLSV\t\t ERROR: Failed to read domain's ghost values" << endl;
               success = false;
               break;
            }
            ptr += components;
         }
      }
      
      delete [] ghostDomains; ghostDomains = NULL;
      delete [] ghostLocalIDs; ghostLocalIDs = NULL;
      
      if (success == false) {
         rv->Delete();
         output = NULL;
         return false;
      } else {
         output = rv;
         return true;
      }
   }
   
} // namespace vlsvplugin


