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

#include <mesh_reader_visit_ucd_multi.h>
#include <mesh_metadata_visit_ucd_multi.h>
#include <duplicate_node_elimination.h>

#include <typeinfo>
#include <unordered_map>
#include <cmath>

#include <DebugStream.h>
#include <avtGhostData.h>

#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

namespace vlsvplugin {

   VisitUCDMultiMeshReader::VisitUCDMultiMeshReader(): MeshReader() { 
      crds_node_x = NULL;
      crds_node_y = NULL;
      crds_node_z = NULL;
      nodeCoordinateArraysRead = false;
   }
   
   VisitUCDMultiMeshReader::~VisitUCDMultiMeshReader() { 
      delete [] crds_node_x; crds_node_x = NULL;
      delete [] crds_node_y; crds_node_y = NULL;
      delete [] crds_node_z; crds_node_z = NULL;
   }
   
   bool VisitUCDMultiMeshReader::readMesh(VLSVReader* vlsv,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDMultiMeshReader::readMesh called, domain: " << domain << endl;
      output = NULL;
      
      // Check that VLSVReader exists:
      if (vlsv == NULL) {
	 debug2 << "VLSV\t\t ERROR: VLSVReader is NULL" << endl;
	 return false;
      }
      
      // Check that metadata is not NULL:
      if (md == NULL) {
	 debug2 << "VLSV\t\t ERROR: MeshMetadata object is NULL" << endl;
	 return false;
      }
      
      // Check that given metadata is of correct type:
      VisitUCDMultiMeshMetadata* const metadata = dynamic_cast<VisitUCDMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
	 debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDMultiMeshMedata" << endl;
	 return false;
      }
      
      debug4 << "VLSV\t\t arraysize:  " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize: " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:   " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:   " << metadata->getDatatype() << endl;
            
      // Get domain offset arrays.
      // NOTE: Offsets are measured in number of (mesh) blocks, not cells:
      const uint64_t* domainOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfo(vlsv,domain,domainOffsets,ghostOffsets,variableOffsets) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
	 return false;
      }
      const uint64_t N_totalBlocks = domainOffsets[domain+1]-domainOffsets[domain];
      const uint64_t N_ghosts      = ghostOffsets[domain+1]-ghostOffsets[domain];
      const uint64_t N_blocks      = N_totalBlocks - N_ghosts;
      debug4 << "VLSV\t\t N_totalBlocks: " << N_totalBlocks << endl;
      debug4 << "VLSV\t\t N_ghosts:      " << N_ghosts << endl;
      debug4 << "VLSV\t\t N_blocks:      " << N_blocks << endl;
            
      // Get mesh bounding box:
      const uint64_t* const bbox = metadata->getMeshBoundingBox();
      if (bbox == NULL) {
	 debug2 << "VLSV\t\t ERROR: Failed to obtain mesh bounding box" << endl;
	 return false;
      }
      N_nodes_x = bbox[0]*bbox[3]+1;
      N_nodes_y = bbox[1]*bbox[4]+1;
      N_nodes_z = bbox[2]*bbox[5]+1;
      const uint64_t blockSize = bbox[3]*bbox[4]*bbox[5];
      debug4 << "VLSV\t\t N_nodes_(x,y,z): " << N_nodes_x << ' ' << N_nodes_y << ' ' << N_nodes_z << endl;
      
      // Read domain's blocks' global indices:
      uint64_t* blockGIDs = NULL;
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",metadata->getName()));
      if (vlsv->read("MESH",attribs,domainOffsets[domain],domainOffsets[domain+1]-domainOffsets[domain],blockGIDs) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read block global indices" << endl;
	 delete [] blockGIDs; blockGIDs = NULL;
	 return false;
      }

      // Read bounding box nodes' (x,y,z) coordinate arrays:
      if (readNodeCoordinateArrays(vlsv,metadata->getName()) == false) {
	 return false;
      }
      
      // For each block, attempt to insert its (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map.
      // The insertion will fail if the node already exists in the unordered_map, in which case 
      // counter is not increased. Counter, i.e. the value of unordered_map for given
      // NodeIndices, tells node's index.
      vtkIdType counter = 0;
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual> nodeIndices;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;

      // NOTE: only inserts inner blocks, not ghosts:
      for (uint64_t block=0; block<N_blocks; ++block) {
	 // Calculate block's (i,j,k) indices in the mesh bounding box:
	 uint64_t i_block = blockGIDs[block];
	 uint64_t k_block = i_block / (bbox[1]*bbox[0]);
	 i_block -= k_block*(bbox[1]*bbox[0]);
	 uint64_t j_block = i_block / bbox[0];
	 i_block -= j_block*bbox[0];

	 debug5 << "VLSV\t\t block GID: " << blockGIDs[block] << " indices: " << i_block << ' ' << j_block << ' ' << k_block << endl;
	 
	 // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
	 for (uint64_t k=0; k<bbox[5]+1; ++k) {
	    const uint64_t k_cell = k_block*bbox[5] + k;
	    for (uint64_t j=0; j<bbox[4]+1; ++j) {
	       const uint64_t j_cell = j_block*bbox[4] + j;
	       for (uint64_t i=0; i<bbox[3]+1; ++i) {
		  const uint64_t i_cell = i_block*bbox[3] + i;
		  result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,k_cell),counter));
		  
		  if (result.second == true)
		    debug5 << "VLSV\t\t\t node inserted: " << i_cell << ' ' << j_cell << ' ' << k_cell << " position: " << counter << endl;
		  
		  if (result.second == true) ++counter;
	       }
	    }
	 }
      }
      
      // unordered_map now contains all unique nodes. Its size is equal to
      // number of nodes in this domain:
      debug4 << "VLSV\t\t domain has " << nodeIndices.size() << " unique nodes" << endl;
      const size_t N_uniqueNodes = nodeIndices.size();
      
      // Create vtkPoints object and copy node coordinates to it:
      vtkPoints* coordinates = vtkPoints::New();
      coordinates->SetNumberOfPoints(N_uniqueNodes);
      float* pointer = reinterpret_cast<float*>(coordinates->GetVoidPointer(0));

      /*
      for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
	   it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
	 const vtkIdType position = it->second;
	 
	 debug5 << "VLSV\t\t Vertex " << position;
	 debug5 << " indices: ";
	 debug5 << it->first.i << ' ' << it->first.j << ' ' << it->first.k;
	 debug5 << " crds: ";
	 debug5 << crds_node_x[it->first.i] << ' ';
	 debug5 << crds_node_y[it->first.j] << ' ';
	 debug5 << crds_node_z[it->first.k] << ' ';
	 debug5 << endl;
	 
	 pointer[3*position+0] = crds_node_x[it->first.i];
	 pointer[3*position+1] = crds_node_y[it->first.j];
	 pointer[3*position+2] = crds_node_z[it->first.k];
      }*/
      
      const string& geometry = metadata->getMeshGeometry();
      if (geometry == VLSV::GEOM_CARTESIAN) {
	 // Cartesian geometry, no coordinate transformation:
	 for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
	      it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
	    const vtkIdType position = it->second;
	    pointer[3*position+0] = crds_node_x[it->first.i];
	    pointer[3*position+1] = crds_node_y[it->first.j];
	    pointer[3*position+2] = crds_node_z[it->first.k];
	 }
      } else if (geometry == VLSV::GEOM_CYLINDRICAL) {
	 // Cylindrical geometry, x' = r cos(phi) y' = r sin(phi) z' = z
	 for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
	      it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
	    const float R   = crds_node_x[it->first.i];
	    const float PHI = crds_node_y[it->first.j];
	    const float Z   = crds_node_z[it->first.k];
	    
	    const vtkIdType position = it->second;
	    pointer[3*position+0] = R*cos(PHI);
	    pointer[3*position+1] = R*sin(PHI);
	    pointer[3*position+2] = Z;
	    
	    //debug5 << "VLSV\t\t (R,PHI,Z): " << R << ' ' << PHI << ' ' << Z << "\t (x,y,z): ";
	    //debug5 << pointer[3*position+0] << ' ' << pointer[3*position+1] << ' ' << pointer[3*position+2] << endl;
	 }
      } else if (geometry == VLSV::GEOM_SPHERICAL) {
	 // Spherical geometry, x' = R sin(theta) cos(phi) y' = R sin(theta) sin(phi) z' = R cos(theta)
	 for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
	      it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
	    const float R     = crds_node_x[it->first.i];
	    const float THETA = crds_node_y[it->first.j];
	    const float PHI   = crds_node_z[it->first.k];
	    
	    const vtkIdType position = it->second;
	    pointer[3*position+0] = R*sin(THETA)*cos(PHI);
	    pointer[3*position+1] = R*sin(THETA)*sin(PHI);
	    pointer[3*position+2] = R*cos(THETA);
	 }
      }
      
      // Create vtkUnstructuredGrid:
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(coordinates);
      coordinates->Delete();
      ugrid->Allocate(N_blocks*blockSize); // FIXME (N_blocks)

      // Add all cells' connectivity information to vtkUnstructuredGrid:
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_blocks; ++block) { // FIXME (N_blocks)
	 // Calculate block's (i,j,k) indices in the mesh bounding box:
	 uint64_t i_block = blockGIDs[block];
	 uint64_t k_block = i_block / (bbox[1]*bbox[0]);
	 i_block -= k_block*(bbox[1]*bbox[0]);
	 uint64_t j_block = i_block / bbox[0];
	 i_block -= j_block*bbox[0];

	 //debug5 << "VLSV\t\t block GID: " << blockGIDs[block] << " indices: " << i_block << ' ' << j_block << ' ' << k_block << endl;
	 
	 // For each cell in the block, find indices of its eight nodes:
	 for (uint64_t k=0; k<bbox[5]; ++k) {
	    for (uint64_t j=0; j<bbox[4]; ++j) {
	       for (uint64_t i=0; i<bbox[3]; ++i) {
		  // Calculate cell's bounding box global indices:
		  const uint64_t i_cell = i_block*bbox[3] + i;
		  const uint64_t j_cell = j_block*bbox[4] + j;
		  const uint64_t k_cell = k_block*bbox[5] + k;
		  
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+1,k_cell+0));
		  vertices[0] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+0,k_cell+0));
		  vertices[1] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+0,k_cell+0));
		  vertices[2] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+1,k_cell+0));
		  vertices[3] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+1,k_cell+1));
		  vertices[4] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+0,k_cell+1));
		  vertices[5] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+0,k_cell+1));
		  vertices[6] = nodeIt->second;
		  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+1,k_cell+1));
		  vertices[7] = nodeIt->second;
		  
		  ugrid->InsertNextCell(cellType,8,vertices);
	       }
	    }
	 }
      }
      delete [] blockGIDs; blockGIDs = NULL;

      output = ugrid;
      return true;
   }

   bool VisitUCDMultiMeshReader::readNodeCoordinateArrays(VLSVReader* vlsv,const std::string& meshName) {
      // Check if coordinate arrays are already cached:
      if (nodeCoordinateArraysRead == true) return nodeCoordinateArraysRead;
      
      nodeCoordinateArraysRead = true;
      delete [] crds_node_x; crds_node_x = NULL;
      delete [] crds_node_y; crds_node_y = NULL;
      delete [] crds_node_z; crds_node_z = NULL;
      
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",meshName));
      
      // Check that node coordinate arrays have correct sizes:
      
      
      // Read node coordinate arrays:
      if (vlsv->read("MESH_NODE_CRDS_X",attribs,0,N_nodes_x,crds_node_x) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read node x-coordinate array" << endl;
	 nodeCoordinateArraysRead = false;
      }
      if (vlsv->read("MESH_NODE_CRDS_Y",attribs,0,N_nodes_y,crds_node_y) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read node y-coordinate array" << endl;
	 nodeCoordinateArraysRead = false;
      }
      if (vlsv->read("MESH_NODE_CRDS_Z",attribs,0,N_nodes_z,crds_node_z) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read node z-coordinate array" << endl;
	 nodeCoordinateArraysRead = false;
      }
      
      // If error occurred while reading data, delete coordinate arrays:
      if (nodeCoordinateArraysRead == false) {
	 delete [] crds_node_x; crds_node_x = NULL;
	 delete [] crds_node_y; crds_node_y = NULL;
	 delete [] crds_node_z; crds_node_z = NULL;
      }
      
      return nodeCoordinateArraysRead;
   }
   
   bool VisitUCDMultiMeshReader::readVariable(VLSVReader* vlsv,MeshMetadata* md,const VariableMetadata& vmd,int domain,float*& output) {
   
      return false;
   }
      
} // namespace vlsvplugin
