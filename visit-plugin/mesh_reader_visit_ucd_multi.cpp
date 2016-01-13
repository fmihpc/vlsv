/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2016 Finnish Meteorological Institute
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

#include <typeinfo>
#include <cmath>

#include <DebugStream.h>
#include <avtGhostData.h>

#include <vtkCellType.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>

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

   void VisitUCDMultiMeshReader::cartesianNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
						       vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_QUAD;
      vtkIdType vertices[4];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // For each cell in the block, find indices of its eight nodes:
         for (uint64_t j=0; j<bbox[4]; ++j) {
            for (uint64_t i=0; i<bbox[3]; ++i) {
               // Calculate cell's bounding box global indices:
               const uint64_t i_cell = i_block*bbox[3] + i;
               const uint64_t j_cell = j_block*bbox[4] + j;
               //const uint64_t k_cell = k_block*bbox[5] + k;
               
               nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+1,0));
               vertices[0] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+0,0));
               vertices[1] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+0,0));
               vertices[2] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+1,0));
               vertices[3] = nodeIt->second;
               
               ugrid->InsertNextCell(cellType,4,vertices);
            }
         }
      }
   }
   
   void VisitUCDMultiMeshReader::cartesianNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						     uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
						     vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
	                
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
   }

   void VisitUCDMultiMeshReader::cylindricalNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
                                                         uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
                                                         vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_QUAD;
      vtkIdType vertices[4];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting cylindrical cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // For each cell in the block, find indices of its four nodes:
         for (uint64_t j=0; j<bbox[4]; ++j) {      
            for (uint64_t i=0; i<bbox[3]; ++i) {
               // Calculate cell's bounding box global indices:
               const uint64_t i_cell = i_block*bbox[3] + i;
               const uint64_t j_cell = j_block*bbox[4] + j;
               
               uint64_t j_cell_plus_1 = j_cell + 1;
               if (yPeriodic == true) {
                  if (j_cell_plus_1 >= (N_nodes_y-1)) j_cell_plus_1 -= (N_nodes_y-1);
               }
               
               nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell_plus_1,0));
               vertices[0] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell       ,0));
               vertices[1] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell       ,0));
               vertices[2] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell_plus_1,0));
               vertices[3] = nodeIt->second;
               
               ugrid->InsertNextCell(cellType,4,vertices);
            }
         }
      }
   }

   void VisitUCDMultiMeshReader::cylindricalNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
                                                         uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
                                                         vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting cylindrical cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // For each cell in the block, find indices of its eight nodes:
         for (uint64_t k=0; k<bbox[5]; ++k) {
            for (uint64_t j=0; j<bbox[4]; ++j) {
               for (uint64_t i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  const uint64_t j_cell = j_block*bbox[4] + j;
                  const uint64_t k_cell = k_block*bbox[5] + k;
                  
                  uint64_t j_cell_plus_1 = j_cell + 1;
                  if (yPeriodic == true) {
                     if (j_cell_plus_1 >= (N_nodes_y-1)) j_cell_plus_1 -= (N_nodes_y-1);
                  }
                  
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell_plus_1,k_cell+0));
                  vertices[0] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell       ,k_cell+0));
                  vertices[1] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell       ,k_cell+0));
                  vertices[2] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell_plus_1,k_cell+0));
                  vertices[3] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell_plus_1,k_cell+1));
                  vertices[4] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell       ,k_cell+1));
                  vertices[5] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell       ,k_cell+1));
                  vertices[6] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell_plus_1,k_cell+1));
                  vertices[7] = nodeIt->second;
                  
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }
      }
   }

   void VisitUCDMultiMeshReader::sphericalNodeLookup(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						     uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
						     vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         for (uint64_t k=0; k<bbox[5]; ++k) {
            for (uint64_t j=0; j<bbox[4]; ++j) {
               for (uint64_t i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  const uint64_t j_cell = j_block*bbox[4] + j;
                  const uint64_t k_cell = k_block*bbox[5] + k;
                  
                  uint64_t k_cell_plus_1 = k_cell + 1;
             
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+1,k_cell       ));
                  vertices[0] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+0,k_cell       ));
                  vertices[1] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+0,k_cell       ));
                  vertices[2] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+1,k_cell       ));
                  vertices[3] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+1,k_cell_plus_1));
                  vertices[4] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+0,j_cell+0,k_cell_plus_1));
                  vertices[5] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+0,k_cell_plus_1));
                  vertices[6] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices(i_cell+1,j_cell+1,k_cell_plus_1));
                  vertices[7] = nodeIt->second;
                  
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }
      }
   }   
   
   void VisitUCDMultiMeshReader::insertCartesianNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         //uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         //i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t j=0; j<bbox[4]+1; ++j) {
            const uint64_t j_cell = j_block*bbox[4] + j;
            for (uint64_t i=0; i<bbox[3]+1; ++i) {
               const uint64_t i_cell = i_block*bbox[3] + i;
               result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,0),counter));
               if (result.second == true) ++counter;
            }
         }
      }
      
   }
   	  
   void VisitUCDMultiMeshReader::insertCartesianNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            const uint64_t k_cell = k_block*bbox[5] + k;
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               const uint64_t j_cell = j_block*bbox[4] + j;
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,k_cell),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
   
   void VisitUCDMultiMeshReader::insertCylindricalNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
                                                          uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t j=0; j<bbox[4]+1; ++j) {
            uint64_t j_cell = j_block*bbox[4] + j;
            
            if (yPeriodic == true) {
               if (j_cell >= (N_nodes_y-1)) j_cell -= (N_nodes_y-1);
            }
            
            for (uint64_t i=0; i<bbox[3]+1; ++i) {
               const uint64_t i_cell = i_block*bbox[3] + i;
               result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,0),counter));
               if (result.second == true) ++counter;
            }
         }
      }
   }
	 
   void VisitUCDMultiMeshReader::insertCylindricalNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            uint64_t k_cell = k_block*bbox[5] + k;
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               uint64_t j_cell = j_block*bbox[4] + j;
               
               if (yPeriodic == true) {
                  if (j_cell >= (N_nodes_y-1)) j_cell -= (N_nodes_y-1);
               }
               
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,k_cell),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
   
   void VisitUCDMultiMeshReader::insertSphericalNodes(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
                                                      uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];

         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            uint64_t k_cell = k_block*bbox[5] + k;
            
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               uint64_t j_cell = j_block*bbox[4] + j;
               
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell,j_cell,k_cell),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
      
   bool VisitUCDMultiMeshReader::readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDMultiMeshReader::readMesh called, domain: " << domain << endl;
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
      VisitUCDMultiMeshMetadata* const metadata = dynamic_cast<VisitUCDMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDMultiMeshMedata" << endl;
         return false;
      }
            
      // Get domain offset arrays.
      // NOTE: Offsets are measured in number of (mesh) blocks, not cells:
      const uint64_t* domainOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfo(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
         return false;
      }
      
      if (domainOffsets == NULL) {
         debug2 << "VLSV\t\t ERROR: domainOffsets is NULL" << endl; return false;
      }
      if (ghostOffsets == NULL) {
         debug2 << "VLSV\t\t ERROR: ghostOffsets is NULL" << endl; return false;
      }
      if (variableOffsets == NULL) {
         debug2 << "VLSV\t\t ERROR: variableOffsets is NULL" << endl; return false;
      }
      
      const uint64_t N_totalBlocks = domainOffsets[domain+1]-domainOffsets[domain];
      const uint64_t N_ghosts      = ghostOffsets[domain+1]-ghostOffsets[domain];
      const uint64_t N_blocks      = N_totalBlocks - N_ghosts;
      debug4 << "VLSV\t\t N_totalBlocks: " << N_totalBlocks << endl;
      debug4 << "VLSV\t\t N_ghosts:      " << N_ghosts << endl;
      debug4 << "VLSV\t\t N_blocks:      " << N_blocks << endl;
            
      // Get mesh bounding box:
      const vector<uint64_t>& bbox = metadata->getMeshBoundingBox();
      if (bbox.size() != 6) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain a valid mesh bounding box" << endl;
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
      if (vlsvReader->read("MESH",attribs,domainOffsets[domain],domainOffsets[domain+1]-domainOffsets[domain],blockGIDs) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read block global indices" << endl;
         delete [] blockGIDs; blockGIDs = NULL;
         return false;
      }

      // Read bounding box nodes' (x,y,z) coordinate arrays:
      if (readNodeCoordinateArrays(vlsvReader,metadata->getName()) == false) {
         return false;
      }
      
      // Get mesh geometry and periodicity:
      const vlsv::geometry::type& geometry = metadata->getMeshGeometry();
      metadata->getMeshPeriodicity(xPeriodic,yPeriodic,zPeriodic);

      // For each block, attempt to insert its (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map.
      // The insertion will fail if the node already exists in the unordered_map, in which case 
      // counter is not increased. Counter, i.e. the value of unordered_map for given
      // NodeIndices, tells node's index.
      //vtkIdType counter = 0;
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual> nodeIndices;

      switch (geometry) {
       case vlsv::geometry::UNKNOWN:
         insertCartesianNodes3D(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         break;
       case vlsv::geometry::CARTESIAN:
         if (metadata->getSpatialDimension() == 2) {
            insertCartesianNodes2D(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         } else {
            insertCartesianNodes3D(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         }
         break;
       case vlsv::geometry::CYLINDRICAL:
         if (metadata->getSpatialDimension() == 2) {
            insertCylindricalNodes2D(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         } else {
            insertCylindricalNodes3D(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         }
         break;
       case vlsv::geometry::SPHERICAL:
         insertSphericalNodes(nodeIndices,N_totalBlocks,blockGIDs,bbox);
         break;
       default:
         break;
      }
      
      // unordered_map now contains all unique nodes. Its size is equal to
      // number of nodes in this domain:
      debug4 << "VLSV\t\t domain has " << nodeIndices.size() << " unique nodes" << endl;
      const size_t N_uniqueNodes = nodeIndices.size();
      
      // Create vtkPoints object and copy node coordinates to it:
      vtkPoints* coordinates = vtkPoints::New();
      coordinates->SetNumberOfPoints(N_uniqueNodes);
      float* pointer = reinterpret_cast<float*>(coordinates->GetVoidPointer(0));

      // Get transformation matrix
      const double* transform = metadata->getTransform();
      float crds[3];
      
      switch (metadata->getMeshGeometry()) {
       case vlsv::geometry::CARTESIAN:
         // Cartesian geometry, no coordinate transformation:
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            const vtkIdType position = it->second;
            crds[0] = crds_node_x[it->first.i];
            crds[1] = crds_node_y[it->first.j];
            crds[2] = crds_node_z[it->first.k];
            
            pointer[3*position+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            pointer[3*position+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            pointer[3*position+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
         }
         break;
       case vlsv::geometry::CYLINDRICAL:
         // Cylindrical geometry, x' = r cos(phi) y' = r sin(phi) z' = z
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            const float R   = crds_node_x[it->first.i];
            const float PHI = crds_node_y[it->first.j];
            const float Z   = crds_node_z[it->first.k];
            
            crds[0] = R*cos(PHI);
            crds[1] = R*sin(PHI);
            crds[2] = Z;

            const vtkIdType position = it->second;
            pointer[3*position+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            pointer[3*position+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            pointer[3*position+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
         }
         break;
       case vlsv::geometry::SPHERICAL:
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            const float R     = crds_node_x[it->first.i];
            const float THETA = crds_node_y[it->first.j];
            const float PHI   = crds_node_z[it->first.k];

            const vtkIdType position = it->second;
            crds[0] = R*sin(THETA)*cos(PHI);
            crds[1] = R*sin(THETA)*sin(PHI);
            crds[2] = R*cos(THETA);

            pointer[3*position+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            pointer[3*position+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            pointer[3*position+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
         }
         break;
       default:
         return false;
         break;
      }
      
      // Create vtkUnstructuredGrid:
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(coordinates);
      coordinates->Delete();
      ugrid->Allocate(N_totalBlocks*blockSize);
      
      // Add all cells' connectivity information to vtkUnstructuredGrid:
      switch (geometry) {
       case vlsv::geometry::UNKNOWN:
         cartesianNodeLookup3D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         break;
       case vlsv::geometry::CARTESIAN:
         if (metadata->getSpatialDimension() == 2) {
            cartesianNodeLookup2D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         } else {
            cartesianNodeLookup3D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         }
         break;
       case vlsv::geometry::CYLINDRICAL:
         if (metadata->getSpatialDimension() == 2) {
            cylindricalNodeLookup2D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         } else {
            cylindricalNodeLookup3D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         }
         break;
       case vlsv::geometry::SPHERICAL:
         sphericalNodeLookup(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         break;
       default:
         cartesianNodeLookup3D(nodeIndices,N_totalBlocks,blockGIDs,bbox,ugrid);
         break;
      }
      delete [] blockGIDs; blockGIDs = NULL;

      // Determine correct values for real and ghost (internal to problem) zones:
      unsigned char cellIsReal = 0;
      unsigned char cellIsGhost = 0;
      avtGhostData::AddGhostZoneType(cellIsGhost,DUPLICATED_ZONE_INTERNAL_TO_PROBLEM);
      
      // Create an array that flags each zone either as real or internal ghost:
      vtkUnsignedCharArray* ghostZones = vtkUnsignedCharArray::New();
      ghostZones->SetName("avtGhostZones");
      ghostZones->Allocate(N_totalBlocks*blockSize);
      for (uint64_t i=0; i<N_blocks*blockSize; ++i) ghostZones->InsertNextValue(cellIsReal);
      for (uint64_t i=N_blocks*blockSize; i<N_totalBlocks*blockSize; ++i) ghostZones->InsertNextValue(cellIsGhost);

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

   bool VisitUCDMultiMeshReader::readNodeCoordinateArrays(vlsv::Reader* vlsvReader,const std::string& meshName) {
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
      if (vlsvReader->read("MESH_NODE_CRDS_X",attribs,0,N_nodes_x,crds_node_x) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node x-coordinate array" << endl;
         nodeCoordinateArraysRead = false;
      }
      if (vlsvReader->read("MESH_NODE_CRDS_Y",attribs,0,N_nodes_y,crds_node_y) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node y-coordinate array" << endl;
         nodeCoordinateArraysRead = false;
      }
      if (vlsvReader->read("MESH_NODE_CRDS_Z",attribs,0,N_nodes_z,crds_node_z) == false) {
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
   
   bool VisitUCDMultiMeshReader::readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDMultiMeshReader::readVariable called, domain: " << domain << endl;
      output = NULL;
      
      // Check that VLSVReader exists:
      if (vlsvReader == NULL) {
         debug3 << "VLSV\t\t ERROR: VLSVReader is NULL" << endl;
         return false;
      }
                
      // Check that metadata is not NULL:
      if (md == NULL) {
         debug3 << "VLSV\t\t ERROR: MeshMetadata object is NULL" << endl;
         return false;
      }
       
      // Check that given mesh metadata is of correct type: FIXME
      VisitUCDMultiMeshMetadata* const metadata = dynamic_cast<VisitUCDMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug3 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDMultiMeshMedata" << endl;
         return false;
      }
      
      // Get mesh bounding box:
      const std::vector<uint64_t>& bbox = metadata->getMeshBoundingBox();
      if (bbox.size() != 6) {
         debug3 << "VLSV\t\t ERROR: Mesh bounding box array is invalid" << endl;
         return false;
      }
      const uint64_t blockSize = bbox[3]*bbox[4]*bbox[5];
      
      // Get domain offset arrays:
      const uint64_t* blockOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfo(vlsvReader,domain,blockOffsets,ghostOffsets,variableOffsets) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
         return false;
      }
      const uint64_t N_totalBlocks = blockOffsets[domain+1] - blockOffsets[domain];
      const uint64_t N_ghosts      = ghostOffsets[domain+1] - ghostOffsets[domain];
      const uint64_t N_blocks      = N_totalBlocks - N_ghosts;
      const uint64_t components    = vmd.vectorSize;

      // Create vtkDataArray for variable data:
      bool success = true;
      //vtkFloatArray* rv = vtkFloatArray::New();
      vtkDoubleArray* rv = vtkDoubleArray::New();
      rv->SetNumberOfComponents(vmd.vectorSize);
      rv->SetNumberOfTuples(N_totalBlocks*blockSize);
      //float* variableData = rv->GetPointer(0);
      double* variableData = rv->GetPointer(0);
      
      // Read variable values from domain's real cells:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",vmd.name));
      attribs.push_back(make_pair("mesh",md->getName()));
      if (vlsvReader->read("VARIABLE",attribs,variableOffsets[domain]*blockSize,N_blocks*blockSize,variableData,false) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's real cell variable data" << endl;
         success = false;
      }
      
      // Read array that tell which domains contain ghost block data:
      list<pair<string,string> > meshAttribs;
      meshAttribs.push_back(make_pair("mesh",md->getName()));
      uint64_t* ghostDomains = NULL;
      if (vlsvReader->read("MESH_GHOST_DOMAINS",meshAttribs,ghostOffsets[domain],N_ghosts,ghostDomains,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_DOMAINS array" << endl;
         delete [] ghostDomains; ghostDomains = NULL;
         success = false;
      }
      
      // Read array that tells local IDs of ghost cell data in each domain:
      uint64_t* ghostLocalIDs = NULL;
      if (vlsvReader->read("MESH_GHOST_LOCALIDS",meshAttribs,ghostOffsets[domain],N_ghosts,ghostLocalIDs,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_LOCALIDS array" << endl;
         delete [] ghostDomains; ghostDomains = NULL;
         delete [] ghostLocalIDs; ghostLocalIDs = NULL;
         success = false;
      }
      
      // Read variable values for domain's ghost cells:
      if (success == true) {
         //float* ptr = variableData + N_blocks*blockSize*components;
         double* ptr = variableData + N_blocks*blockSize*components;
         for (uint64_t i=0; i<N_ghosts; ++i) {
            const uint64_t ghostDomainID    = ghostDomains[i];
            const uint64_t ghostValueOffset = (variableOffsets[ghostDomainID] + ghostLocalIDs[i])*blockSize;
            if (vlsvReader->read("VARIABLE",attribs,ghostValueOffset,blockSize,ptr,false) == false) {
               debug2 << "VLSV\t\t ERROR: Failed to read domain's ghost values" << endl;
               success = false;
               break;
            }
            ptr += blockSize*components;
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
