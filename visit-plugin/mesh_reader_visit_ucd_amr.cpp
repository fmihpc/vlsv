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

#include <mesh_reader_visit_ucd_amr.h>
#include <mesh_metadata_visit_ucd_amr.h>

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

#include <vlsv_amr.h>

using namespace std;

namespace vlsvplugin {

   VisitUCDAMRReader::VisitUCDAMRReader(): MeshReader() { 
      crds_node_x = NULL;
      crds_node_y = NULL;
      crds_node_z = NULL;
      nodeCoordinateArraysRead = false;
   }
   
   VisitUCDAMRReader::~VisitUCDAMRReader() { 
      delete [] crds_node_x; crds_node_x = NULL;
      delete [] crds_node_y; crds_node_y = NULL;
      delete [] crds_node_z; crds_node_z = NULL;
   }

   void VisitUCDAMRReader::cartesianNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
						       vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_QUAD;
      vtkIdType vertices[4];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t mul = pow(2,maxRefinementLevel-refLevel);
         
         // For each cell in the block, find indices of its eight nodes:
         for (uint64_t j=0; j<bbox[4]; ++j) {
            for (uint64_t i=0; i<bbox[3]; ++i) {
               // Calculate cell's bounding box global indices:
               const uint64_t i_cell = i_block*bbox[3] + i;
               const uint64_t j_cell = j_block*bbox[4] + j;
               
               nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+1)*mul,0));
               vertices[0] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+0)*mul,0));
               vertices[1] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+0)*mul,0));
               vertices[2] = nodeIt->second;
               nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+1)*mul,0));
               vertices[3] = nodeIt->second;
               
               ugrid->InsertNextCell(cellType,4,vertices);
            }
         }
      }
   }
   
   void VisitUCDAMRReader::cartesianNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						     uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
						     vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t mul = pow(2,maxRefinementLevel-refLevel);
         
         // For each cell in the block, find indices of its eight nodes:
         for (uint64_t k=0; k<bbox[5]; ++k) {
            for (uint64_t j=0; j<bbox[4]; ++j) {
               for (uint64_t i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  const uint64_t j_cell = j_block*bbox[4] + j;
                  const uint64_t k_cell = k_block*bbox[5] + k;
                  
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+1)*mul,(k_cell+0)*mul));
                  vertices[0] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+0)*mul,(k_cell+0)*mul));
                  vertices[1] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+0)*mul,(k_cell+0)*mul));
                  vertices[2] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+1)*mul,(k_cell+0)*mul));
                  vertices[3] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+1)*mul,(k_cell+1)*mul));
                  vertices[4] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+0)*mul,(k_cell+1)*mul));
                  vertices[5] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+0)*mul,(k_cell+1)*mul));
                  vertices[6] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+1)*mul,(k_cell+1)*mul));
                  vertices[7] = nodeIt->second;
                  
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }
      }
   }

   void VisitUCDAMRReader::cylindricalNodeLookup2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							 uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
							 vtkUnstructuredGrid* ugrid) {
//      #warning FIXME Cylindrical node lookup 2D
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
               
               // If mesh is periodic in phi coordinate, set a
               // node with phi=(2*pi) to be equal to a node with phi=0:
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
   
   void VisitUCDAMRReader::cylindricalNodeLookup3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							 uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
							 vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting cylindrical cells to unstructured AMR mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t mul = pow(2,maxRefinementLevel-refLevel);
         
         // For each cell in the block, find indices of its eight nodes:
         for (uint64_t k=0; k<bbox[5]; ++k) {
            for (uint64_t j=0; j<bbox[4]; ++j) {
               for (uint64_t i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  const uint64_t j_cell = j_block*bbox[4] + j;
                  const uint64_t k_cell = k_block*bbox[5] + k;
                  
                  // If mesh is periodic in phi coordinate, set a
                  // node with phi=(2*pi) to be equal to a node with phi=0:
                  uint64_t j_cell_plus_1 = j_cell + 1;
                  if (yPeriodic == true) {
                     if (j_cell_plus_1 >= (N_nodes_y-1)*(refLevel+1)) j_cell_plus_1 -= (N_nodes_y-1)*(refLevel+1);
                  }
                  
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell_plus_1)*mul,(k_cell+0)*mul));
                  vertices[0] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell       )*mul,(k_cell+0)*mul));
                  vertices[1] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell       )*mul,(k_cell+0)*mul));
                  vertices[2] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell_plus_1)*mul,(k_cell+0)*mul));
                  vertices[3] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell_plus_1)*mul,(k_cell+1)*mul));
                  vertices[4] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell       )*mul,(k_cell+1)*mul));
                  vertices[5] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell       )*mul,(k_cell+1)*mul));
                  vertices[6] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell_plus_1)*mul,(k_cell+1)*mul));
                  vertices[7] = nodeIt->second;
                  
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }
      }
   }
   
   void VisitUCDAMRReader::sphericalNodeLookup(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						     uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
						     vtkUnstructuredGrid* ugrid) {
      const int cellType = VTK_HEXAHEDRON;
      vtkIdType vertices[8];
      unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator nodeIt;
      debug5 << "VLSV\t\t Inserting Cartesian cells to unstructured mesh:" << endl;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t mul = pow(2,maxRefinementLevel-refLevel);
         
         for (uint64_t k=0; k<bbox[5]; ++k) {
            const uint64_t k_cell = k_block*bbox[5] + k;
            uint64_t k_cell_plus_1 = k_cell + 1;
	    
            // If mesh is periodic in phi coordinate, set a
            // node with phi=(2*pi) to be equal to a node with phi=0:
            if (zPeriodic == true) {
               if (k_cell_plus_1 >= (N_nodes_z-1)*(refLevel+1)) {
                  k_cell_plus_1 -= (N_nodes_z-1)*(refLevel+1);
               }
            }
            
            for (uint64_t j=0; j<bbox[4]; ++j) {
               for (uint64_t i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  const uint64_t j_cell = j_block*bbox[4] + j;
                  
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+1)*mul,(k_cell       )*mul));
                  vertices[0] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+0)*mul,(k_cell       )*mul));
                  vertices[1] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+0)*mul,(k_cell       )*mul));
                  vertices[2] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+1)*mul,(k_cell       )*mul));
                  vertices[3] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+1)*mul,(k_cell_plus_1)*mul));
                  vertices[4] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+0)*mul,(j_cell+0)*mul,(k_cell_plus_1)*mul));
                  vertices[5] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+0)*mul,(k_cell_plus_1)*mul));
                  vertices[6] = nodeIt->second;
                  nodeIt = nodeIndices.find(NodeIndices((i_cell+1)*mul,(j_cell+1)*mul,(k_cell_plus_1)*mul));
                  vertices[7] = nodeIt->second;
                  
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }
      }
   }   
   
   void VisitUCDAMRReader::insertCartesianNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t multiplier = pow(2,maxRefinementLevel-refLevel);
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t j=0; j<bbox[4]+1; ++j) {
            const uint64_t j_cell = j_block*bbox[4] + j;
            for (uint64_t i=0; i<bbox[3]+1; ++i) {
               const uint64_t i_cell = i_block*bbox[3] + i;
               result = nodeIndices.insert(make_pair(NodeIndices(i_cell*multiplier,j_cell*multiplier,0),counter));
               if (result.second == true) ++counter;
            }
         }
      }
      
   }
   	  
   void VisitUCDAMRReader::insertCartesianNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;

      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t multiplier = pow(2,maxRefinementLevel-refLevel);
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            const uint64_t k_cell = k_block*bbox[5] + k;
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               const uint64_t j_cell = j_block*bbox[4] + j;
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell*multiplier,j_cell*multiplier,k_cell*multiplier),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
   
   void VisitUCDAMRReader::insertCylindricalNodes2D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox) {
//      #warning FIXME Cylindrical AMR node insert 2D
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         uint64_t j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t j=0; j<bbox[4]+1; ++j) {
            uint64_t j_cell = j_block*bbox[4] + j;
            
            // If mesh is periodic in phi coordinate, set a
            // node with phi=(2*pi) to be equal to a node with phi=0:
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
   
   void VisitUCDAMRReader::insertCylindricalNodes3D(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
							  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t multiplier = pow(2,maxRefinementLevel-refLevel);
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            uint64_t k_cell = k_block*bbox[5] + k;
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               uint64_t j_cell = j_block*bbox[4] + j;
               
               // If mesh is periodic in phi coordinate, set a
               // node with phi=(2*pi) to be equal to a node with phi=0:
               if (yPeriodic == true) {
                  if (j_cell >= (N_nodes_y-1)*(refLevel+1)) j_cell -= (N_nodes_y-1)*(refLevel+1);
               }
               
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell*multiplier,j_cell*multiplier,k_cell*multiplier),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
   
   void VisitUCDAMRReader::insertSphericalNodes(std::unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>& nodeIndices,
						      uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox) {
      vtkIdType counter = 0;
      pair<unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::iterator,bool> result;
      
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint32_t refLevel,i_block,j_block,k_block;
         vlsv::calculateCellIndices(blockGIDs[block],refLevel,i_block,j_block,k_block);
         const uint32_t multiplier = pow(2,maxRefinementLevel-refLevel);
         
         // Attempt to insert all (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map:
         for (uint64_t k=0; k<bbox[5]+1; ++k) {
            uint64_t k_cell = k_block*bbox[5] + k;
            
            // If mesh is periodic in phi coordinate, set a
            // node with phi=(2*pi) to be equal to a node with phi=0:
            if (zPeriodic == true) {
               if (k_cell >= (N_nodes_z-1)*(refLevel+1)) {
                  k_cell -= (N_nodes_z-1)*(refLevel+1);
               }
            }
            
            for (uint64_t j=0; j<bbox[4]+1; ++j) {
               uint64_t j_cell = j_block*bbox[4] + j;
               
               for (uint64_t i=0; i<bbox[3]+1; ++i) {
                  const uint64_t i_cell = i_block*bbox[3] + i;
                  result = nodeIndices.insert(make_pair(NodeIndices(i_cell*multiplier,j_cell*multiplier,k_cell*multiplier),counter));
                  if (result.second == true) ++counter;
               }
            }
         }
      }
      
   }
      
   bool VisitUCDAMRReader::readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDAMRReader::readMesh called, domain: " << domain << endl;
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
      VisitUCDAMRMetadata* const metadata = dynamic_cast<VisitUCDAMRMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDAMRMedata" << endl;
         return false;
      }
      
      maxRefinementLevel = metadata->getMaximumRefinementLevel();
      debug4 << "VLSV\t\t arraysize:            " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize:           " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:             " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:             " << metadata->getDatatype() << endl;
      debug4 << "VLSV\t\t max refinement level: " << maxRefinementLevel << endl;

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
      
      // Init AMR mesh:
      if (vlsv::initMesh(bbox[0],bbox[1],bbox[2],maxRefinementLevel) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to init AMR mesh" << endl;
         return false;
      }

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

      const uint32_t multiplier = pow(2,maxRefinementLevel);
      const double dx_node = (crds_node_x[1] - crds_node_x[0])/multiplier;
      const double dy_node = (crds_node_y[1] - crds_node_y[0])/multiplier;
      const double dz_node = (crds_node_z[1] - crds_node_z[0])/multiplier;
      uint64_t i_node_base,j_node_base,k_node_base;

      const double* transform = metadata->getTransform();

      switch (metadata->getMeshGeometry()) {
       case vlsv::geometry::CARTESIAN:
         // Cartesian geometry, no coordinate transformation:
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            const vtkIdType position = it->second;
            
            i_node_base = it->first.i / multiplier;
            j_node_base = it->first.j / multiplier;
            k_node_base = it->first.k / multiplier;

            // Node xyz coordinates
            float crds[3];
            crds[0] = crds_node_x[i_node_base] + (it->first.i-i_node_base*multiplier)*dx_node;
            crds[1] = crds_node_y[j_node_base] + (it->first.j-j_node_base*multiplier)*dy_node;
            crds[2] = crds_node_z[k_node_base] + (it->first.k-k_node_base*multiplier)*dz_node;

            pointer[3*position+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            pointer[3*position+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            pointer[3*position+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
         }
         break;
       case vlsv::geometry::CYLINDRICAL:
         // Cylindrical geometry, x' = r cos(phi) y' = r sin(phi) z' = z
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            i_node_base = it->first.i / multiplier;
            j_node_base = it->first.j / multiplier;
            k_node_base = it->first.k / multiplier;
            
            const float R   = crds_node_x[i_node_base] + (it->first.i-i_node_base*multiplier)*dx_node;
            const float PHI = crds_node_y[j_node_base] + (it->first.j-j_node_base*multiplier)*dy_node;
            const float Z   = crds_node_z[k_node_base] + (it->first.k-k_node_base*multiplier)*dz_node;

            float crds[3];
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
         // Spherical geometry, x' = r sin(theta) cos(phi), y' = r sin(theta) sin(phi), z' = r cos(theta)
         for (unordered_map<NodeIndices,vtkIdType,NodeHash,NodesAreEqual>::const_iterator
              it=nodeIndices.begin(); it!=nodeIndices.end(); ++it) {
            i_node_base = it->first.i / multiplier;
            j_node_base = it->first.j / multiplier;
            k_node_base = it->first.k / multiplier;
            
            const float R     = crds_node_x[i_node_base] + (it->first.i-i_node_base*multiplier)*dx_node;
            const float THETA = crds_node_y[j_node_base] + (it->first.j-j_node_base*multiplier)*dy_node;
            const float PHI   = crds_node_z[k_node_base] + (it->first.k-k_node_base*multiplier)*dz_node;

            float crds[3];
            crds[0] = R*sin(THETA)*cos(PHI);
            crds[1] = R*sin(THETA)*sin(PHI);
            crds[2] = R*cos(THETA);
            
            const vtkIdType position = it->second;
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

   bool VisitUCDAMRReader::readNodeCoordinateArrays(vlsv::Reader* vlsvReader,const std::string& meshName) {
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
   
   bool VisitUCDAMRReader::readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDAMRReader::readVariable called, domain: " << domain << endl;
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
      VisitUCDAMRMetadata* const metadata = dynamic_cast<VisitUCDAMRMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug3 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDAMRMedata" << endl;
         return false;
      }
      
      // Get mesh bounding box:
      const uint64_t* bbox = metadata->getMeshBoundingBox();
      if (bbox == NULL) {
         debug3 << "VLSV\t\t ERROR: Mesh bounding box array is NULL" << endl;
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
      vtkDoubleArray* rv = vtkDoubleArray::New();
      rv->SetNumberOfComponents(vmd.vectorSize);
      rv->SetNumberOfTuples(N_totalBlocks*blockSize);
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
