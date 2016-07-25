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

#include <mesh_reader_visit_ucd_generic_multi.h>
#include <mesh_metadata_visit_ucd_generic_multi.h>

#include <typeinfo>
#include <cmath>

#include <DebugStream.h>
#include <avtGhostData.h>

#include <vtkCellType.h>
#include <vtkCell.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>

using namespace std;

namespace vlsvplugin {

   VisitUCDGenericMultiMeshReader::VisitUCDGenericMultiMeshReader(): MeshReader() {
      nodeCoordinateArraysRead = false;
   }
   
   VisitUCDGenericMultiMeshReader::~VisitUCDGenericMultiMeshReader() { }

   bool VisitUCDGenericMultiMeshReader::readCellVariable(vlsv::Reader* vlsvReader,VisitUCDGenericMultiMeshMetadata*  metadata,
							 const VariableMetadata& vmd,int domain,void*& output) {
      // Get mesh bounding box:
      const uint64_t* bbox = metadata->getMeshBoundingBox();
      if (bbox == NULL) {
         debug3 << "VLSV\t\t ERROR: Mesh bounding box array is NULL" << endl;
         return false;
      }
      const uint64_t blockSize =
        bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
        * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
        * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      
      // Get domain offset arrays:
      const uint64_t* blockOffsets = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (metadata->getDomainInfoZones(vlsvReader,domain,blockOffsets,ghostOffsets,variableOffsets) == false) {
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
      attribs.push_back(make_pair("mesh",metadata->getName()));
      if (vlsvReader->read("VARIABLE",attribs,variableOffsets[domain]*blockSize,N_blocks*blockSize,variableData,false) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read domain's real cell variable data" << endl;
         success = false;
      }
      
      // Read array that tells which domains contain ghost block data:
      list<pair<string,string> > meshAttribs;
      meshAttribs.push_back(make_pair("mesh",metadata->getName()));
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

      // Deallocate memory and exit:
      delete [] ghostDomains; ghostDomains = NULL;
      delete [] ghostLocalIDs; ghostLocalIDs = NULL;
      
      if (success == false) {
         rv->Delete();
         output = NULL;
      } else {
         output = rv;
      }
      
      return success;
   }
   
   bool VisitUCDGenericMultiMeshReader::readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDGenericMultiMeshReader::readMesh called, domain: " << domain << endl;
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
      VisitUCDGenericMultiMeshMetadata* const metadata = dynamic_cast<VisitUCDGenericMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDGenericMultiMeshMedata" << endl;
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
      if (metadata->getDomainInfoZones(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets) == false) {
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
      const uint64_t blockSize =
          bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
        * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
        * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      
      // Create vtkPoints object and read node coordinates to it:
      const size_t N_uniqueNodes = metadata->getNumberOfNodes(domain);
      debug4 << "VLSV\t\t Mesh has " << N_uniqueNodes << " unique nodes" << endl;
      debug5 << "VLSV\t\t Node coordinates datatype: " << metadata->getNodeDatatype() << " datasize: " << metadata->getNodeDataSize() << endl;
      
      const int vtkDatatype = getVtkDatatype(metadata->getNodeDatatype(),metadata->getNodeDataSize());
      debug5 << "VLSV\t\t VTK datatype is: " << vtkDatatype << endl;

      vtkPoints* coordinates = vtkPoints::New(vtkDatatype);
      if (readNodeCoordinateArrays(vlsvReader,metadata->getName(),metadata->getNodeOffset(domain),metadata->getNumberOfNodes(domain),coordinates) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node coordinates" << endl;
      }
      debug5 << "VLSV\t Read node coordinates, offset: " << metadata->getNodeOffset(domain) << " nodes: " << metadata->getNumberOfNodes(domain) << endl;

      // Create vtkUnstructuredGrid:
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(coordinates);
      coordinates->Delete();
      ugrid->Allocate(N_totalBlocks*blockSize);
      debug5 << "VLSV\t\t Allocated a VTK unstructured mesh with " << N_totalBlocks*blockSize << " connectivity array size" << endl;

      // Read cell connectivity from VLSV file:
      list<pair<string,string> > xmlAttributes;
      xmlAttributes.push_back(make_pair("name",metadata->getName()));
      const uint32_t connectivityOffset = metadata->getZoneOffset(domain);
      const uint32_t connectivitySize = metadata->getCellConnectivitySize(domain);
      debug4 << "VLSV\t\t offset to connectivity array: " << connectivityOffset << " size: " << connectivitySize << endl;

      vtkIdType* cellConnectivity = NULL;
      if (vlsvReader->read("MESH",xmlAttributes,connectivityOffset,connectivitySize,cellConnectivity,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read cell connectivity" << endl;
      }

      size_t counter = 0;
      while (counter < connectivitySize) {
         // First value indices VTK cell type:
         const uint32_t cellType = cellConnectivity[counter];
         const int vtkCelltype = vlsvplugin::getVtkCelltype(cellType);
         ++counter;
	 
         // Size of entry is determined by cell type:
         const int entrySize = cellConnectivity[counter];
         ++counter;

         ugrid->InsertNextCell(vtkCelltype,entrySize,cellConnectivity+counter);
         
         debug5 << "VLSV\t\t inserted cell of type: " << vtkCelltype << " entry size: " << entrySize << " nodes: ";
         for (int i=0; i<entrySize; ++i) debug5 << cellConnectivity[counter+i] << ' ';
         debug5 << endl;

         counter += entrySize;
      }
      delete [] cellConnectivity; cellConnectivity = NULL;

      
      debug5 << "VLSV\t\t ugrid reports " << ugrid->GetNumberOfCells() << " cells and " << ugrid->GetNumberOfPoints() << " nodes" << endl;

      double x[3];
      for (vtkIdType p=0; p<ugrid->GetNumberOfPoints(); ++p) {
         ugrid->GetPoint(p,x);
         debug5 << "VLSV\t\t point " << p << ":\t " << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
      }
      
      vtkCell* cellPtr;
      for (vtkIdType c=0; c<ugrid->GetNumberOfCells(); ++c) {
         cellPtr = ugrid->GetCell(c);
         debug5 << "VLSV\t\t cell " << c << " points: " << cellPtr->GetNumberOfPoints() << " ids: ";         
         for (int i=0; i<cellPtr->GetNumberOfPoints(); ++i) {
            debug5 << cellPtr->GetPointId(i) << ' ';
         }
         debug5 << endl;
      }

      output = ugrid;
      return true;
   }

   bool VisitUCDGenericMultiMeshReader::readNodeCoordinateArrays(vlsv::Reader* vlsvReader,const std::string& meshName,
								 uint64_t offset,uint64_t amount,vtkPoints* points) {
      // Check if coordinate arrays are already cached:
      if (nodeCoordinateArraysRead == true) return nodeCoordinateArraysRead;
      nodeCoordinateArraysRead = true;
      
      debug4 << "VLSV\t\t attempting to read " << amount << " node coordinates starting from offset: " << offset << endl;

      double* buffer = NULL;
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",meshName));
      if (vlsvReader->read("MESH_NODE_CRDS",attribs,offset,amount,buffer,true) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node coordinates" << endl;
         nodeCoordinateArraysRead = false;
         delete [] buffer; buffer = NULL;
         return nodeCoordinateArraysRead;
      }

      for (uint64_t i=0; i<amount; ++i) {
         //points->SetPoint(i,buffer[3*i+0],buffer[3*i+1],buffer[3*i+2]);
         points->InsertNextPoint(buffer[3*i+0],buffer[3*i+1],buffer[3*i+2]);
      }
      
      double* ptr = reinterpret_cast<double*>(buffer);
      for (uint64_t i=0; i<amount; ++i) {
         debug5 << "VLSV\t\t node " << i << " coords: " << ptr[i*3+0] << ' ' << ptr[i*3+1] << ' ' << ptr[i*3+2] << endl;
      }

      delete [] buffer; buffer = NULL;
      return nodeCoordinateArraysRead;
   }
   
   bool VisitUCDGenericMultiMeshReader::readNodeVariable(vlsv::Reader* vlsvReader,VisitUCDGenericMultiMeshMetadata*  metadata,
                                                         const VariableMetadata& vmd,int domain,void*& output) {
      return false;      
   }
   
   bool VisitUCDGenericMultiMeshReader::readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output) {
      debug2 << "VLSV\t VisitUCDGenericMultiMeshReader::readVariable called, domain: " << domain << endl;
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
      VisitUCDGenericMultiMeshMetadata* const metadata = dynamic_cast<VisitUCDGenericMultiMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
         debug3 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitUCDGenericMultiMeshMedata" << endl;
         return false;
      }

      /*
      if (vmd.centering == vlsvplugin::ZONE_CENTERED) {
	 return readCellVariable(vlsvReader,metadata,vmd,domain,output);
      } else {
	 output = NULL;
	 return false;
      }*/
      
      // Get domain offset arrays:
      const uint64_t* dummy = NULL;
      const uint64_t* ghostOffsets  = NULL;
      const uint64_t* variableOffsets = NULL;
      if (vmd.centering == vlsvplugin::NODE_CENTERED) {
         if (metadata->getDomainInfoNodes(vlsvReader,domain,dummy,ghostOffsets,variableOffsets) == false) {
            debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
            return false;
         }
      } else {
         if (metadata->getDomainInfoZones(vlsvReader,domain,dummy,ghostOffsets,variableOffsets) == false) {
            debug2 << "VLSV\t\t ERROR: Failed to obtain domain metadata" << endl;
            return false;
         }
      }

      // Figure out total number of variable values to read:
      uint32_t blockSize = 1;
      uint64_t N_totalValues = 0;
      uint64_t N_ghostValues = 0;
      uint64_t N_localValues = 0;
      const uint64_t N_components  = vmd.vectorSize;

      if (vmd.centering == vlsvplugin::NODE_CENTERED) {
         N_totalValues = metadata->getNumberOfTotalNodes(domain);
         N_ghostValues = metadata->getNumberOfGhostNodes(domain);
         N_localValues = metadata->getNumberOfLocalNodes(domain);
      } else {
         N_totalValues = metadata->getNumberOfTotalZones(domain);
         N_ghostValues = metadata->getNumberOfGhostZones(domain);
         N_localValues = metadata->getNumberOfLocalZones(domain);
         
         // Get mesh bounding box:
         const uint64_t* bbox = metadata->getMeshBoundingBox();
         if (bbox == NULL) {
            debug3 << "VLSV\t\t ERROR: Mesh bounding box array is NULL" << endl;
            return false;
         }
         blockSize =
           bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
           * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
           * bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      }
      
      /*
      // Get mesh bounding box:
      const uint64_t* bbox = metadata->getMeshBoundingBox();
      if (bbox == NULL) {
	 debug3 << "VLSV\t\t ERROR: Mesh bounding box array is NULL" << endl;
	 return false;
      }
      const uint64_t blockSize =
	  bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
	* bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
	* bbox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      
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
      */
      // Create vtkDataArray for variable data:
      bool success = true;
      //vtkFloatArray* rv = vtkFloatArray::New();
      vtkDoubleArray* rv = vtkDoubleArray::New();
      rv->SetNumberOfComponents(vmd.vectorSize);
      rv->SetNumberOfTuples(N_totalValues*blockSize);
      //float* variableData = rv->GetPointer(0);
      double* variableData = rv->GetPointer(0);
      
      // Read variable values from domain's real cells:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",vmd.name));
      attribs.push_back(make_pair("mesh",md->getName()));
      if (vlsvReader->read("VARIABLE",attribs,variableOffsets[domain]*blockSize,N_localValues*blockSize,variableData,false) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read domain's real cell variable data" << endl;
	 success = false;
      }
      
      // Read array that tell which domains contain ghost block data:
      list<pair<string,string> > meshAttribs;
      meshAttribs.push_back(make_pair("mesh",md->getName()));
      uint64_t* ghostDomains = NULL;
      if (vlsvReader->read("MESH_GHOST_DOMAINS",meshAttribs,ghostOffsets[domain],N_ghostValues,ghostDomains,true) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_DOMAINS array" << endl;
	 delete [] ghostDomains; ghostDomains = NULL;
	 success = false;
      }
      
      // Read array that tells local IDs of ghost cell data in each domain:
      uint64_t* ghostLocalIDs = NULL;
      if (vlsvReader->read("MESH_GHOST_LOCALIDS",meshAttribs,ghostOffsets[domain],N_ghostValues,ghostLocalIDs,true) == false) {
	 debug2 << "VLSV\t\t ERROR: Failed to read domain's MESH_GHOST_LOCALIDS array" << endl;
	 delete [] ghostDomains; ghostDomains = NULL;
	 delete [] ghostLocalIDs; ghostLocalIDs = NULL;
	 success = false;
      }
      
      // Read variable values for domain's ghost cells:
      if (success == true) {
	 //float* ptr = variableData + N_blocks*blockSize*components;
	 double* ptr = variableData + N_localValues*blockSize*N_components;
	 for (uint64_t i=0; i<N_ghostValues; ++i) {
	    const uint64_t ghostDomainID    = ghostDomains[i];
	    const uint64_t ghostValueOffset = (variableOffsets[ghostDomainID] + ghostLocalIDs[i])*blockSize;
	    if (vlsvReader->read("VARIABLE",attribs,ghostValueOffset,blockSize,ptr,false) == false) {
	       debug2 << "VLSV\t\t ERROR: Failed to read domain's ghost values" << endl;
	       success = false;
	       break;
	    }
	    ptr += blockSize*N_components;
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
