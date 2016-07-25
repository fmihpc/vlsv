/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2015 Finnish Meteorological Institute
 *  Copyright 2016 Arto Sandroos
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
#include <vtkIdList.h>
#include <vtkSmartPointer.h>

using namespace std;

namespace vlsvplugin {

   typedef unordered_map<vtkIdType,vtkIdType> MyCont;
   typedef unordered_map<tuple<vtkIdType,vtkIdType,vtkIdType>,vtkIdType,MyHash<vtkIdType> > NodeMap;

   typedef vtkIdType (*nodeFunction2D)(const double* transform,
                                      const float* crds_node_x,const float* crds_node_y,
                                      const bool& yPeriodic,const uint64_t& N_nodes_y,
					                  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
                                      vtkUnstructuredGrid* ugrid,vtkPoints* coordinates);

   typedef vtkIdType (*nodeFunction3D)(const double* transform,
                                      const float* crds_node_x,const float* crds_node_y,const float* crds_node_z,
                                      const bool& yPeriodic,const bool& zPeriodic,const uint64_t& N_nodes_y,const uint64_t& N_nodes_z,
					                  uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
                                      vtkUnstructuredGrid* ugrid,vtkPoints* coordinates);

   map<vlsv::geometry::type,nodeFunction2D> nodeFunctions2D;
   map<vlsv::geometry::type,nodeFunction3D> nodeFunctions3D;

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

   template<class CONT,class GEOMETRY> inline
   vtkIdType insertNodes2D(const double* transform,
                           const float* crds_node_x,const float* crds_node_y,
                           const bool& yPeriodic,const uint64_t& N_nodes_y,
					       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
                           vtkUnstructuredGrid* ugrid,vtkPoints* coordinates) {
      
      debug2 << "VLSV\t\t insertNodes2D called" << endl;
      debug3 << "VLSV\t\t bbox: " << bbox[0] << ' ' << bbox[1] << ' ' << bbox[2] << ' ';
      debug3 << bbox[3] << ' ' << bbox[4] << ' ' << bbox[5] << endl;

      // Cartesian cells are technically VTK_PIXELs, but 
      // VTK_QUAD happens to work all basic geometries
      const int cellType = VTK_QUAD;
      const int N_VERTICES = 4;

      vtkIdType counter = 0;

      GEOMETRY geometry;

      const vtkIdType SX  = bbox[0]+1;

      // Reserving enough capacity in nodeMapping gives about 20% performance
      // boost for large Cartesian grids. Preallocating vtkUnstructuredGrid or 
      // vtkPoints doesn't seem to have any performance effects.
      const vtkIdType initialCapacity = static_cast<vtkIdType>(ceil(1.1*N_totalBlocks));
      CONT nodeMapping;
      nodeMapping.reserve(initialCapacity);
      ugrid->Allocate(initialCapacity,1000);
      coordinates->Allocate(initialCapacity,1000);

      auto nodeIndex = [&](const vtkIdType& i,const vtkIdType& j) {
         return j*SX + i;
      };

      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         vtkIdType j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];

         for (vtkIdType j=0; j<bbox[4]; ++j) {
            const vtkIdType j_cell = j_block*bbox[4] + j;
            for (vtkIdType i=0; i<bbox[3]; ++i) {
               const vtkIdType i_cell = i_block*bbox[3] + i;

               vtkIdType vertices[N_VERTICES];
               int count=0;

               for (int n=0; n<4; ++n) {
                  const vtkIdType in = min(1,n%3);
                  vtkIdType jn = n/2;

                  // Take (possible) periodicity into account:
                  if (yPeriodic == true) if (j_cell+jn >= N_nodes_y-1) jn = -j_cell;

                  // Global ID of the node we're looking for:
                  const vtkIdType node_index = nodeIndex(i_cell+in,j_cell+jn);

                  // If the node already exists, store it's position to vertices.
                  // Otherwise it needs to be created. 
                  // Note: old version of this code attempted to insert all nodes, 
                  // the new version here is 100% faster!
                  const auto result = nodeMapping.find(node_index);
                  if (result != nodeMapping.end()) {
                     vertices[count] = result->second;
                  } else {
                     nodeMapping.emplace(make_pair(node_index,counter));
                     vertices[count] = counter;
                     ++counter;
                     float crds[3];
                     geometry(crds_node_x,crds_node_y,crds,i_cell+in,j_cell+jn);
                     vlsvplugin::applyTransform(crds,transform);
                     coordinates->InsertNextPoint(crds);                           
                  }

                  ++count;
               }
               ugrid->InsertNextCell(cellType,N_VERTICES,vertices);
            }
         }
      }
      debug3 << "VLSV\t\t mesh has " << nodeMapping.size() << " unique nodes" << endl;
      debug3 << "VLSV\t\t coordinates has " << coordinates->GetNumberOfPoints() << " points" << endl;
      debug3 << "VLSV\t\t ugrid has " << ugrid->GetNumberOfCells() << " cells" << endl;

      return counter;
   }

   template<class CONT,class GEOMETRY> inline
   vtkIdType insertNodes3D(const double* transform,
                           const float* crds_node_x,const float* crds_node_y,const float* crds_node_z,
                           const bool& yPeriodic,const bool& zPeriodic,const uint64_t& N_nodes_y,const uint64_t& N_nodes_z,
					       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const uint64_t* bbox,
                           vtkUnstructuredGrid* ugrid,vtkPoints* coordinates) {

      // Cartesian cells are technically VTK_VOXELs, but 
      // VTK_HEXAHEDRON happens to work all basic geometries
      const int cellType = VTK_HEXAHEDRON;
      
      GEOMETRY geometry;      

      // Reserving enough capacity in nodeMapping gives about 15% performance
      // boost for large Cartesian grids. Preallocating vtkUnstructuredGrid or 
      // vtkPoints doesn't seem to have any performance effects.
      const vtkIdType initialCapacity = static_cast<vtkIdType>(ceil(1.1*N_totalBlocks));
      CONT nodeMapping;
      nodeMapping.reserve(initialCapacity);
      ugrid->Allocate(initialCapacity,1000);
      coordinates->Allocate(initialCapacity,1000);

      const vtkIdType SX  = bbox[0]+1;
      const vtkIdType SY  = bbox[1]+1;
      const vtkIdType SXY = SX*SY;
      auto nodeIndex = [&](const vtkIdType& i,const vtkIdType& j,const vtkIdType& k) {
         return k*SXY + j*SX + i;
      };

      // The counter 'counter' counts the number of nodes inserted to nodeMapping
      vtkIdType counter = 0;
      for (uint64_t block=0; block<N_totalBlocks; ++block) {
         uint64_t i_block = blockGIDs[block];
         vtkIdType k_block = i_block / (bbox[1]*bbox[0]);
         i_block -= k_block*(bbox[1]*bbox[0]);
         vtkIdType j_block = i_block / bbox[0];
         i_block -= j_block*bbox[0];

         for (vtkIdType k=0; k<bbox[5]; ++k) {
            const vtkIdType k_cell = k_block*bbox[5] + k;
            for (vtkIdType j=0; j<bbox[4]; ++j) {
               const vtkIdType j_cell = j_block*bbox[4] + j;
               for (vtkIdType i=0; i<bbox[3]; ++i) {
                  // Calculate cell's bounding box global indices:
                  const vtkIdType i_cell = i_block*bbox[3] + i;

                  vtkIdType vertices[8];
                  int count=0;
                  for (vtkIdType z=0; z<2; ++z) {
                     vtkIdType kn = z;
                     //for (vtkIdType jn=0; jn<2; ++jn) for (vtkIdType in=0; in<2; ++in) {
                     for (int n=0; n<4; ++n) {
                        const vtkIdType in = min(1,n%3);
                        vtkIdType jn = n/2;

                        // Take (possible) periodicity into account:
                        if (yPeriodic == true) if (j_cell+jn >= N_nodes_y-1) jn = -j_cell;
                        if (zPeriodic == true) if (k_cell+kn >= N_nodes_z-1) kn = -k_cell;
                       
                        // Global ID of the node we're looking for:
                        const vtkIdType node_index = nodeIndex(i_cell+in,j_cell+jn,k_cell+kn);

                        // If the node already exists, store it's position to vertices.
                        // Otherwise it needs to be created. 
                        // Note: old version of this code attempted to insert all nodes, 
                        // the new version here is 100% faster!
                        const auto result = nodeMapping.find(node_index);
                        if (result != nodeMapping.end()) {
                           vertices[count] = result->second;
                        } else {
                           nodeMapping.emplace(make_pair(node_index,counter));
                           vertices[count] = counter;
                           ++counter;
                           float crds[3];
                           geometry(crds_node_x,crds_node_y,crds_node_z,crds,i_cell+in,j_cell+jn,k_cell+kn);
                           vlsvplugin::applyTransform(crds,transform);
                           coordinates->InsertNextPoint(crds);                           
                        }
                        ++count;
                     }
                  }
                  ugrid->InsertNextCell(cellType,8,vertices);
               }
            }
         }

      }
      return counter;
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

      nodeFunctions2D[vlsv::geometry::UNKNOWN]     = insertNodes2D<MyCont,CartesianGeometry2D<float,float,vtkIdType> >;
      nodeFunctions2D[vlsv::geometry::CARTESIAN]   = insertNodes2D<MyCont,CartesianGeometry2D<float,float,vtkIdType> >;
      nodeFunctions2D[vlsv::geometry::CYLINDRICAL] = insertNodes2D<MyCont,CylindricalGeometry2D<float,float,vtkIdType> >;
      
      nodeFunctions3D[vlsv::geometry::UNKNOWN]     = insertNodes3D<MyCont,CartesianGeometry3D<float,float,vtkIdType> >;
      nodeFunctions3D[vlsv::geometry::CARTESIAN]   = insertNodes3D<MyCont,CartesianGeometry3D<float,float,vtkIdType> >;
      nodeFunctions3D[vlsv::geometry::CYLINDRICAL] = insertNodes3D<MyCont,CylindricalGeometry3D<float,float,vtkIdType> >;
      nodeFunctions3D[vlsv::geometry::SPHERICAL]   = insertNodes3D<MyCont,SphericalGeometry3D<float,float,vtkIdType> >;

      debug4 << "VLSV\t\t arraysize:  " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize: " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:   " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:   " << metadata->getDatatype() << endl;
            
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
      if (metadata->getSpatialDimension() == 3) {
         N_nodes_z = bbox[2]*bbox[5]+1;
      } else {
         N_nodes_z = 1;
      }
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
      if (readNodeCoordinateArrays(vlsvReader,metadata) == false) {
         return false;
      }
      
      // Get mesh geometry and periodicity:
      const vlsv::geometry::type& geometry = metadata->getMeshGeometry();
      metadata->getMeshPeriodicity(xPeriodic,yPeriodic,zPeriodic);

      // Force Cartesian meshes to be non-periodic for correct visualization:
      if (geometry == vlsv::geometry::CARTESIAN) {
         xPeriodic = false; yPeriodic = false; zPeriodic = false;
      }

      // For each block, attempt to insert its (WX+1)*(WY+1)*(WZ+1) nodes into unordered_map.
      // The insertion will fail if the node already exists in the unordered_map, in which case 
      // counter is not increased. Counter, i.e. the value of unordered_map for given
      // NodeIndices, tells node's index.
      NodeMap nodeIndices;
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->Allocate(N_totalBlocks*blockSize);
      vtkPoints* coordinates = vtkPoints::New();
      size_t N_uniqueNodes = 0;

      // Get transformation matrix
      const double* transform = metadata->getTransform();

      if (metadata->getSpatialDimension() == 2) {
         if (nodeFunctions2D.find(geometry) != nodeFunctions2D.end()) {
            N_uniqueNodes = nodeFunctions2D[geometry](transform,crds_node_x,crds_node_y,
                                                      yPeriodic,N_nodes_y,
                                                      N_totalBlocks,blockGIDs,bbox,ugrid,coordinates);
         }
      } else {
         if (nodeFunctions3D.find(geometry) != nodeFunctions3D.end()) {
            N_uniqueNodes = nodeFunctions3D[geometry](transform,crds_node_x,crds_node_y,crds_node_z,
                                                      yPeriodic,zPeriodic,N_nodes_y,N_nodes_z,
                                                      N_totalBlocks,blockGIDs,bbox,ugrid,coordinates);
         }
      }
      
      // unordered_map now contains all unique nodes. Its size is equal to
      // number of nodes in this domain:
      debug4 << "VLSV\t\t domain has " << N_uniqueNodes << " unique nodes" << endl;

      // Create vtkUnstructuredGrid:
      ugrid->SetPoints(coordinates);
      coordinates->Delete();
      delete [] blockGIDs; blockGIDs = NULL;

      // Determine correct values for real and ghost (internal to problem) zones:
      unsigned char cellIsReal = 0;
      unsigned char cellIsGhost = 0;
      avtGhostData::AddGhostZoneType(cellIsGhost,DUPLICATED_ZONE_INTERNAL_TO_PROBLEM);
      
      // Create an array that flags each zone either as real or internal ghost:
      vtkUnsignedCharArray* ghostZones = vtkUnsignedCharArray::New();
      ghostZones->SetName("avtGhostZones");
      ghostZones->Allocate(N_totalBlocks*blockSize);
      debug4 << "VLSV\t\t Flagging cells 0-" << N_blocks*blockSize-1 << " as existing" << endl;
      debug4 << "VLSV\t\t Flagging cells " << N_blocks*blockSize << "-" << N_totalBlocks*blockSize-1 << " as ghosts" << endl;
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

   bool VisitUCDMultiMeshReader::readNodeCoordinateArrays(vlsv::Reader* vlsvReader,vlsvplugin::VisitMeshMetadata* md) {
      debug1 << "VLSV\t\t readNodeCoordinateArrays called" << endl;
      // Check if coordinate arrays are already cached:
      if (nodeCoordinateArraysRead == true) return nodeCoordinateArraysRead;
      
      nodeCoordinateArraysRead = true;
      delete [] crds_node_x; crds_node_x = NULL;
      delete [] crds_node_y; crds_node_y = NULL;
      delete [] crds_node_z; crds_node_z = NULL;
      
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",md->getName()));
      
      // Check that node coordinate arrays have correct sizes:
      
      
      // Read node coordinate arrays:
      if (vlsvReader->read("MESH_NODE_CRDS_X",attribs,0,N_nodes_x,crds_node_x) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node x-coordinate array" << endl;
         nodeCoordinateArraysRead = false;
      } else {
         debug2 << "VLSV\t\t Read " << N_nodes_x << " node x-coordinates" << endl;
      }
      if (vlsvReader->read("MESH_NODE_CRDS_Y",attribs,0,N_nodes_y,crds_node_y) == false) {
         debug2 << "VLSV\t\t ERROR: Failed to read node y-coordinate array" << endl;
         nodeCoordinateArraysRead = false;
      } else {
         debug2 << "VLSV\t\t Read " << N_nodes_y << " node y-coordinates" << endl;
      }
      if (md->getSpatialDimension() == 3) {
         if (vlsvReader->read("MESH_NODE_CRDS_Z",attribs,0,N_nodes_z,crds_node_z) == false) {
            debug2 << "VLSV\t\t ERROR: Failed to read node z-coordinate array" << endl;
            nodeCoordinateArraysRead = false;
         } else {
            debug2 << "VLSV\t\t Read " << N_nodes_z << " node z-coordinates" << endl;
         }
      } else {
         delete [] crds_node_z; crds_node_z = NULL;
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
