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

#ifndef DUPLICATE_NODE_ELIMINATION_H
#define DUPLICATE_NODE_ELIMINATION_H

#include <stdint.h>
#include <functional>
#include <tuple>
#include <cmath>
#include <vector>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

namespace vlsvplugin {
   
   template<typename T,typename U> inline
   void applyTransform(T* crds,const U* const transform) {
      T tmp[3];
      tmp[0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
      tmp[1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
      tmp[2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
      std::copy(tmp,tmp+3,crds);
   }

   template<typename S,typename T,typename U>
   class CartesianGeometry2D {
   public:
      void operator()(const S* crds_node_x,const S* crds_node_y,T* crds,const U& i,const U& j) {
         crds[0] = crds_node_x[i];
         crds[1] = crds_node_y[j];
         crds[2] = 0.0;
      }
   };

   template<typename S,typename T,typename U>
   class CartesianGeometry3D {
   public:
      void operator()(const S* crds_node_x,const S* crds_node_y,const S* crds_node_z,T* crds,const U& i,const U& j,const U& k) {
         crds[0] = crds_node_x[i];
         crds[1] = crds_node_y[j];
         crds[2] = crds_node_z[k];
      }
   };

   template<typename S,typename T,typename U>
   class CylindricalGeometry2D {
   public:
      void operator()(const S* crds_node_x,const S* crds_node_y,T* crds,const U& i,const U& j) {
         const T R   = crds_node_x[i];
         const T PHI = crds_node_y[j];
            
         crds[0] = R*cos(PHI);
         crds[1] = R*sin(PHI);
         crds[2] = 0.0;
      }
   };

   template<typename S,typename T,typename U>
   class CylindricalGeometry3D {
   public:
      void operator()(const S* crds_node_x,const S* crds_node_y,const S* crds_node_z,T* crds,const U& i,const U& j,const U& k) {
         const T R   = crds_node_x[i];
         const T PHI = crds_node_y[j];
         const T Z   = crds_node_z[k];
            
         crds[0] = R*cos(PHI);
         crds[1] = R*sin(PHI);
         crds[2] = Z;
      }
   };

   template<typename S,typename T,typename U>
   class SphericalGeometry3D {
   public:
      void operator()(const S* crds_node_x,const S* crds_node_y,const S* crds_node_z,T* crds,const U& i,const U& j,const U& k) {
         const T R     = crds_node_x[i];
         const T THETA = crds_node_y[j];
         const T PHI   = crds_node_z[k];

         crds[0] = R*sin(THETA)*cos(PHI);
         crds[1] = R*sin(THETA)*sin(PHI);
         crds[2] = R*cos(THETA);
      }
   };

   /** Insert all nodes for the given cells to a VTK unstructured grid while simultaneously
    *  eliminating duplicate node.
    *
    *  @tparam CONT Associative container that is used in duplicate node elimination.
    *  @tparam GEOMETRY Class that calculates Cartesian xyz coordinates from given coordinates 
    *          specified in another (Cartesian, Cylindrical, Spherical) coordinate system.
    */
   template<class CONT,class GEOMETRY> inline
   vtkIdType insertNodes2D(const double* transform,
                           const float* crds_node_x,const float* crds_node_y,
                           const bool& yPeriodic,const uint64_t& N_nodes_y,
					       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
                           vtkUnstructuredGrid* ugrid,vtkPoints* coordinates) {

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

      return counter;
   }

   /** Insert all nodes for the given cells to a VTK unstructured grid while simultaneously
    *  eliminating duplicate node.
    *
    *  @tparam CONT Associative container that is used in duplicate node elimination.
    *  @tparam GEOMETRY Class that calculates Cartesian xyz coordinates from given coordinates 
    *          specified in another (Cartesian, Cylindrical, Spherical) coordinate system.
    */
   template<class CONT,class GEOMETRY> inline
   vtkIdType insertNodes3D(const double* transform,
                           const float* crds_node_x,const float* crds_node_y,const float* crds_node_z,
                           const bool& yPeriodic,const bool& zPeriodic,const uint64_t& N_nodes_y,const uint64_t& N_nodes_z,
					       uint64_t N_totalBlocks,const uint64_t* blockGIDs,const std::vector<uint64_t>& bbox,
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


   template<class T> inline 
   void hash_combine(std::size_t& seed,const T& v) {
      seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
   }
   
   template<class T>
   class MyHash {
   public:
      inline size_t operator()(const std::tuple<T,T,T>& t) const {
         size_t seed = 0;
         hash_combine<T>(seed,std::get<0>(t));
         hash_combine<T>(seed,std::get<1>(t));
         hash_combine<T>(seed,std::get<2>(t));
         return seed;
      }
   };

   /** Maximum value for node index, equals to 2^21, node (i,j,k) indices are packed into a 64-bit integer (3*21=63).*/
   const int32_t maxNodeIndex = 2097152-1;
   
   /** Struct for storing (i,j,k) indices of nodes in the mesh. This 
    * struct is used to eliminate duplicate nodes from multi-domain
    * mesh pieces written to VLSV file.*/
   struct NodeIndices {
      int32_t i;            /**< i-index of the node.*/
      int32_t j;            /**< j-index of the node.*/
      int32_t k;            /**< k-index of the node.*/
      
      /** Constructor.
       * @param i i-index of the node.
       * @param j j-index of the node.
       * @param k k-index of the node.*/
      NodeIndices(int32_t i,int32_t j,int32_t k);
   };
   
   /** Comparator object for struct NodeIndices, required for being
    * able to store NodeIndices to unordered_map.*/
   struct NodesAreEqual {
      bool operator()(const NodeIndices& first,const NodeIndices& second) const;
   };
   
   /** Hash function implementation for object NodeIndices, required for 
    * being able to store NodeIndices to unordered_map.*/
   struct NodeHash {
      uint64_t operator()(const NodeIndices& node) const;
   };
   
} // namespace vlsvplugin

#endif
