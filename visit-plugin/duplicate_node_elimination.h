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
