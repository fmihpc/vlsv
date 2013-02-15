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

#ifndef DUPLICATE_NODE_ELIMINATION_H
#define DUPLICATE_NODE_ELIMINATION_H

#include <stdint.h>

namespace vlsvplugin {
   
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
