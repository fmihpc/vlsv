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

#include <duplicate_node_elimination.h>

using namespace std;

namespace vlsvplugin {

   /** Constructor.
    * @param i i-index of the node.
    * @param j j-index of the node.
    * @param k k-index of the node.*/
   NodeIndices::NodeIndices(int32_t i,int32_t j,int32_t k): i(i),j(j),k(k) { }

   /** Compare given NodeIndices objects for equality. Objects 
    * are equal if their (i,j,k) indices are equal.
    * @param first First NodeIndices object to compare.
    * @param second Second NodeIndices object to compare.
    * @return If true, given objects are equal.*/
   bool NodesAreEqual::operator()(const NodeIndices& first,const NodeIndices& second) const {
      if (first.i != second.i) return false;
      if (first.j != second.j) return false;
      if (first.k != second.k) return false;
      return true;
   }
   
   /** Calculate hash value for given NodeIndices object.
    * @param Node NodeIndices object whose hash value is to be calculated.
    * @return Calculated hash value.*/
   uint64_t NodeHash::operator()(const NodeIndices& node) const {
      uint64_t result = 0;
      uint64_t tmp = node.i % maxNodeIndex;
      result = (result | tmp);
      
      tmp = node.j % maxNodeIndex;
      result = (result | (tmp << 21));
      
      tmp = node.k % maxNodeIndex;
      result = (result | (tmp << 42));
      return result;
   }   
   
} // namespace vlsvplugin
