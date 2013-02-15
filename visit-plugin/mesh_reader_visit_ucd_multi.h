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

#ifndef MESH_READER_VISIT_UCD_MULTI_H
#define MESH_READER_VISIT_UCD_MULTI_H

#include <mesh_reader.h>

namespace vlsvplugin {

   class VisitUCDMultiMeshReader: public MeshReader {
    public:
      VisitUCDMultiMeshReader();
      virtual ~VisitUCDMultiMeshReader();
      
      virtual bool readMesh(VLSVReader* vlsv,MeshMetadata* md,int domain,void*& output);
      virtual bool readVariable(VLSVReader* vlsv,MeshMetadata* md,const VariableMetadata& vmd,int domain,float*& output);
      
    protected:
      float* crds_node_x;
      float* crds_node_y;
      float* crds_node_z;
      bool nodeCoordinateArraysRead;
      uint64_t N_nodes_x;
      uint64_t N_nodes_y;
      uint64_t N_nodes_z;
      
      virtual bool readNodeCoordinateArrays(VLSVReader* vlsv,const std::string& meshName);
   };

} // namespace plugin

#endif
