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

#ifndef MESH_READER_VISIT_POINT_H
#define MESH_READER_VISIT_POINT_H

#include <mesh_reader.h>

namespace vlsvplugin {
   
   class VisitPointMeshReader: public MeshReader {
    public:
      VisitPointMeshReader();
      ~VisitPointMeshReader();
      
      bool readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output);
      bool readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output);
      
    protected:
	      
   };
   
} // namespace vlsvplugin

#endif
