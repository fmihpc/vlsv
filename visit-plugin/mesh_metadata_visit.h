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

//#pragma once

#ifndef MESH_METADATA_VISIT_H
#define MESH_METADATA_VISIT_H

#include <vector>

#include <avtSTMDFileFormat.h>

#include <mesh_vtk.h>
#include <mesh_metadata.h>

namespace vlsvplugin {
   class VisitMeshMetadata: public MeshMetadata {
    public:
      VisitMeshMetadata();
      virtual ~VisitMeshMetadata();
      
      virtual int getBlockOrigin() const;
      virtual avtMeshType getMeshType() const;
      virtual std::string getMeshTypeString() const;

      virtual bool read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);
      
    protected:

      virtual bool readDomainMetadata(vlsv::Reader* vlsvReader);

      int blockOrigin;
      avtMeshType meshType;       /**< Mesh type, one of AVT_POINT_MESH, ... */
      std::string meshTypeString; /**< String representation of mesh type.*/      
   };
} // namespace vlsvplugin

#endif
