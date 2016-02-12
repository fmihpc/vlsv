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

//#pragma once

#ifndef MESH_METADATA_VISIT_H
#define MESH_METADATA_VISIT_H

#include <vector>

#include <avtSTMDFileFormat.h>

#include <mesh_vtk.h>
#include <mesh_metadata.h>

/** Class VisitMeshMetadata is an interface class which provide
 *  additional information about VLSV mesh(es) to VisIt. 
 *  Member functions are typically called in avtVlsvFileFormat.C file.
 *  
 *  VisIt plugin only uses mesh metadata pointers that are of type VisitMeshMetadata.
 *  Wrapper classes for all existing metadata classes are provided 
 *  in mesh_metadata_visit_classes.* files.
 */
namespace vlsvplugin {
   class VisitMeshMetadata: public virtual MeshMetadata {
    public:
      VisitMeshMetadata();
      virtual ~VisitMeshMetadata();
      
      virtual int getBlockOrigin() const;
      virtual avtMeshType getAvtMeshType() const =0;
      virtual std::string getAvtMeshTypeString() const =0;
      
    protected:
      virtual bool checkVlsvMeshType(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs);
      virtual bool readDomainMetadata(vlsv::Reader* vlsvReader);
      virtual bool readVariables(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs);

      int blockOrigin;
   };
} // namespace vlsvplugin

#endif
