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

#include <mesh_metadata_visit_point.h>

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   VisitPointMeshMetadata::VisitPointMeshMetadata(): VisitMeshMetadata() { }
   
   VisitPointMeshMetadata::~VisitPointMeshMetadata() { }
   
   const std::string& VisitPointMeshMetadata::getCorrectVlsvMeshType() const {
      return vlsv::mesh::STRING_POINT;
   }

   bool VisitPointMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      if (meshMetadataRead == true) return meshMetadataRead;

      // Call superclass read function:
      if (VisitMeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;
      bool success = true;

      meshType = AVT_POINT_MESH;
      meshTypeString = "AVT_POINT_MESH";
      spatialDimension = 3;
      topologicalDimension = 0;

      auto it = attribs.find("arraysize");
      if (it == attribs.end()) success = false;
      N_totalZones = atoi(it->second.c_str());
      N_localZones = N_totalZones;
      N_ghostZones = 0;

      if (success == true) meshMetadataRead = true;
      return meshMetadataRead;
   }
   
} // namespace vlsvplugin
