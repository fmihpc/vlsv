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

#include <mesh_metadata_visit_point.h>

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   VisitPointMeshMetadata::VisitPointMeshMetadata(): VisitMeshMetadata() { }
   
   VisitPointMeshMetadata::~VisitPointMeshMetadata() { }

   const vlsv::geometry::type& VisitPointMeshMetadata::getMeshGeometry() const {return geometry;}

   uint64_t VisitPointMeshMetadata::getNumberOfGhostNodes(int domain) const {
      return MeshMetadata::getNumberOfGhostZones();
   }
   
   uint64_t VisitPointMeshMetadata::getNumberOfGhostZones(int domain) const {
      return MeshMetadata::getNumberOfGhostZones();
   }

   uint64_t VisitPointMeshMetadata::getNumberOfLocalNodes(int domain) const {
      return MeshMetadata::getNumberOfLocalZones();
   }
   
   uint64_t VisitPointMeshMetadata::getNumberOfLocalZones(int domain) const {
      return MeshMetadata::getNumberOfLocalZones();
   }

   uint64_t VisitPointMeshMetadata::getNumberOfTotalNodes(int domain) const {
      return MeshMetadata::getNumberOfTotalZones();
   }
   
   uint64_t VisitPointMeshMetadata::getNumberOfTotalZones(int domain) const {
      return MeshMetadata::getNumberOfTotalZones();
   }
   
   bool VisitPointMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Call superclass read function:
      if (VisitMeshMetadata::read(vlsvReader,attribs) == false) return false;
      
      // Check that we are reading point mesh metadata:
      bool success = true;
      map<string,string>::const_iterator it = attribs.find("type");
      if (it == attribs.end()) {
	 debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
	 success = false;
      }
      
      if (it->second != vlsv::mesh::STRING_POINT) {
	 debug3 << "VLSV\t\t ERROR: Mesh type is '" << it->second << "', should be '" << vlsv::mesh::STRING_POINT << "'" << endl;
	 success = false;
      } else {
	 meshType = AVT_POINT_MESH;
	 meshTypeString = "AVT_POINT_MESH";
	 spatialDimension = 3;
	 topologicalDimension = 0;
      }

      it = attribs.find("arraysize");
      if (it == attribs.end()) success = false;
      N_totalZones = atoi(it->second.c_str());
      N_localZones = N_totalZones;
      N_ghostZones = 0;

      // Get coordinate system:
      it = attribs.find("geometry");
      if (it == attribs.end()) geometry = vlsv::geometry::CARTESIAN;
      else {
	 geometry = vlsv::getMeshGeometry(it->second);
      }
      
      return success;
   }
   
} // namespace vlsvplugin
