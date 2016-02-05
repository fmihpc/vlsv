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

#include <mesh_metadata_visit_quad_multi.h>

#include <list>

using namespace std;

namespace vlsvplugin {
   
   QuadMultiMeshMetadata::QuadMultiMeshMetadata(): MeshMetadata() { }
   
   QuadMultiMeshMetadata::~QuadMultiMeshMetadata() { }
   
   const std::string& QuadMultiMeshMetadata::getCorrectVlsvMeshType() const {
      return vlsv::mesh::STRING_QUAD_MULTI;
   }

   bool QuadMultiMeshMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
						  const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      // Check that VLSVReader exists:
      if (vlsvReader == NULL) {
         stringstream ss;
         ss << "ERROR: VLSVReader is NULL";
         return exitWithError(ss);
      }

      // Read domain info:
      if (readDomains(vlsvReader) == false) return false;

      domainOffsets = MeshMetadata::zoneDomainOffsets.data();
      ghostOffsets = MeshMetadata::zoneGhostOffsets.data();
      variableOffsets = MeshMetadata::zoneVariableOffsets.data();
      return true;
   }

   const float* QuadMultiMeshMetadata::getMeshBoundingBox() {return meshCoordinates.data();}
   
   bool QuadMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function. If it fails, meshMetadataRead has value 'false'.
      if (MeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;
      
      spatialDimension = 3;
      topologicalDimension = 3;

      // Figure out total number of cells in the mesh:
      auto it = attribs.find("arraysize");
      if (it == attribs.end()) return false;
      MeshMetadata::N_totalZones = atoi(it->second.c_str());
      
      // Read XML tag 'MESH_ZONES' to figure out how many domains the mesh has:
      map<string,string> attribsOut;
      list<pair<string,string> > attribsIn;
      attribsIn.push_back(make_pair("mesh",getName()));
      if (vlsvReader->getArrayAttributes("MESH_ZONES",attribsIn,attribsOut) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_ZONES' attributes for mesh '" << getName() << "' from VLSV file";
         return exitWithError(ss);
      }      
      it = attribsOut.find("arraysize");
      if (it == attribsOut.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_ZONES' XML tag does not have attribute 'arraysize'";
         return exitWithError(ss);
      } else {         
         MeshMetadata::N_domains = atoi(it->second.c_str());
      }

      meshMetadataRead = true;
      return meshMetadataRead;
   }
   
   bool QuadMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      return MeshMetadata::readDomainMetadata(vlsvReader);
   }

} // namespace vlsvplugin
