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

#include <mesh_metadata_visit_ucd_multi.h>

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   UCDMultiMeshMetadata::UCDMultiMeshMetadata() { }
   
   UCDMultiMeshMetadata::~UCDMultiMeshMetadata() { }
   
   const std::string& UCDMultiMeshMetadata::getCorrectVlsvMeshType() const {
      return vlsv::mesh::STRING_UCD_MULTI;
   }

   bool UCDMultiMeshMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
						                         const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t UCDMultiMeshMetadata::getDomainInfo called, domain: " << domain << endl;
      
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         debug3 << "VLSV\t\t ERROR: vlsv::Reader is NULL" << endl;
         return false;
      }
	     
      // Read domain info:
      if (readDomains(vlsvReader) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read quad multimesh domains" << endl;
         return false;
      }
      
      debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
      debug4 << (this->zoneDomainOffsets)[domain] << ' ' << (this->zoneDomainOffsets)[domain+1];
      debug4 << " ghost offsets: " << (this->zoneGhostOffsets)[domain] << ' ' << (this->zoneGhostOffsets)[domain+1];
      debug4 << " variable offsets: " << (this->zoneVariableOffsets)[domain] << ' ' << (this->zoneVariableOffsets)[domain+1];
      debug4 << endl;

      domainOffsets = this->zoneDomainOffsets.data();
      ghostOffsets = this->zoneGhostOffsets.data();
      variableOffsets = this->zoneVariableOffsets.data();
      return true;
   }

   bool UCDMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t UCDMultiMeshMetadata::read called" << endl;
      
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function:
      if (MeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;
      
      spatialDimension = 3;
      topologicalDimension = 3;
      
      // Check for definition of spatial dimensions (could be two):
      auto it = attribs.find("spatial_dimension");
      if (it != attribs.end()) {
         spatialDimension = atoi(it->second.c_str());
         if (spatialDimension < 2 || spatialDimension > 3) {
            debug3 << "VLSV\t\t ERROR: Spatial dimension must be 2 or 3, value '" << it->second << "' given in XML attributes" << endl;
            return false;
         }
         topologicalDimension = spatialDimension;
      }

      // Figure out total number of cells in the mesh:
      it = attribs.find("arraysize");
      if (it == attribs.end()) return false;
      MeshMetadata::N_totalZones = atoi(it->second.c_str());

      // Read XML tag 'MESH_ZONES' to figure out how many domains the mesh has:
      map<string,string> attribsOut;
      list<pair<string,string> > attribsIn;
      attribsIn.push_back(make_pair("mesh",getName()));
      if (vlsvReader->getArrayAttributes("MESH_DOMAIN_SIZES",attribsIn,attribsOut) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_DOMAIN_SIZES' attributes for mesh '" << getName() << "' from VLSV file" << endl;
         return false;
      }
      it = attribsOut.find("arraysize");
      if (it == attribsOut.end()) {
         debug3 << "VLSV\t\t ERROR: Array 'MESH_DOMAIN_SIZES' XML tag does not have attribute 'arraysize'" << endl;
         return false;
      } else {
         debug3 << "VLSV\t\t Mesh has " << it->second << " domains" << endl;
         MeshMetadata::N_domains = atoi(it->second.c_str());
      }
      
      meshMetadataRead = true;
      return true;
   }

   bool UCDMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      return MeshMetadata::readDomainMetadata(vlsvReader);
   }
      
} // namespace vlsvplugin
