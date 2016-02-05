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

#include <mesh_metadata_visit_ucd_amr.h>

#include <list>

using namespace std;

namespace vlsvplugin {

   UCDAMRMetadata::UCDAMRMetadata() { }
   
   UCDAMRMetadata::~UCDAMRMetadata() { }
   
   const std::string& UCDAMRMetadata::getCorrectVlsvMeshType() const {
      return vlsv::mesh::STRING_UCD_AMR;
   }

   bool UCDAMRMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
						  const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         stringstream ss;
         ss << "ERROR: vlsv::Reader is NULL" << endl;
         return MeshMetadata::exitWithError(ss);
      }
	     
      // Read domain info:
      if (readDomains(vlsvReader) == false) return false;
      
      domainOffsets   = MeshMetadata::zoneDomainOffsets.data();
      ghostOffsets    = MeshMetadata::zoneGhostOffsets.data();
      variableOffsets = MeshMetadata::zoneVariableOffsets.data();
      return true;
   }

   bool UCDAMRMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
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
            stringstream ss;
            ss << "ERROR: Spatial dimension must be 2 or 3, value '" << it->second << "' given in XML attributes";
            return MeshMetadata::exitWithError(ss);
         }
         topologicalDimension = spatialDimension;
      }

      // Figure out total number of cells in the mesh:
      it = attribs.find("arraysize");
      if (it == attribs.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH' does not have an XML attribute called 'arraysize'";
         return MeshMetadata::exitWithError(ss);
      }
      MeshMetadata::N_totalZones = atoi(it->second.c_str());

      // Read XML tag 'MESH_ZONES' to figure out how many domains the mesh has:
      map<string,string> attribsOut;
      list<pair<string,string> > attribsIn;
      attribsIn.push_back(make_pair("mesh",getName()));
      if (vlsvReader->getArrayAttributes("MESH_DOMAIN_SIZES",attribsIn,attribsOut) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_DOMAIN_SIZES' attributes for mesh '" << getName() << "'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'.";
         return MeshMetadata::exitWithError(ss);
      }
      it = attribsOut.find("arraysize");
      if (it == attribsOut.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_DOMAIN_SIZES' XML tag does not have attribute 'arraysize'";
         return MeshMetadata::exitWithError(ss);
      } else {         
         MeshMetadata::N_domains = atoi(it->second.c_str());
      }

      meshMetadataRead = true;
      return true;
   }

   bool UCDAMRMetadata::readDomains(vlsv::Reader* vlsvReader) {
      // Exit if domain metadata has already been read:
      if (domainMetadataRead == true) return true;
      domainMetadataRead = false;      
      
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         stringstream ss;
         ss << "ERROR: vlsv::Reader is NULL";
         return MeshMetadata::exitWithError(ss);
      }
      
      // Read mesh bounding box:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",name));
      meshBoundingBox.resize(6);
      uint64_t* ptr = meshBoundingBox.data();
      if (vlsvReader->read("MESH_BBOX",attribs,0,6,ptr,false) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_BBOX'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'.";
         return MeshMetadata::exitWithError(ss);
      }
      blockSize = meshBoundingBox[3]*meshBoundingBox[4]*meshBoundingBox[5];
      
      // Read domain info:
      int64_t* domainInfo = NULL;
      if (vlsvReader->read("MESH_DOMAIN_SIZES",attribs,0,MeshMetadata::N_domains,domainInfo) == false) {
         delete [] domainInfo; domainInfo = NULL;
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_ZONES'. " << endl;
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'.";         
         return MeshMetadata::exitWithError(ss);
      }
      
      // Calculate offsets where data for each domain begins:
      zoneDomainOffsets.resize(MeshMetadata::N_domains+1); zoneDomainOffsets.shrink_to_fit();
      zoneGhostOffsets.resize(MeshMetadata::N_domains+1); zoneGhostOffsets.shrink_to_fit();
      zoneVariableOffsets.resize(MeshMetadata::N_domains+1); zoneVariableOffsets.shrink_to_fit();
      zoneDomainOffsets[0]   = 0;
      zoneGhostOffsets[0]    = 0;
      zoneVariableOffsets[0] = 0;
      for (uint64_t i=0; i<MeshMetadata::N_domains; ++i) {
         zoneDomainOffsets[i+1]   = zoneDomainOffsets[i] + domainInfo[i*2];
         zoneGhostOffsets[i+1]    = zoneGhostOffsets[i] + domainInfo[i*2+1];
         zoneVariableOffsets[i+1] = zoneVariableOffsets[i] + domainInfo[i*2+0]-domainInfo[i*2+1];
      }
      delete [] domainInfo; domainInfo = NULL;
      
      // Compute total number of real and ghost cells:
      MeshMetadata::N_ghostZones = zoneGhostOffsets[MeshMetadata::N_domains];
      MeshMetadata::N_localZones  = zoneVariableOffsets[MeshMetadata::N_domains];
   
      domainMetadataRead = true;
      return true;
   }
      
} // namespace vlsvplugin
