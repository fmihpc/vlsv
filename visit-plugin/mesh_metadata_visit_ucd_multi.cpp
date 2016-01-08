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

#include <mesh_metadata_visit_ucd_multi.h>

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   VisitUCDMultiMeshMetadata::VisitUCDMultiMeshMetadata() { }
   
   VisitUCDMultiMeshMetadata::~VisitUCDMultiMeshMetadata() { }
   
   bool VisitUCDMultiMeshMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
						                         const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitUCDMultiMeshMetadata::getDomainInfo called, domain: " << domain << endl;
      
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

   bool VisitUCDMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitUCDMultiMeshMetadata::read called" << endl;
      
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function:
      if (VisitMeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;

      // Check that we are reading multi-domain unstructured mesh metadata:
      auto it = attribs.find("type");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
         return false;
      }
      if (it->second != vlsv::mesh::STRING_UCD_MULTI) {
         debug3 << "VLSV\t\t ERROR: Mesh type is '" << it->second << "', should be '" << vlsv::mesh::STRING_UCD_MULTI << "'" << endl;
         return false;
      }
      
      meshType = AVT_UNSTRUCTURED_MESH;
      meshTypeString = "AVT_UNSTRUCTURED_MESH";
      spatialDimension = 3;
      topologicalDimension = 3;
      
      // Check for definition of spatial dimensions (could be two):
      it = attribs.find("spatial_dimension");
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

      // Get mesh geometry:
      it = attribs.find("geometry");
      if (it == attribs.end()) {
         geometry = vlsv::geometry::CARTESIAN;
      } else {
         geometry = vlsv::getMeshGeometry(it->second);
         if (geometry == vlsv::geometry::UNKNOWN) geometry = vlsv::geometry::CARTESIAN;
      }

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
         VisitMeshMetadata::N_domains = atoi(it->second.c_str());
      }
      
      // Read 'VARIABLE' XML tags to parse names of variables in this mesh:
      set<string> variableNames;
      if (vlsvReader->getUniqueAttributeValues("VARIABLE","name",variableNames) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read variable names" << endl;
         return false;
      }

      // Remove variables that do not belong to this mesh:
      debug4 << "VLSV\t\t Found variables:" << endl;
      for (auto& it : variableNames) {
         map<string,string> attribsOut;
         map<string,string>::const_iterator mapIt;
         list<pair<string,string> > attribsIn;
         attribsIn.push_back(make_pair("mesh",getName()));
         attribsIn.push_back(make_pair("name",it));
	 
         // Skip variables belonging to other meshes:
         if (vlsvReader->getArrayAttributes("VARIABLE",attribsIn,attribsOut) == false) continue;
         debug4 << "VLSV\t\t\t '" << it << "' vectorsize: " << attribsOut["vectorsize"] << endl;
	 
         bool success = true;
	 
         // Parse variable vector size:
         uint64_t vectorSize = 1;
         mapIt = attribsOut.find("vectorsize");
         if (mapIt == attribsOut.end()) {
            debug3 << "VLSV\t\t ERROR: Variable '" << it << "' XML tag does not contain attribute 'vectorsize'" << endl;
            success = false;
         } else {
            vectorSize = atoi(mapIt->second.c_str());
         }
	 
         // By default assume that variable is zone-centered:
         vlsvplugin::VariableCentering centering = vlsvplugin::ZONE_CENTERED;
         mapIt = attribsOut.find("centering");
         if (mapIt != attribsOut.end()) {
            if (mapIt->second == "zone") centering = vlsvplugin::ZONE_CENTERED;
            else if (mapIt->second == "node") centering = vlsvplugin::NODE_CENTERED;
            else {
               debug3 << "VLSV\t\t ERROR: Variable '" << it << "' has unsupported centering!" << endl;
               success = false;
            }
         }
      
         if (success == false) continue;
	 
         MeshMetadata::variableMetadata.push_back(vlsvplugin::VariableMetadata(centering,it,vectorSize));
      }

      meshMetadataRead = true;
      return true;
   }

   bool VisitUCDMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      return VisitMeshMetadata::readDomainMetadata(vlsvReader);
   }
      
} // namespace vlsvplugin
