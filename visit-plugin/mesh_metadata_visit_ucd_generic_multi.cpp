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

#include <mesh_metadata_visit_ucd_generic_multi.h>

#include <list>

using namespace std;

namespace vlsvplugin {

   UCDGenericMultiMeshMetadata::UCDGenericMultiMeshMetadata(): MeshMetadata() { }
   
   UCDGenericMultiMeshMetadata::~UCDGenericMultiMeshMetadata() { }
   
   const std::string& UCDGenericMultiMeshMetadata::getCorrectVlsvMeshType() const {
      return vlsv::mesh::STRING_UCD_GENERIC_MULTI;
   }

   bool UCDGenericMultiMeshMetadata::getDomainInfoNodes(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
							                            const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         stringstream ss;
         ss << "ERROR: vlsv::Reader is NULL";
         return MeshMetadata::exitWithError(ss);
      }
      
      // Read domain info:
      if (readDomains(vlsvReader) == false) return false;

      domainOffsets   = this->nodeDomainOffsets.data();
      ghostOffsets    = this->nodeGhostOffsets.data();
      variableOffsets = this->nodeVariableOffsets.data();
      return true;
   }
   
   bool UCDGenericMultiMeshMetadata::getDomainInfoZones(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
							                            const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         stringstream ss;
         ss << "ERROR: vlsv::Reader is NULL";
         return MeshMetadata::exitWithError(ss);
      }
	     
      // Read domain info:
      if (readDomains(vlsvReader) == false) return false;
      
      domainOffsets   = this->zoneDomainOffsets.data();
      ghostOffsets    = this->zoneGhostOffsets.data();
      variableOffsets = this->zoneVariableOffsets.data();
      return true;
   }

   const uint64_t UCDGenericMultiMeshMetadata::getCellConnectivitySize(int domain) const {
      return zoneConnectivityOffsets[domain+1]-zoneConnectivityOffsets[domain];
   }

   vlsv::datatype::type UCDGenericMultiMeshMetadata::getNodeDatatype() const {      
      return nodeDatatype;
   }
   
   int UCDGenericMultiMeshMetadata::getNodeDataSize() const {      
      return nodeDataSize;
   }
   
   bool UCDGenericMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function:
      if (MeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;
      
      // Parse values from XML tag 'MESH' associated with this mesh.
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
      }
      
      // Figure out total number of cells in the mesh:
      it = attribs.find("arraysize");
      if (it == attribs.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH' does not have an XML attribute called 'arraysize'.";
         return MeshMetadata::exitWithError(ss);
      }
      MeshMetadata::N_totalZones = atoi(it->second.c_str());

      // Read XML tag 'MESH_DOMAIN_SIZES' to figure out how many domains the mesh has:
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
      
      // Read XML tag 'MESH_NODE_CRDS' to figure out node coordinate datatype and size:
      attribsOut.clear();
      if (vlsvReader->getArrayAttributes("MESH_NODE_CRDS",attribsIn,attribsOut) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_NODE_CRDS' attributes for mesh '" << getName() << "'. " << endl;
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'.";
         return MeshMetadata::exitWithError(ss);
      }
      it = attribsOut.find("datatype");
      if (it == attribsOut.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_NODE_CRDS' XML tag does not have attribute 'datatype'";
         return MeshMetadata::exitWithError(ss);
      } else {
         nodeDatatype = vlsv::getVLSVDatatype(it->second);
      }      
      
      it = attribsOut.find("datasize");
      if (it == attribsOut.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_NODE_CRDS' XML tag does not have attribute 'datasize'";
         return MeshMetadata::exitWithError(ss);
      } else {
         nodeDataSize = atoi(it->second.c_str());
      }      

      meshMetadataRead = true;
      return true;
   }

   bool UCDGenericMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
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
      meshBoundingBox.resize(vlsv::ucdgenericmulti::bbox::SIZE);
      uint64_t* ptr = meshBoundingBox.data();
      if (vlsvReader->read("MESH_BBOX",attribs,0,vlsv::ucdgenericmulti::bbox::SIZE,ptr,false) == false) {
         stringstream ss;
         ss << "VLSV\t\t ERROR: Failed to read array 'MESH_BBOX'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'";
         return MeshMetadata::exitWithError(ss);
      }
      blockSize =
        meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
        * meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
        * meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      
      // Check that MESH_DOMAIN_SIZES has correct size:
      map<string,string> attribsOut;
      if (vlsvReader->getArrayAttributes("MESH_DOMAIN_SIZES",attribs,attribsOut) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_DOMAIN_SIZES' attributes for mesh '" << getName() << "'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'";
         return MeshMetadata::exitWithError(ss);
      }
      auto it = attribsOut.find("vectorsize");
      if (it == attribsOut.end()) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_DOMAIN_SIZES' XML tag does not have attribute 'vectorsize'";
         return MeshMetadata::exitWithError(ss);
      }
      if (atol(it->second.c_str()) != vlsv::ucdgenericmulti::domainsizes::SIZE) {
         stringstream ss;
         ss << "ERROR: Array 'MESH_DOMAIN_SIZES' has incorrect vectorsize '" << it->second
            << "' when it should be '" << vlsv::ucdgenericmulti::domainsizes::SIZE << "'";
         return MeshMetadata::exitWithError(ss);
      }
      
      // Read domain info:
      int64_t* domainInfo = NULL;
      if (vlsvReader->read("MESH_DOMAIN_SIZES",attribs,0,MeshMetadata::N_domains,domainInfo) == false) {
         delete [] domainInfo; domainInfo = NULL;
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_DOMAIN_SIZES'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'";
         return MeshMetadata::exitWithError(ss);
      }
      
      // Calculate offsets where data for each domain begins:
      nodeDomainOffsets.resize(MeshMetadata::N_domains+1); nodeDomainOffsets.shrink_to_fit();
      nodeGhostOffsets.resize(MeshMetadata::N_domains+1); nodeGhostOffsets.shrink_to_fit();
      nodeVariableOffsets.resize(MeshMetadata::N_domains+1); nodeVariableOffsets.shrink_to_fit();
      zoneGhostOffsets.resize(MeshMetadata::N_domains+1); zoneGhostOffsets.shrink_to_fit();
      zoneDomainOffsets.resize(MeshMetadata::N_domains+1); zoneDomainOffsets.shrink_to_fit();
      zoneVariableOffsets.resize(MeshMetadata::N_domains+1); zoneVariableOffsets.shrink_to_fit();
      nodeDomainOffsets[0]   = 0;
      zoneDomainOffsets[0]   = 0;
      nodeGhostOffsets[0]    = 0;
      zoneGhostOffsets[0]    = 0;
      zoneVariableOffsets[0] = 0;
      nodeVariableOffsets[0] = 0;
      MeshMetadata::N_localZones = 0;
      for (uint64_t i=0; i<MeshMetadata::N_domains; ++i) {
         // Count number of local zones in domain i:
         const uint64_t localZones = 
           domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS]
           - domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
         MeshMetadata::N_localZones += localZones;
         
         zoneGhostOffsets[i+1] = zoneGhostOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
         zoneDomainOffsets[i+1] = zoneDomainOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS]
           - domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
         nodeDomainOffsets[i+1] = nodeDomainOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_NODES];
         nodeGhostOffsets[i+1] = nodeGhostOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_NODES];
         
         // Offset to domain i's cell variable values is the sum of local blocks, i.e., 
         // total cells - ghost cells values, in all previous domains:
         zoneVariableOffsets[i+1] = zoneVariableOffsets[i] 
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS]
           - domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
         
         // Offset to domain i's node variable values is the sum of local nodes, i.e., 
         // total nodes - ghost nodes values, in all previous domains:
         nodeVariableOffsets[i+1] = nodeVariableOffsets[i] 
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_NODES]
           - domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_NODES];
      }      
      delete [] domainInfo; domainInfo = NULL;
      MeshMetadata::N_ghostZones = zoneGhostOffsets[MeshMetadata::N_domains];

      // Read offsets into cell connectivity array:
      if (vlsvReader->read("MESH_OFFSETS",attribs,0,MeshMetadata::N_domains,domainInfo) == false) {
         delete [] domainInfo; domainInfo = NULL;
         stringstream ss;
         ss << "ERROR: Failed to read array 'MESH_CELL_OFFSETS'. ";
         ss << "vlsv::Reader says '" << vlsvReader->getLastError() << "'";
         return MeshMetadata::exitWithError(ss);
      }
      
      zoneConnectivityOffsets.resize(MeshMetadata::N_domains+1); zoneConnectivityOffsets.shrink_to_fit();
      for (uint64_t i=0; i<MeshMetadata::N_domains; ++i) {
         zoneConnectivityOffsets[i+1] = zoneConnectivityOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::offsets::SIZE+vlsv::ucdgenericmulti::offsets::ZONE_ENTRIES];
      }
      delete [] domainInfo; domainInfo = NULL;
      domainMetadataRead = true;
      return true;
   }

} // namespace vlsvplugin
