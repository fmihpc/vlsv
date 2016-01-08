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

#include <mesh_metadata_visit_ucd_generic_multi.h>

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   VisitUCDGenericMultiMeshMetadata::VisitUCDGenericMultiMeshMetadata(): VisitMeshMetadata() { }
   
   VisitUCDGenericMultiMeshMetadata::~VisitUCDGenericMultiMeshMetadata() { }
   
   bool VisitUCDGenericMultiMeshMetadata::getDomainInfoNodes(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
							     const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::getDomainInfoNodes called, domain: " << domain << endl;
       
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
      
      debug4 << "VLSV\t\t Domain: " << domain << " node offsets: ";
      debug4 << nodeDomainOffsets[domain] << ' ' << nodeDomainOffsets[domain+1];
      debug4 << " ghost offsets: " << nodeGhostOffsets[domain] << ' ' << nodeGhostOffsets[domain+1];
      debug4 << " variable offsets: " << nodeVariableOffsets[domain] << ' ' << nodeVariableOffsets[domain+1];
      debug4 << endl;

      domainOffsets   = this->nodeDomainOffsets.data();
      ghostOffsets    = this->nodeGhostOffsets.data();
      variableOffsets = this->nodeVariableOffsets.data();
      return true;
   }
   
   bool VisitUCDGenericMultiMeshMetadata::getDomainInfoZones(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
							     const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::getDomainInfoZones called, domain: " << domain << endl;
      
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
      debug4 << zoneDomainOffsets[domain] << ' ' << zoneDomainOffsets[domain+1];
      debug4 << " ghost offsets: " << zoneGhostOffsets[domain] << ' ' << zoneGhostOffsets[domain+1];
      debug4 << " variable offsets: " << (this->zoneVariableOffsets)[domain] << ' ' << (this->zoneVariableOffsets)[domain+1];
      debug4 << endl;
      
      domainOffsets   = this->zoneDomainOffsets.data();
      ghostOffsets    = this->zoneGhostOffsets.data();
      variableOffsets = this->zoneVariableOffsets.data();
      return true;
   }

   const uint64_t VisitUCDGenericMultiMeshMetadata::getCellConnectivitySize(int domain) const {
      return zoneDomainOffsets[domain+1]-zoneDomainOffsets[domain];
   }

   vlsv::datatype::type VisitUCDGenericMultiMeshMetadata::getNodeDatatype() const {
      return nodeDatatype;
   }
   
   int VisitUCDGenericMultiMeshMetadata::getNodeDataSize() const {return nodeDataSize;}
   
   bool VisitUCDGenericMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::read called" << endl;
      
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function:
      if (VisitMeshMetadata::read(vlsvReader,attribs) == false) return false;
      meshMetadataRead = false;
      
      // Parse values from XML tag 'MESH' associated with this mesh.
      // Check that we are reading multi-domain unstructured mesh metadata:
      auto it = attribs.find("type");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
         return false;
      }
      if (it->second != vlsv::mesh::STRING_UCD_GENERIC_MULTI) {
         debug3 << "VLSV\t\t ERROR: Mesh type is '" << it->second << "', should be '" << vlsv::mesh::STRING_UCD_GENERIC_MULTI << "'" << endl;
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
      }
      
      // Figure out total number of cells in the mesh:
      it = attribs.find("arraysize");
      if (it == attribs.end()) return false;
      MeshMetadata::N_totalZones = atoi(it->second.c_str());

      // Get mesh geometry:
      it = attribs.find("geometry");
      if (it == attribs.end()) geometry = vlsv::geometry::CARTESIAN;
      else {
         geometry = vlsv::getMeshGeometry(it->second);
         if (geometry == vlsv::geometry::UNKNOWN) geometry = vlsv::geometry::CARTESIAN;
      }

      // Read XML tag 'MESH_DOMAIN_SIZES' to figure out how many domains the mesh has:
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
      
      // Read XML tag 'MESH_NODE_CRDS' to figure out node coordinate datatype and size:
      attribsOut.clear();
      if (vlsvReader->getArrayAttributes("MESH_NODE_CRDS",attribsIn,attribsOut) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_NODE_CRDS' attributes for mesh '" << getName() << "' from VLSV file" << endl;
         return false;
      }
      it = attribsOut.find("datatype");
      if (it == attribsOut.end()) {
         debug3 << "VLSV\t\t ERROR: Array 'MESH_NODE_CRDS' XML tag does not have attribute 'datatype'" << endl;
         return false;
      } else {
         nodeDatatype = vlsv::getVLSVDatatype(it->second);
      }
      debug4 << "VLSV\t\t Node coordinate datatype: " << it->second << endl;
      
      it = attribsOut.find("datasize");
      if (it == attribsOut.end()) {
         debug3 << "VLSV\t\t ERROR: Array 'MESH_NODE_CRDS' XML tag does not have attribute 'datasize'" << endl;
         return false;
      } else {
         nodeDataSize = atoi(it->second.c_str());
      }
      debug4 << "VLSV\t\t Node coordinate datasize: " << nodeDataSize << endl;
      
      // Read 'VARIABLE' XML tags to parse names of variables in this mesh:
      set<string> variableNames;
      if (vlsvReader->getUniqueAttributeValues("VARIABLE","name",variableNames) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read variable names" << endl;
         return false;
      }
      
      // Remove variables that do not belong to this mesh:
      debug4 << "VLSV\t\t Found variables:" << endl;
      for (const auto& it : variableNames) {
         map<string,string> attribsOut;
         list<pair<string,string> > attribsIn;
         attribsIn.push_back(make_pair("mesh",getName()));
         attribsIn.push_back(make_pair("name",it));
         
         // Skip variables belonging to other meshes:
         if (vlsvReader->getArrayAttributes("VARIABLE",attribsIn,attribsOut) == false) continue;
         debug4 << "VLSV\t\t\t '" << it << "' vectorsize: " << attribsOut["vectorsize"] << endl;
         
         bool success = true;
         
         // Parse variable vector size:
         uint64_t vectorSize = 1;
         auto mapIt = attribsOut.find("vectorsize");
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

   bool VisitUCDGenericMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      // Exit if domain metadata has already been read:
      if (domainMetadataRead == true) return true;
      domainMetadataRead = false;
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::readDomains called" << endl;

      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
         debug3 << "VLSV\t\t ERROR: vlsv::Reader is NULL" << endl;
         return false;
      }
      
      // Read mesh bounding box:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",name));
      meshBoundingBox.resize(vlsv::ucdgenericmulti::bbox::SIZE);
      uint64_t* ptr = meshBoundingBox.data();
      if (vlsvReader->read("MESH_BBOX",attribs,0,vlsv::ucdgenericmulti::bbox::SIZE,ptr,false) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_BBOX'" << endl;
         return false;
      }
      blockSize =
        meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X]
        * meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y]
        * meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z];
      
      debug4 << "VLSV\t\t Mesh size in blocks: " 
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::X_BLOCKS] << ' ' 
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::Y_BLOCKS] << ' ' 
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::Z_BLOCKS];
      debug4 << " block sizes: " 
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X] << ' '
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y] << ' ' 
        << meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z] << endl;
      
      // Check that MESH_DOMAIN_SIZES has correct vectorsize (=4):
      map<string,string> attribsOut;
      if (vlsvReader->getArrayAttributes("MESH_DOMAIN_SIZES",attribs,attribsOut) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_DOMAIN_SIZES' attributes for mesh '" << getName() << "' from VLSV file" << endl;
         return false;
      }
      auto it = attribsOut.find("vectorsize");
      if (it == attribsOut.end()) {
         debug3 << "VLSV\t\t ERROR: Array 'MESH_DOMAIN_SIZES' XML tag does not have attribute 'vectorsize'" << endl;
         return false;
      }
      if (atol(it->second.c_str()) != vlsv::ucdgenericmulti::domainsizes::SIZE) {
         debug3 << "VLSV\t\t ERROR: Array 'MESH_DOMAIN_SIZES' has incorrect vectorsize '" << it->second
           << "' when it should be '" << vlsv::ucdgenericmulti::domainsizes::SIZE << "'" << endl;
         return false;
      }
      
      // Read domain info:
      int64_t* domainInfo = NULL;
      if (vlsvReader->read("MESH_DOMAIN_SIZES",attribs,0,VisitMeshMetadata::N_domains,domainInfo) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_DOMAIN_SIZES'" << endl;
         delete [] domainInfo; domainInfo = NULL;
         return false;
      }
      
      // Calculate offsets where data for each domain begins:
      nodeDomainOffsets.resize(VisitMeshMetadata::N_domains+1); nodeDomainOffsets.shrink_to_fit();
      nodeGhostOffsets.resize(VisitMeshMetadata::N_domains+1); nodeGhostOffsets.shrink_to_fit();
      nodeVariableOffsets.resize(VisitMeshMetadata::N_domains+1); nodeVariableOffsets.shrink_to_fit();
      zoneGhostOffsets.resize(VisitMeshMetadata::N_domains+1); zoneGhostOffsets.shrink_to_fit();
      zoneVariableOffsets.resize(VisitMeshMetadata::N_domains+1); zoneVariableOffsets.shrink_to_fit();
      nodeDomainOffsets[0]   = 0;
      nodeGhostOffsets[0]    = 0;
      zoneGhostOffsets[0]    = 0;
      zoneVariableOffsets[0] = 0;
      nodeVariableOffsets[0] = 0;
      MeshMetadata::N_localZones = 0;
      for (auto i=0; i<VisitMeshMetadata::N_domains; ++i) {
         // Count number of local zones in domain i:
         const uint64_t localZones = 
           domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::TOTAL_BLOCKS]
           - domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
         MeshMetadata::N_localZones += localZones;
         
         zoneGhostOffsets[i+1] = zoneGhostOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::domainsizes::SIZE+vlsv::ucdgenericmulti::domainsizes::GHOST_BLOCKS];
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
      MeshMetadata::N_ghostZones = zoneGhostOffsets[VisitMeshMetadata::N_domains];
      
      // Read offsets into cell connectivity array:
      if (vlsvReader->read("MESH_OFFSETS",attribs,0,VisitMeshMetadata::N_domains,domainInfo) == false) {
         debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_CELL_OFFSETS'" << endl;
         delete [] domainInfo; domainInfo = NULL;
         return false;
      }
      
      zoneDomainOffsets.resize(VisitMeshMetadata::N_domains+1); zoneDomainOffsets.shrink_to_fit();
      zoneDomainOffsets[0] = 0;
      for (auto i=0; i<VisitMeshMetadata::N_domains; ++i) {
         zoneDomainOffsets[i+1] = zoneDomainOffsets[i]
           + domainInfo[i*vlsv::ucdgenericmulti::offsets::SIZE+vlsv::ucdgenericmulti::offsets::ZONE_ENTRIES];
      }
      delete [] domainInfo; domainInfo = NULL;
      domainMetadataRead = true;
      return true;
   }

} // namespace vlsvplugin
