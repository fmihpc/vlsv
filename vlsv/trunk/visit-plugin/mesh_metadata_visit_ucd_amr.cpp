/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#include <DebugStream.h>
#include <list>

using namespace std;

namespace vlsvplugin {

   VisitUCDAMRMetadata::VisitUCDAMRMetadata() {
      domainMetadataRead = false;
      meshMetadataRead = false;
      domainOffsets = NULL;
      ghostOffsets = NULL;
      meshBoundingBox = NULL;
      variableOffsets = NULL;
   }
   
   VisitUCDAMRMetadata::~VisitUCDAMRMetadata() {
      delete [] domainOffsets; domainOffsets = NULL;
      delete [] ghostOffsets; ghostOffsets = NULL;
      delete [] meshBoundingBox; meshBoundingBox = NULL;
      delete [] variableOffsets; variableOffsets = NULL;
   }
   
   bool VisitUCDAMRMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
						  const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitUCDAMRMetadata::getDomainInfo called, domain: " << domain << endl;
      
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
	 debug3 << "VLSV\t\t ERROR: vlsv::Reader is NULL" << endl;
	 return false;
      }
	     
      // Read domain info:
      if (readDomains(vlsvReader) == false) {
	 debug3 << "VLSV\t\t ERROR: Failed to read ucd amr mesh domains" << endl;
	 return false;
      }
      
      debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
      debug4 << (this->domainOffsets)[domain] << ' ' << (this->domainOffsets)[domain+1];
      debug4 << " ghost offsets: " << (this->ghostOffsets)[domain] << ' ' << (this->ghostOffsets)[domain+1];
      debug4 << " variable offsets: " << (this->variableOffsets)[domain] << ' ' << (this->variableOffsets)[domain+1];
      debug4 << endl;
      
      domainOffsets   = this->domainOffsets;
      ghostOffsets    = this->ghostOffsets;
      variableOffsets = this->variableOffsets;
      return true;
   }

   uint64_t VisitUCDAMRMetadata::getBlockSize() const {return blockSize;}
   
   const uint64_t* VisitUCDAMRMetadata::getDomainOffsets() {return domainOffsets;}
   
   const uint64_t* VisitUCDAMRMetadata::getGhostOffsets() {return ghostOffsets;}
   
   const uint64_t* VisitUCDAMRMetadata::getMeshBoundingBox() {return meshBoundingBox;}

   const vlsv::geometry::type& VisitUCDAMRMetadata::getMeshGeometry() const {return geometry;}
   
   uint64_t VisitUCDAMRMetadata::getNumberOfGhostNodes(int domain) const {
      return 0;
   }
   
   uint64_t VisitUCDAMRMetadata::getNumberOfGhostZones(int domain) const {
      return blockSize*(ghostOffsets[domain+1]-ghostOffsets[domain]);
   }
   
   uint64_t VisitUCDAMRMetadata::getNumberOfLocalNodes(int domain) const {
      return 0;
   }
   
   uint64_t VisitUCDAMRMetadata::getNumberOfLocalZones(int domain) const {
      return getNumberOfTotalZones(domain)-getNumberOfGhostZones(domain);
   }
   
   uint64_t VisitUCDAMRMetadata::getNumberOfTotalNodes(int domain) const {
      return 0;
   }
   
   uint64_t VisitUCDAMRMetadata::getNumberOfTotalZones(int domain) const {
      return blockSize*(domainOffsets[domain+1]-domainOffsets[domain]);
   }
   
   const uint64_t* VisitUCDAMRMetadata::getVariableOffsets() {return variableOffsets;}

   bool VisitUCDAMRMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitUCDAMRMetadata::read called" << endl;
      
      // Exit if mesh metadata has already been read:
      if (meshMetadataRead == true) return true;
      
      // Call superclass read function:
      if (VisitMeshMetadata::read(vlsvReader,attribs) == false) return false;
      
      // Check that we are reading multi-domain unstructured mesh metadata:
      map<string,string>::const_iterator it = attribs.find("type");
      if (it == attribs.end()) {
	 debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
	 return false;
      }
      if (it->second != vlsv::mesh::STRING_UCD_AMR) {
	 debug3 << "VLSV\t\t ERROR: Mesh type is '" << it->second << "', should be '" << vlsv::mesh::STRING_UCD_AMR << "'" << endl;
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
      if (it == attribs.end()) geometry = vlsv::geometry::CARTESIAN;
      else {
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
      for (set<string>::const_iterator it=variableNames.begin(); it!=variableNames.end(); ++it) {
	 map<string,string> attribsOut;
	 map<string,string>::const_iterator mapIt;
	 list<pair<string,string> > attribsIn;
	 attribsIn.push_back(make_pair("mesh",getName()));
	 attribsIn.push_back(make_pair("name",*it));
	 
	 // Skip variables belonging to other meshes:
	 if (vlsvReader->getArrayAttributes("VARIABLE",attribsIn,attribsOut) == false) continue;
	 debug4 << "VLSV\t\t\t '" << *it << "' vectorsize: " << attribsOut["vectorsize"] << endl;
	 
	 bool success = true;
	 
	 // Parse variable vector size:
	 uint64_t vectorSize = 1;
	 mapIt = attribsOut.find("vectorsize");
	 if (mapIt == attribsOut.end()) {
	    debug3 << "VLSV\t\t ERROR: Variable '" << *it << "' XML tag does not contain attribute 'vectorsize'" << endl;
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
	       debug3 << "VLSV\t\t ERROR: Variable '" << *it << "' has unsupported centering!" << endl;
	       success = false;
	    }
	 }
      
	 if (success == false) continue;
	 
	 MeshMetadata::variableMetadata.push_back(vlsvplugin::VariableMetadata(centering,*it,vectorSize));
      }

      meshMetadataRead = true;
      return true;
   }

   bool VisitUCDAMRMetadata::readDomains(vlsv::Reader* vlsvReader) {
      // Exit if domain metadata has already been read:
      if (domainMetadataRead == true) return true;
      debug2 << "VLSV\t\t VisitUCDAMRMetadata::readDomains called" << endl;
      
      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
	 debug3 << "VLSV\t\t ERROR: vlsv::Reader is NULL" << endl;
	 return false;
      }
      
      // Read mesh bounding box:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",name));
      delete [] meshBoundingBox; meshBoundingBox = NULL;
      if (vlsvReader->read("MESH_BBOX",attribs,0,6,meshBoundingBox) == false) {
	 debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_BBOX'" << endl;
	 return false;
      }
      blockSize = meshBoundingBox[3]*meshBoundingBox[4]*meshBoundingBox[5];
      debug4 << "VLSV\t\t Mesh size in blocks: " << meshBoundingBox[0] << ' ' << meshBoundingBox[1] << ' ' << meshBoundingBox[2];
      debug4 << " block sizes: " << meshBoundingBox[3] << ' ' << meshBoundingBox[4] << ' ' << meshBoundingBox[5] << endl;
      
      // Read domain info:
      int64_t* domainInfo = NULL;
      if (vlsvReader->read("MESH_DOMAIN_SIZES",attribs,0,VisitMeshMetadata::N_domains,domainInfo) == false) {
	 debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_ZONES'" << endl;
	 delete [] domainInfo; domainInfo = NULL;
	 return false;
      }
      
      // Calculate offsets where data for each domain begins:
      delete [] domainOffsets;   domainOffsets   = new uint64_t[VisitMeshMetadata::N_domains+1];
      delete [] ghostOffsets;    ghostOffsets    = new uint64_t[VisitMeshMetadata::N_domains+1];
      delete [] variableOffsets; variableOffsets = new uint64_t[VisitMeshMetadata::N_domains+1];
      domainOffsets[0]   = 0;
      ghostOffsets[0]    = 0;
      variableOffsets[0] = 0;
      for (int64_t i=0; i<VisitMeshMetadata::N_domains; ++i) {
	 domainOffsets[i+1]   = domainOffsets[i] + domainInfo[i*2];
	 ghostOffsets[i+1]    = ghostOffsets[i] + domainInfo[i*2+1];
	 variableOffsets[i+1] = variableOffsets[i] + domainInfo[i*2+0]-domainInfo[i*2+1];
      }
      delete [] domainInfo; domainInfo = NULL;
      
      // Compute total number of real and ghost cells:
      MeshMetadata::N_ghostZones = ghostOffsets[VisitMeshMetadata::N_domains];
      MeshMetadata::N_localZones  = variableOffsets[VisitMeshMetadata::N_domains];
   
      domainMetadataRead = true;
      return true;
   }
      
} // namespace vlsvplugin
