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

#include <mesh_metadata.h>

#include <DebugStream.h>

using namespace std;

namespace vlsvplugin {
   
   VariableMetadata::VariableMetadata(VariableCentering centering,const std::string& name,uint64_t vectorSize):
     centering(centering),name(name),vectorSize(vectorSize) { }
   
   MeshMetadata::MeshMetadata() {
      blockSize = 1;
      domainMetadataRead = false;
      meshMetadataRead = false;
      geometry = vlsv::geometry::UNKNOWN;
      name   = "";
      N_domains    = 1;
      N_ghostNodes = 0;
      N_ghostZones = 0;
      N_localNodes = 0;
      N_localZones = 0;
      N_totalNodes = 0;
      N_totalZones = 0;
      xLabel = "x-coordinate";
      yLabel = "y-coordinate";
      zLabel = "z-coordinate";
      xPeriodic = false;
      yPeriodic = false;
      zPeriodic = false;
      xUnits = "";
      yUnits = "";
      zUnits = "";
      maxRefinementLevel = 0;

      transformName = "";
      for (int i=0; i<16; ++i) transform[i] = 0;
      transform[0 ] = 1;
      transform[5 ] = 1;
      transform[10] = 1;
      transform[15] = 1;
   }
   
   MeshMetadata::~MeshMetadata() { }
   
   bool MeshMetadata::checkVlsvMeshType(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      // Check that we are reading multi-domain unstructured mesh metadata:
      auto it = attribs.find("type");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
         return false;
      }
      if (it->second != getCorrectVlsvMeshType()) {
         debug3 << "VLSV\t\t ERROR: Mesh type is '" << it->second << "', should be '" << getCorrectVlsvMeshType() << "'" << endl;
         return false;
      }
      return true;
   }

   void MeshMetadata::getBlockWidths(uint64_t& blockWidthX,uint64_t& blockWidthY,uint64_t& blockWidthZ) const {
      if (meshBoundingBox.size() < 6) {
         blockWidthX = 1;
         blockWidthY = 1;
         blockWidthZ = 1;
      } else {
         blockWidthX = meshBoundingBox[3];
         blockWidthY = meshBoundingBox[4];
         blockWidthZ = meshBoundingBox[5];
      }
   }

   uint64_t MeshMetadata::getMaximumRefinementLevel() const {return maxRefinementLevel;}

   const std::vector<uint64_t>& MeshMetadata::getMeshBoundingBox() const {return meshBoundingBox;}

   void MeshMetadata::getMeshPeriodicity(bool& xPeriodic,bool& yPeriodic,bool& zPeriodic) const {
      xPeriodic = this->xPeriodic;
      yPeriodic = this->yPeriodic;
      zPeriodic = this->zPeriodic;
   }
   
   const vlsv::geometry::type& MeshMetadata::getMeshGeometry() const {return geometry;}

   std::string MeshMetadata::getName() const {return name;}

   uint64_t MeshMetadata::getNodeDomainOffset(uint64_t domain) const {return nodeDomainOffsets[domain];}

   uint64_t MeshMetadata::getZoneDomainOffset(uint64_t domain) const {return zoneDomainOffsets[domain];}

   uint64_t MeshMetadata::getNumberOfDomains() const {return N_domains;}

   //uint64_t MeshMetadata::getNumberOfGhostNodes() const {return N_ghostNodes;}
   
   uint64_t MeshMetadata::getNumberOfGhostNodes(uint64_t domain) const {
      return nodeGhostOffsets[domain+1]-nodeGhostOffsets[domain];
   }

   //uint64_t MeshMetadata::getNumberOfGhostZones() const {return N_ghostZones;}

   uint64_t MeshMetadata::getNumberOfGhostZones(uint64_t domain) const {
      return blockSize*(zoneGhostOffsets[domain+1]-zoneGhostOffsets[domain]);
   }

   //uint64_t MeshMetadata::getNumberOfLocalNodes() const {return N_localNodes;}
   
   uint64_t MeshMetadata::getNumberOfLocalNodes(uint64_t domain) const {
      return getNumberOfTotalNodes(domain)-getNumberOfGhostNodes(domain);
   }

   //uint64_t MeshMetadata::getNumberOfLocalZones() const {return N_localZones;}

   uint64_t MeshMetadata::getNumberOfLocalZones(uint64_t domain) const {
      return getNumberOfTotalZones(domain)-getNumberOfGhostZones(domain);
   }

   //uint64_t MeshMetadata::getNumberOfTotalNodes() const {return N_totalNodes;}
   
   uint64_t MeshMetadata::getNumberOfTotalNodes(uint64_t domain) const {
      return nodeVariableOffsets[domain+1]-nodeVariableOffsets[domain];
   }

   uint64_t MeshMetadata::getNumberOfTotalZones() const {return N_totalZones;}

   uint64_t MeshMetadata::getNumberOfTotalZones(uint64_t domain) const {
      return blockSize*(zoneDomainOffsets[domain+1]-zoneDomainOffsets[domain]);
   }

   int MeshMetadata::getSpatialDimension() const {return spatialDimension;}
   
   int MeshMetadata::getTopologicalDimension() const {return topologicalDimension;}

   const double* MeshMetadata::getTransform() const {return transform;}

   const std::vector<VariableMetadata>& MeshMetadata::getVariables() const {return variableMetadata;}
   
   void MeshMetadata::getAxisLabels(std::string& xLabel,std::string& yLabel,std::string& zLabel) const {
      xLabel = this->xLabel;
      yLabel = this->yLabel;
      zLabel = this->zLabel;
   }
   
   void MeshMetadata::getAxisUnits(std::string& xUnits,std::string& yUnits,std::string& zUnits) const {
      xUnits = this->xUnits;
      yUnits = this->yUnits;
      zUnits = this->zUnits;
   }

   bool MeshMetadata::hasTransform() const {
      if (transformName.size() > 0) return true;
      return false;
   }

   bool MeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Check that reader exists:
      meshMetadataRead = false;
      if (vlsvReader == NULL) return meshMetadataRead;

      bool success = true;
      
      // Check that we are reading correct mesh type:
      if (checkVlsvMeshType(vlsvReader,attribs) == false) return false;

      // Get mandatory mesh parameters
      // Get mesh name:
      auto it = attribs.find("name"); 
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'name'" << endl;
         success = false; 
      } else {
         name = it->second;
      }
      // Get mesh type:
      it = attribs.find("type");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does not have attribute 'type'" << endl;
         success = false;
      } else {
         vlsvMeshType = vlsv::getMeshType(it->second);
      }
      // Get mesh geometry:
      if (readMeshGeometry(vlsvReader,attribs) == false) success = false;

      // Get variables that belong to this mesh:
      if (readVariables(vlsvReader,attribs) == false) success = false;

      // Get optional mesh parameters:
      transformName = "";
      it = attribs.find("xlabel"); if (it != attribs.end()) xLabel = it->second;
      it = attribs.find("ylabel"); if (it != attribs.end()) yLabel = it->second;
      it = attribs.find("zlabel"); if (it != attribs.end()) zLabel = it->second;
      it = attribs.find("xunits"); if (it != attribs.end()) xUnits = it->second;
      it = attribs.find("yunits"); if (it != attribs.end()) yUnits = it->second;
      it = attribs.find("zunits"); if (it != attribs.end()) zUnits = it->second;
      it = attribs.find("xperiodic"); if (it != attribs.end()) if (it->second == "yes") xPeriodic = true;
      it = attribs.find("yperiodic"); if (it != attribs.end()) if (it->second == "yes") yPeriodic = true;
      it = attribs.find("zperiodic"); if (it != attribs.end()) if (it->second == "yes") zPeriodic = true;
      it = attribs.find("max_refinement_level"); if (it != attribs.end()) maxRefinementLevel = atoi(it->second.c_str());
      it = attribs.find("transform"); if (it != attribs.end()) transformName = it->second;

      // Read the transform matrix, if defined:
      if (hasTransform() == true) {
         // Attempt to read transform matrix tag:
         list<pair<string,string> > attribsIn;
         map<string,string> attribsOut;

         attribsIn.push_back(make_pair("name",transformName));
         if (vlsvReader->getArrayAttributes("TRANSFORM",attribsIn,attribsOut) == false) {
            debug3 << "VLSV\t\t ERROR: Failed to read array 'TRANSFORM' attributes for mesh '" << getName() << "' from VLSV file" << endl;
            transformName = ""; return false;
         }

         // Check arraysize and vectorsize for correctness:
         auto xml = attribsOut.find("arraysize");
         if (atoi(xml->second.c_str()) != 16) {
            debug3 << "VLSV\t\t ERROR: Transform matrix has wrong arraysize '" << atoi(xml->second.c_str()) << "', should be 16" << endl;
            transformName = ""; return false;
         }
         
         xml = attribsOut.find("vectorsize");
         if (atoi(xml->second.c_str()) != 1) {
            debug3 << "VLSV\t\t ERROR: Transform matrix has wrong vectorsize '" << atoi(xml->second.c_str()) << "', should be 1" << endl;
            transformName = ""; return false;
         }

         // Read the matrix components:
         double* ptr = transform;
         if (vlsvReader->read("TRANSFORM",attribsIn,0,16,ptr,false) == false) {
            debug3 << "VLSV\t\t ERROR: Failed to read the transform matrix for mesh '" << getName() << "'" << endl;
            transformName = ""; return false;
         }
      }

      if (success == true) meshMetadataRead = true;
      return meshMetadataRead;
   }

   /** Read mesh geometry from file and store it to MeshMetadata::geometry.
    * @param vlsvReader vlsv::Reader that has input file open.
    * @param attribs XML attributes for the array 'MESH'.
    * @return If true, mesh geometry was successfully read.*/
   bool MeshMetadata::readMeshGeometry(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Get mesh geometry:
      auto it = attribs.find("geometry");
      if (it == attribs.end()) {
         geometry = vlsv::geometry::CARTESIAN;
      } else {
         geometry = vlsv::getMeshGeometry(it->second);
         if (geometry == vlsv::geometry::UNKNOWN) geometry = vlsv::geometry::CARTESIAN;
      }
      return true;
   }

   /** Read variable metadata from input file and store to MeshMetadata::variableMetadata.
    * @param vlsvReader vlsv::Reader that has input file open.
    * @param attribs XML attributes for the array 'MESH'.
    * @return If true, variables were read successfully.*/
   bool MeshMetadata::readVariables(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
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
	 
         variableMetadata.push_back(vlsvplugin::VariableMetadata(centering,it,vectorSize));
      }

      return true;
   }
}


