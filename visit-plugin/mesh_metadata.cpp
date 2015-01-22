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

#include <mesh_metadata.h>

#include <DebugStream.h>

using namespace std;

namespace vlsvplugin {
   
   VariableMetadata::VariableMetadata(VariableCentering centering,const std::string& name,int vectorSize): 
     centering(centering),name(name),vectorSize(vectorSize) { }
   
   MeshMetadata::MeshMetadata() {
      name   = "";
      arraySize = 0;
      dataSize = 0;
      vectorSize = 0;
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

   uint64_t MeshMetadata::getArraySize() const {return arraySize;}
   
   uint64_t MeshMetadata::getDataSize() const {return dataSize;}
   
   vlsv::datatype::type MeshMetadata::getDatatype() const {return datatype;}
   
   void MeshMetadata::getMeshPeriodicity(bool& xPeriodic,bool& yPeriodic,bool& zPeriodic) const {
      xPeriodic = this->xPeriodic;
      yPeriodic = this->yPeriodic;
      zPeriodic = this->zPeriodic;
   }
   
   uint64_t MeshMetadata::getMaximumRefinementLevel() const {return maxRefinementLevel;}
   
   std::string MeshMetadata::getName() const {return name;}

   uint64_t MeshMetadata::getNumberOfGhostNodes() const {return N_ghostNodes;}
   
   uint64_t MeshMetadata::getNumberOfGhostZones() const {return N_ghostZones;}

   uint64_t MeshMetadata::getNumberOfLocalNodes() const {return N_localNodes;}
   
   uint64_t MeshMetadata::getNumberOfLocalZones() const {return N_localZones;}

   uint64_t MeshMetadata::getNumberOfTotalNodes() const {return N_totalNodes;}
   
   uint64_t MeshMetadata::getNumberOfTotalZones() const {return N_totalZones;}

   const double* MeshMetadata::getTransform() const {return transform;}

   const std::vector<VariableMetadata>& MeshMetadata::getVariables() const {return variableMetadata;}
   
   uint64_t MeshMetadata::getVectorSize() const {return vectorSize;}
   
   std::string MeshMetadata::getXLabel() const {return xLabel;}
   
   std::string MeshMetadata::getYLabel() const {return yLabel;}
   
   std::string MeshMetadata::getZLabel() const {return zLabel;}
   
   std::string MeshMetadata::getXUnits() const {return xUnits;}
   
   std::string MeshMetadata::getYUnits() const {return yUnits;}   

   std::string MeshMetadata::getZUnits() const {return zUnits;}

   bool MeshMetadata::hasTransform() const {
      if (transformName.size() > 0) return true;
      return false;
   }

   bool MeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      // Check that reader exists:
      if (vlsvReader == NULL) return false;
      
      bool success = true;
      map<string,string>::const_iterator it;
      
      // Get mandatory mesh parameters
      // Get mesh name:
      it = attribs.find("name"); 
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'name'" << endl;
         success = false; 
      } else {
         name = it->second;
      }
      
      // Get arraysize, vectorsize, datasize, datatype:
      it = attribs.find("arraysize"); 
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'arraysize'" << endl;
         success = false;
      } else {
         arraySize = atoi(it->second.c_str());
      }
      
      it = attribs.find("vectorsize");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'vectorsize'" << endl;
         success = false;
      } else {
         vectorSize = atoi(it->second.c_str());
      }
      
      it = attribs.find("datasize");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'datasize'" << endl;
         success = false;
      } else {
         dataSize = atoi(it->second.c_str());
      }
      
      it = attribs.find("datatype");
      if (it == attribs.end()) {
         debug3 << "VLSV\t\t ERROR: XML tag does have attribute 'datatype'" << endl;
         success = false;
      } else {
         datatype = vlsv::getVLSVDatatype(it->second);
      }

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
         map<string,string>::const_iterator xml = attribsOut.find("arraysize");
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

      return success;
   }
}


