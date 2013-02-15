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
      N_ghostCells = 0;
      N_realCells  = 0;
      N_totalCells = 0;
      xLabel = "x-coordinate";
      yLabel = "y-coordinate";
      zLabel = "z-coordinate";
      xUnits = "";
      yUnits = "";
      zUnits = "";
   }
   
   MeshMetadata::~MeshMetadata() { }

   uint64_t MeshMetadata::getArraySize() const {return arraySize;}
   
   uint64_t MeshMetadata::getDataSize() const {return dataSize;}
   
   vlsv::datatype::type MeshMetadata::getDatatype() const {return datatype;}
   
   std::string MeshMetadata::getName() const {return name;}

   uint64_t MeshMetadata::getNumberOfGhostCells() const {return N_ghostCells;}

   uint64_t MeshMetadata::getNumberOfRealCells() const {return N_realCells;}

   uint64_t MeshMetadata::getNumberOfTotalCells() const {return N_totalCells;}
   
   const std::vector<VariableMetadata>& MeshMetadata::getVariables() const {return variableMetadata;}
   
   uint64_t MeshMetadata::getVectorSize() const {return vectorSize;}
   
   std::string MeshMetadata::getXLabel() const {return xLabel;}
   
   std::string MeshMetadata::getYLabel() const {return xLabel;}
   
   std::string MeshMetadata::getZLabel() const {return xLabel;}
   
   std::string MeshMetadata::getXUnits() const {return xUnits;}
   
   std::string MeshMetadata::getYUnits() const {return yUnits;}   

   std::string MeshMetadata::getZUnits() const {return zUnits;}
   
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
      it = attribs.find("xlabel"); if (it != attribs.end()) xLabel = it->second;
      it = attribs.find("ylabel"); if (it != attribs.end()) yLabel = it->second;
      it = attribs.find("zlabel"); if (it != attribs.end()) zLabel = it->second;
      it = attribs.find("xunit"); if (it != attribs.end()) xUnits = it->second;
      it = attribs.find("yunit"); if (it != attribs.end()) yUnits = it->second;
      it = attribs.find("zunit"); if (it != attribs.end()) zUnits = it->second;
      
      return success;
   }
}


