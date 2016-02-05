/** This file is part of VLSV file format.
 * 
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

#include <mesh_metadata_visit_classes.h>

#include <DebugStream.h>

using namespace std;

namespace vlsvplugin {
   VisitPointMeshMetadata::VisitPointMeshMetadata(): PointMeshMetadata(),VisitMeshMetadata() { }

   VisitPointMeshMetadata::~VisitPointMeshMetadata() { }

   avtMeshType VisitPointMeshMetadata::getAvtMeshType() const {return AVT_POINT_MESH;}
   
   std::string VisitPointMeshMetadata::getAvtMeshTypeString() const {return "AVT_POINT_MESH";}
   
   const std::string& VisitPointMeshMetadata::getCorrectVlsvMeshType() const {return PointMeshMetadata::getCorrectVlsvMeshType();}

   bool VisitPointMeshMetadata::readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return VisitMeshMetadata::readVariables(vlsv,attribs);
   }



   VisitUCDMultiMeshMetadata::VisitUCDMultiMeshMetadata(): UCDMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitUCDMultiMeshMetadata::~VisitUCDMultiMeshMetadata() { }

   avtMeshType VisitUCDMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitUCDMultiMeshMetadata::getCorrectVlsvMeshType() const {return UCDMultiMeshMetadata::getCorrectVlsvMeshType();}
   
   bool VisitUCDMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      return VisitMeshMetadata::readDomainMetadata(vlsvReader);
   }

   bool VisitUCDMultiMeshMetadata::readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return VisitMeshMetadata::readVariables(vlsv,attribs);
   }



   VisitQuadMultiMeshMetadata::VisitQuadMultiMeshMetadata(): QuadMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitQuadMultiMeshMetadata::~VisitQuadMultiMeshMetadata() { }

   avtMeshType VisitQuadMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitQuadMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitQuadMultiMeshMetadata::getCorrectVlsvMeshType() const {return QuadMultiMeshMetadata::getCorrectVlsvMeshType();}

   bool VisitQuadMultiMeshMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			                                      const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitQuadMultiMeshMetadata::getDomainInfo called, domain: " << domain << endl;

      const bool rvalue = QuadMultiMeshMetadata::getDomainInfo(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
         debug4 << MeshMetadata::zoneDomainOffsets[domain] << ' ' << MeshMetadata::zoneDomainOffsets[domain+1];
         debug4 << " ghost offsets: " << MeshMetadata::zoneGhostOffsets[domain] << ' ' << MeshMetadata::zoneGhostOffsets[domain+1];
         debug4 << " variable offsets: " << MeshMetadata::zoneVariableOffsets[domain] << ' ' << MeshMetadata::zoneVariableOffsets[domain+1];
         debug4 << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitQuadMultiMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitQuadMultiMeshMetadata::read called" << endl;

      const bool rvalue = QuadMultiMeshMetadata::read(vlsvReader,attribs);

      if (rvalue == true) {
         debug3 << "VLSV\t\t Mesh has " << MeshMetadata::N_domains << " domains" << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitQuadMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      return VisitMeshMetadata::readDomainMetadata(vlsvReader);
   }

   bool VisitQuadMultiMeshMetadata::readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return VisitMeshMetadata::readVariables(vlsv,attribs);
   }



   VisitUCDAMRMetadata::VisitUCDAMRMetadata(): UCDAMRMetadata(),VisitMeshMetadata() { }

   VisitUCDAMRMetadata::~VisitUCDAMRMetadata() { }

   avtMeshType VisitUCDAMRMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDAMRMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitUCDAMRMetadata::getCorrectVlsvMeshType() const {return UCDAMRMetadata::getCorrectVlsvMeshType();}
   
   bool VisitUCDAMRMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			                               const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t VisitUCDAMRMetadata::getDomainInfo called, domain: " << domain << endl;

      const bool rvalue = UCDAMRMetadata::getDomainInfo(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
         debug4 << MeshMetadata::zoneDomainOffsets[domain] << ' ' << MeshMetadata::zoneDomainOffsets[domain+1];
         debug4 << " ghost offsets: " << MeshMetadata::zoneGhostOffsets[domain] << ' ' << MeshMetadata::zoneGhostOffsets[domain+1];
         debug4 << " variable offsets: " << MeshMetadata::zoneVariableOffsets[domain] << ' ' << MeshMetadata::zoneVariableOffsets[domain+1];
         debug4 << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitUCDAMRMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitUCDAMRMetadata::read called" << endl;

      const bool rvalue = UCDAMRMetadata::read(vlsvReader,attribs);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Mesh has " << MeshMetadata::N_domains << " domains" << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitUCDAMRMetadata::readDomains(vlsv::Reader* vlsvReader) {
      debug2 << "VLSV\t\t UCDAMRMetadata::readDomains called" << endl;

      const bool rvalue = UCDAMRMetadata::readDomains(vlsvReader);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Mesh size in blocks: " << meshBoundingBox[0] << ' ' << meshBoundingBox[1] << ' ' << meshBoundingBox[2];
         debug4 << " block sizes: " << meshBoundingBox[3] << ' ' << meshBoundingBox[4] << ' ' << meshBoundingBox[5] << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitUCDAMRMetadata::readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return VisitMeshMetadata::readVariables(vlsv,attribs);
   }



   VisitUCDGenericMultiMeshMetadata::VisitUCDGenericMultiMeshMetadata(): UCDGenericMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitUCDGenericMultiMeshMetadata::~VisitUCDGenericMultiMeshMetadata() { }

   avtMeshType VisitUCDGenericMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDGenericMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}
   
   const std::string& VisitUCDGenericMultiMeshMetadata::getCorrectVlsvMeshType() const {
      return UCDGenericMultiMeshMetadata::getCorrectVlsvMeshType();
   }

   const uint64_t VisitUCDGenericMultiMeshMetadata::getCellConnectivitySize(int domain) const {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::getCellConnectivitySize called" << endl;
      return UCDGenericMultiMeshMetadata::getCellConnectivitySize(domain);
   }

   bool VisitUCDGenericMultiMeshMetadata::getDomainInfo(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
					                                    const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t UCDMultiMeshMetadata::getDomainInfo called, domain: " << domain << endl;
      
      const bool rvalue = UCDGenericMultiMeshMetadata::getDomainInfoNodes(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets);

      if (rvalue == false) {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      } else {
         debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
         debug4 << MeshMetadata::zoneDomainOffsets[domain] << ' ' << MeshMetadata::zoneDomainOffsets[domain+1];
         debug4 << " ghost offsets: " << MeshMetadata::zoneGhostOffsets[domain] << ' ' << MeshMetadata::zoneGhostOffsets[domain+1];
         debug4 << " variable offsets: " << MeshMetadata::zoneVariableOffsets[domain] << ' ' << MeshMetadata::zoneVariableOffsets[domain+1];
         debug4 << endl;
      }

      return rvalue;
   }

   bool VisitUCDGenericMultiMeshMetadata::getDomainInfoNodes(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			                                                 const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      debug2 << "VLSV\t\t UCDGenericMultiMeshMetadata::getDomainInfoNodes called, domain: " << domain << endl;

      const bool rvalue = UCDGenericMultiMeshMetadata::getDomainInfoNodes(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Domain: " << domain << " node offsets: ";
         debug4 << MeshMetadata::nodeDomainOffsets[domain] << ' ' << MeshMetadata::nodeDomainOffsets[domain+1];
         debug4 << " ghost offsets: " << MeshMetadata::nodeGhostOffsets[domain] << ' ' << MeshMetadata::nodeGhostOffsets[domain+1];
         debug4 << " variable offsets: " << MeshMetadata::nodeVariableOffsets[domain] << ' ' << MeshMetadata::nodeVariableOffsets[domain+1];
         debug4 << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      } 

      return rvalue;
   }

   bool VisitUCDGenericMultiMeshMetadata::getDomainInfoZones(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			                                                 const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets) {
      const bool rvalue = UCDGenericMultiMeshMetadata::getDomainInfoZones(vlsvReader,domain,domainOffsets,ghostOffsets,variableOffsets);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Domain: " << domain << " cell offsets: ";
         debug4 << MeshMetadata::zoneDomainOffsets[domain] << ' ' << MeshMetadata::zoneDomainOffsets[domain+1];
         debug4 << " ghost offsets: " << MeshMetadata::zoneGhostOffsets[domain] << ' ' << MeshMetadata::zoneGhostOffsets[domain+1];
         debug4 << " variable offsets: " << MeshMetadata::zoneVariableOffsets[domain] << ' ' << MeshMetadata::zoneVariableOffsets[domain+1];
         debug4 << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      } 

      return rvalue;
   }

   int VisitUCDGenericMultiMeshMetadata::getNodeDataSize() const {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::getNodeDataSize called" << endl;
      return UCDGenericMultiMeshMetadata::getNodeDataSize();
   }

   vlsv::datatype::type VisitUCDGenericMultiMeshMetadata::getNodeDatatype() const {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::getNodeDatatype called" << endl;
      return UCDGenericMultiMeshMetadata::getNodeDatatype();
   }

   bool VisitUCDGenericMultiMeshMetadata::read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::read called" << endl;
      const bool rvalue = UCDGenericMultiMeshMetadata::read(vlsv,attribs);

      if (rvalue == false) {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      } else {
         debug4 << "VLSV\t\t Mesh has " << MeshMetadata::N_domains << " domains" << endl;
         debug4 << "VLSV\t\t Mesh has " << MeshMetadata::N_totalZones << " zones" << endl;
         debug4 << "VLSV\t\t Node coordinate datatype: " << UCDGenericMultiMeshMetadata::nodeDatatype << endl;
         debug4 << "VLSV\t\t Node coordinate datasize: " << UCDGenericMultiMeshMetadata::nodeDataSize << endl;
      }

      return rvalue;
   }

   bool VisitUCDGenericMultiMeshMetadata::readDomains(vlsv::Reader* vlsvReader) {
      debug2 << "VLSV\t\t VisitUCDGenericMultiMeshMetadata::readDomains called" << endl;

      const bool rvalue = UCDGenericMultiMeshMetadata::readDomains(vlsvReader);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Mesh size in blocks: " 
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::X_BLOCKS] << ' ' 
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::Y_BLOCKS] << ' ' 
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::Z_BLOCKS];
         debug4 << " block sizes: " 
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_X] << ' '
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Y] << ' ' 
            << MeshMetadata::meshBoundingBox[vlsv::ucdgenericmulti::bbox::BLOCK_WIDTH_Z] << endl;
         debug4 << "VLSV\t\t Mesh has " << MeshMetadata::N_localZones << " local zones in total" << endl;
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   bool VisitUCDGenericMultiMeshMetadata::readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return VisitMeshMetadata::readVariables(vlsv,attribs);
   }
}