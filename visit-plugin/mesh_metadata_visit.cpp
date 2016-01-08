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

#include <DebugStream.h>

#include "mesh_metadata_visit.h"

using namespace std;

namespace vlsvplugin {
   VisitMeshMetadata::VisitMeshMetadata(): MeshMetadata() {
      domainMetadataRead = false;
      blockOrigin = 0;
      blockSize = 1;
      N_domains = 1;
   }
   
   VisitMeshMetadata::~VisitMeshMetadata() { }
   
   uint64_t VisitMeshMetadata::getBlockSize() const {return blockSize;}
   
   uint64_t VisitMeshMetadata::getNodeDomainOffset(uint64_t domain) const {return nodeDomainOffsets[domain];}

   int VisitMeshMetadata::getBlockOrigin() const {return blockOrigin;}
   
   const std::vector<uint64_t>& VisitMeshMetadata::getMeshBoundingBox() const {return meshBoundingBox;}

   avtMeshType VisitMeshMetadata::getMeshType() const {return meshType;}
   
   std::string VisitMeshMetadata::getMeshTypeString() const {return meshTypeString;}

   uint64_t VisitMeshMetadata::getNumberOfDomains() const {return N_domains;}
   
   uint64_t VisitMeshMetadata::getNumberOfGhostNodes(uint64_t domain) const {
      return nodeGhostOffsets[domain+1]-nodeGhostOffsets[domain];
   }

   uint64_t VisitMeshMetadata::getNumberOfGhostZones(uint64_t domain) const {
      return blockSize*(zoneGhostOffsets[domain+1]-zoneGhostOffsets[domain]);
   }

   uint64_t VisitMeshMetadata::getNumberOfLocalNodes(uint64_t domain) const {
      return getNumberOfTotalNodes(domain)-getNumberOfGhostNodes(domain);
   }

   uint64_t VisitMeshMetadata::getNumberOfLocalZones(uint64_t domain) const {
      return getNumberOfTotalZones(domain)-getNumberOfGhostZones(domain);
   }

   uint64_t VisitMeshMetadata::getNumberOfTotalNodes(uint64_t domain) const {
      return nodeVariableOffsets[domain+1]-nodeVariableOffsets[domain];
   }

   uint64_t VisitMeshMetadata::getNumberOfTotalZones() const {return N_totalZones;}

   uint64_t VisitMeshMetadata::getNumberOfTotalZones(uint64_t domain) const {
      return blockSize*(zoneDomainOffsets[domain+1]-zoneDomainOffsets[domain]);
   }

   int VisitMeshMetadata::getSpatialDimension() const {return spatialDimension;}
   
   int VisitMeshMetadata::getTopologicalDimension() const {return topologicalDimension;}

   uint64_t VisitMeshMetadata::getZoneDomainOffset(uint64_t domain) const {return zoneDomainOffsets[domain];}

   bool VisitMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      return MeshMetadata::read(vlsvReader,attribs);
   }

   bool VisitMeshMetadata::readDomainMetadata(vlsv::Reader* vlsvReader) {
      // Exit if domain metadata has already been read:
      if (domainMetadataRead == true) return true;
      domainMetadataRead = false;
      debug2 << "VLSV\t\t VisitMeshMetadata::readDomains called" << endl;

      // Check that vlsv::Reader exists:
      if (vlsvReader == NULL) {
	     debug3 << "VLSV\t\t ERROR: vlsv::Reader is NULL" << endl;
	     return false;
      }

      // Read mesh bounding box:
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("mesh",name));
      if (vlsvMeshType == vlsv::mesh::UCD_MULTI) {
         meshBoundingBox.resize(6);   
         uint64_t* ptr = meshBoundingBox.data();
         if (vlsvReader->read("MESH_BBOX",attribs,0,6,ptr,false) == false) {
		    debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_BBOX'" << endl;
		    return false;
         }
         blockSize = meshBoundingBox[3]*meshBoundingBox[4]*meshBoundingBox[5];
         debug4 << "VLSV\t\t Mesh size in blocks: " << meshBoundingBox[0] << ' ' << meshBoundingBox[1] << ' ' << meshBoundingBox[2];
         debug4 << " block sizes: " << meshBoundingBox[3] << ' ' << meshBoundingBox[4] << ' ' << meshBoundingBox[5] << endl;
      } else if (vlsvMeshType == vlsv::mesh::QUAD_MULTI) {
         meshCoordinates.resize(6);
         float* ptr = meshCoordinates.data();
         if (vlsvReader->read("MESH_BBOX",attribs,0,6,ptr,false) == false) {
	        debug3 << "VLSV\t\t ERROR: Failed to read array 'MESH_BBOX'" << endl;
	        return false;
         }
         debug4 << "VLSV\t\t Mesh corner coordinates: " << meshCoordinates[0] << ' ' << meshCoordinates[1] << ' ' << meshCoordinates[2];
         debug4 << " cell sizes: " << meshCoordinates[3] << ' ' << meshCoordinates[4] << ' ' << meshCoordinates[5] << endl;
      }

      // Read domain sizes:
      uint64_t* domainInfo = NULL;
      string domainArrayName;
      if (vlsvMeshType == vlsv::mesh::QUAD_MULTI) domainArrayName = "MESH_ZONES";
      else if (vlsvMeshType == vlsv::mesh::UCD_MULTI) domainArrayName = "MESH_DOMAIN_SIZES";
      else {
         debug3 << "VLSV\t\t ERROR: Unsupported mesh type '" << vlsvMeshType << "' in " << __FILE__ << ":" << __LINE__ << endl;
         return false;
      }
      if (vlsvReader->read(domainArrayName,attribs,0,VisitMeshMetadata::N_domains,domainInfo) == false) {
		debug3 << "VLSV\t\t ERROR: Failed to read array '" << domainArrayName << "' in " << __FILE__ << ":" << __LINE__ << endl;
		debug3 << "VLSV\t\t Message from VLSV Reader '" << vlsvReader->getLastError() << "'" << endl;
		delete [] domainInfo; domainInfo = NULL;
		return false;
      }

      zoneDomainOffsets.resize(VisitMeshMetadata::N_domains+1); zoneDomainOffsets.shrink_to_fit();
      zoneGhostOffsets.resize(VisitMeshMetadata::N_domains+1); zoneGhostOffsets.shrink_to_fit();
      zoneVariableOffsets.resize(VisitMeshMetadata::N_domains+1); zoneVariableOffsets.shrink_to_fit();

      zoneDomainOffsets[0]   = 0;
      zoneGhostOffsets[0]    = 0;
      zoneVariableOffsets[0] = 0;
      for (auto i=0; i<VisitMeshMetadata::N_domains; ++i) {
         zoneDomainOffsets[i+1]   = zoneDomainOffsets[i] + domainInfo[i*2];
         zoneGhostOffsets[i+1]    = zoneGhostOffsets[i] + domainInfo[i*2+1];
         zoneVariableOffsets[i+1] = zoneVariableOffsets[i] + domainInfo[i*2+0]-domainInfo[i*2+1];
      }
      delete [] domainInfo; domainInfo = NULL;

      // Compute total number of real and ghost cells:
      MeshMetadata::N_ghostZones = zoneGhostOffsets[VisitMeshMetadata::N_domains];
      MeshMetadata::N_localZones  = zoneVariableOffsets[VisitMeshMetadata::N_domains];

      domainMetadataRead = true;
      return true;
   }
   
} // namespace vlsvplugin


