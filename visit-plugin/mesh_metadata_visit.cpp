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

#include <DebugStream.h>

#include "mesh_metadata_visit.h"

using namespace std;

namespace vlsvplugin {
   VisitMeshMetadata::VisitMeshMetadata(): MeshMetadata() {
      blockOrigin = 0;
   }
   
   VisitMeshMetadata::~VisitMeshMetadata() { }

   int VisitMeshMetadata::getBlockOrigin() const {return blockOrigin;}

   /** Get VisIt / VTK mesh type that corresponds to the VLSV mesh.
    * Type of VisIt / VTK mesh stored in the VLSV mesh.*/
   avtMeshType VisitMeshMetadata::getAvtMeshType() const {return meshType;}
   
   /** Get string representation of VisIt / VTK mesh type that corresponds to the VLSV mesh.
    * @return String representation of VisIt / VTK mesh.*/
   std::string VisitMeshMetadata::getAvtMeshTypeString() const {return meshTypeString;}

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
      for (uint64_t i=0; i<VisitMeshMetadata::N_domains; ++i) {
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


