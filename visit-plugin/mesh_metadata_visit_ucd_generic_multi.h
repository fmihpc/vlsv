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

#ifndef MESH_METADATA_VISIT_UCD_GENERIC_MULTI_H
#define MESH_METADATA_VISIT_UCD_GENERIC_MULTI_H

#include <stdint.h>

#include <mesh_metadata_visit.h>

namespace vlsvplugin {
   
   class VisitUCDGenericMultiMeshMetadata: public VisitMeshMetadata {
    public:
      VisitUCDGenericMultiMeshMetadata();
      virtual ~VisitUCDGenericMultiMeshMetadata();
      
      uint64_t getBlockSize() const;
      bool getDomainInfoNodes(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			      const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets);
      bool getDomainInfoZones(vlsv::Reader* vlsvReader,int domain,const uint64_t*& domainOffsets,
			      const uint64_t*& ghostOffsets,const uint64_t*& variableOffsets);
      const uint64_t getCellConnectivitySize(int domain) const;
      const uint64_t* getDomainOffsets();
      const uint64_t* getGhostOffsets();
      const uint64_t* getMeshBoundingBox();
      int getNodeDataSize() const;
      vlsv::datatype::type getNodeDatatype() const;
      const uint64_t getNodeOffset(int domain) const;
      const uint64_t getNumberOfNodes(int domain) const;
      const vlsv::geometry::type& getMeshGeometry() const;
      uint64_t getNumberOfGhostNodes(int domain) const;
      uint64_t getNumberOfGhostZones(int domain) const;
      uint64_t getNumberOfLocalNodes(int domain) const;
      uint64_t getNumberOfLocalZones(int domain) const;
      uint64_t getNumberOfTotalNodes(int domain) const;
      uint64_t getNumberOfTotalZones(int domain) const;
      const uint64_t* getVariableOffsets();
      const uint64_t getZoneOffset(int domain) const;
      
      bool read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs);
      
    protected:

      bool domainMetadataRead;        /**< If true, domain metadata has been read.*/
      bool meshMetadataRead;          /**< If true, mesh metadata has been read.*/

      vlsv::datatype::type nodeDatatype;
      int nodeDataSize;
      
      vlsv::geometry::type geometry;  /**< Mesh geometry (Cartesian, Cylindrical, etc.).*/
      uint64_t blockSize;             /**< Number of cells in one block/patch. This will 
				       * probably be deprecated in the future.*/
      uint64_t* ghostNodeOffsets;
      uint64_t* ghostZoneOffsets;
      uint64_t* meshBoundingBox;      /**< Mesh bounding box information. Contains elements 
				       * defined in enumeration vlsv::ucdgenericmulti::bbox::elements.*/
      uint64_t* nodeOffsets;          /**< Offsets into array containing node coordinates. Tells where 
				       * data for each domain begins. Size of array is number of domains+1.*/
      uint64_t* nodeVariableOffsets;  /**< Offsets into arrays containing node-centered variables.
				       * Size of array is number of domains+1.*/
      uint64_t* zoneOffsets;          /**< Offsets into array containing zone connectivity entries. Tells
				       * where data for each domain begins. Size of array is number of domains+1.*/
      uint64_t* zoneVariableOffsets;  /**< Offsets into arrays containing zone-centered variables.
				       * Size of array is number of domains+1.*/

      bool readDomains(vlsv::Reader* vlsvReader);
   };
   
} // namespace vlsvplugin

#endif
