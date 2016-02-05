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

//#pragma once

#ifndef MESH_METADATA_H
#define MESH_METADATA_H

#include <stdint.h>
#include <string>
#include <vector>
#include <vlsv_reader.h>

namespace vlsvplugin {
   
   enum VariableCentering {
      NODE_CENTERED,           /**< Indicates node variable.*/
      ZONE_CENTERED            /**< Indicates cell variable.*/
   };
   
   struct VariableMetadata {
      VariableMetadata(VariableCentering centering,const std::string& name,uint64_t vectorSize);
      
      VariableCentering centering;     /**< Variable centering, either NODE_CENTERED or ZONE_CENTERED.*/
      std::string name;                /**< Name of the variable.*/
      uint64_t vectorSize;             /**< Size of the data vector, i.e., how many values per cell.*/
   };
   
   /** Virtual base class for storing mesh metadata.*/
   class MeshMetadata {
    public:
      MeshMetadata();      
      virtual ~MeshMetadata();

      const std::string& getErrorString() const;

      virtual void getBlockWidths(uint64_t& blockWidthX,uint64_t& blockWidthY,uint64_t& blockWidthZ) const;
      virtual uint64_t getMaximumRefinementLevel() const;
      virtual const vlsv::geometry::type& getMeshGeometry() const;
      virtual const std::vector<uint64_t>& getMeshBoundingBox() const;
      virtual void getMeshPeriodicity(bool& xPeriodic,bool& yPeriodic,bool& zPeriodic) const;
      virtual std::string getName() const;

      virtual uint64_t getNodeDomainOffset(uint64_t domain) const;
      virtual uint64_t getZoneDomainOffset(uint64_t domain) const;

      virtual uint64_t getNumberOfDomains() const;
      virtual uint64_t getNumberOfGhostNodes(uint64_t domain) const;
      virtual uint64_t getNumberOfLocalNodes(uint64_t domain) const;
      virtual uint64_t getNumberOfTotalNodes(uint64_t domain) const;

      virtual uint64_t getNumberOfGhostZones(uint64_t domain) const;
      virtual uint64_t getNumberOfLocalZones(uint64_t domain) const;      
      virtual uint64_t getNumberOfTotalZones() const;
      virtual uint64_t getNumberOfTotalZones(uint64_t domain) const;

      virtual int getSpatialDimension() const;
      virtual int getTopologicalDimension() const;

      virtual const double* getTransform() const;
      virtual const std::vector<VariableMetadata>& getVariables() const;

      virtual void getAxisLabels(std::string& xLabel,std::string& yLabel,std::string& zLabel) const;
      virtual void getAxisUnits(std::string& xUnits,std::string& yUnits,std::string& zUnits) const;
      virtual bool hasTransform() const;
      
      virtual bool read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);

    protected:

      bool exitWithError(const std::string& s) const;
      bool exitWithError(const std::stringstream& ss) const;

      virtual const std::string& getCorrectVlsvMeshType() const = 0;
      virtual bool checkVlsvMeshType(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);
      virtual bool readDomainMetadata(vlsv::Reader* vlsvReader);
      virtual bool readMeshGeometry(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);
      virtual bool readVariables(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);

      mutable std::string errorString;        /**< Description of the most recent error that has occurred, if any.*/

      bool domainMetadataRead;                /**< If true, domain metadata has been read.*/
      bool meshMetadataRead;                  /**< If true, mesh metadata has been read.*/
      std::vector<VariableMetadata> variableMetadata;

      int spatialDimension;       /**< Spatial dimension of mesh.*/
      int topologicalDimension;   /**< Topological dimension of mesh.*/

      int32_t maxRefinementLevel; /**< Maximum refinement level in mesh, equals zero if
                                   * the mesh type does not support refinement, or if the 
                                   * mesh is not refined.*/
      
      uint64_t blockSize;
      uint64_t N_domains;       /**< Number of domains in the mesh.*/
      uint64_t N_ghostNodes;    /**< Total number of ghost nodes in the mesh, summed over all domains.*/
      uint64_t N_ghostZones;    /**< Total number of ghost zones in the mesh, summed over all domains.*/
      uint64_t N_localNodes;    /**< Total number of local nodes in the mesh, summer over all domains.*/
      uint64_t N_localZones;    /**< Total number of local zones in the mesh, summed over all domains.*/
      uint64_t N_totalNodes;    /**< Total number of zones (local+ghost) in the mesh, summer over all domains.*/
      uint64_t N_totalZones;    /**< Total number of zones (local+ghost) in the mesh, summed over all domains.*/
      
      vlsv::geometry::type geometry;  /**< Mesh geometry (Cartesian, Cylindrical, etc.).*/
      vlsv::mesh::type vlsvMeshType;
      std::string name;         /**< Name of the mesh.*/
      std::string xLabel;       /**< x-coordinate axis label.*/
      std::string yLabel;       /**< y-coordinate axis label.*/
      std::string zLabel;       /**< z-coordinate axis label.*/
      bool xPeriodic;           /**< If true, mesh is periodic in x.*/
      bool yPeriodic;           /**< If true, mesh is periodic in y.*/
      bool zPeriodic;           /**< If true, mesh is periodic in z.*/
      std::string xUnits;       /**< Unit for x-coordinate.*/
      std::string yUnits;       /**< Unit for y-coordinate.*/
      std::string zUnits;       /**< Unit for z-coordinate.*/

      std::string transformName; /**< Name of the transformation matrix that needs to be applied 
                                  * to the node coordinates. Zero length string indicates no matrix.*/
      double transform[16];      /**< Components of the transform matrix, defaults to identity matrix.*/

      std::vector<uint64_t> meshBoundingBox;
      std::vector<float> meshCoordinates;

      std::vector<uint64_t> nodeDomainOffsets;
      std::vector<uint64_t> nodeGhostOffsets;
      std::vector<uint64_t> nodeVariableOffsets;

      std::vector<uint64_t> zoneConnectivityOffsets; /**< For each domain, an offset into zone connectivity array
                                                      * that tells where to start to read data from.*/
      std::vector<uint64_t> zoneDomainOffsets;
      std::vector<uint64_t> zoneGhostOffsets;
      std::vector<uint64_t> zoneVariableOffsets;
   };
}

#endif
