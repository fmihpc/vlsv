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

#ifndef MESH_METADATA_H
#define MESH_METADATA_H

#include <map>
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
      VariableMetadata(VariableCentering centering,const std::string& name,int vectorSize);
      
      VariableCentering centering;
      std::string name;
      int vectorSize;
   };
   
   /** Virtual base class for storing mesh metadata.*/
   class MeshMetadata {
    public:
      MeshMetadata();      
      virtual ~MeshMetadata();

      virtual uint64_t getArraySize() const;
      virtual uint64_t getDataSize() const;
      virtual vlsv::datatype::type getDatatype() const;
      virtual uint64_t getMaximumRefinementLevel() const;
      virtual void getMeshPeriodicity(bool& xPeriodic,bool& yPeriodic,bool& zPeriodic) const;
      virtual std::string getName() const;
      virtual uint64_t getNumberOfGhostNodes() const;
      virtual uint64_t getNumberOfGhostNodes(int domain) const =0;
      virtual uint64_t getNumberOfGhostZones() const;
      virtual uint64_t getNumberOfGhostZones(int domain) const =0;
      virtual uint64_t getNumberOfLocalNodes() const;
      virtual uint64_t getNumberOfLocalNodes(int domain) const =0;
      virtual uint64_t getNumberOfLocalZones() const;
      virtual uint64_t getNumberOfLocalZones(int domain) const =0;
      virtual uint64_t getNumberOfTotalNodes() const;
      virtual uint64_t getNumberOfTotalNodes(int domain) const =0;
      virtual uint64_t getNumberOfTotalZones() const;
      virtual uint64_t getNumberOfTotalZones(int domain) const =0;
      virtual const double* getTransform() const;
      virtual const std::vector<VariableMetadata>& getVariables() const;
      virtual uint64_t getVectorSize() const;
      virtual std::string getXLabel() const;
      virtual std::string getYLabel() const;
      virtual std::string getZLabel() const;
      virtual std::string getXUnits() const;
      virtual std::string getYUnits() const;
      virtual std::string getZUnits() const;
      virtual bool hasTransform() const;
      
      virtual bool read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);
      
    protected:
      uint64_t arraySize;
      uint64_t vectorSize;
      uint64_t dataSize;
      vlsv::datatype::type datatype;
      std::vector<VariableMetadata> variableMetadata;

      int32_t maxRefinementLevel; /**< Maximum refinement level in mesh, equals zero if
                                   * the mesh type does not support refinement, or if the 
                                   * mesh is not refined.*/
      
      uint64_t N_ghostNodes;    /**< Total number of ghost nodes in the mesh, summed over all domains.*/
      uint64_t N_ghostZones;    /**< Total number of ghost zones in the mesh, summed over all domains.*/
      uint64_t N_localNodes;    /**< Total number of local nodes in the mesh, summer over all domains.*/
      uint64_t N_localZones;    /**< Total number of local zones in the mesh, summed over all domains.*/
      uint64_t N_totalNodes;    /**< Total number of zones (local+ghost) in the mesh, summer over all domains.*/
      uint64_t N_totalZones;    /**< Total number of zones (local+ghost) in the mesh, summed over all domains.*/
      
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
   };
}

#endif
