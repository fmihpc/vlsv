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
      virtual void getMeshPeriodicity(bool& xPeriodic,bool& yPeriodic,bool& zPeriodic) const;
      virtual std::string getName() const;
      virtual uint64_t getNumberOfGhostCells() const;
      virtual uint64_t getNumberOfGhostCells(int domain) const =0;
      virtual uint64_t getNumberOfRealCells() const;
      virtual uint64_t getNumberOfRealCells(int domain) const =0;
      virtual uint64_t getNumberOfTotalCells() const;
      virtual uint64_t getNumberOfTotalCells(int domain) const =0;
      virtual const std::vector<VariableMetadata>& getVariables() const;
      virtual uint64_t getVectorSize() const;
      virtual std::string getXLabel() const;
      virtual std::string getYLabel() const;
      virtual std::string getZLabel() const;
      virtual std::string getXUnits() const;
      virtual std::string getYUnits() const;
      virtual std::string getZUnits() const;
      
      virtual bool read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs);
      
    protected:
      uint64_t arraySize;
      uint64_t vectorSize;
      uint64_t dataSize;
      vlsv::datatype::type datatype;
      std::vector<VariableMetadata> variableMetadata;

      uint64_t N_ghostCells;    /**< Total number of ghost cells in the mesh, summed over all domains.*/
      uint64_t N_realCells;     /**< Total number of real cells in the mesh, summed over all domains.*/
      uint64_t N_totalCells;    /**< Total number of cells (real+ghost) in the mesh, summed over all domains.*/
      
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
   };
}

#endif
