/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2014 Finnish Meteorological Institute
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

#ifndef AMR_MESH_H
#define AMR_MESH_H

#include <vector>
#include <stdint.h>
#include <unordered_map>

class AmrMesh {
 public:
   AmrMesh(const uint32_t& Nx0,const uint32_t& Ny0,const uint32_t& Nz0,
	   const uint32_t& xCells,const uint32_t& yCells,const uint32_t& zCells,
	   const uint8_t& maxRefLevel);
   ~AmrMesh();

   const std::unordered_map<uint64_t,uint8_t>::iterator begin();
   const std::unordered_map<uint64_t,uint8_t>::iterator end();
   
   bool initialize(const double& xmin,const double& xmax,const double& ymin,const double& ymax,
		   const double& zmin,const double& zmax,const uint8_t refLevel=0);
   bool finalize();
   
   void getChildren(const uint64_t& globalID,std::vector<uint64_t>& children);
   uint64_t getGlobalID(const uint32_t& refLevel,const uint32_t& i,const uint32_t& j,const uint32_t& k);
   void getNeighbors(const uint64_t& globalID,std::vector<uint64_t>& neighborIDs);
   void getIndices(const uint64_t& globalID,uint32_t& refLevel,uint32_t& i,uint32_t& j,uint32_t& k);
   uint64_t getParent(const uint64_t& globalID);
   void getSiblingNeighbors(const uint64_t& globalID,std::vector<uint64_t>& nbrs);
   void getSiblings(const uint64_t& globalID,std::vector<uint64_t>& siblings);
   size_t size() const;

   bool coarsen(const uint64_t& globalID);
   bool refine(const uint64_t& globalID);

   bool write(const std::string fileName);

 private:
   AmrMesh();

   uint32_t bbox[6];
   std::unordered_map<uint64_t,uint8_t> globalIDs;
   bool initialized;
   std::vector<uint64_t> offsets;
   uint64_t N_blocks0;
   uint8_t refLevelMaxAllowed;
   
   double meshLimits[6];
};

#endif
