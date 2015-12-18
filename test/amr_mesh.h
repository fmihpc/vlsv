/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2015 Finnish Meteorological Institute
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

#pragma once

#ifndef AMR_MESH_H
#define AMR_MESH_H

#include <vector>
#include <stdint.h>
#include <unordered_map>
#include <limits>

#warning TODO Add 2D mode
#warning TOTO Add periodic boundaries

namespace amr {
   
   typedef uint64_t GlobalID;
   typedef uint32_t LocalID;

   const GlobalID INVALID_GLOBALID = std::numeric_limits<GlobalID>::max();
   const LocalID INVALID_LOCALID = std::numeric_limits<LocalID>::max();

   typedef int (*CallbackCoarsenBlock)(const GlobalID siblingIDs[8],const LocalID siblingIndices[8],
				       const GlobalID globalID,LocalID& index);
   typedef int (*CallbackCreateBlock)(const GlobalID& globalID,LocalID& index);
   typedef int (*CallbackDeleteBlock)(const GlobalID& globalID,const LocalID& index);
   typedef int (*CallbackRefineBlock)(const GlobalID& globalID,const LocalID& index,
				      const GlobalID childrenIDs[8],LocalID childrenIndices[8]); 
   
   class AmrMesh {
    public:
      AmrMesh(const uint32_t& Nx0,const uint32_t& Ny0,const uint32_t& Nz0,
	      const uint32_t& xCells,const uint32_t& yCells,const uint32_t& zCells,
	      const uint8_t& maxRefLevel);
      ~AmrMesh();
      
      bool initialize(const double& xmin,const double& xmax,const double& ymin,const double& ymax,
		      const double& zmin,const double& zmax,const uint8_t refLevel=0);
      bool finalize();
      bool registerCallbacks(CallbackCoarsenBlock coarsenBlock,CallbackCreateBlock createBlock,
			     CallbackDeleteBlock deleteBlock,CallbackRefineBlock refineBlock);

      std::unordered_map<GlobalID,LocalID>::iterator begin();
      std::unordered_map<GlobalID,LocalID>::iterator end();

      LocalID  get(const GlobalID& globalID) const;
      bool     getBlockCoordinates(const GlobalID& globalID,double coords[3]) const;
      bool     getBlockSize(const GlobalID& globalID,double size[3]) const;
      bool     getCellSize(const GlobalID& globalID,double size[3]) const;
      void     getChildren(const GlobalID& globalID,std::vector<GlobalID>& children);
      GlobalID getGlobalID(const double& x,const double& y,const double& z);
      GlobalID getGlobalID(const uint32_t& refLevel,const uint32_t& i,const uint32_t& j,const uint32_t& k);
      void     getNeighbors(const GlobalID& globalID,std::vector<GlobalID>& neighborIDs);
      void     getIndices(const GlobalID& globalID,uint32_t& refLevel,uint32_t& i,uint32_t& j,uint32_t& k) const;
      GlobalID getParent(const GlobalID& globalID);
      void     getSiblingNeighbors(const GlobalID& globalID,std::vector<GlobalID>& nbrs);
      void     getSiblings(const GlobalID& globalID,GlobalID siblings[8]);
      void     getSiblings(const GlobalID& globalID,std::vector<GlobalID>& siblings);
      bool     set(const GlobalID& globalID,const LocalID& localID);
      size_t   size() const;

      bool checkBlock(const GlobalID& globalID);
      bool checkMesh();
      bool coarsen(const GlobalID& globalID);
      bool refine(const GlobalID& globalID);

      bool write(const std::string fileName);

    private:
      AmrMesh();
      
      uint32_t bbox[6];                                /**< Mesh bounding box. First three elements contain the 
							* number of blocks in x,y,z directions at the base grid level.
							* Last three indices contain the number of cells in each 
							* block in x,y,z directions.*/
      double meshLimits[6];                            /**< Physical coordinates of the mesh boundaries. Elements contain
							* 0 = Minimum x-coordinate.
							* 1 = Maximum x-coordinate.
							* 2 = Minimum y-coordinate.
							* 3 = Maximum y-coordinate.
							* 4 = Minimum z-coordinate.
							* 5 = Maximum z-coordinate.*/
      std::vector<uint64_t> offsets;                   /**< Block global ID offsets for each refinement level.*/
      uint64_t N_blocks0;                              /**< Number of blocks at base grid level, equal to bbox[0]*bbox[1]*bbox[2].*/
      uint8_t refLevelMaxAllowed;                      /**< Maximum allowed refinement level.*/
      
      std::unordered_map<GlobalID,LocalID> globalIDs;  /**< Hash table containing all existing blocks (global ID,refinement level).*/
      bool initialized;                                /**< If true, mesh was initialized successfully.*/

      CallbackCoarsenBlock callbackCoarsenBlock;
      CallbackCreateBlock callbackCreateBlock;
      CallbackDeleteBlock callbackDeleteBlock;
      CallbackRefineBlock callbackRefineBlock;
   };

} // namespace amr

#endif
