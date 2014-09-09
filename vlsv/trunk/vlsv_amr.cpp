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

#include <algorithm>
#include <cmath>
#include <vector>
#include "vlsv_amr.h"

using namespace std;

static vector<uint64_t> offsets;
static uint64_t maxRefinementLevel = 0;
static uint64_t Nx0 = 0;
static uint64_t Ny0 = 0;
static uint64_t Nz0 = 0;
static uint64_t N_cells0 = 0;

namespace vlsv {

   /** Calculate i,j,k indices, and the refinement level, of a block from the given global ID.
    * @param globalID Global ID of the block.
    * @param refLevel Refinement level of the block, zero value indicates no refinement.
    * @param i Block i-index is written here.
    * @param j Block j-index is written here.
    * @param k Block k-index is written here.*/
   void calculateCellIndices(const uint64_t& globalID,uint32_t& refLevel,uint32_t& i,uint32_t& j,uint32_t& k) {
      refLevel   = upper_bound(offsets.begin(),offsets.end(),globalID)-offsets.begin()-1;
      const uint64_t cellOffset = offsets[refLevel];

      const uint64_t multiplier = pow(2,refLevel);
      const uint64_t Nx = Nx0 * multiplier;
      const uint64_t Ny = Ny0 * multiplier;
      //const uint64_t Nz = Nz0 * multiplier;
      
      uint64_t index = globalID - cellOffset;
      k = index / (Ny*Nx);
      index -= k*Ny*Nx;
      j = index / Nx;
      i = index - j*Nx;
   }

   /** Calculate the global ID of a block from the given i,j,k indices.
    * @param refLevel Refinement level of the block.
    * @param i Block i-index.
    * @param j Block j-index.
    * @param k Block k-index.
    * @return The unique global ID of the block.*/
   uint64_t calculateGlobalID(const uint32_t& refLevel,const uint32_t& i,const uint32_t& j,const uint32_t& k) {
      const uint32_t multiplier = pow(2,refLevel);
      return offsets[refLevel] + k*Ny0*Nx0*multiplier*multiplier + j*Nx0*multiplier + i;
   }

   /** Initialize internal variables so that functions calculating block global IDs etc. work correctly.
    * @param Nx0 Number of blocks in x-direction in unrefined mesh.
    * @param Ny0 Number of blocks in y-direction in unrefined mesh.
    * @param Nz0 Number of blocks in z-direction in unrefined mesh.
    * @param maxRefinementLevel Maximum mesh refinement level, zero value indicates that 
    * adaptive mesh refinement is not used.
    * @return If true, initialization was successful.*/
   bool initMesh(const uint64_t& Nx0,const uint64_t& Ny0,const uint64_t Nz0,const uint64_t& maxRefinementLevel) {
      // Copy values to global variables:
      ::maxRefinementLevel = maxRefinementLevel;
      ::Nx0 = Nx0;
      ::Ny0 = Ny0;
      ::Nz0 = Nz0;
      ::N_cells0 = Nx0*Ny0*Nz0;
      
      // Clear offsets vector:
	{
	   vector<uint64_t> dummy;
	   offsets.swap(dummy);
	}
      offsets.resize(maxRefinementLevel+1);
      
      offsets[0] = 0;
      for (size_t i=1; i<maxRefinementLevel+1; ++i) {
	 offsets[i] = offsets[i-1] + N_cells0 * pow(8,i-1);
      }
      
      return true;
   }
   
} // namespace vlsv

