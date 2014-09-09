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

#ifndef VLSV_AMR_H
#define VLSV_AMR_H

#include <stdint.h>

namespace vlsv {

   /** This header file contains definitions for functions that are used to 
    * calculate cell indices in meshes that support adaptive refinement.
    * Function initMesh must be called prior to calling any other functions.*/
   
   void calculateCellIndices(const uint64_t& globalID,uint32_t& refLevel,uint32_t& i,uint32_t& j,uint32_t& k);
   uint64_t calculateGlobalID(const uint32_t& refLevel,const uint32_t& i,const uint32_t& j,const uint32_t& k);
   bool initMesh(const uint64_t& Nx0,const uint64_t& Ny0,const uint64_t Nz0,const uint64_t& maxRefinementLevel);

} // namespace vlsv

#endif
