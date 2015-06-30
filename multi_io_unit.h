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

#ifndef MULTI_IO_UNIT_H
#define MULTI_IO_UNIT_H

#include <stdint.h>
#include <mpi.h>

namespace vlsv {
   /** Definition of a parallel file I/O unit. Processes can 
    * define zero or more file action units for a collective file I/O
    * operation in VLSVWriter and VLSVParReader. Those classes will then 
    * arrange all defined I/O operations to be performed with a single MPI call
    * to optimize performance.
    */
   struct Multi_IO_Unit {
    public:
      Multi_IO_Unit(char* array,const MPI_Datatype& mpiType,const uint64_t& amount);
      
      char* array;              /**< Pointer to data to be written.*/
      MPI_Datatype mpiType;     /**< MPI datatype of data that is written.*/
      uint64_t amount;          /**< How many elements of type mpiType are to be written.*/

    private:
      
      /** Private default constructor to prevent creation of empty multiwrite units.*/
      Multi_IO_Unit();
   };
}

#endif