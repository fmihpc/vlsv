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

#include <cstdlib>
#include <iostream>

#include "multi_io_unit.h"

using namespace std;

namespace vlsv {

   /** Constructor for struct Multi_IO_Unit.
    * NOTE: MPI datatypes passed to Multi_IO_Unit are not committed or freed, 
    * thus one should only use native MPI datatypes, e.g. as returned by 
    * getMPIDatatype function.
    * @param array In vlsv::Writer this is a pointer to array whose contents are
    * written to file. In vlsv::ParallelReader this is a pointer to array where data 
    * from file is to be read.
    * @param mpiType MPI datatype defining the I/O operation.
    * @param amount Amount of data to be written or read, in units of mpiType.*/
   Multi_IO_Unit::Multi_IO_Unit(char* array,const MPI_Datatype& mpiType,const uint64_t& amount): array(array),mpiType(mpiType),amount(amount) { }
   
} // namespace vlsv
