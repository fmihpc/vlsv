/* This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013,2015 Finnish Meteorological Institute
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
#include <limits>

#include "mpiconversion.h"
#include "vlsv_common.h"
#include "vlsv_common_mpi.h"

using namespace std;

namespace vlsv {

   /** Get maximum number of bytes that can be read from a file using a single collective MPI routine.
    * This should equal to numeric_limits<MPI_Count>::max() in MPI 3.0+, but in older 
    * library versions the 'count' parameter is simply an integer.
    * @return Maximum number of bytes read using a single MPI collective routine.*/
   size_t getMaxBytesPerRead() {
      return numeric_limits<int>::max();
      //return numeric_limits<MPI_Count>::max();
   }
   
   /** Get maximum number of bytes that can be written to a file using a single collective MPI routine.
    * This should equal to numeric_limits<MPI_Count>::max() in MPI 3.0+, but in older 
    * library versions the 'count' parameter is simply an integer.
    * @return Maximum number of bytes written using a single MPI collective routine.*/
   size_t getMaxBytesPerWrite() {
      return 2000000000;
      return numeric_limits<int>::max();
      //return numeric_limits<MPI_Count>::max();
   }

   MPI_Datatype getMPIDatatype(datatype::type dt,uint64_t dataSize) {
      switch (dt) {
       case datatype::UNKNOWN:
         // TEST
         cerr << "(VLSV) ERROR: VLSV::getMPIDatatype called with datatype::UNKNOWN datatype, returning MPI_DATATYPE_NULL!" << endl;
         // END TEST
         return MPI_DATATYPE_NULL;
         break;
       case datatype::INT:
         switch (dataSize) {
          case (sizeof(int8_t)):
            return MPI_Type<int8_t>();
            break;
          case (sizeof(int16_t)):
            return MPI_Type<int16_t>();
            break;
          case (sizeof(int32_t)):
            return MPI_Type<int32_t>();
            break;
          case (sizeof(int64_t)):
            return MPI_Type<int64_t>();
            break;
          default:
            cerr << "(VLSV) ERROR: VLSV::getMPIDatatype called with datatype::INT datatype and unsupported datasize of " << dataSize << "!" << endl;
            return MPI_DATATYPE_NULL;
            break;
         }
       case datatype::UINT:
         switch (dataSize) {
          case (sizeof(uint8_t)):
            return MPI_Type<uint8_t>();
            break;
          case (sizeof(uint16_t)):
            return MPI_Type<uint16_t>();
            break;
          case (sizeof(uint32_t)):
            return MPI_Type<uint32_t>();
            break;
          case (sizeof(uint64_t)):
            return MPI_Type<uint64_t>();
            break;
          default:
            cerr << "(VLSV) ERROR: VLSV::getMPIDatatype called with datatype::UINT datatype and unsupported datasize of " << dataSize << "!" << endl;
            return MPI_DATATYPE_NULL;
            break;
         }
       case datatype::FLOAT:
         switch (dataSize) {
          case (sizeof(float)):
            return MPI_Type<float>();
            break;
          case (sizeof(double)):
            return MPI_Type<double>();
            break;
          case (sizeof(long double)):
            return MPI_Type<long double>();
            break;
          default:
            cerr << "(VLSV) ERROR: VLSV::getMPIDatatype called with datatype:FLOAT datatype and unsupported datasize of " << dataSize << "!" << endl;
            return MPI_DATATYPE_NULL;
            break;
         }
       default:
         return MPI_DATATYPE_NULL;
         break;
      }
   }

} // namespace vlsv
