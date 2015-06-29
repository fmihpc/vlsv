/* This file is part of VLSV file format.
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

#include "mpiconversion.h"
#include "vlsv_common.h"
#include "vlsv_common_mpi.h"

using namespace std;

namespace vlsv {

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
