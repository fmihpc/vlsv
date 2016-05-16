/* This file is part of VLSV file format.
 * 
 *  Copyright 2011-2016 Finnish Meteorological Institute
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
#include <cstring>

#include "mpiconversion.h"
#include "vlsv_common.h"
#include "vlsv_common_mpi.h"

using namespace std;

static const uint64_t MAX_MPI_FILE_IO_BYTES = 2147460000;

namespace vlsv {

   /** Get maximum number of bytes that can be read from a file using a single collective MPI routine.
    * This should equal to numeric_limits<MPI_Count>::max() in MPI 3.0+, but in older 
    * library versions the 'count' parameter is simply an integer.
    * @return Maximum number of bytes read using a single MPI collective routine.*/
   uint64_t getMaxBytesPerRead() {
      // For some obscure reason OpenMPI can only write 2147479552 bytes with a single 
      // collective call. So I'm manually setting the max bytes to bit less than that.
      return MAX_MPI_FILE_IO_BYTES;
   }
   
   /** Get maximum number of bytes that can be written to a file using a single collective MPI routine.
    * This should equal to numeric_limits<MPI_Count>::max() in MPI 3.0+, but in older 
    * library versions the 'count' parameter is simply an integer.
    * @return Maximum number of bytes written using a single MPI collective routine.*/
   uint64_t getMaxBytesPerWrite() {
      return MAX_MPI_FILE_IO_BYTES;
   }


   bool broadcast(const std::string& input,std::string& output,MPI_Comm comm,const int& masterRank) {
      bool success = true;
      size_t length = input.size();
      if (MPI_Bcast(&length,1,MPI_Type<size_t>(),masterRank,comm) != MPI_SUCCESS) success = false;
      if (checkSuccess(success,comm) == false) return false;

      int myRank;
      MPI_Comm_rank(comm,&myRank);

      // Copy string from input to tmp on master process only:
      char* tmp = new char[length+1];
      if (myRank == masterRank) {
         #ifdef WINDOWS
            strncpy_s(tmp,length+1,input.c_str(),length+1);
         #else
            strncpy(tmp,input.c_str(),length+1);
         #endif
      }

      // Broadcast input string to all processes:
      if (MPI_Bcast(tmp,length+1,MPI_Type<char>(),masterRank,comm) != MPI_SUCCESS) success = false;
      output = tmp;
      delete [] tmp; tmp = nullptr;
      return checkSuccess(success,comm);
   }

   /** Check that all processes in communicator comm called this function 
    * with myStatus=true.
    * @param myStatus The status flag of this process.
    * @param comm MPI Communicator.
    * @return If true, all processes called this function with myStatus=true.*/
   bool checkSuccess(const bool& myStatus,MPI_Comm comm) {
      // If myStatus is false, set mySuccess to value 1:
      int32_t mySuccess = 0;
      int32_t globalSuccess = 0;
      if (myStatus == false) mySuccess = 1;

      // Sum mySuccess values to master process and broadcast to all processes:
      MPI_Reduce(&mySuccess,&globalSuccess,1,MPI_Type<int32_t>(),MPI_SUM,0,comm);
      MPI_Bcast(&globalSuccess,1,MPI_Type<int32_t>(),0,comm);

      // If globalSuccess equals zero all processes called this function with myStatus set to 'true':
      if (globalSuccess == 0) return true;
      return false;
   }

   /** Get the MPI primitive datatype corresponding to the given VLSV datatype.
    * Note that the VLSV 'datatype' here is either INT, UINT, or FLOAT.
    * @param dt VLSV datatype.
    * @param dataSize Byte size of the VLSV datatype.
    * @return MPI Datatype that corresponds to the given VLSV datatype.*/
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
          #ifndef WINDOWS
          case (sizeof(long double)):
            return MPI_Type<long double>();
            break;
	      #endif
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
