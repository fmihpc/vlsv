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

#include "mpiconversion.h"
#include "vlsv_common.h"

using namespace std;

namespace vlsv {
   
   /** Check that all processes have the same success status.
    * @param myStatus If true, this process has successfully executed all code.
    * @return If true, all processes called checkSuccess with myStatus set to 'true'.*/
   bool checkSuccess(const bool& myStatus,MPI_Comm comm) {
      // If myStatus is false, set mySuccess to value 1:
      int32_t mySuccess = 0;
      int32_t globalSuccess = 0;
      if (myStatus == false) mySuccess = 1;

      // Sum mySuccess values to master process and broadcast to all processes:
      MPI_Reduce(&mySuccess,&globalSuccess,1,MPI_Type<int32_t>(),MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&globalSuccess,1,MPI_Type<int32_t>(),0,MPI_COMM_WORLD);

      // If globalSuccess equals zero all processes called this function with myStatus set to 'true':
      if (globalSuccess == 0) return true;
      return false;
   }

   uint64_t convUInt64(const char* const ptr,const bool& swapEndian) {
      if (swapEndian == false) return *(reinterpret_cast<const uint64_t*>(ptr));
      int index = 0;
      uint64_t tmp = 0;
      char* const ptrtmp = reinterpret_cast<char*>(&tmp);
      for (int i=sizeof(uint64_t)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
   
   unsigned char detectEndianness() {
      const int number = 1;
      const char* const ptr = reinterpret_cast<const char*>(&number);
      if (ptr[0] == 1) return datatype::ENDIANNESS_LITTLE;
      else return datatype::ENDIANNESS_BIG;
   }

   const std::string& getMeshGeometry(geometry::type geom) {
      switch (geom) {
       case geometry::UNKNOWN:
	 return geometry::STRING_UNKNOWN;
	 break;
       case geometry::CARTESIAN:
	 return geometry::STRING_CARTESIAN;
	 break;
       case geometry::CYLINDRICAL:
	 return geometry::STRING_CYLINDRICAL;
	 break;
       case geometry::SPHERICAL:
	 return geometry::STRING_SPHERICAL;
	 break;
       case geometry::UNSTRUCTURED:
	 return geometry::STRING_UNSTRUCTURED;
	 break;
       default:
	 return geometry::STRING_UNKNOWN;
	 break;
      }
   }

   geometry::type getMeshGeometry(const std::string& s) {
      if (s == geometry::STRING_UNKNOWN) return geometry::UNKNOWN;
      else if (s == geometry::STRING_CARTESIAN) return geometry::CARTESIAN;
      else if (s == geometry::STRING_CYLINDRICAL) return geometry::CYLINDRICAL;
      else if (s == geometry::STRING_SPHERICAL) return geometry::SPHERICAL;
      else if (s == geometry::STRING_UNSTRUCTURED) return geometry::UNSTRUCTURED;
      return geometry::UNKNOWN;
   }

   datatype::type getVLSVDatatype(const std::string& s) {
      if (s == "int") return datatype::INT;
      else if (s == "uint") return datatype::UINT;
      else if (s == "float") return datatype::FLOAT;
      else return datatype::UNKNOWN;
   }

} // namespace vlsv
