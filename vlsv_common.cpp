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

#include <cstdlib>
#include <iostream>
#include <sstream>

#include "vlsv_common.h"

using namespace std;

namespace vlsv {

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

   std::string getStringDatatype(const vlsv::datatype::type& dt) {
      switch (dt) {
       case datatype::UNKNOWN:
         return "unknown";
         break;
       case datatype::INT:
         return "int";
         break;
       case datatype::UINT:
         return "uint";
         break;
       case datatype::FLOAT:
         return "float";
         break;
       default:
         return "unknown";
         break;
      }
   }
   
   datatype::type getVLSVDatatype(const std::string& s) {
      if (s == "int") return datatype::INT;
      else if (s == "uint") return datatype::UINT;
      else if (s == "float") return datatype::FLOAT;
      else return datatype::UNKNOWN;
   }

   /** Print the data rate corresponding to given number of bytes and time in human readable format.
    * @param bytes Number of bytes written or read.
    * @param t Time spent in reading or writing bytes in seconds.
    * @return String representation of the data rate.*/
   std::string printDataRate(const uint64_t& bytes,const double& t) {
      stringstream s;
      if (bytes/t > 1e12) s << bytes/t/1e12 << " TB/s";
      else if (bytes/t > 1e9) s << bytes/t/1e9 << " GB/s";
      else if (bytes/t > 1e6) s << bytes/t/1e6 << " MB/s";
      else if (bytes/t > 1e3) s << bytes/t/1e3 << " kB/s";
      else s << bytes/t << " B/s";
      return s.str();
   }

   const std::string getErrorString(const vlsv::error::type& errorCode) {
      switch (errorCode) {
       case error::NONE:
         return "No errors";
         break;
       case error::UNKNOWN:
         return "Unknown or unsupported error code";
         break;
       case error::READ_CWD_FAIL:
         return "Failed to cwd to input file dir";
         break;
       case error::READ_FILE_BAD:
         return "Failed to open input file";
         break;
       case error::READ_FILE_ALREADY_OPEN:
         return "vlsv::Reader already has an open input file";
         break;
       case error::READ_FILE_ENDIANNESS:
         return "Failed to read file endianness";
         break;
       case error::READ_NO_FOOTER:
         return "Input file broken, footer not found";
         break;
       case error::READ_FOOTER_OFFSET:
         return "Failed to read footer position";
         break;
       case error::READ_FOOTER:
         return "Failed to read footer";
         break;
       default:
         return "Unknown or unsupported error code";
         break;
      }
   }
   
} // namespace vlsv
