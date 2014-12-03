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

} // namespace vlsv
