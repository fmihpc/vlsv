/** This file is part of VLSV file format.
 * 
 *  Copyright 2011, 2012 Finnish Meteorological Institute
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

unsigned char detectEndianness() {
   const int number = 1;
   const char* const ptr = reinterpret_cast<const char*>(&number);
   if (ptr[0] == 1) return VLSV::LITTLE_END;
   else return VLSV::BIG_END;
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
