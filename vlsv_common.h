/** This file is part of VLSV file format.
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

#ifndef VLSV_COMMON_H
#define VLSV_COMMON_H

#include <stdint.h>

namespace VLSV {
   const unsigned char LITTLE_END = 0;
   const unsigned char BIG_END    = 1;
   
   enum datatype {UNKNOWN,INT,UINT,FLOAT};
   
   const std::string MESH_POINT = "point";
   const std::string MESH_QUAD = "quad";
}

unsigned char detectEndianness();

int8_t convInt8(const char* const ptr,const bool& swapEndian=false);
int16_t convInt16(const char* const ptr,const bool& swapEndian=false);
int32_t convInt32(const char* const ptr,const bool& swapEndian=false);
int64_t convInt64(const char* const ptr,const bool& swapEndian=false);
uint8_t convUInt8(const char* const ptr,const bool& swapEndian=false);
uint16_t convUInt16(const char* const ptr,const bool& swapEndian=false);
uint32_t convUInt32(const char* const ptr,const bool& swapEndian=false);
uint64_t convUInt64(const char* const ptr,const bool& swapEndian=false);
float convReal4(const char* const ptr,const bool& swapEndian=false);
double convReal8(const char* const ptr,const bool& swapEndian=false);

#endif
