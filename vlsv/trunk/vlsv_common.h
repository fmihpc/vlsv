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
