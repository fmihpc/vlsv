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
