#include <cstdlib>
#include <iostream>
#include <time.h>
#include <cmath>
#include <vector>
#include <algorithm>

#include "../vlsv_amr.h"

using namespace std;

/*int calcLevel(int GID,int size) {
   int value = GID/size;
   if (value == 0) return 0;

   double a = log10(1.0*GID/size) / log10(8.0);
   cout << GID << '\t' << value << '\t' << a << '\t' << 1+floor(a) << endl;
   return floor(a);
}*/

int main(int argn,char* args[]) {
   int size = 64;
   
   int maxRefLevel = 6;
   vlsv::initMesh(4,4,4,maxRefLevel);
   srand(time(NULL));

   int N = 30;
   for (int i=0; i<N; ++i) {
      uint64_t gid = static_cast<int>(2000*(1.0*rand())/RAND_MAX);
      
      uint32_t i_ind,j_ind,k_ind;
      uint32_t refLevel;
      vlsv::calculateCellIndices(gid,refLevel,i_ind,j_ind,k_ind);
      
      cout << gid << '\t' << refLevel << "\t indices: " << i_ind << ' ' << j_ind << ' ' << k_ind << endl;
   }

   return 0;
}



