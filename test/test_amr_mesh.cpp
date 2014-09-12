#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <map>
#include <cmath>
#include <mpi.h>
#include <time.h>

#include "amr_mesh.h"

using namespace std;

int main(int argn,char* args[]) {
   bool success = true;

   // Init MPI:
   MPI_Init(&argn,&args);
   
   int Nx0 = 5;
   int Ny0 = 5;
   int Nz0 = 1;
   int xCells = 1;
   int yCells = 1;
   int zCells = 1;
   int maxRefLevel = 8;

   AmrMesh mesh(Nx0,Ny0,Nz0,xCells,yCells,zCells,maxRefLevel);
   if (mesh.initialize(0,Nx0,0,Ny0,0,Nz0,0) == false) {
      cerr << "mesh failed to init" << endl;
      return 1;
   }

   srand(time(NULL));
   for (int i=0; i<100; ++i) {
      const uint64_t index = mesh.size()*1.0*rand()/RAND_MAX;
      
      unordered_map<uint64_t,uint8_t>::iterator it = mesh.begin();
      for (int j=0; j<index; ++j) ++it;
      if (mesh.refine(it->first) == false) {
	 cerr << "refine failed" << endl;
	 return 1;
      }
   }
   
   if (mesh.checkMesh() == true) {
      cout << "Mesh checks OK after refines" << endl;
   } else {
      cout << "Mesh IS NOT OK after refines" << endl;
   }
   
   if (mesh.write("amr_refined.vlsv") == false) {
      cerr << "failed to write mesh" << endl;
      return 1;
   }

   for (int i=0; i<150; ++i) {
      const uint64_t index = mesh.size()*1.0*rand()/RAND_MAX;
      
      unordered_map<uint64_t,uint8_t>::iterator it = mesh.begin();
      for (int j=0; j<index; ++j) ++it;
      mesh.coarsen(it->first);
   }
   
   if (mesh.checkMesh() == true) {
      cout << "Mesh checks OK after coarsenings" << endl;
   } else {
      cout << "Mesh IS NOT OK after coarsenings" << endl;
   }
   
   if (mesh.write("amr_coarsen.vlsv") == false) {
      cerr << "failed to write mesh" << endl;
      return 1;
   }

   MPI_Finalize();
   if (success == false) return 1;
   return 0;
}




