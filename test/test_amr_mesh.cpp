#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <map>
#include <cmath>
#include <mpi.h>
#include <time.h>

#include "amr_mesh.h"

using namespace std;

int CallbackCoarsen(const amr::GlobalID siblingIDs[8],const amr::LocalID siblingIndices[8],
                    const amr::GlobalID globalID,amr::LocalID& index) {
   cout << "Coarsening following blocks" << endl;
   for (int i=0; i<8; ++i) {
      cout << '\t' << siblingIDs[i] << '\t' << siblingIndices[i] << endl;
   }
   
   index = amr::INVALID_LOCALID;
   return true;
}

int CallbackCreate(const amr::GlobalID& globalID,amr::LocalID& index) {
   cout << "Block " << globalID << " created" << endl;
   return true;
}

int CallbackDelete(const amr::GlobalID& globalID,const amr::LocalID& index) {
   cout << "Block " << globalID << " deleted, data stored in " << index << endl;
   return true;
}

int CallbackRefine(const amr::GlobalID& globalID,const amr::LocalID& index,
                   const amr::GlobalID childrenIDs[8],amr::LocalID childrenIndices[8]) {
   cout << "Block " << globalID << " refined, data stored in " << index << "\t. Creating following blocks" << endl;
   for (int i=0; i<8; ++i) {
      childrenIndices[i] = amr::INVALID_LOCALID;
      cout << '\t' << childrenIDs[i] << '\t' << childrenIndices[i] << endl;
   }
   return true;
}

int main(int argn,char* args[]) {
   bool success = true;

   // Init MPI:
   MPI_Init(&argn,&args);
   
   int Nx0 = 5;
   int Ny0 = 5;
   int Nz0 = 5;
   int xCells = 4;
   int yCells = 1;
   int zCells = 1;
   int maxRefLevel = 4;

   amr::AmrMesh mesh(Nx0,Ny0,Nz0,xCells,yCells,zCells,maxRefLevel);
   if (mesh.registerCallbacks(CallbackCoarsen,CallbackCreate,CallbackDelete,CallbackRefine) == false) {
      cerr << "failed to register callbacks" << endl;
      return 1;
   }
   if (mesh.initialize(0,Nx0,0,Ny0,0,Nz0,0) == false) {
      cerr << "mesh failed to init" << endl;
      return 1;
   }

   srand(time(NULL));
   const int N_refines = 200;
   double t_ref_total = 0.0;
   for (int i=0; i<N_refines; ++i) {
      const uint64_t index = mesh.size()*1.0*rand()/RAND_MAX;

      unordered_map<amr::GlobalID,amr::LocalID>::iterator it = mesh.begin();
      for (int j=0; j<index; ++j) ++it;

      cout << "Refinement " << i+1 << "/" << N_refines << endl;
      const double t_ref_start = MPI_Wtime();
      if (mesh.refine(it->first) == false) {
         cerr << "refine failed" << endl;
         //return 1;
      }
      t_ref_total += (MPI_Wtime() - t_ref_start);
   }

   if (mesh.checkMesh() == true) {
      cout << "Mesh checks OK after refines, time \t" << t_ref_total/N_refines << "\t per refine" << endl;
   } else {
      cout << "Mesh IS NOT OK after refines" << endl;
   }
   
   if (mesh.write("amr_refined.vlsv") == false) {
      cerr << "failed to write mesh" << endl;
      return 1;
   }

   /*
   const int N_coarsens = 2000;
   double t_coa_total = 0.0;
   for (int i=0; i<N_coarsens; ++i) {
      const uint64_t index = mesh.size()*1.0*rand()/RAND_MAX;
      
      unordered_map<amr::GlobalID,amr::LocalID>::iterator it = mesh.begin();
      for (int j=0; j<index; ++j) ++it;
      
      const double t_coa_start = MPI_Wtime();
      mesh.coarsen(it->first);
      t_coa_total += (MPI_Wtime() - t_coa_start);
   }
   
   if (mesh.checkMesh() == true) {
      cout << "Mesh checks OK after coarsenings, time \t" << t_coa_total/N_coarsens << "\t per coarsen" << endl;
   } else {
      cout << "Mesh IS NOT OK after coarsenings" << endl;
   }
   
   if (mesh.write("amr_coarsen.vlsv") == false) {
      cerr << "failed to write mesh" << endl;
      return 1;
   }*/

   MPI_Finalize();
   if (success == false) return 1;
   return 0;
}
