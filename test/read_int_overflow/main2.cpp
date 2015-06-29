#include <cstdlib>
#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argn,char* args[]) {
   if (argn != 2) {
      cerr << endl;
      cerr << "USAGE: ./main2 <entries>" << endl;
      cerr << "Size of entry is 400 MB, so using entries=5 should work but entries=6 fails" << endl;
      cerr << endl;
      return 1;
   }
   const int entries = atoi(args[1]);

   int myrank,processes;
   MPI_Init(&argn,&args);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   
   // Open output file
   int accessMode = (MPI_MODE_WRONLY | MPI_MODE_CREATE);
   string filename = "test_file.dat";
   MPI_File file_ptr;
   MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.c_str()),accessMode,MPI_INFO_NULL,&file_ptr);

   // One struct entry is 50e6*sizeof(size_t) bytes = 400 MB
   // so entries=5 works but entries=6 fails
   size_t elements = 50000000; // 50e6

   // Create array containing all output data
   size_t* array = new size_t[entries*elements];
   for (size_t i=0; i<entries*elements; ++i) array[i] = i;

   // Create MPI_Struct that should write all data at once
   int* blockLengths       = new int[entries];
   MPI_Aint* displacements = new MPI_Aint[entries];
   MPI_Datatype* datatypes = new MPI_Datatype[entries];
   for (int i=0; i<entries; ++i) {
      blockLengths[i]  = elements;
      displacements[i] = (MPI_Aint)(i*entries*elements);
      datatypes[i]     = MPI_UNSIGNED_LONG;
   }
   MPI_Datatype outputType;
   MPI_Type_create_struct(entries,blockLengths,displacements,datatypes,&outputType);
   MPI_Type_commit(&outputType);
   
   // Write the data and close file
   MPI_File_write_at_all(file_ptr,0,array,1,outputType,MPI_STATUS_IGNORE);
   MPI_Type_free(&outputType);
   MPI_File_close(&file_ptr);

   delete [] array;
   delete [] blockLengths;
   delete [] displacements;
   delete [] datatypes;
   MPI_Finalize();
   return 0;
}
