#include <cstdlib>
#include <iostream>
#include <sstream>
#include <mpi.h>
#include <limits>

using namespace std;

const string filename = "test_file.dat";

void write_contiguous(const int& myrank,const int& processes,const size_t& entries,const size_t& elements,
                      const bool& useBottomAsBase) {
   // Open output file
   MPI_File file_ptr;
   MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&file_ptr);

   // Open output file
   size_t* array = new size_t[entries*elements];
   for (size_t i=0; i<entries*elements; ++i) array[i] = i;

   // File offset (in bytes) for this process
   MPI_Offset myFileOffset = myrank*entries*elements*sizeof(size_t);

   // Create MPI_Struct that should read or write all data at once
   int* blockLengths       = new int[entries];
   MPI_Aint* displacements = new MPI_Aint[entries];
   MPI_Datatype* datatypes = new MPI_Datatype[entries];
   MPI_Aint baseAddress;
   if (useBottomAsBase == true) {
      baseAddress = 0;
   } else {
      MPI_Get_address(array,&baseAddress);
   }

   for (int i=0; i<entries; ++i) {
      MPI_Aint address;
      MPI_Get_address(array+i*elements,&address);
      
      blockLengths[i]  = elements;
      displacements[i] = address - baseAddress;
      datatypes[i]     = MPI_UNSIGNED_LONG;
   }
   MPI_Datatype outputType;
   MPI_Type_create_struct(entries,blockLengths,displacements,datatypes,&outputType);
   MPI_Type_commit(&outputType);

   // Write the data
   MPI_Barrier(MPI_COMM_WORLD);
   double t_start_write = MPI_Wtime();
   if (useBottomAsBase == true) 
     MPI_File_write_at_all(file_ptr,myFileOffset,MPI_BOTTOM,1,outputType,MPI_STATUS_IGNORE);
   else
     MPI_File_write_at_all(file_ptr,myFileOffset,array,1,outputType,MPI_STATUS_IGNORE);
   MPI_Barrier(MPI_COMM_WORLD);
   double t_end = MPI_Wtime();
   double t_end_write;
   MPI_Reduce(&t_end,&t_end_write,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (myrank == 0) {
      stringstream ss;
      ss << "Wrote " << processes*entries*elements*sizeof(size_t)/1.0e9 << " GB data in " << (t_end_write-t_start_write) << " seconds, ";
      ss << "approx. data rate " << processes*entries*elements*sizeof(size_t)/1.0e9/(t_end_write-t_start_write) << " GB/s" << endl;
      cout << ss.str() << flush;
   }

   // Reallocate the MPI datatype
   MPI_Type_free(&outputType);
   MPI_File_close(&file_ptr);   
   delete [] array;
   delete [] blockLengths;
   delete [] displacements;
   delete [] datatypes;
}

void write_random_chunks(const int& myrank,const int& processes,const size_t& entries,const size_t& elements,
                         const bool& useBottomAsBase) {
   // Open output file
   MPI_File file_ptr;
   MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.c_str()),MPI_MODE_WRONLY | MPI_MODE_CREATE,MPI_INFO_NULL,&file_ptr);
   
   // Open output file
   size_t** container = new size_t* [entries];
   for (size_t i=0; i<entries; ++i) container[i] = new size_t[elements];
   for (size_t i=0; i<entries; ++i) for (size_t j=0; j<elements; ++j) {
      container[i][j] = i*elements+j;
   }
   
   // File offset (in bytes) for this process
   MPI_Offset myFileOffset = myrank*entries*elements*sizeof(size_t);
   
   // Create MPI_Struct that should read or write all data at once
   int* blockLengths       = new int[entries];
   MPI_Aint* displacements = new MPI_Aint[entries];
   MPI_Datatype* datatypes = new MPI_Datatype[entries];
   MPI_Aint baseAddress;
   if (useBottomAsBase == true) {
      baseAddress = 0;
   } else {
      MPI_Get_address(container[0],&baseAddress);
   }

   for (int i=0; i<entries; ++i) {
      MPI_Aint address;
      MPI_Get_address(container[i],&address);

      blockLengths[i]  = elements;
      displacements[i] = address - baseAddress;
      datatypes[i]     = MPI_UNSIGNED_LONG;
   }
   
   // Randomize the chunk order:
   srand(time(NULL));
   for (int n=0; n<500; ++n) {
      int i = static_cast<int>(entries*1.0*rand()/RAND_MAX);
      int j = static_cast<int>(entries*1.0*rand()/RAND_MAX);
      if (i == j) continue;
      if (i < 0 || i >= entries) continue;
      if (j < 0 || j >= entries) continue;
      const MPI_Aint dummyDisp = displacements[i];
      displacements[i] = displacements[j];
      displacements[j] = dummyDisp;
   }

   MPI_Datatype outputType;
   MPI_Type_create_struct(entries,blockLengths,displacements,datatypes,&outputType);
   MPI_Type_commit(&outputType);

   // Write the data
   MPI_Barrier(MPI_COMM_WORLD);
   double t_start_write = MPI_Wtime();
   if (useBottomAsBase == true)
     MPI_File_write_at_all(file_ptr,myFileOffset,MPI_BOTTOM,1,outputType,MPI_STATUS_IGNORE);
   else
     MPI_File_write_at_all(file_ptr,myFileOffset,container[0],1,outputType,MPI_STATUS_IGNORE);
   MPI_Barrier(MPI_COMM_WORLD);
   double t_end = MPI_Wtime();
   double t_end_write;
   MPI_Reduce(&t_end,&t_end_write,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   
   if (myrank == 0) {
      stringstream ss;
      ss << "Wrote " << processes*entries*elements*sizeof(size_t)/1.0e9 << " GB data in " << (t_end_write-t_start_write) << " seconds, ";
      ss << "approx. data rate " << processes*entries*elements*sizeof(size_t)/1.0e9/(t_end_write-t_start_write) << " GB/s" << endl;
      cout << ss.str() << flush;
   }
   
   // Reallocate the MPI datatype
   MPI_Type_free(&outputType);
   MPI_File_close(&file_ptr);
   for (size_t i=0; i<entries; ++i) {delete [] container[i]; container[i] = NULL;}
   delete [] container;
   delete [] blockLengths;
   delete [] displacements;
   delete [] datatypes;
}

void read_simple(const int& myrank,const int& processes,const size_t& entries,const size_t& elements) {
   MPI_File file_ptr;
   MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.c_str()),MPI_MODE_RDONLY,MPI_INFO_NULL,&file_ptr);
   
   // File offset (in bytes) for this process
   MPI_Offset myFileOffset = myrank*entries*elements*sizeof(size_t);
   
   // Create array for data
   size_t* array = new size_t[entries*elements];
   for (size_t i=0; i<entries*elements; ++i) array[i] = numeric_limits<size_t>::max();

   // Read data and measure how long it took
   bool success = true;
   MPI_Barrier(MPI_COMM_WORLD);
   double t_start_read = MPI_Wtime();
   MPI_File_read_at_all(file_ptr,myFileOffset,array,entries*elements,MPI_UNSIGNED_LONG,MPI_STATUS_IGNORE);
   MPI_Barrier(MPI_COMM_WORLD);
   double t_end = MPI_Wtime();
   double t_end_read;
   MPI_Reduce(&t_end,&t_end_read,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   
   if (myrank == 0) {
      stringstream ss;
      ss << "Read  " << processes*entries*elements*sizeof(size_t)/1.0e9 << " GB data in " << (t_end_read-t_start_read) << " seconds, ";
      ss << "approx. data rate " << processes*entries*elements*sizeof(size_t)/1.0e9/(t_end_read-t_start_read) << " GB/s" << endl;
      cout << ss.str() << flush;
   }

   // Check that data was read correctly
   size_t incorrectEntries = 0;
   for (size_t i=0; i<entries*elements; ++i) {
      if (array[i] != i) {
         ++incorrectEntries;
         success = false;
      }
   }
   if (incorrectEntries > 0) {
      cerr << "Found " << incorrectEntries << " ERROR(s) in read data" << endl;
   } else {
      //cerr << "Input data correct" << endl;
   }
   
   delete [] array;
   MPI_File_close(&file_ptr);
}

void read_contiguous(const int& myrank,const int& processes,const size_t& entries,const size_t& elements,void* basePtr,
                     MPI_Aint& baseAddress,const bool& useArrayAsBase) {
   MPI_File file_ptr;
   MPI_File_open(MPI_COMM_WORLD,const_cast<char*>(filename.c_str()),MPI_MODE_RDONLY,MPI_INFO_NULL,&file_ptr);

   // File offset (in bytes) for this process
   MPI_Offset myFileOffset = myrank*entries*elements*sizeof(size_t);

   // Create array for data
   size_t* array = new size_t[entries*elements];
   for (size_t i=0; i<entries*elements; ++i) array[i] = numeric_limits<size_t>::max();

   if (useArrayAsBase == true) {
      basePtr = array;
      MPI_Get_address(basePtr,&baseAddress);
   }

   // Create MPI_Struct that should read or write all data at once
   int* blockLengths       = new int[entries];
   MPI_Aint* displacements = new MPI_Aint[entries];
   MPI_Datatype* datatypes = new MPI_Datatype[entries];
   for (int i=0; i<entries; ++i) {
      MPI_Aint address;
      MPI_Get_address(array+i*elements,&address);
      
      blockLengths[i]  = elements;
      displacements[i] = address - baseAddress;
      datatypes[i]     = MPI_UNSIGNED_LONG;
   }
   MPI_Datatype outputType;
   MPI_Type_create_struct(entries,blockLengths,displacements,datatypes,&outputType);
   MPI_Type_commit(&outputType);

   // Read data and measure how long it took
   bool success = true;
   MPI_Barrier(MPI_COMM_WORLD);
   double t_start_read = MPI_Wtime();
   MPI_File_read_at_all(file_ptr,myFileOffset,basePtr,1,outputType,MPI_STATUS_IGNORE);
   MPI_Barrier(MPI_COMM_WORLD);
   double t_end = MPI_Wtime();
   double t_end_read;
   MPI_Reduce(&t_end,&t_end_read,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (myrank == 0) {
      stringstream ss;
      ss << "Read  " << processes*entries*elements*sizeof(size_t)/1.0e9 << " GB data in " << (t_end_read-t_start_read) << " seconds, ";
      ss << "approx. data rate " << processes*entries*elements*sizeof(size_t)/1.0e9/(t_end_read-t_start_read) << " GB/s" << endl;
      cout << ss.str() << flush;
   }

   // Check that data was read correctly
   size_t incorrectEntries = 0;
   for (size_t i=0; i<entries*elements; ++i) {
      if (array[i] != i) {
         ++incorrectEntries;
         success = false;
      }  
   }
   if (incorrectEntries > 0) {
      cerr << "Found " << incorrectEntries << " ERROR(s) in read data" << endl;
   } else {
      //cerr << "Input data correct" << endl;
   }

   MPI_Type_free(&outputType);
   delete [] datatypes;
   delete [] displacements;
   delete [] blockLengths;
   delete [] array;
   MPI_File_close(&file_ptr);
}

int main(int argn,char* args[]) {
   int myrank,processes;
   MPI_Init(&argn,&args);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);

   if (argn != 3) {
      if (myrank == 0) {
         cerr << endl;
         cerr << "USAGE: ./main2 <entries> <elements>" << endl;
         cerr << "Size of entry is 400 MB, so using entries=5 should work but entries=6 fails " << endl;
         cerr << "if you use elements=50000000." << endl;
         cerr << endl;
         cerr << "Each entry consists of <elements> 64-bit integers" << endl;
         cerr << endl;
      }

      MPI_Finalize();
      return 1;
   }
   double t_begin = MPI_Wtime();
   const int entries = atoi(args[1]);
   const size_t elements = atol(args[2]);

   MPI_Aint baseAddress = 0;
   if (myrank == 0) cout << endl << "Write contiguous data, either array pointer as base or MPI_BOTTOM:" << endl;
   write_contiguous(myrank,processes,entries,elements,false);
   //write_contiguous(myrank,processes,entries,elements,true);

   if (myrank == 0) cout << endl << "Reading contiguous data, basic MPI_Datatype only:" << endl;
   read_simple(myrank,processes,entries,elements);
   
   if (myrank == 0) cout << endl << "Reading contiguous data, either array pointer as base or MPI_BOTTOM:" << endl;
   read_contiguous(myrank,processes,entries,elements,NULL,baseAddress,true);
   read_contiguous(myrank,processes,entries,elements,NULL,baseAddress,false);

   if (myrank == 0) cout << endl << "Write chunks in random order, either array pointer as base or MPI_BOTTOM:" << endl;
   write_random_chunks(myrank,processes,entries,elements,false);
   write_random_chunks(myrank,processes,entries,elements,true);

   MPI_Finalize();
   return 0;
}
