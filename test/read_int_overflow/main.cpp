#include <cstdlib>
#include <iostream>
#include <list>
#include <map>

#include "../../vlsv_common.h"
#include "../../vlsv_writer.h"
#include "../../vlsv_reader_parallel.h"

using namespace std;
using namespace vlsv;

bool read(const int& myrank,const size_t& elements,double& dataRate) {
   bool success = true;

   // Calculate offset from array start
   const size_t offset = myrank*elements;

   size_t* array = new size_t[elements];
   for (size_t i=0; i<elements; ++i) array[i] = numeric_limits<size_t>::max();
   char* ptr = reinterpret_cast<char*>(array);

   MPI_Aint address;
   MPI_Get_address(array,&address);
   //cerr << "P#" << myrank << " array address is " << address << endl;

   // Read data from file
   ParallelReader vlsvReader;
   if (vlsvReader.open("test_file.vlsv",MPI_COMM_WORLD,0) == false) {
      cerr << "failed to open input file!" << endl;
      delete [] array; return false;
   }

   list<pair<string,string> > attribs;
   if (vlsvReader.readArray("ARRAY",attribs,offset,elements,ptr) == false) {
      cerr << "readArray failed!" << endl;
      success = false;
   }
   dataRate = vlsvReader.getBytesRead()/vlsvReader.getReadTime();

   // Verify that the data is correct:
   size_t errors = 0;
   for (size_t i=0; i<elements; ++i) {      
      if (array[i] != i) {
         if (errors == 0) {
            stringstream ss;
            ss << "P#" << myrank;// << " array contents invalid!" << endl;
            ss << " element " << i << " should equal " << i << " but is has value " << array[i] << endl;
            cerr << ss.str();
         }
         success = false;
         ++errors;
      }
   }
   if (errors > 0) {
      cerr << "I found " << errors << " errors in total" << endl;
   }

   if (vlsvReader.close() == false) success = false;
   delete [] array; array = NULL;

   // Return the same value on all processes
   int globalResult;
   int result = 0;
   if (success == false) result = 1;
   MPI_Allreduce(&result,&globalResult,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   if (globalResult == 0) return true;
   else return false;
}

bool verifyWriteMaster(const size_t& elements) {
   bool success = true;
   
   ifstream in;
   in.open("test_file.vlsv",ifstream::in | ifstream::binary);
   if (in.good() == false) {
      cerr << "Failed to open file for verifying" << endl;
      return false;
   }

   // Go to the start of data
   in.seekg(16,ios_base::beg);

   // Make array for data
   size_t* array = new size_t[elements];
   
   double t_read = 0;
   int processes;
   MPI_Comm_size(MPI_COMM_WORLD,&processes);
   for (int p=0; p<processes; ++p) {
      // Invalidate array contents
      for (size_t i=0; i<elements; ++i) array[i] = numeric_limits<size_t>::max();
      char* ptr = reinterpret_cast<char*>(array);
   
      // Read data
      double t_start = MPI_Wtime();
      in.read(ptr,elements*sizeof(size_t));
      t_read += (MPI_Wtime() - t_start);

      // Verify
      size_t errors = 0;
      size_t firstErrorPos = 0;
      for (size_t i=0; i<elements; ++i) {
         if (errors == 0) firstErrorPos = i;
         if (array[i] != i) {
            ++errors;
            success = false;
         }
      }
      if (errors > 0) {
         cout << "Found " << errors << " errors in written data from process " << p << " first erroneous index " << firstErrorPos << " data rate ";
         cout << elements*sizeof(size_t)/t_read/1.0e9 << endl;
      } else {
         //cout << "Process " << p << " data correct, read data rate " << elements*sizeof(size_t)/t_read/1.0e9 << " GB/s" << endl;
      }
   }

   delete [] array; array = NULL;
   return success;
}

bool write(const size_t& elements,double& dataRate) {
   bool success = true;
   
   // Create data that is written to file
   size_t* array = new size_t[elements];
   for (size_t i=0; i<elements; ++i) array[i] = i;
   char* ptr = reinterpret_cast<char*>(array);

   Writer vlsvWriter;
   if (vlsvWriter.open("test_file.vlsv",MPI_COMM_WORLD,0) == false) {delete [] array; return false;}

   map<string,string> attribs;
   if (vlsvWriter.writeArray("ARRAY",attribs,getStringDatatype<size_t>(),elements,1,sizeof(size_t),ptr) == false) {
      success = false;
   }
   dataRate = vlsvWriter.getBytesWritten()/vlsvWriter.getWriteTime();
   vlsvWriter.close();
   delete [] array; array = NULL;
   
   // Check that data on file is correct:
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   if (myrank == 0) {
      if (verifyWriteMaster(elements) == true) {
         //cout << "Write test result: PASS" << endl;
      } else {
         //cout << "Write test result: FAILURE" << endl;
         success = false;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
   
   int myResult = 0;
   if (success == false) myResult = 1;
   int globalResult;
   MPI_Allreduce(&myResult,&globalResult,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
   if (globalResult == 0) return true;
   else return false;
}

int main(int argn,char* args[]) {
   int rvalue = 0;
   int myrank,processes;
   MPI_Init(&argn,&args);
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   MPI_Comm_size(MPI_COMM_WORLD,&processes);

   if (argn != 2) {
      cerr << endl;
      cerr << "USAGE: ./main <elements>" << endl;
      cerr << "<elements> Number of size_t values written per process." << endl;
      cerr << endl;
      return 1;
   }

   size_t elements = atol(args[1]);
   //cout << "elements is " << elements << endl;

   if (myrank == 0) {
      string perProcessUnits = "B";
      double bytesPerProcess = elements*sizeof(size_t);
      if (bytesPerProcess >= 1e9) {bytesPerProcess /= 1e9; perProcessUnits = "GB";}
      else if (bytesPerProcess >= 1e6) {bytesPerProcess /= 1e6; perProcessUnits = "MB";}

      string totalBytesUnits = "B";
      double totalBytes = elements*sizeof(size_t)*processes;
      if (totalBytes >= 1e9) {totalBytes /= 1e9; totalBytesUnits = "GB";}
      else if (totalBytes >= 1e6) {totalBytes /= 1e6; totalBytesUnits = "MB";}

      cout << "Writing " << bytesPerProcess << ' ' << perProcessUnits << " per process using ";
      cout << processes << " MPI processes, " << totalBytes << ' ' << totalBytesUnits << " in total" << endl;
   }

   // Write data to file
   double writeDataRate;
   if (write(elements,writeDataRate) == true) {
      if (myrank == 0) {
         string dataRateUnits = "B/s";
         if (writeDataRate > 1e9) {writeDataRate /= 1e9; dataRateUnits = "GB/s";}
         else if (writeDataRate > 1e6) {writeDataRate /= 1e6; dataRateUnits = "MB/s";}
         cout << "Write test: SUCCESS \t\t approximate datarate " << writeDataRate << ' ' << dataRateUnits << endl;
      }
   } else {
      cout << "Write test: FAILED" << endl;
      rvalue = 1;
   }
   if (rvalue != 0) {
      MPI_Finalize();
      return rvalue;
   }
   
   // Read data from file
   double readDataRate;
   if (read(myrank,elements,readDataRate) == true) {
      if (myrank == 0) {
         string dataRateUnits = "B/s";
         if (readDataRate > 1e9) {readDataRate /= 1e9; dataRateUnits = "GB/s";}
         else if (readDataRate > 1e6) {readDataRate /= 1e6; dataRateUnits = "MB/s";}
         cout << "Read test : SUCCESS \t\t approximate datarate " << readDataRate << ' ' << dataRateUnits << endl;
      }
   } else {
      cout << "Read test : FAILED" << endl;
      rvalue = 1;
   }
   
   MPI_Finalize();
   return rvalue;
}
