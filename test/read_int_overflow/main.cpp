#include <cstdlib>
#include <iostream>
#include <map>

#include "../../vlsv_common.h"
#include "../../vlsv_writer.h"

using namespace std;
using namespace vlsv;

int main(int argn,char* args[]) {
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
   cout << "elements is " << elements << endl;

   if (myrank == 0) {
      string perProcessUnits = "B";
      double bytesPerProcess = elements*sizeof(size_t);
      if (bytesPerProcess >= 1e9) {bytesPerProcess /= 1e9; perProcessUnits = "GB";}
      else if (bytesPerProcess >= 1e6) {bytesPerProcess /= 1e6; perProcessUnits = "MB";}

      string totalBytesUnits = "B";
      double totalBytes = elements*sizeof(size_t)*processes;
      if (totalBytes >= 1e9) {totalBytes /= 1e9; totalBytesUnits = "GB";}
      else if (totalBytes >= 1e6) {totalBytes /= 1e6; totalBytesUnits = "MB";}

      cout << "Writing " << bytesPerProcess << ' ' << perProcessUnits << " per process, " << totalBytes << ' ' << totalBytesUnits << " in total" << endl;
   }

   size_t* array = new size_t[elements];
   for (size_t i=0; i<elements; ++i) array[i] = i;

   char* ptr = reinterpret_cast<char*>(array);
   Writer vlsvWriter;
   vlsvWriter.open("test_file.vlsv",MPI_COMM_WORLD,0);

   map<string,string> attribs;
   double t_start = MPI_Wtime();
   if (vlsvWriter.writeArray("ARRAY",attribs,getStringDatatype<size_t>(),elements,1,sizeof(size_t),ptr) == false) {
      cerr << "Failed to write array" << endl;
   }
   cout << "Approximate output datarate was " << (1.0*elements*sizeof(size_t)*processes)/1e9/(MPI_Wtime()-t_start) << " GB/s" << endl;

   vlsvWriter.close();
   delete [] array; array = NULL;
   MPI_Finalize();
   return 0;
}
