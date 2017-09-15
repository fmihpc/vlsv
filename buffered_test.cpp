#include <cstdlib>
#include <iostream>
#include <iomanip> 
#include <cmath>
#include <sstream>
#include <ctime>
#include <array>
#include <algorithm>
#include <limits>
#include <numeric>

#include "vlsv_writer.h"
#include "vlsv_amr.h"

using namespace vlsv;

// use the buffered API or not
#define BUFFERED

std::vector<int> sizes = {540,554,540,21,65,56,53,21,64,24};
std::vector<int> sizes2 = {1,2,10,100,1000,10000};
std::vector<int> parts = {1,3,2,1,1,1,10,1,1,1};
std::vector<int> parts2 = {1,2,5,10};


int *intData;
int myrank;


void writeNint(vlsv::Writer &vlsv, int N, int writeParts, int vectorSize)
{
	std::cout << myrank << " added " << std::endl;
	std::map<std::string, std::string> xmlAttributes;
	xmlAttributes["Arr"] = "name";

	vlsv.startMultiwrite(getStringDatatype<int>(),N*writeParts,vectorSize,sizeof(int));
	for (size_t i = 0; i < writeParts; i++)
	{
		vlsv.addMultiwriteUnit((char*)intData,N);
	}
	vlsv.endMultiwrite("arrayName"+std::to_string(N),xmlAttributes);
}

void setBuffers(vlsv::Writer &vlsv, int bufferSize)
{
	bool same = false;
#ifdef BUFFERED
	if(!same)
	{
		if(myrank == 1)	   
			vlsv.setBuffer(bufferSize);
		else
			vlsv.setBuffer(bufferSize/10);
	}
	else
		vlsv.setBuffer(bufferSize);
#endif
}

int main(int argc,char* argv[]) {
   	bool success = true;

	if (argc == 1)
	{
		std::cout << "usage : a.out buffer_size"
		std::cout << "for MPI with up to 10 processes: mpiexec -n 10 a.out buffer_size"
	}

   // Init MPI:
   	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int bufferSize = std::atoi(argv[1]);

	int totSize = 10000;
	intData = new int[totSize];

	for (size_t i = 0; i < totSize; i++)
	{
		intData[i] = myrank;
	}

	
	{
		// all writing the same size messages
		vlsv::Writer vlsv;
		setBuffers(vlsv, bufferSize);
	   	if (vlsv.open("allTheSame.out"+std::to_string(bufferSize),MPI_COMM_WORLD,0) == false) {
			  success = false;
		 	 MPI_Finalize();	
		 	 return 1;
		}
		for(auto size: sizes2)
		for(auto part: parts)
			writeNint(vlsv, size, part, 1);

		if (vlsv.close() == false) {
			success = false;
	   	}
	}
	
	
	{
		// all writing different size messages
		vlsv::Writer vlsv;

		setBuffers(vlsv, bufferSize);
	   	if (vlsv.open("allDifferent.out"+std::to_string(bufferSize),MPI_COMM_WORLD,0) == false) {
			  success = false;
		 	 MPI_Finalize();	
		 	 return 1;
		}

		writeNint(vlsv, sizes[myrank], parts[myrank], 1);

		if (vlsv.close() == false) {
			success = false;
	   	}
	}
	
	{
		// all but one writing
		vlsv::Writer vlsv;

		setBuffers(vlsv, bufferSize);
	   	if (vlsv.open("allBut1.out"+std::to_string(bufferSize),MPI_COMM_WORLD,0) == false) {
			  success = false;
		 	 MPI_Finalize();	
		 	 return 1;
		}

		if(myrank == 4)
			writeNint(vlsv, 0, 1, 1);	
		else
			writeNint(vlsv, sizes[myrank], parts[myrank], 1);


		if (vlsv.close() == false) {
			success = false;
	   	}
	}
	
	{
		// only one writing 
		vlsv::Writer vlsv;

		setBuffers(vlsv, bufferSize);
	   	if (vlsv.open("only1.out"+std::to_string(bufferSize),MPI_COMM_WORLD,0) == false) {
			  success = false;
		 	 MPI_Finalize();	
		 	 return 1;
		}

		if(myrank == 5)
			writeNint(vlsv, sizes[myrank], 1, 1);
		else
			writeNint(vlsv, 0, 0, 1);	
		
		if (vlsv.close() == false) {
			success = false;
	   	}
	}

	   MPI_Finalize();
 	  if (success == false) return 1;
	return 0;
}
