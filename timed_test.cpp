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
#include <random>

#include "vlsv_writer.h"
#include "vlsv_amr.h"

using namespace vlsv;

std::vector<int> sizes;
std::vector<int> parts;


int *intData;
int myRank, nRanks;


void writeNint(vlsv::Writer &vlsv, int N, int writeParts, int vectorSize)
{
	std::map<std::string, std::string> xmlAttributes;
	xmlAttributes["Arr"] = "name";

	vlsv.startMultiwrite(getStringDatatype<int>(),N*writeParts,vectorSize,sizeof(int));
	for (int i = 0; i < writeParts; i++)
	{
		vlsv.addMultiwriteUnit((char*)intData,N);
	}
	vlsv.endMultiwrite("arrayName"+std::to_string(N),xmlAttributes);
}


int main(int argc,char* argv[]) {
   	bool success = true;

	if (argc < 6 || (argc >= 6 && argc % 2 != 0))
	{
		std::cout << "usage : srun timed_test buffer_size min_rank_part_size max_rank_part_size min_rank_part_count max_rank_part_count [mpiio_hint mpiio_value]" << std::endl;
	}

        // Init MPI:
   	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
	
	const double GiB = 1024*1024*1024;
	
	int bufferSize = std::atoi(argv[1]);
        int minPartSize = std::atoi(argv[2]);
        int maxPartSize = std::atoi(argv[3]);
        int minPartCount = std::atoi(argv[4]);
        int maxPartCount = std::atoi(argv[5]);

	MPI_Info MPIinfo;
	if (argc == 6) {
 		MPIinfo = MPI_INFO_NULL;
	} else {
		MPI_Info_create(&MPIinfo);
		for (int i=6; i<argc; i+=2)
		{
			MPI_Info_set(MPIinfo, argv[i], argv[i+1]);
		}
	}
	if(myRank == 0) {
		std::cout << "INFO: This will write int data, int has size " << sizeof(int) << std::endl;
		std::cout << "INFO: This will write between " << minPartCount << " and " << maxPartCount << " parts of between " << minPartSize*sizeof(int) / GiB << " and " << maxPartSize*sizeof(int) / GiB << " GiB on each of the " << nRanks << " tasks." << std::endl;
	}
	
	intData = new int[maxPartSize];

	for (int i = 0; i < maxPartSize; i++)
	{
		intData[i] = myRank;
	}

	std::default_random_engine generatorSizes(42), generatorParts(43);
	std::uniform_int_distribution<int> randomSizes(minPartSize, maxPartSize), randomParts(minPartCount, maxPartCount);
	for (int i=0; i<nRanks; i++) {
		sizes.push_back(randomSizes(generatorSizes));
		parts.push_back(randomParts(generatorParts));
	}
	
	if(myRank == 0) {
		std::stringstream stream;
		stream << "# Rank \t parts \t size (GiB) \t total (GiB) \t time (s) \t speed (GiB/s)" << std::endl;
		std::cerr << stream.str();
	}
	
	vlsv::Writer vlsv;
	vlsv.setBuffer(bufferSize);
	if (vlsv.open("file.out",MPI_COMM_WORLD,0, MPIinfo) == false) {
		success = false;
		MPI_Finalize();
		return 1;
	}
	
	const double partSize = sizes[myRank]*sizeof(int) / GiB;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	const double tStart = MPI_Wtime();
	writeNint(vlsv, sizes[myRank], parts[myRank], 1);
	const double tTime = MPI_Wtime() - tStart;
	
	std::stringstream stream;
        stream << myRank << "\t" << parts[myRank] << "\t" << partSize << "\t" << parts[myRank]*partSize << "\t" << tTime << "\t" << parts[myRank]*partSize / tTime << std::endl;
        std::cerr << stream.str();
	
	if (vlsv.close() == false) {
		success = false;
	}
	
	// Do this here as the close call does the actual writing when buffering
	MPI_Barrier(MPI_COMM_WORLD);
	if(myRank == 0) {
		const double totalTime = MPI_Wtime() - tStart;
		uint64_t totalParts = 0;
		double totalSize = 0;
		for (int i=0; i<nRanks; i++) {
			totalParts += parts[i];
			totalSize += parts[i]*sizes[i];
		}
		totalSize *= sizeof(int) / GiB;
		std::stringstream stream;
		stream << nRanks << "\t" << totalParts << "\t" << totalSize / totalParts << "\t" << totalSize << "\t" << totalTime << "\t" << totalSize / totalTime << std::endl;
		std::cerr << stream.str();
	}
	
	MPI_Finalize();
 	if (success == false) return 1;
	return 0;
}
