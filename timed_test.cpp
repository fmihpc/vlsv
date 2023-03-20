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
#include "vlsv_reader_parallel.h"
#include "vlsv_amr.h"

using namespace vlsv;

std::vector<int> chunkSizes;
std::vector<int> chunkCounts;


std::vector<int> intData;
int myRank, nRanks;


void writeNint(vlsv::Writer &vlsv, int chunkSize, int chunkCount)
{
	std::map<std::string, std::string> xmlAttributes;
	xmlAttributes["Arr"] = "name";

	vlsv.startMultiwrite(getStringDatatype<int>(),chunkSize*chunkCount,1 /*vectorSize */,sizeof(int));
	for (int i = 0; i < chunkCount; i++)
	{
		vlsv.addMultiwriteUnit(intData.data(),chunkSize);
	}
	vlsv.endMultiwrite("arrayName"+std::to_string(chunkSize),xmlAttributes);
}

void readNint(vlsv::ParallelReader &vlsv, int chunkSize, int chunkCount, uint64_t fileOffset)
{
	std::list<std::pair<std::string,std::string> > xmlAttributes;
	xmlAttributes.push_back({"Arr",  "name"});

	vlsv.startMultiread("Arr", xmlAttributes);
	for (int i = 0; i < chunkCount; i++)
	{
		vlsv.addMultireadUnit((char*)intData.data() + chunkSize*i*sizeof(int), chunkSize);
	}
	vlsv.endMultiread(fileOffset);
}


int main(int argc,char* argv[]) {
   	bool success = true;

	if (argc < 6 || (argc >= 6 && argc % 2 != 0))
	{
		std::cout << "usage : srun timed_test buffer_size min_rank_chunk_size max_rank_chunk_size min_rank_chunk_count max_rank_chunk_count [mpiio_hint mpiio_value]" << std::endl;
	}

        // Init MPI:
   	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
	
	const double GiB = 1024*1024*1024;
	
	int bufferSize = std::atoi(argv[1]);
        int minChunkSize = std::atoi(argv[2]);
        int maxChunkSize = std::atoi(argv[3]);
        int minChunkCount = std::atoi(argv[4]);
        int maxChunkCount = std::atoi(argv[5]);

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
		std::cout << "INFO: This will write between " << minChunkCount << " and " << maxChunkCount << " chunks of between " << minChunkSize*sizeof(int) / GiB << " and " << maxChunkSize*sizeof(int) / GiB << " GiB on each of the " << nRanks << " tasks." << std::endl;
	}
	
	std::default_random_engine generatorChunkSizes(42), generatorChunkCounts(43);
	std::uniform_int_distribution<int> randomChunkSizes(minChunkSize, maxChunkSize), randomChunkCounts(minChunkCount, maxChunkCount);
	for (int i=0; i<nRanks; i++) {
		chunkSizes.push_back(randomChunkSizes(generatorChunkSizes));
		chunkCounts.push_back(randomChunkCounts(generatorChunkCounts));
	}
	
	if(myRank == 0) {
		std::stringstream stream;
		stream << "# RD/WR \t Rank \t chunks \t size (GiB) \t total (GiB) \t time (s) \t speed (GiB/s)" << std::endl;
		std::cerr << stream.str();
	}
	
	
	// WRITE
	
	// For writing, we need only one chunk that'll be written a number of times
	intData.assign(chunkSizes[myRank], myRank);
	
	vlsv::Writer vlsvWriter;
	vlsvWriter.setBuffer(bufferSize);
	if (vlsvWriter.open("file.out",MPI_COMM_WORLD,0, MPIinfo) == false) {
		success = false;
		MPI_Finalize();
		return 1;
	}
	
	const double chunkSize = chunkSizes[myRank]*sizeof(int) / GiB;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double tStart = MPI_Wtime();
	writeNint(vlsvWriter, chunkSizes[myRank], chunkCounts[myRank]);
	double tTime = MPI_Wtime() - tStart;
	
	std::stringstream stream;
        stream << "WR\t" << myRank << "\t" << chunkCounts[myRank] << "\t" << chunkSize << "\t" << chunkCounts[myRank]*chunkSize << "\t" << tTime << "\t" << chunkCounts[myRank]*chunkSize / tTime << std::endl;
        std::cerr << stream.str();
	
	if (vlsvWriter.close() == false) {
		success = false;
	}
	
	// Do this here as the close call does the actual writing when buffering
	MPI_Barrier(MPI_COMM_WORLD);
	uint64_t totalChunkCount = 0;
	double totalSize = 0;
	if(myRank == 0) {
		const double totalTime = MPI_Wtime() - tStart;
		for (int i=0; i<nRanks; i++) {
			totalChunkCount += chunkCounts[i];
			totalSize += chunkCounts[i]*chunkSizes[i];
		}
		totalSize *= sizeof(int) / GiB;
		std::stringstream stream;
		stream << "WR\t" << nRanks << "\t" << totalChunkCount << "\t" << totalSize / totalChunkCount << "\t" << totalSize << "\t" << totalTime << "\t" << totalSize / totalTime << std::endl;
		std::cerr << stream.str();
	}
	
	
	// READ
	
	// For reading we want the buffer to have the full size
	intData.resize(chunkSizes[myRank]*chunkCounts[myRank]);
	
	vlsv::ParallelReader vlsvReader;
	if (vlsvReader.open("file.out",MPI_COMM_WORLD,0, MPIinfo) == false) {
		success = false;
		MPI_Finalize();
		return 1;
	}
	
	uint64_t myFileOffset = 0;
	for(int t=0; t<myRank; t++) {
		myFileOffset += chunkSizes[t]*chunkCounts[t]*sizeof(int);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	tStart = MPI_Wtime();
	readNint(vlsvReader, chunkSizes[myRank], chunkCounts[myRank], myFileOffset);
	tTime = MPI_Wtime() - tStart;
	
	stream.clear();
	stream.str(std::string());
        stream << "RD\t" << myRank << "\t" << chunkCounts[myRank] << "\t" << chunkSize << "\t" << chunkCounts[myRank]*chunkSize << "\t" << tTime << "\t" << chunkCounts[myRank]*chunkSize / tTime << std::endl;
        std::cerr << stream.str();
	
	if (vlsvReader.close() == false) {
		success = false;
	}
	stream.clear();
        stream.str(std::string());
	for (auto i : intData) {
		stream << myRank << " " << intData.size() << " " << intData[i] << std::endl;
		std::cerr << stream.str();
		stream.clear();
		stream.str(std::string());
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(myRank == 0) {
		const double totalTime = MPI_Wtime() - tStart;
		std::stringstream stream;
		stream << "RD\t" << nRanks << "\t" << totalChunkCount << "\t" << totalSize / totalChunkCount << "\t" << totalSize << "\t" << totalTime << "\t" << totalSize / totalTime << std::endl;
		std::cerr << stream.str();
	}
	
	MPI_Finalize();
 	if (success == false) return 1;
	return 0;
}
