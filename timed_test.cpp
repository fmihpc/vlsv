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
#include <omp.h>

#include "vlsv_writer.h"
#include "vlsv_reader_parallel.h"
#include "vlsv_amr.h"

using namespace vlsv;

std::vector<int> chunkSizes;
std::vector<int> chunkCounts;


std::vector<int> intData_w, intData_r;
int myRank, nRanks;


void writeNint(vlsv::Writer &vlsv, int chunkSize, int chunkCount)
{
	std::map<std::string, std::string> xmlAttributes;
	xmlAttributes["Arr"] = "name";

	vlsv.startMultiwrite(getStringDatatype<int>(),chunkSize*chunkCount,1 /*vectorSize */,sizeof(int));
	for (int i = 0; i < chunkCount; i++)
	{
		vlsv.addMultiwriteUnit(intData_w.data(),chunkSize);
	}
	vlsv.endMultiwrite("arrayName", xmlAttributes);
}

void readNint(vlsv::ParallelReader &vlsv, int chunkSize, int chunkCount, uint64_t fileOffset)
{
	std::list<std::pair<std::string,std::string> > xmlAttributes;
	xmlAttributes.push_back({"Arr",  "name"});

	vlsv.startMultiread("arrayName", xmlAttributes);
	for (int i = 0; i < chunkCount; i++)
	{
		vlsv.addMultireadUnit((char*)(intData_r.data() + chunkSize*i), chunkSize);
	}
	vlsv.endMultiread(fileOffset);
}

// Every element of intData_r should equal myRank, since writeNint() fills
// intData_w with myRank before writing it out.
bool checkNint(int myRank)
{
	size_t errorCount = 0;
	const size_t n = intData_r.size();
	#pragma omp parallel for reduction(+:errorCount)
	for (size_t i = 0; i < n; i++) {
		errorCount += (intData_r[i] != myRank);
	}
	if (errorCount > 0) {
		std::cerr << "ERROR: Rank " << myRank << " read " << errorCount << "/" << intData_r.size()
			<< " incorrect values (expected all == " << myRank << ")" << std::endl;
	}
	return errorCount == 0;
}


int main(int argc,char* argv[]) {
   	bool success = true;

	if (argc < 7)
	{
		std::cout << "usage : srun timed_test R|W|RW buffer_size min_rank_chunk_size max_rank_chunk_size min_rank_chunk_count max_rank_chunk_count [WR mpiio_hint mpiio_value] [RD mpiio_hint mpiio_value] [CHECK]" << std::endl;
		exit(42);
	}

        // Init MPI:
   	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

	// If present, CHECK must be the last argument; it enables a correctness
	// check of the data read back in (reduces argc so the rest of the
	// positional/MPI-hint parsing below is unaffected).
	bool checkData = false;
	if (argv[argc-1] == std::string("CHECK")) {
		checkData = true;
		argc--;
	}
	
	const double GiB = 1024*1024*1024;
	
	enum mode {R, W, WR};
	enum mode rwmode;
	if(argv[1] == std::string("R")) {
		rwmode = R;
	} else if (argv[1] == std::string("W")) {
		rwmode = W;
	} else if (argv[1] == std::string("WR") || argv[1] == std::string("RW")) {
		rwmode = WR;
	} else {
		std::cout << "The first argument must be R, W, or RW to only read, only write, or read and write in the same pass." << std::endl;
		exit(42);
	}
	int bufferSize = std::atoi(argv[2]);
        int minChunkSize = std::atoi(argv[3]);
        int maxChunkSize = std::atoi(argv[4]);
        int minChunkCount = std::atoi(argv[5]);
        int maxChunkCount = std::atoi(argv[6]);

	MPI_Info MPIinfo_wr, MPIinfo_rd;
	if (argc == 7) {
 		MPIinfo_wr = MPIinfo_rd = MPI_INFO_NULL;
	} else {
		int i=7, wr_count=0, rd_count=0;
		MPI_Info_create(&MPIinfo_wr);
		MPI_Info_create(&MPIinfo_rd);
		MPI_Info * MPIinfo_ptr;
		while(true) {
			if(i >= argc) {
				break;
			}
			if(argv[i] == std::string("WR")) {
				MPIinfo_ptr = &MPIinfo_wr;
				i++;
				wr_count++;
			}
			if(argv[i] == std::string("RD")) {
				MPIinfo_ptr = &MPIinfo_rd;
				i++;
				rd_count++;
			}
			MPI_Info_set(*MPIinfo_ptr, argv[i], argv[i+1]);
			i += 2;
		}
		if(wr_count == 0) {
			MPIinfo_wr = MPI_INFO_NULL;
		}
		if(rd_count == 0) {
			MPIinfo_rd = MPI_INFO_NULL;
		}
	}
	if(myRank == 0) {
		std::string action;
		switch(rwmode) {
			case R:
				action = std::string("only read");
				break;
			case W:
				action = std::string("only write");
				break;
			case WR:
				action = std::string("write, then read");
				break;
			default:
				break;
		}
		std::cout << "INFO: This will " << action << " int data, int has size " << sizeof(int) << std::endl;
		std::cout << "INFO: This will " << action << " between " << minChunkCount << " and " << maxChunkCount << " chunks of between " << minChunkSize*sizeof(int) / GiB << " and " << maxChunkSize*sizeof(int) / GiB << " GiB on each of the " << nRanks << " tasks." << std::endl;
	}
	
	std::default_random_engine generatorChunkSizes(42), generatorChunkCounts(43);
	std::uniform_int_distribution<int> randomChunkSizes(minChunkSize, maxChunkSize), randomChunkCounts(minChunkCount, maxChunkCount);
	for (int i=0; i<nRanks; i++) {
		chunkSizes.push_back(randomChunkSizes(generatorChunkSizes));
		chunkCounts.push_back(randomChunkCounts(generatorChunkCounts));
	}

	const double chunkSize = chunkSizes[myRank]*sizeof(int) / GiB;
	uint64_t totalChunkCount = 0;
	double totalSize = 0;
	
	std::stringstream stream;

	if(myRank == 0) {
		for (int i=0; i<nRanks; i++) {
			totalChunkCount += chunkCounts[i];
			totalSize += chunkCounts[i]*chunkSizes[i];
		}
		totalSize *= sizeof(int) / GiB;
		
		stream << "# RD/WR \t Rank \t chunks \t size (GiB) \t total (GiB) \t time (s) \t speed (GiB/s)" << std::endl;
		std::cerr << stream.str();
		stream.clear();
		stream.str(std::string());
	}
	
	// WRITE
	if(rwmode != R) {
		// For writing, we need only one chunk that'll be written a number of times
		intData_w.assign(chunkSizes[myRank], myRank);
		
		vlsv::Writer vlsvWriter;
		vlsvWriter.setBuffer(bufferSize);
		if (vlsvWriter.open("file.out",MPI_COMM_WORLD,0, MPIinfo_wr) == false) {
			success = false;
			MPI_Finalize();
			return 1;
		}
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		double tStart = MPI_Wtime();
		writeNint(vlsvWriter, chunkSizes[myRank], chunkCounts[myRank]);
		double tTime = MPI_Wtime() - tStart;
		
		stream << "WR\t" << myRank << "\t" << chunkCounts[myRank] << "\t" << chunkSize << "\t" << chunkCounts[myRank]*chunkSize << "\t" << tTime << "\t" << chunkCounts[myRank]*chunkSize / tTime << std::endl;
		std::cerr << stream.str();
		stream.clear();
		stream.str(std::string());
		
		if (vlsvWriter.close() == false) {
			success = false;
		}
		
		// Do this here as the close call does the actual writing when buffering
		MPI_Barrier(MPI_COMM_WORLD);
		if(myRank == 0) {
			const double totalTime = MPI_Wtime() - tStart;
			stream << "WR\t" << nRanks << "\t" << totalChunkCount << "\t" << totalSize / totalChunkCount << "\t" << totalSize << "\t" << totalTime << "\t" << totalSize / totalTime << std::endl;
			std::cerr << stream.str();
			stream.clear();
			stream.str(std::string());
		}
	}
	
	// READ
	if(rwmode != W) {
		// For reading we want the buffer to have the full size
		intData_r.resize(chunkSizes[myRank]*chunkCounts[myRank]);

		vlsv::ParallelReader vlsvReader;
		if (vlsvReader.open("file.out",MPI_COMM_WORLD,0, MPIinfo_rd) == false) {
			success = false;
			MPI_Finalize();
			return 1;
		}
		
		uint64_t myFileOffset = 0;
		for(int t=0; t<myRank; t++) {
			myFileOffset += chunkSizes[t]*chunkCounts[t];
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		double tStart = MPI_Wtime();
		readNint(vlsvReader, chunkSizes[myRank], chunkCounts[myRank], myFileOffset);
		double tTime = MPI_Wtime() - tStart;

		stream << "RD\t" << myRank << "\t" << chunkCounts[myRank] << "\t" << chunkSize << "\t" << chunkCounts[myRank]*chunkSize << "\t" << tTime << "\t" << chunkCounts[myRank]*chunkSize / tTime << std::endl;
		std::cerr << stream.str();
		stream.clear();
		stream.str(std::string());

		if (checkData) {
			bool localOk = checkNint(myRank);
			int localOkInt = localOk ? 1 : 0;
			int globalOkInt = 0;
			MPI_Reduce(&localOkInt, &globalOkInt, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
			if (!localOk) {
				success = false;
			}
			if (myRank == 0) {
				std::cout << "INFO: Data correctness check " << (globalOkInt ? "PASSED" : "FAILED") << std::endl;
			}
		}

		if (vlsvReader.close() == false) {
			success = false;
		}

		MPI_Barrier(MPI_COMM_WORLD);
		if(myRank == 0) {
			const double totalTime = MPI_Wtime() - tStart;
			stream << "RD\t" << nRanks << "\t" << totalChunkCount << "\t" << totalSize / totalChunkCount << "\t" << totalSize << "\t" << totalTime << "\t" << totalSize / totalTime << std::endl;
			std::cerr << stream.str();
			stream.clear();
			stream.str(std::string());
		}
	}

	MPI_Finalize();
 	if (success == false) return 1;
	return 0;
}
