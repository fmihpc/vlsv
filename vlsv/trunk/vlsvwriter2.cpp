#include <cstdlib>
#include <iostream>
#include <fstream>

#include "mpiconversion.h"
#include "vlsvwriter.h"

using namespace std;

/** Constructor for WriteUnit. Simply copies parameters to internal variables.
 * @param array Pointer to the data to be written.
 * @param mpiType MPI datatype describing the type of data to be written.
 * @param amount Number of MPI datatypes in array.
 */
VLSV::WriteUnit::WriteUnit(char* array,const MPI_Datatype& mpiType,const uint64_t& amount): array(array),mpiType(mpiType),amount(amount) { }

/** Constructor for VLSVWriter.
 * In multithreaded mode the mutexes and condition variables are 
 * initialized here.
 * In multithreaded mode only the master thread should be allowed to 
 * call the constructor.
 */
VLSVWriter::VLSVWriter() {
   blockLengths = NULL;
   bytesPerProcess = NULL;
   displacements = NULL;
   endMultiwriteCounter = 0;
   fileOpen = false;
   initialized = false;
   multiwriteFinalized = false;
   multiwriteInitialized = false;
   N_multiwriteUnits = 0;
   offset = 0;
   offsets = NULL;
   types = NULL;
   xmlWriter = NULL;
}

/** Destructor for VLSVWriter. Deallocates XML writer.
 * In multithreaded mode mutexes and condition variables are 
 * destroyed here.
 * In multithreaded mode only the master thread should be allowed to 
 * call the destructor.
 */
VLSVWriter::~VLSVWriter() {
   if (blockLengths != NULL) {cerr << "bl " << blockLengths;}
   if (bytesPerProcess != NULL) {cerr << "bl " << bytesPerProcess;}
   if (displacements != NULL) {cerr << "bl " << displacements;}
   if (offsets != NULL) {cerr << "bl " << offsets;}
   if (types != NULL) {cerr << "bl " << types;}
   if (xmlWriter != NULL) {cerr << "bl " << xmlWriter;}
   
   delete [] blockLengths; blockLengths = NULL;
   delete [] bytesPerProcess; bytesPerProcess = NULL;
   delete [] displacements; displacements = NULL;
   delete [] offsets; offsets = NULL;
   delete [] types; types = NULL;
   delete xmlWriter; xmlWriter = NULL;
}

bool VLSVWriter::addMultiwriteUnit(char* array,const uint64_t& arrayElements,const int& threadID) {
   // Check that startMultiwrite has initialized correctly:
   if (multiwriteInitialized == false) return false;

   // Each thread records their multiwrite units to per-thread storage,
   // so there is no need to synchronize access to vector multiwriteUnits:
   multiwriteUnits[threadID].push_back(VLSV::WriteUnit(array,getMPIDatatype(vlsvType,dataSize),arrayElements*vectorSize));
   return true;
}

/** Close a file that has been previously opened by calling VLSVWriter::open.
 * After the file has been closed the MPI master process appends an XML footer 
 * to the end of the file, and writes an offset to the footer to the start of 
 * the file.
 * In multithreaded mode only the master thread is allowed to 
 * participate in the file closing process, all other threads will 
 * block until the master thread has finished.
 * In multithreaded mode it is safe to call this function from all processes.
 * @param threadID Thread ID of the thread calling this function. If multithreaded 
 * mode is not used the master thread ID should be used here. Defaults to 
 * value zero.
 * @return If true, the file was closed successfully. If false, a file may not 
 * have been opened successfully by VLSVWriter::open.
 */
bool VLSVWriter::close(const int& threadID) {
   // If a file was never opened, exit immediately:
   if (fileOpen == false) {
      return false;
   }
   
   // Close MPI file:
   MPI_Barrier(comm);
   MPI_File_close(&fileptr);
   
   // Master process appends footer to the end of binary file:
   if (myrank == masterRank) {
      fstream footer;
      footer.open(fileName.c_str(),fstream::out | fstream::app);
      footer.seekp(0,std::ios_base::end);
      
      // Get put position
      uint64_t footerOffset = footer.tellp();
        
      // Write footer:
      xmlWriter->print(footer);
      delete xmlWriter; xmlWriter = NULL;
      footer.close();
      
      // Write header position to the beginning of binary file:
      footer.open(fileName.c_str(),fstream::in | fstream::out | fstream::binary | fstream::ate);
      char* ptr = reinterpret_cast<char*>(&footerOffset);
      footer.seekp(sizeof(uint64_t));
      footer.write(ptr,sizeof(uint64_t));
      footer.close();
   }

   delete [] offsets; offsets = NULL;
   delete [] bytesPerProcess; bytesPerProcess = NULL;

   // Wait for master process to finish:
   MPI_Barrier(comm);
   
   // Wake up other threads (multithreaded mode only):
   fileOpen = false;
   return true;
}

MPI_Datatype VLSVWriter::getMPIDatatype(VLSV::datatype dt,uint64_t dataSize) {
   switch (dt) {
    case VLSV::INT:
      switch (dataSize) {
       case (sizeof(int8_t)):
	 return MPI_Type<int8_t>();
	 break;
       case (sizeof(int16_t)):
	 return MPI_Type<int16_t>();
	 break;
       case (sizeof(int32_t)):
	 return MPI_Type<int32_t>();
	 break;
       case (sizeof(int64_t)):
	 return MPI_Type<int64_t>();
	 break;
      }
    case VLSV::UINT:
      switch (dataSize) {
       case (sizeof(uint8_t)):
	 return MPI_Type<uint8_t>();
	 break;
       case (sizeof(uint16_t)):
	 return MPI_Type<uint16_t>();
	 break;
       case (sizeof(uint32_t)):
	 return MPI_Type<uint32_t>();
	 break;
       case (sizeof(uint64_t)):
	 return MPI_Type<uint64_t>();
	 break;
      }
    case VLSV::FLOAT:
      switch (dataSize) {
       case (sizeof(float)):
	 return MPI_Type<float>();
	 break;
       case (sizeof(double)):
	 return MPI_Type<double>();
	 break;
       case (sizeof(long double)):
	 return MPI_Type<long double>();
	 break;
      }
   }
}

VLSV::datatype VLSVWriter::getVLSVDatatype(const string& s) {
   if (s == "int") return VLSV::INT;
   else if (s == "uint") return VLSV::UINT;
   else if (s == "float") return VLSV::FLOAT;
   else return VLSV::UNKNOWN;
}

/** Open a VLSV file for parallel output. The file is opened on all processes 
 * in the given communicator. Additionally, master MPI process writes a 
 * header into the output file and caches a footer which will be written 
 * in VLSVWriter::close.
 * In multithreaded mode only the master thread will participate in the file 
 * opening process, all other threads will simply block until master thread 
 * has finished.
 * In multithreaded mode it is safe to call this function with all threads.
 * @param fname The name of the output file.
 * @param comm MPI communicator used in writing.
 * @param masterProcessID ID of the MPI master process.
 * @param mpiThreadingLevel Threading level of MPI, obtained from MPI_Init_thread. 
 * If multithreaded mode is not used a value MPI_THREAD_SINGLE should be used here.
 * Defaults to MPI_THREAD_SINGLE.
 * @param N_threads Total number of threads using VLSV writer, this value is 
 * required for allocation of per-thread variables.
 * @param threadID Thread ID of the thread calling this function. If multithreaded 
 * mode is not used a value zero should be used here. Defaults to value zero.
 * @param masterThreadID Thread ID of the master thread. If multithreaded mode is 
 * not used the master thread ID should be used here. Defalts to value zero.
 * @return If true, a file was opened successfully.
 */
bool VLSVWriter::open(const std::string& fname,MPI_Comm comm,const int& masterProcessID,
		      const int& mpiThreadingLevel,const int& N_threads,const int& threadID,const int& masterThreadID) {
   MPI_Comm_dup(comm,&(this->comm));
   masterRank              = masterProcessID;
   this->masterThreadID    = masterThreadID;
   this->mpiThreadingLevel = mpiThreadingLevel;
   this->N_threads         = N_threads;
   MPI_Comm_rank(this->comm,&myrank);
   MPI_Comm_size(this->comm,&N_processes);

   // Allocate per-thread storage:
   multiwriteOffsets.resize(N_threads);
   multiwriteUnits.resize(N_threads);
   
   // All processes in communicator comm open the same file. If a file with the 
   // given name already exists it is deleted:
   int accessMode = (MPI_MODE_WRONLY | MPI_MODE_CREATE);
   MPI_Info MPIinfo = MPI_INFO_NULL;
   fileName = fname;
   MPI_File_delete(const_cast<char*>(fname.c_str()),MPI_INFO_NULL);
   if (MPI_File_open(comm,const_cast<char*>(fileName.c_str()),accessMode,MPIinfo,&fileptr) != MPI_SUCCESS) {
      return fileOpen;
   }   
   
   offset = 0;           //offset set to 0 when opening a new file
   MPI_File_set_view(fileptr,0,MPI_BYTE,MPI_BYTE,const_cast<char*>("native"),MPI_INFO_NULL);
       
   // Only master rank needs these arrays:
   if (myrank == masterRank) {
      offsets = new MPI_Offset[N_processes];
      bytesPerProcess = new uint64_t[N_processes];    
   }
    
   // Master opens an XML tree for storing the footer:
   if (myrank == masterRank) {
      xmlWriter     = new MuXML();
      XMLNode* root = xmlWriter->getRoot();
      xmlWriter->addNode(root,"VLSV","");
   }
   
   // Master writes 2 64bit integers to the start of file. 
   // Second value will be overwritten in close() function to tell 
   // the position of footer:
   bool success = true;
   if (myrank == masterRank) {
      // Write file endianness to the first byte:
      uint64_t endianness = 0;
      unsigned char* ptr = reinterpret_cast<unsigned char*>(&endianness);
      ptr[0] = detectEndianness();
      if (MPI_File_write_at(fileptr,0,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
      if (MPI_File_write_at(fileptr,8,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;	
      offset += 2*sizeof(uint64_t); //only master rank keeps a running count
   }

   // Check that everything is OK, if not then we need to close the file here. 
   // Master process needs to broadcast the status to all other processes:
   MPI_Bcast(&success,sizeof(bool),MPI_BYTE,masterRank,comm);   
   if (success == false) {
      MPI_File_close(&fileptr);
      MPI_File_delete(const_cast<char*>(fileName.c_str()),MPI_INFO_NULL);
      delete [] offsets; offsets = NULL;
      delete [] bytesPerProcess; bytesPerProcess = NULL;
   }

   initialized = true;
   fileOpen    = success;
   return fileOpen;
}

bool VLSVWriter::startMultiwrite(const string& datatype,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,const int& threadID) {
   // Clear per-thread storage:
   multiwriteUnits[threadID].clear();
   multiwriteOffsets[threadID] = std::numeric_limits<unsigned int>::max();
   
   // Array datatype and byte size of each vector element are determined 
   // from the template parameter, other values are copied from parameters:
   this->dataType   = datatype;
   this->vlsvType   = getVLSVDatatype(datatype);
   this->arraySize  = arraySize;
   this->vectorSize = vectorSize;
   this->dataSize   = dataSize;
   multiwriteFinalized = false;
   N_multiwriteUnits = 0;
   endMultiwriteCounter = 0;
   
   // Gather the number of bytes written by every process to MPI master process:
   myBytes = arraySize * vectorSize * dataSize;
   MPI_Gather(&myBytes,1,MPI_Type<uint64_t>(),bytesPerProcess,1,MPI_Type<uint64_t>(),masterRank,comm);
   
   // MPI master process calculates an offset to the output file for all processes:
   if (myrank == masterRank) {
      offsets[0] = offset;
      for (int i=1; i<N_processes; ++i) offsets[i] = offsets[i-1] + bytesPerProcess[i-1];
   }
   
   // MPI master scatters offsets:
   MPI_Scatter(offsets,1,MPI_Type<uint64_t>(),&offset,1,MPI_Type<uint64_t>(),masterRank,comm);

   multiwriteInitialized = true;
   return multiwriteInitialized;
}

/** Write multiwrite units to file.
 * In multithreaded mode it is safe to call this function with every thread. This function contains several 
 * thread synchronization points, in particular MPI calls have been restricted to the master thread.
 * TODO: In MPI_THREAD_SERIALIZED and MPI_THREAD_MULTIPLE modes it might be more efficient to call 
 * MPI functions with all threads.
 * @param tagName Name of the XML tag for this array.
 * @param attribs Attributes for the XML tag.
 * @param threadID Thread ID of the thread calling this function. If multithreading mode is not 
 * used, masterThreadID value should be used here. Defaults to zero value.
 * @return If true, array was successfully written to file. In multithreaded mode all threads 
 * calling this function will return the same value.
 */
bool VLSVWriter::endMultiwrite(const std::string& tagName,const std::map<std::string,std::string>& attribs,const int& threadID) {
   // Master thread allocates memory for an MPI_Struct that is used to 
   // write all multiwrite units with a single collective call. In 
   // multithreaded mode other threads wait until master is finished:
   if (threadID == masterThreadID) {
      // Count the total number of multiwrite units:
      N_multiwriteUnits = 0;
      for (size_t i=0; i<multiwriteUnits.size(); ++i) {
	 N_multiwriteUnits += multiwriteUnits[i].size();
      }
      // Calculate offset for each thread:
      multiwriteOffsets[0] = 0;
      for (size_t i=1; i<multiwriteUnits.size(); ++i) {
	 multiwriteOffsets[i] = multiwriteOffsets[i-1] + multiwriteUnits[i].size();
      }
      // Allocate memory for an MPI struct:
      blockLengths  = new int[N_multiwriteUnits];
      displacements = new MPI_Aint[N_multiwriteUnits];
      types         = new MPI_Datatype[N_multiwriteUnits];
      // Calculate a global offset pointer for MPI struct, i.e. an 
      // offset which is used to calculate the displacements:
      bool found = false;
      multiwriteOffsetPointer = NULL;
      for (size_t i=0; i<multiwriteUnits.size(); ++i) {
	 if (multiwriteUnits[i].size() == 0) continue;
	 multiwriteOffsetPointer = multiwriteUnits[i].begin()->array;
      }
   }
   
   // Every thread copies its multiwrite units into a global MPI struct:
   if (multiwriteUnits[threadID].size() > 0) {
      unsigned int offset = multiwriteOffsets[threadID];
      unsigned int i = 0;
      for (list<VLSV::WriteUnit>::const_iterator it=multiwriteUnits[threadID].begin(); it!=multiwriteUnits[threadID].end(); ++it) {
	 //cerr << "TID#" << threadID << " copying to offset: " << offset+i << endl;
	 blockLengths[offset+i]  = (*it).amount;
	 displacements[offset+i] = (*it).array - multiwriteOffsetPointer;
	 types[offset+i]         = (*it).mpiType;
	 ++i;
      }
   }
   
   // Write data to file:
   if (threadID == masterThreadID) {
      if (N_multiwriteUnits > 0) {
	 // Create an MPI struct containing the multiwrite units:
	 MPI_Datatype outputType;
	 MPI_Type_create_struct(N_multiwriteUnits,blockLengths,displacements,types,&outputType);
	 MPI_Type_commit(&outputType);

	 // Synchronize MPI processes. This is to make sure that file writing is not started 
	 // until every process has finished adding their multiwrite units:
	 MPI_Barrier(comm);
	 
	 // Write data to output file with a single collective call:
	 MPI_File_write_at_all(fileptr,offset,multiwriteOffsetPointer,1,outputType,MPI_STATUS_IGNORE);
	 MPI_Type_free(&outputType);
      } else {
	 // We have no data to write but we need to participate in the collective call anyway:
	 MPI_File_write_at_all(fileptr,offset,NULL,0,MPI_BYTE,MPI_STATUS_IGNORE);
      }
      
      // Deallocate memory:
      delete [] blockLengths; blockLengths = NULL;
      delete [] displacements; displacements = NULL;
      delete [] types; types = NULL;
   }

   // MPI master process writes footer tag:
   if (myrank == masterRank) {
      // Only the master thread is allowed to update XML tree:
      if (threadID == masterThreadID) {
	 // Count total number of bytes written to file:
	 uint64_t totalBytes = 0;
	 for(int i=0;i<N_processes;i++) totalBytes+=bytesPerProcess[i];
	 
	 XMLNode* root = xmlWriter->getRoot();
	 XMLNode* xmlnode = xmlWriter->find("VLSV",root);
	 XMLNode* node = xmlWriter->addNode(xmlnode,tagName,offset);
	 for (map<string,string>::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
	    xmlWriter->addAttribute(node,it->first,it->second);
	 }
	 xmlWriter->addAttribute(node,"vectorsize",vectorSize);
	 xmlWriter->addAttribute(node,"arraysize",totalBytes/dataSize/vectorSize);
	 xmlWriter->addAttribute(node,"datatype",dataType);
	 xmlWriter->addAttribute(node,"datasize",dataSize);
	 
	 // Update global file offset:
	 offset +=totalBytes;
      }
   }
   
   multiwriteInitialized = false;
   return true;;
}

bool VLSVWriter::writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,const std::string& dataType,
			    const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,char* array,const int& threadID) {
   bool success = true;
   if (startMultiwrite(dataType,arraySize,vectorSize,dataSize,threadID) == false) {
      success = false; return success;
   }
   if (addMultiwriteUnit(array,arraySize,threadID) == false) {
      success = false; return success;
   }
   
   if (endMultiwrite(arrayName,attribs,threadID) == false) {
      success = false; return success;
   }
   return success;
}

