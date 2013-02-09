/** This file is part of VLSV file format.
 * 
 *  Copyright 2011, 2012 Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <iostream>
#include <string.h>

#include "vlsv_reader_parallel.h"

using namespace std;

VLSVParReader::VLSVParReader(): VLSVReader() {
   multireadStarted = false;
}
   
VLSVParReader::~VLSVParReader() {
   close();
}

bool VLSVParReader::close() {
   multireadStarted = false;
   if (fileOpen == false) return true;
   MPI_File_close(&filePtr);
   
   if (myRank == masterRank) filein.close();
   return true;
}

bool VLSVParReader::getArrayAttributes(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribsIn,
				       std::map<std::string,std::string>& attribsOut) const {
   bool success = true;
   
   // Master process reads footer:
   if (myRank == masterRank) {
      success = VLSVReader::getArrayAttributes(tagName,attribsIn,attribsOut);
   }
   
   // Check that master process read array attributes correctly:
   int globalSuccess = 0;
   if (success == true) globalSuccess = 1;
   MPI_Bcast(&globalSuccess,1,MPI_Type<int>(),masterRank,comm);
   if (globalSuccess == 0) return false;
   
   // Master process distributes contents of map attribsOut to all processes:
   // Broadcast number of entries in attribsOut:
   size_t N_attribs = 0;
   if (myRank == masterRank) N_attribs = attribsOut.size();
   MPI_Bcast(&N_attribs,1,MPI_Type<size_t>(),masterRank,comm);
   
   const unsigned int maxLength = 512;
   char attribName[maxLength];
   char attribValue[maxLength];
   if (myRank == masterRank) {
      for (map<string,string>::const_iterator it=attribsOut.begin(); it!=attribsOut.end(); ++it) {
	 strncpy(attribName,it->first.c_str(),maxLength-1);
	 strncpy(attribValue,it->second.c_str(),maxLength-1);
	 MPI_Bcast(attribName,maxLength,MPI_Type<char>(),masterRank,comm);
	 MPI_Bcast(attribValue,maxLength,MPI_Type<char>(),masterRank,comm);
      }
   } else {
      for (size_t i=0; i<N_attribs; ++i) {
	 MPI_Bcast(attribName,maxLength,MPI_Type<char>(),masterRank,comm);
	 MPI_Bcast(attribValue,maxLength,MPI_Type<char>(),masterRank,comm);
	 attribsOut[attribName] = attribValue;
      }
   }
   
   return success;
}

bool VLSVParReader::getArrayInfoMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				       uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& dataSize) {
   if (myRank != masterRank) {
      cerr << "(VLSVPARREADER): getArrayInfoMaster called on process #" << myRank << endl;
      exit(1);
   }
   return VLSVReader::getArrayInfo(tagName,attribs,arraySize,vectorSize,dataType,dataSize);
}

bool VLSVParReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
   bool success = true;
   if (myRank == masterRank) {
      success = VLSVReader::loadArray(tagName,attribs);
   }
   
   // Check that master read array info correctly:
   int globalSuccess = 0;
   if (success == true) globalSuccess = 1;
   MPI_Bcast(&globalSuccess,1,MPI_Type<int>(),masterRank,comm);
   if (globalSuccess == 0) success = false;
   if (success == false) return false;
   
   // Master broadcasts array info to all processes:
   MPI_Bcast(&arrayOpen.offset,    1,MPI_Type<streamoff>(),masterRank,comm);
   MPI_Bcast(&arrayOpen.arraySize, 1,MPI_Type<uint64_t>(), masterRank,comm);
   MPI_Bcast(&arrayOpen.vectorSize,1,MPI_Type<uint64_t>(), masterRank,comm);
   MPI_Bcast(&arrayOpen.dataType,  1,MPI_Type<int>(),      masterRank,comm);
   MPI_Bcast(&arrayOpen.dataSize,  1,MPI_Type<uint64_t>(), masterRank,comm);
   return success;
}

bool VLSVParReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				 uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& byteSize) {
   if (getArrayInfo(tagName,attribs) == false) return false;
   
   // Copy values to output variables:
   arraySize  = arrayOpen.arraySize;
   vectorSize = arrayOpen.vectorSize;
   dataType   = arrayOpen.dataType;
   byteSize   = arrayOpen.dataSize;
   return true;
}

bool VLSVParReader::getUniqueAttributeValues(const std::string& tagName,const std::string& attribName,std::set<std::string>& output) const {
   bool success = true;
   
   if (myRank == masterRank) {
      success = VLSVReader::getUniqueAttributeValues(tagName,attribName,output);
   }
   
   // Check that master read array info correctly:
   uint8_t masterSuccess = 0;
   if (success == true) masterSuccess = 1;
   MPI_Bcast(&masterSuccess,1,MPI_Type<uint8_t>(),masterRank,comm);
   if (masterSuccess == 0) success = false;
   if (success == false) return false;
   
   // Master broadcasts number of entries in set 'output':
   size_t N_entries = output.size();
   MPI_Bcast(&N_entries,1,MPI_Type<size_t>(),masterRank,comm);
   
   // Master broadcasts values in set 'output' to all other processes
   // who insert the received values to set 'output':
   const unsigned int maxLength = 512;
   char attribValue[maxLength];
   if (myRank == masterRank) {
      for (set<string>::const_iterator it=output.begin(); it!=output.end(); ++it) {
	 strncpy(attribValue,it->c_str(),maxLength-1);
	 MPI_Bcast(attribValue,maxLength,MPI_Type<char>(),masterRank,comm);
      }
   } else {
      for (size_t i=0; i<N_entries; ++i) {
	 MPI_Bcast(attribValue,maxLength,MPI_Type<char>(),masterRank,comm);
	 output.insert(attribValue);
      }
   }
   
   return success;
}
 
/** Add a file read unit. Note that multiReadStart function must have been called 
 * to initialize multiread mode before calling this function.
 * @param amount Number of array elements to read.
 * @param buffer Pointer to memory location where data from file is placed.
 * @return If true, multiread unit was added successfully.*/
bool VLSVParReader::multiReadAddUnit(const uint64_t& amount,char* buffer) {
   bool success = true;
   if (multireadStarted == false) return false;
   multiReadUnits.push_back(VLSV::Multi_IO_Unit(buffer,multiReadVectorType,amount));
   return success;
}

/** End multiread mode and read all data from file.
 * @param offset Offset into input array for this process, in units of array elements.
 * @return If true, data was read successfully.*/
bool VLSVParReader::multiReadEnd(const uint64_t& offset) {
   bool success = true;
   if (multireadStarted == false) return false;
   
   // Create an MPI datatype for reading all data at once:
   const uint64_t N_reads = multiReadUnits.size();
   int* blockLengths  = new int[N_reads];
   MPI_Aint* displacements = new MPI_Aint[N_reads];
   MPI_Datatype* datatypes = new MPI_Datatype[N_reads];
   
   // Copy multiread units to arrays defined above so that 
   // we can create an MPI struct that is used to read all data 
   // with a single collective call:
   MPI_Aint address;
   uint64_t counter = 0;
   for (list<VLSV::Multi_IO_Unit>::const_iterator it=multiReadUnits.begin(); it!=multiReadUnits.end(); ++it) {
      if (it->amount == 0) {
	 // MPI works with empty reads but datatype cannot be MPI_DATATYPE_NULL:
	 blockLengths[counter] = 0;
	 displacements[counter] = 0;
	 datatypes[counter] = MPI_Type<char>();
      } else {
	 MPI_Get_address(it->array,&address);
	 blockLengths[counter] = it->amount;
	 displacements[counter] = address;
	 datatypes[counter] = it->mpiType;
      }
      ++counter;
   }
   
   // Create MPI datatype containing all reads:
   MPI_Datatype readType;
   MPI_Type_create_struct(N_reads,blockLengths,displacements,datatypes,&readType);
   delete [] blockLengths;
   delete [] displacements;
   delete [] datatypes;
   multiReadUnits.clear();
   
   // Commit datatype and read everything in parallel:
   MPI_Type_commit(&readType);
   const uint64_t byteOffset = arrayOpen.offset + offset*arrayOpen.vectorSize*arrayOpen.dataSize;
   if (MPI_File_read_at_all(filePtr,byteOffset,MPI_BOTTOM,1,readType,MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      success = false;
   }
   MPI_Type_free(&readType);
   MPI_Type_free(&multiReadVectorType);
   multireadStarted = false;
   return success;
}

/** Start multiread mode. In multiread mode processes add zero or more file I/O units 
 * that define the data read from an array in VLSV file and where it is placed in memory.
 * File I/O units are defined by calling multiReadAddUnit. Data is not actually read until 
 * multiReadEnd is called. XML tag name and contents of list 'attribs' need to uniquely 
 * define the array.
 * @param tagName Array XML tag name in VLSV file.
 * @param attribs Additional attributes limiting array search.
 * @return If true, multiread mode was started successfully.
 * @see multiReadAddUnit.
 * @see multiReadEnd.*/
bool VLSVParReader::multiReadStart(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
   if (fileOpen == false) return false;
   bool success = true;
   multiReadUnits.clear();
   if (getArrayInfo(tagName,attribs) == false) return false;
   if (MPI_Type_contiguous(arrayOpen.vectorSize*arrayOpen.dataSize,MPI_Type<char>(),&multiReadVectorType) != MPI_SUCCESS) success = false;
   if (success == true) multireadStarted = true;
   return success;
}

/** Open a VLSV file for parallel reading.
 * @param fname Name of the VLSV file.
 * @param comm MPI communicator used in collective MPI operations.
 * @param masterRank MPI rank of master process.
 * @param mpiInfo Additional MPI info for optimizing file I/O.
 * @return If true, VLSV file was opened successfully.*/
bool VLSVParReader::open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo) {
   bool success = true;
   this->comm = comm;
   this->masterRank = masterRank;
   MPI_Comm_rank(comm,&myRank);
   MPI_Comm_size(comm,&processes);
   multireadStarted = false;
   
   // Attempt to open the given input file using MPI:
   fileName = fname;
   int accessMode = MPI_MODE_RDONLY;
   if (MPI_File_open(comm,const_cast<char*>(fileName.c_str()),accessMode,mpiInfo,&filePtr) != MPI_SUCCESS) success = false;
   else fileOpen = true;
   
   // Only master process reads file footer and endianness. This is done using 
   // VLSVReader open() member function:
   if (myRank == masterRank) {
      VLSVReader::open(fname);
   }
   
   // Check that all processes have opened the file successfully:
   unsigned char globalSuccess = 0;
   if (success == true) globalSuccess = 1;
   unsigned char* results = new unsigned char[processes];
   MPI_Allgather(&globalSuccess,1,MPI_Type<unsigned char>(),results,1,MPI_Type<unsigned char>(),comm);
   for (int i=0; i<processes; ++i) if (results[i] == 0) success = false;
   delete [] results; results = NULL;
   if (success == false) return success;
   
   // Broadcast file endianness to all processes:
   MPI_Bcast(&endiannessFile,1,MPI_Type<unsigned char>(),masterRank,comm);
   
   return success;
}
   
bool VLSVParReader::readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
				    const uint64_t& begin,const uint64_t& amount,char* buffer) {
   if (myRank != masterRank) {
      cerr << "(VLSVPARREADER) readArrayMaster erroneously called on process #" << myRank << endl;
      exit(1);
   }
   // VLSVreadArray reads offset from XML tree into master only
   return VLSVReader::readArray(tagName,attribs,begin,amount,buffer);
}

/** Read data from an array in VLSV file using collective MPI file I/O operations.
 * XML tag name and contents of list 'attribs' need to uniquely define the array.
 * @param tagName Array XML tag name.
 * @param attribs Additional attributes limiting array search.
 * @param begin First array element read, i.e. this process' offset into the array.
 * @param amount Number of array elements to read.
 * @param buffer Buffer in which data is read from VLSV file.
 * @return If true, this process read its data successfully.*/
bool VLSVParReader::readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			      const uint64_t& begin,const uint64_t& amount,char* buffer) {
   if (fileOpen == false) return false;
   bool success = true;
   
   // Fetch array info to all processes:
   if (getArrayInfo(tagName,attribs) == false) return false;
   const MPI_Offset start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   const int readBytes    = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   
   // Read data on all processes in parallel:
   if (MPI_File_read_at_all(filePtr,start,buffer,readBytes,MPI_Type<char>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      success = false;
   }
   return success;
}
