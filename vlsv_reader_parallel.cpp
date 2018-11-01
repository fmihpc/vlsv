/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2015 Finnish Meteorological Institute
 *  Copyright 2016 Arto Sandroos
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

#include "vlsv_common_mpi.h"
#include "vlsv_reader_parallel.h"

using namespace std;

namespace vlsv {

   /** Default constructor for class ParallelReader.*/
   ParallelReader::ParallelReader(): Reader() {
      multireadStarted = false;
   }
   
   /** Destructor for class ParallelReader. It closes the input file (if still open).*/
   ParallelReader::~ParallelReader() {
      close();
   }

   /** Add a multi-read unit. Note that startMmultiread must have been called 
    * to initialize multi-read mode before calling this function. The data is not 
    * actually read until endMultiread is called.
    * @param buffer Pointer to allocated memory location where data from file is placed.
    * @param arrayElements Number of array elements to read.
    * @return If true, multiread unit was added successfully.
    * @see startMultiread.
    * @see endMultiread.*/
   bool ParallelReader::addMultireadUnit(char* buffer,const uint64_t& arrayElements) {
      bool success = true;
      if (multireadStarted == false) return false;
      
      // Ignore zero size reads:
      if (arrayElements == 0) return true;
      
      // Get the byte size of the MPI primitive datatype (MPI_INT etc.) used here:
      int datatypeBytesize;
      MPI_Type_size(getMPIDatatype(arrayOpen.dataType,arrayOpen.dataSize),&datatypeBytesize);

      // Calculate the maximum number of array elements written using a single multi-write.
      // Note: element = vector of size vectorSize, each vector element has byte size of datatypeBytesize.
      const size_t maxElementsPerRead = getMaxBytesPerRead() / (datatypeBytesize*arrayOpen.vectorSize);

      // Split the multi-read if the array has more elements than what we can 
      // read from output file using a single MPI collective:
      if (arrayElements > maxElementsPerRead) {
         // Calculate how many collectives this process needs:
         size_t N_reads = arrayElements / maxElementsPerRead;
         if (arrayElements % maxElementsPerRead != 0) ++N_reads;

         // Add N_reads multi-read units:
         for (size_t i=0; i<N_reads; ++i) {
            auto elements = maxElementsPerRead;
            if ((i+1)*maxElementsPerRead >= arrayElements) elements = arrayElements - i*maxElementsPerRead;

            const size_t byteOffset = maxElementsPerRead*arrayOpen.vectorSize*datatypeBytesize;
            multiReadUnits.push_back(
                Multi_IO_Unit(buffer+i*byteOffset,
                              getMPIDatatype(arrayOpen.dataType,arrayOpen.dataSize),
                              elements*arrayOpen.vectorSize));
         }
      } else {
         multiReadUnits.push_back(
             Multi_IO_Unit(buffer,
                           getMPIDatatype(arrayOpen.dataType,arrayOpen.dataSize),
                           arrayElements*arrayOpen.vectorSize));
      }
      
      return success;
   }

   /** End multi-read mode and read data from file on all processes.
    * Note that startMultiread and one or more addMultireadUnit calls must have 
    * been made before calling this function.
    * @param arrayOffset Offset into input array relative to array start on file. 
    * Determines where this process starts to read its data. Offset is given in 
    * units of array elements.
    * @return If true, all processes read their data successfully.
    * @see startMultiread.
    * @see addMultireadUnit.*/
   bool ParallelReader::endMultiread(const uint64_t& arrayOffset) {
      bool success = true;
      if (multireadStarted == false) success = false;
      if (checkSuccess(success,comm) == false) return false;

      // Calculate how many collective MPI calls are needed to 
      // read all the data from input file:
      size_t inputBytesize    = 0;
      size_t myCollectiveCalls = 0;
      if (multiReadUnits.size() > 0) myCollectiveCalls = 1;

      vector<pair<list<Multi_IO_Unit>::iterator,list<Multi_IO_Unit>::iterator> > multireadList;
      auto first = multiReadUnits.begin();
      auto last  = multiReadUnits.begin();
      for (auto it=multiReadUnits.begin(); it!=multiReadUnits.end(); ++it) {
         if (inputBytesize + it->amount*arrayOpen.dataSize > getMaxBytesPerRead()) {
            multireadList.push_back(make_pair(first,last));
            first = it; last = it;

            inputBytesize = 0;
            ++myCollectiveCalls;
         }
         inputBytesize += it->amount*arrayOpen.dataSize;
         ++last;
      }
      multireadList.push_back(make_pair(first,last));
      
      // Reduce the maximum number of needed collective reads to all processes:
      size_t N_collectiveCalls;
      MPI_Allreduce(&myCollectiveCalls,&N_collectiveCalls,1,MPI_Type<size_t>(),MPI_MAX,comm);

      // If more collective calls are made than what this process needs, 
      // insert dummy reads to the end of multireadList:
      if (N_collectiveCalls > multireadList.size()) {
         const size_t N_dummyCalls = N_collectiveCalls-multireadList.size();
         for (size_t i=0; i<N_dummyCalls; ++i) {
            multireadList.push_back(make_pair(multiReadUnits.end(),multiReadUnits.end()));
         }
      }

      MPI_Offset unitOffset = 0;
      unitOffset  = arrayOpen.offset;                                    // Offset of array start relative to file start
      unitOffset += arrayOffset*arrayOpen.vectorSize*arrayOpen.dataSize; // Byte offset relative to array start where 
                                                                         // this process starts to read data from.

      for (size_t i=0; i<multireadList.size(); ++i) {
         if (flushMultiread(i,unitOffset,multireadList[i].first,multireadList[i].second) == false) success = false;
         for (auto it=multireadList[i].first; it!=multireadList[i].second; ++it) {
            unitOffset += it->amount*arrayOpen.dataSize;
         }
      }

      multireadStarted = false;
      return checkSuccess(success,comm);
   }

   /** Close the input file.
    * @return If true, the input file was closed successfully.
    * @see vlsv::ParallelReader::open().*/
   bool ParallelReader::close() {
      multireadStarted = false;
      if (parallelFileOpen == true) {
         MPI_File_close(&filePtr);
         parallelFileOpen = false;
      }

      if (myRank == masterRank) filein.close();
      return true;
   }

   /** Get the XML attributes for the given array on all processes. The file 
    * footer is read by master process only, who then broadcasts the contents 
    * to all processes. This function returns the same value on all processes.
    * @param tagName Name of the array's XML tag. Only significant on master process.
    * @param attribsIn XML tag attributes that uniquely define the array. Only significant on master process.
    * @param attribsOut XML tag attributes read from the input file. The contents of attribsOut will be the same on all processes.
    * @return If true, array attributes were read successfully. All processes return the same value.*/
   bool ParallelReader::getArrayAttributes(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribsIn,
                                           std::map<std::string,std::string>& attribsOut) const {
      bool success = true;
   
      // Master process reads footer:
      if (myRank == masterRank) {
         success = Reader::getArrayAttributes(tagName,attribsIn,attribsOut);
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
         for (auto& it: attribsOut) {
            #ifdef WINDOWS
               strncpy_s(attribName ,maxLength,it.first.c_str() ,it.first.size() );
               strncpy_s(attribValue,maxLength,it.second.c_str(),it.second.size());
            #else
               strncpy(attribName ,it.first.c_str() ,maxLength-1);
               strncpy(attribValue,it.second.c_str(),maxLength-1);
            #endif
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

   bool ParallelReader::getArrayInfoMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                                           uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& dataSize) {
      if (myRank != masterRank) {
         cerr << "(PARALLEL READER): getArrayInfoMaster called on process #" << myRank << endl;
         exit(1);
      }
      return Reader::getArrayInfo(tagName,attribs,arraySize,vectorSize,dataType,dataSize);
   }

   /** Read array metadata to all processes.
    * @param tagName Name of the XML tag corresponding to the array. Only significant on master process.
    * @param attribs A list of attribute,value pairs that uniquely identify the array. Only significant on master process.
    * @return If true, array metadata was read successfully. All processes return the same value.*/
   bool ParallelReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
      bool success = true;
      if (myRank == masterRank) {
         success = Reader::loadArray(tagName,attribs);
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

   /** Read array metadata to all processes.
    * @param Name of the array's XML tag. Only significant on master process.
    * @param attribs A list of attribute,value pairs that uniquely identify the array. Only significant on master process.
    * @param arraySize Total number of elements in array. Same value on all processes upon successful exit.
    * @param vectorSize Size of the data vector in array. Same value on all processes upon successful exit.
    * @param dataType Datatype of each vector element. Same value on all processes upon successful exit.
    * @param byteSize Byte size of the datatype. Same value on all processes upon successful exit.
    * @return If true, array metadata was successfully read and output variables contain meaningful values.
    * Same return value is returned on all processes.*/
   bool ParallelReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                                     uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& byteSize) {
      if (getArrayInfo(tagName,attribs) == false) return false;
   
      // Copy values to output variables:
      arraySize  = arrayOpen.arraySize;
      vectorSize = arrayOpen.vectorSize;
      dataType   = arrayOpen.dataType;
      byteSize   = arrayOpen.dataSize;
      return true;
   }

   /** Get the amount of bytes read from input file so far. On master process this 
    * function returns the total number of bytes read by all processes. On other 
    * processes the returned value is the number of bytes read by that process.
    * This function must be called simultaneously by all processes in the communicator.
    * @return Number of bytes read.*/
   uint64_t ParallelReader::getBytesRead() {
      uint64_t totalBytesRead;
      MPI_Reduce(&bytesRead,&totalBytesRead,1,MPI_Type<uint64_t>(),MPI_SUM,masterRank,comm);
      return totalBytesRead;
   }
   
   /** Get the time spent by this process in file I/O. The approximate data rate is 
    * getBytesRead() / getReadTime() when calculated on master process in the communicator 
    * used to read the file.
    * @return Time in seconds spent in file I/O.*/
   double ParallelReader::getReadTime() const {
      return readTime;
   }

   /** Get unique XML attribute values for given tag name. This function 
    * must be called by all processes simultaneously.
    * @param tagName Name of the XML tag. Only significant on master process.
    * @param attribName Name of the attribute. Only significant on master process.
    * @param output Unique attribute values are inserted here. The content will be the same on all processes.
    * @return If true, attribute values were read successfully. All processes return the same value.*/
   bool ParallelReader::getUniqueAttributeValues(const std::string& tagName,const std::string& attribName,
                                                 std::set<std::string>& output) const {
      bool success = true;

      // First the master process reads unique attribute values:                                              
      if (myRank == masterRank) {
         success = Reader::getUniqueAttributeValues(tagName,attribName,output);
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
         for (auto& it: output) {
            #ifdef WINDOWS
               strncpy_s(attribValue,maxLength,it.c_str(),it.size());
            #else
               strncpy(attribValue,it.c_str(),maxLength-1);
            #endif
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

   bool ParallelReader::flushMultiread(const size_t& unit,const MPI_Offset& fileOffset,
                                       std::list<Multi_IO_Unit>::iterator& start,std::list<Multi_IO_Unit>::iterator& stop) {
      bool success = true;

      // Count the number of multi-read units read:
      uint64_t N_multiReadUnits = 0;
      for (auto it=start; it!=stop; ++it) {
         ++N_multiReadUnits;
      }

      // Create an MPI datatype for reading all units with a single collective call:
      int* blockLengths       = new int[N_multiReadUnits];
      MPI_Aint* displacements = new MPI_Aint[N_multiReadUnits];
      MPI_Datatype* datatypes = new MPI_Datatype[N_multiReadUnits];

      char* multireadOffsetPointer = NULL;
      if (N_multiReadUnits > 0) multireadOffsetPointer = start->array;
      
      // Copy pointers etc. to MPI struct:
      size_t i=0;
      size_t amount = 0;
      for (auto it=start; it!=stop; ++it) {
         blockLengths[i]  = it->amount;
         displacements[i] = it->array - multireadOffsetPointer;
         datatypes[i]     = it->mpiType;

         int datatypeBytesize;
         MPI_Type_size(it->mpiType,&datatypeBytesize);
         amount += it->amount * datatypeBytesize;

         ++i;
      }

      // Read data from file:
      if (N_multiReadUnits > 0) {
         // Create an MPI struct containing the multiread units:
         MPI_Datatype inputType;
         MPI_Type_create_struct(N_multiReadUnits,blockLengths,displacements,datatypes,&inputType);
         MPI_Type_commit(&inputType);

         // Read data from output file with a single collective call:
         const auto t_start = MPI_Wtime();
         MPI_File_read_at_all(filePtr,fileOffset,multireadOffsetPointer,1,inputType,MPI_STATUS_IGNORE);
         readTime += (MPI_Wtime() - t_start);
         MPI_Type_free(&inputType);
         
         bytesRead += amount;
      } else {
         // Process has no data to read but needs to participate in the collective call to prevent deadlock:
         const auto t_start = MPI_Wtime();
         MPI_File_read_at_all(filePtr,fileOffset,NULL,0,MPI_BYTE,MPI_STATUS_IGNORE);
         readTime += (MPI_Wtime() - t_start);
      }

      delete [] blockLengths;
      delete [] displacements;
      delete [] datatypes;
      return success;
   }

   /** Open a VLSV file for parallel reading.
    * @param fname Name of the VLSV file. Only significant on master process.
    * @param comm MPI communicator used in collective MPI operations.
    * @param masterRank MPI rank of master process. Must be the same value on all processes.
    * @param mpiInfo Additional MPI info for optimizing file I/O. Must be the same value on all processes.
    * @return If true, VLSV file was opened successfully. All processes return the same value.*/
   bool ParallelReader::open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo) {
      bool success = true;
      this->comm = comm;
      this->masterRank = masterRank;
      MPI_Comm_rank(comm,&myRank);
      MPI_Comm_size(comm,&processes);
      multireadStarted = false;

      // Broadcast input file name to all processes:
      if (broadcast(fname,fileName,this->comm,masterRank) == false) return false;

      // Attempt to open the given input file using MPI:
      fileName = fname;
      int accessMode = MPI_MODE_RDONLY;
      if (MPI_File_open(comm,const_cast<char*>(fileName.c_str()),accessMode,mpiInfo,&filePtr) != MPI_SUCCESS) success = false;
      else parallelFileOpen = true;
      
      if (success == false) cerr << "Failed to open parallel file" << endl;
      
      // Only master process reads file footer and endianness. This is done using 
      // VLSVReader open() member function:
      if (myRank == masterRank) {
         if (Reader::open(fname) == false) success = false;
      }
      if (success == false) cerr << "MASTER failed to open VLSV file" << endl;

      // Broadcast file endianness to all processes:
      MPI_Bcast(&endiannessFile,1,MPI_Type<unsigned char>(),masterRank,comm);

      bytesRead = 0;
      return checkSuccess(success,this->comm);
   }

   bool ParallelReader::readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
					                    const uint64_t& begin,const uint64_t& amount,char* buffer) {
      if (myRank != masterRank) {
         cerr << "(PARALLEL READER) readArrayMaster erroneously called on process #" << myRank << endl;
         exit(1);
      }
      // readArray reads offset from XML tree into master only
      return Reader::readArray(tagName,attribs,begin,amount,buffer);
   }

   /** Read data from an array in VLSV file using collective MPI file I/O operations.
    * XML tag name and contents of list 'attribs' need to uniquely define the array.
    * @param tagName Array XML tag name. Only significant on master process.
    * @param attribs Additional attributes limiting array search. Only significant on master process.
    * @param begin First array element read, i.e. this process' offset into the array.
    * @param amount Number of array elements to read.
    * @param buffer Buffer in which data is read from VLSV file.
    * @return If true, array contents were successfully read. All processes return the same value.*/
   bool ParallelReader::readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                                  const uint64_t& begin,const uint64_t& amount,char* buffer) {
      bool success = true;

      // Fetch array info to all processes:
      if (getArrayInfo(tagName,attribs) == false) return false;
      const MPI_Offset start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
      const uint64_t readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;

      // If readBytes is larger than getMaxBytesPerRead() this process needs 
      // more than one collective call to read in all the data.
      const uint64_t maxBytes = getMaxBytesPerRead();
      uint64_t myExtraCollectiveReads = readBytes / maxBytes;

      // There's always at least one read:
      ++myExtraCollectiveReads;

      // Reduce the max number of required collective calls to all 
      // processes to prevent deadlock:
      uint64_t globalExtraCollectiveReads;
      MPI_Allreduce(&myExtraCollectiveReads,&globalExtraCollectiveReads,1,MPI_Type<uint64_t>(),MPI_MAX,comm);

      // Read data:
      const auto t_start = MPI_Wtime();
      uint64_t offset = 0;
      for (uint64_t counter=0; counter<globalExtraCollectiveReads; ++counter) {
         char*  pos;
         uint64_t readSize;

         if (counter < (myExtraCollectiveReads-1)) {
            pos      = buffer + counter*maxBytes;
            readSize = maxBytes;
         } else if (counter == (myExtraCollectiveReads-1)) {
            // This is the last read operation for this process
            pos      = buffer + counter*maxBytes;
            readSize = readBytes - counter*maxBytes;
         } else {
            // Nothing left to read
            pos      = NULL;
            readSize = 0;
         }

         MPI_Status status;
         if (MPI_File_read_at_all(filePtr,start+counter*maxBytes,pos,readSize,MPI_BYTE,&status) != MPI_SUCCESS) {
            success = false;
         }

         // Check that we got everything we requested:
         if(readSize>0) {
            int bytesReceived;
            MPI_Get_count(&status,MPI_BYTE,&bytesReceived);
            if (bytesReceived != readSize) {
               stringstream ss;
               ss << "ERROR in vlsv::ParallelReader! I only got " << bytesReceived << "/" << readSize;
               ss << " bytes in " << __FILE__ << ":" << __LINE__ << endl;
               cerr << ss.str();
               success = false;
            }
         }

         offset += readSize;
      }
      readTime  += (MPI_Wtime() - t_start);
      bytesRead += amount*arrayOpen.vectorSize*arrayOpen.dataSize;

      return checkSuccess(success,comm);
   }

   /** Start multi-read mode. In multi-read mode processes add zero or more file I/O units 
    * that define the data that is read from an array in VLSV file, and where it is placed in memory.
    * File I/O units are defined by calling addMultireadUnit. Data is not actually read until 
    * endMultiread is called. XML tag name and contents of list 'attribs' need to uniquely 
    * define the array.
    * @param tagName Array XML tag name in VLSV file. Only significant on master process.
    * @param attribs Additional attributes that uniquely define the array in file. Only significant on master process.
    * @return If true, multi-read mode was started successfully. All processes return the same value.
    * @see addMultireadUnit.
    * @see endMultiread.*/
   bool ParallelReader::startMultiread(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
      if (parallelFileOpen == false) return false;
      bool success = true;
      multiReadUnits.clear();
      if (getArrayInfo(tagName,attribs) == false) {
         return false;
      }
      if (success == true) multireadStarted = true;
      return success;
   }

} // namespace vlsv
