/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013 Finnish Meteorological Institute
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
#include <fstream>

#ifdef PROFILE
   #include <profiler.h>
   static int fileOpenID = -1;
   static int startMWgatherID = -1;
   static int startMWscatterID = -1;
   static int endMWstructID = -1;
   static int endMWwriteID = -1;
   static int addMWunitID = -1;
   static int startMWID = -1;
   static int addMWID = -1;
   static int endMWID = -1;
#endif

#include "mpiconversion.h"
#include "vlsv_common_mpi.h"
#include "vlsv_writer.h"

using namespace std;

namespace vlsv {
   /** Check that all processes have the same success status.
    * @param myStatus If true, this process has successfully executed all code.
    * @return If true, all processes called checkSuccess with myStatus set to 'true'.*/
   bool checkSuccess(const bool& myStatus,MPI_Comm comm) {
      // If myStatus is false, set mySuccess to value 1:
      int32_t mySuccess = 0;
      int32_t globalSuccess = 0;
      if (myStatus == false) mySuccess = 1;
      
      // Sum mySuccess values to master process and broadcast to all processes:
      MPI_Reduce(&mySuccess,&globalSuccess,1,MPI_Type<int32_t>(),MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Bcast(&globalSuccess,1,MPI_Type<int32_t>(),0,MPI_COMM_WORLD);
      
      // If globalSuccess equals zero all processes called this function with myStatus set to 'true':
      if (globalSuccess == 0) return true;
      return false;
   }
       
   /** Constructor for Writer.*/
   Writer::Writer() {
      blockLengths = NULL;
      bytesPerProcess = NULL;
      displacements = NULL;
      endMultiwriteCounter = 0;
      fileOpen = false;
      initialized = false;
      multiwriteFinalized = false;
      multiwriteInitialized = false;
      multiwriteOffsetPointer = NULL;
      N_multiwriteUnits = 0;
      offset = 0;
      offsets = NULL;
      types = NULL;
      xmlWriter = NULL;
      comm = MPI_COMM_NULL;
   }

   /** Destructor for Writer. Deallocates XML writer.*/
   Writer::~Writer() {
      if (fileOpen == true) close();
      if (comm != MPI_COMM_NULL) MPI_Comm_free(&comm);
      delete [] blockLengths; blockLengths = NULL;
      delete [] bytesPerProcess; bytesPerProcess = NULL;
      delete [] displacements; displacements = NULL;
      delete [] offsets; offsets = NULL;
      delete [] types; types = NULL;
      delete xmlWriter; xmlWriter = NULL;
   }

   bool Writer::addMultiwriteUnit(char* array,const uint64_t& arrayElements) {
   // Check that startMultiwrite has initialized correctly:
      #ifdef PROFILE
         profile::start("add mw unit",addMWunitID);
      #endif

      // addMultiwriteUnit does not call any collective MPI functions so 
      // it is safe to exit immediately here if error(s) have occurred:
      if (initialized == false) return false;
      if (multiwriteInitialized == false) return false;
      
      // Do not add zero-size arrays to multiWriteUnits:
      if (arrayElements == 0) return true;
      
      multiwriteUnits[0].push_back(Multi_IO_Unit(array,getMPIDatatype(vlsvType,dataSize),arrayElements*vectorSize));

      #ifdef PROFILE
         profile::stop();
      #endif
      return true;
   }

   /** Close a file that has been previously opened by calling Writer::open.
    * After the file has been closed the MPI master process appends an XML footer 
    * to the end of the file, and writes an offset to the footer to the start of 
    * the file.
    * @return If true, the file was closed successfully. If false, a file may not 
    * have been opened successfully by Writer::open.*/
   bool Writer::close() {
      // If a file was never opened, exit immediately:
      if (fileOpen == false) return false;
      
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
         const double t_start = MPI_Wtime();
         xmlWriter->print(footer);
         writeTime += (MPI_Wtime() - t_start);
         footer.close();

         // Write header position to the beginning of binary file:
         footer.open(fileName.c_str(),fstream::in | fstream::out | fstream::binary | fstream::ate);
         char* ptr = reinterpret_cast<char*>(&footerOffset);
         footer.seekp(sizeof(uint64_t));
         footer.write(ptr,sizeof(uint64_t));
         footer.close();
      }

      initialized = false;
      delete [] blockLengths; blockLengths = NULL;
      delete [] bytesPerProcess; bytesPerProcess = NULL;
      delete [] displacements; displacements = NULL;
      delete [] offsets; offsets = NULL;
      delete [] types; types = NULL;
      delete xmlWriter; xmlWriter = NULL;

      // Wait for master process to finish:
      MPI_Barrier(comm);
      fileOpen = false;
      MPI_Comm_free(&comm);
      return true;
   }
   
   /** Get the total amount of bytes written to VLSV file. This function 
    * returns a meaningful value at master process only.
    * @return Total number of bytes written to output files by all processes.*/
   uint64_t Writer::getBytesWritten() const {return bytesWritten;}

   /** Open a VLSV file for parallel output. The file is opened on all processes 
    * in the given communicator. Additionally, master MPI process writes a 
    * header into the output file and caches a footer which will be written 
    * in Writer::close. If a file has already been opened and Writer::open 
    * is called again, the currently open file is closed before the new file is opened.
    * @param fname The name of the output file.
    * @param comm MPI communicator used in writing.
    * @param masterProcessID ID of the MPI master process.
    * @return If true, a file was opened successfully.*/
   bool Writer::open(const std::string& fname,MPI_Comm comm,const int& masterProcessID,MPI_Info mpiInfo) {
      #ifdef PROFILE
         profile::start("open",fileOpenID);
      #endif
      bool success = true;
   
      // If a file with the same name has already been opened, return immediately.
      // Otherwise close the currently open file before opening the new file.
      if (fileOpen == true) {
         if (fname == fileName) return true;
         close();
      }
   
      MPI_Comm_dup(comm,&(this->comm));
      masterRank = masterProcessID;
      MPI_Comm_rank(this->comm,&myrank);
      MPI_Comm_size(this->comm,&N_processes);
      bytesWritten = 0;
      writeTime = 0;

      // Allocate per-thread storage:
      multiwriteOffsets.resize(1);
      multiwriteUnits.resize(1);
   
      // All processes in communicator comm open the same file. If a file with the 
      // given name already exists it is deleted. Note: We found out that MPI_File_open 
      // failed quite often in meteo, at least when writing many small files. It was 
      // possibly caused by MPI_File_delete call, that's the reason for the barrier.
      int accessMode = (MPI_MODE_WRONLY | MPI_MODE_CREATE);
      fileName = fname;
      if (myrank == masterRank) MPI_File_delete(const_cast<char*>(fname.c_str()),mpiInfo);
      MPI_Barrier(comm);
      if (MPI_File_open(comm,const_cast<char*>(fileName.c_str()),accessMode,mpiInfo,&fileptr) != MPI_SUCCESS) {
         fileOpen = false;
         return fileOpen;
      }
      #ifdef PROFILE
         profile::stop();
      #endif
   
      offset = 0;           //offset set to 0 when opening a new file
      MPI_File_set_view(fileptr,0,MPI_BYTE,MPI_BYTE,const_cast<char*>("native"),mpiInfo);
       
      // Only master process needs these arrays:
      if (myrank == masterRank) {
         offsets = new MPI_Offset[N_processes];
         bytesPerProcess = new uint64_t[N_processes];    
      }
      
      // Master process opens an XML tree for storing the footer:
      if (myrank == masterRank) {
         xmlWriter     = new muxml::MuXML();
         muxml::XMLNode* root = xmlWriter->getRoot();
         xmlWriter->addNode(root,"VLSV","");
      }
      
      // Master writes 2 64bit integers to the start of file. 
      // Second value will be overwritten in close() function to tell 
      // the position of footer:
      if (myrank == masterRank) {
         // Write file endianness to the first byte:
         uint64_t endianness = 0;
         unsigned char* ptr = reinterpret_cast<unsigned char*>(&endianness);
         ptr[0] = detectEndianness();
         const double t_start = MPI_Wtime();
         if (MPI_File_write_at(fileptr,0,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
         if (MPI_File_write_at(fileptr,8,&endianness,1,MPI_Type<uint64_t>(),MPI_STATUS_IGNORE) != MPI_SUCCESS) success = false;
         writeTime += (MPI_Wtime() - t_start);
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
         delete xmlWriter; xmlWriter = NULL;
      }

      initialized = true;
      fileOpen    = success;
      return fileOpen;
   }

   bool Writer::startMultiwrite(const string& datatype,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize) {
      // Check that all processes have made it this far without error(s):
      bool success = true;
      if (fileOpen == false) success = false;
      if (initialized == false) success = false;
      if (checkSuccess(success,comm) == false) return false;

      // Clear per-thread storage:
        {
           vector<list<Multi_IO_Unit> > dummy(1);
           multiwriteUnits.swap(dummy);
        }
      multiwriteOffsets[0] = numeric_limits<unsigned int>::max();

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
      #ifdef PROFILE
         profile::start("gather",startMWgatherID);
      #endif
      myBytes = arraySize * vectorSize * dataSize;
      MPI_Gather(&myBytes,1,MPI_Type<uint64_t>(),bytesPerProcess,1,MPI_Type<uint64_t>(),masterRank,comm);
      #ifdef PROFILE
         profile::stop();
      #endif

      // MPI master process calculates an offset to the output file for all processes:
      if (myrank == masterRank) {
         offsets[0] = offset;
         for (int i=1; i<N_processes; ++i) offsets[i] = offsets[i-1] + bytesPerProcess[i-1];
      }
   
      // MPI master scatters offsets:
      #ifdef PROFILE
         profile::start("scatter",startMWscatterID);
      #endif
      MPI_Scatter(offsets,1,MPI_Type<uint64_t>(),&offset,1,MPI_Type<uint64_t>(),masterRank,comm);
      #ifdef PROFILE
         profile::stop();
      #endif

      multiwriteInitialized = true;
      return multiwriteInitialized;
   }

   /** Write multiwrite units to file.
    * @param tagName Name of the XML tag for this array.
    * @param attribs Attributes for the XML tag.
    * @return If true, array was successfully written to file.*/
   bool Writer::endMultiwrite(const std::string& tagName,const std::map<std::string,std::string>& attribs) {
      // Check that multiwrite mode has started successfully on all processes:
      bool success = true;
      if (initialized == false) success = false;
      if (multiwriteInitialized == false) success = false;
      if (checkSuccess(success,comm) == false) return false;

      // Allocate memory for an MPI_Struct that is used to 
      // write all multiwrite units with a single collective call.

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
      multiwriteOffsetPointer = NULL;
      for (size_t i=0; i<multiwriteUnits.size(); ++i) {
         if (multiwriteUnits[i].size() == 0) continue;
         multiwriteOffsetPointer = multiwriteUnits[i].begin()->array;
      }
   
      // Every thread copies its multiwrite units into a global MPI struct:
      // DEPRECATED
      if (multiwriteUnits[0].size() > 0) {
         unsigned int offset = multiwriteOffsets[0];
         unsigned int i = 0;
         for (list<Multi_IO_Unit>::const_iterator it=multiwriteUnits[0].begin(); it!=multiwriteUnits[0].end(); ++it) {
            blockLengths[offset+i]  = (*it).amount;
            displacements[offset+i] = (*it).array - multiwriteOffsetPointer;
            types[offset+i]         = (*it).mpiType;
            ++i;
         }
      }

      // Write data to file:
      if (N_multiwriteUnits > 0) {
         // Create an MPI struct containing the multiwrite units:
         #ifdef PROFILE
            profile::start("struct creation",endMWstructID);
         #endif
         MPI_Datatype outputType;
         MPI_Type_create_struct(N_multiwriteUnits,blockLengths,displacements,types,&outputType);
         MPI_Type_commit(&outputType);
         #ifdef PROFILE
            profile::stop();
         #endif
      
         // Write data to output file with a single collective call:
         #ifdef PROFILE
            profile::start("write",endMWwriteID);
         #endif
         const double t_start = MPI_Wtime();
         MPI_File_write_at_all(fileptr,offset,multiwriteOffsetPointer,1,outputType,MPI_STATUS_IGNORE);
         writeTime += (MPI_Wtime() - t_start);
         MPI_Type_free(&outputType);
         #ifdef PROFILE
            profile::stop();
         #endif
      } else {
         // Process has no data to write but needs to participate in the collective call to prevent deadlock:
         #ifdef PROFILE
            profile::start("write",endMWwriteID);
         #endif
         const double t_start = MPI_Wtime();
         MPI_File_write_at_all(fileptr,offset,NULL,0,MPI_BYTE,MPI_STATUS_IGNORE);
         writeTime += (MPI_Wtime() - t_start);
         #ifdef PROFILE
            profile::stop();
         #endif
      }

      // Deallocate memory:
      delete [] blockLengths; blockLengths = NULL;
      delete [] displacements; displacements = NULL;
      delete [] types; types = NULL;

      // MPI master process writes footer tag:
      if (myrank == masterRank) {
         // Count total number of bytes written to file:
         uint64_t totalBytes = 0;
         for (int i=0; i<N_processes; ++i) totalBytes += bytesPerProcess[i];
         
         muxml::XMLNode* root = xmlWriter->getRoot();
         muxml::XMLNode* xmlnode = xmlWriter->find("VLSV",root);
         muxml::XMLNode* node = xmlWriter->addNode(xmlnode,tagName,offset);
         for (map<string,string>::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
            xmlWriter->addAttribute(node,it->first,it->second);
         }
         xmlWriter->addAttribute(node,"vectorsize",vectorSize);
         xmlWriter->addAttribute(node,"arraysize",totalBytes/dataSize/vectorSize);
         xmlWriter->addAttribute(node,"datatype",dataType);
         xmlWriter->addAttribute(node,"datasize",dataSize);
         
         // Update global file offset:
         offset += totalBytes;
         bytesWritten += totalBytes;
      }

      multiwriteInitialized = false;
      return success;
   }

   bool Writer::writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,const std::string& dataType,
                           const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,const char* array) {
      // Check that everything is OK before continuing:
      if (initialized == false) return false;
      if (fileOpen == false) return false;
   
      bool success = true;
      #ifdef PROFILE
         profile::start("start multiwrite",startMWID);
      #endif
      if (startMultiwrite(dataType,arraySize,vectorSize,dataSize) == false) {
         success = false; return success;
      }
   
      #ifdef PROFILE
         profile::stop();
         profile::start("add multiwrite",addMWID);
      #endif
      char* arrayPtr = const_cast<char*>(array);
      if (addMultiwriteUnit(arrayPtr,arraySize) == false) {
         success = false; return success;
      }
   
      #ifdef PROFILE
         profile::stop();
         profile::start("end multiwrite",endMWID);
      #endif
      if (endMultiwrite(arrayName,attribs) == false) {
         success = false; return success;
      }
      #ifdef PROFILE
         profile::stop();
      #endif
   
      return success;
   }

} // namespace vlsv
