/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2016 Finnish Meteorological Institute
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

#ifndef VLSV_READER_PARALLEL_H
#define VLSV_READER_PARALLEL_H

#include <mpi.h>

#include "vlsv_reader.h"
#include "mpiconversion.h"
#include "multi_io_unit.h"

namespace vlsv {

   class ParallelReader: public Reader {
    public:
      ParallelReader();
      ~ParallelReader();
   
      bool close();
      bool getArrayAttributes(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribsIn,
                              std::map<std::string,std::string>& attribsOut) const;
      bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                        uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& byteSize);
      bool getArrayInfoMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                              uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& dataSize);
      uint64_t getBytesRead();
      double getReadTime() const;
      bool getUniqueAttributeValues(const std::string& tagName,const std::string& attribName,std::set<std::string>& output) const;
      bool open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo=MPI_INFO_NULL);
      bool readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                           const uint64_t& begin,const uint64_t& amount,char* buffer);
      bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                     const uint64_t& begin,const uint64_t& amount,char* buffer);

      bool addMultireadUnit(char* buffer,const uint64_t& amount);
      bool endMultiread(const uint64_t& arrayOffset);
      bool startMultiread(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);

      template<typename T>
      bool read(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                const uint64_t& begin,const uint64_t& amount,T*& buffer,bool allocateMemory=true);
      template<typename T>
      bool readParameter(const std::string& parameterName,T& value);

    private:
      uint64_t bytesRead;             /**< Number of bytes read by this process.*/
      MPI_Comm comm;                  /**< MPI communicator used to read the file.*/
      MPI_File filePtr;               /**< MPI file pointer to input file.*/
      int masterRank;                 /**< MPI rank of master process.*/
      bool multireadStarted;          /**< If true, multiread mode has been initialized successfully.*/
      int myRank;                     /**< MPI rank of this process in communicator comm.*/
      bool parallelFileOpen;          /**< If true, all processes have opened input file successfully.*/
      int processes;                  /**< Number of MPI processes in communicator comm.*/
      double readTime;                /**< Time spent in seconds to read bytesRead bytes by this process.*/

      std::list<Multi_IO_Unit> multiReadUnits;

      bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
      bool flushMultiread(const size_t& unit,const MPI_Offset& currentOffset,std::list<Multi_IO_Unit>::iterator& start,std::list<Multi_IO_Unit>::iterator& stop);
   };

   template<typename T>
   bool ParallelReader::read(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                             const uint64_t& begin,const uint64_t& amount,T*& outBuffer,bool allocateMemory) {
      // Get array info to all processes:
      if (ParallelReader::getArrayInfo(tagName,attribs) == false) {
         return false;
      }

      // Check that requested read is inside the array:
      if (begin > arrayOpen.arraySize || (begin+amount) > arrayOpen.arraySize) return false;

      char* buffer = new char[amount*arrayOpen.vectorSize*arrayOpen.dataSize];
      if (ParallelReader::readArray(tagName,attribs,begin,amount,buffer) == false) {
         delete [] buffer; buffer = NULL;
         return false;
      }

      // Copy data from temporary buffer to output array:
      if (allocateMemory == true) outBuffer = new T[amount*arrayOpen.vectorSize];
      char* ptr = buffer;
      for (uint64_t i=0; i<amount; ++i) {
         for (uint64_t j=0; j<arrayOpen.vectorSize; ++j) {
            convertValue<T>(outBuffer[i*arrayOpen.vectorSize+j],ptr,arrayOpen.dataType,arrayOpen.dataSize,false);
            ptr += arrayOpen.dataSize;
         }
      }
      delete [] buffer; buffer = NULL;
      return true;
   }

   /** Read the value of a parameter. All processes must call this function simultaneously.
    * @param parameterName Name of the parameter. Only significant on master process.
    * @param value Variable where parameter's value will be written. Will be the same value on all processes upon successful exit.
    * @return If true, parameter was successfully read. All processes return the same value.*/
   template<typename T> inline
   bool ParallelReader::readParameter(const std::string& parameterName,T& value) {
      bool success = true;
      // Master process reads parameter value:
      if (myRank == masterRank) {
         if (Reader::readParameter(parameterName,value) == false) success = false;
      }

      // Master process broadcasts parameter value and 
      // value of variable 'status' to all other processes:
      uint8_t masterSuccess = 0;
      if (success == true) masterSuccess = 1;
      MPI_Bcast(&value,1,MPI_Type<T>(),masterRank,comm);
      MPI_Bcast(&masterSuccess,1,MPI_Type<uint8_t>(),masterRank,comm);
      success = (masterSuccess > 0);
      return success;
   }
}
#endif


