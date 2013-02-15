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
      bool getUniqueAttributeValues(const std::string& tagName,const std::string& attribName,std::set<std::string>& output) const;
      bool multiReadStart(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
      bool multiReadAddUnit(const uint64_t& amount,char* buffer);
      bool multiReadEnd(const uint64_t& offset);
      bool open(const std::string& fname,MPI_Comm comm,const int& masterRank,MPI_Info mpiInfo);
      bool readArrayMaster(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			   const uint64_t& begin,const uint64_t& amount,char* buffer);
      bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
		     const uint64_t& begin,const uint64_t& amount,char* buffer);
      template<typename T>
      bool readParameter(const std::string& parameterName,T& value);
   
    private:
      MPI_Comm comm;                  /**< MPI communicator used to read the file.*/
      MPI_File filePtr;               /**< MPI file pointer to input file.*/
      int masterRank;                 /**< MPI rank of master process.*/
      bool multireadStarted;          /**< If true, multiread mode has been initialized successfully.*/
      int myRank;                     /**< MPI rank of this process in communicator comm.*/
      int processes;                  /**< Number of MPI processes in communicator comm.*/
   
      bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
      
      MPI_Datatype multiReadVectorType;
      MPI_Datatype multiReadArrayType;
      std::list<Multi_IO_Unit> multiReadUnits;
   };


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
      if (masterSuccess > 0) success = true;
      return success;
   }

}
   
#endif
