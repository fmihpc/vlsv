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

#ifndef VLSV_WRITER_H
#define VLSV_WRITER_H

#include <stdint.h>
#include <mpi.h>
#include <limits>

#include "muxml.h"
#include "mpiconversion.h"
#include "vlsv_common.h"
#include "multi_io_unit.h"

/** VLSV file format writer.
 * 
 * List of currently supported XML attributes (vlsv2silo knows what to do with these):
 * MESH (Quad,Point):
 * xlabel (string)               Name of x-coordinate, displayed in VisIt coordinate labels.
 * ylabel (string)               Name of y-coordinate, displayed in VisIt coordinate labels.
 * zlabel (string)               Name of z-coordinate, displayed in VisIt coordinate labels.
 * xunit (string)                Unit of x-coordinate, displayed in VisIt coordinate labels.
 * yunit (string)                Unit of y-coordinate, displayed in VisIt coordinate labels.
 * zunit (string)                Unit of z-coordinate, displayed in VisIt coordinate labels.
 * xscaling (float)              Scale factor for x-coordinate. For example value "1.0e-3" converts meters into kilometers.
 * yscaling (float)              Scale factor for y-coordinate. For example value "1.0e-3" converts meters into kilometers.
 * zscaling (float)              Scale factor for z-coordinate. For example value "1.0e-3" converts meters into kilometers.
 * meshinfo (string)             If defined, mesh information (labels,units,scaling) is read from a mesh with given name.
 * 
 * PARAMETER:
 * These are used to pass parameters (=single values) to VisIt. These should only be written by the master process.
 * time                          Simulation time of the VLSV file, VisIt shows this value as 'Time'.
 * timestep                      Current value of time step, VisIt shows this value as 'Timestep'.
 */

namespace vlsv {

   class Writer {
    public:
      Writer();
      ~Writer();

      bool addMultiwriteUnit(char* array,const uint64_t& arrayElements);
      bool close();
      uint64_t getBytesWritten() const;
      bool endMultiwrite(const std::string& tagName,const std::map<std::string,std::string>& attribs);
      bool open(const std::string& fname,MPI_Comm comm,const int& masterProcessID);
      bool startMultiwrite(const std::string& datatype,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize);
      bool writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,const std::string& dataType,
		      const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,const char* array);
   
      // ***** TEMPLATE WRAPPER FUNCTIONS ***** //

      template<typename T> 
      bool addMultiwriteUnit(const T* array,const uint64_t& arrayElements);
      
      template<typename T>
      bool startMultiwrite(const uint64_t& arraySize,const uint64_t& vectorSize);
      
      template<typename T> 
      bool writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,
		      const uint64_t& arraySize,const uint64_t& vectorSize,const T* array);
      
      template<typename T>
      bool writeParameter(const std::string& parameterName,const T* const array);
      
      template<typename T>
      bool writeWithReduction(const std::string& arrayName,const std::map<std::string,std::string>& attribs,
			      const uint64_t& arraySize,T* array,MPI_Op operation);
   
    private:

      uint64_t arraySize;                     /**< Number of array elements this process will write.*/
      int* blockLengths;                      /**< Used in creation of an MPI_Struct in endMultiwrite.*/
      uint64_t* bytesPerProcess;              /**< Array with N_processes elements. Used to gather myBytes.*/
      uint64_t bytesWritten;                  /**< Total amount of bytes written to output file,
					       * significant at master process only.*/
      MPI_Comm comm;                          /**< MPI communicator used in I/O.*/
      uint64_t dataSize;                      /**< Byte size of each element in data vector, must have
					       * the same value on all participating processes.*/
      std::string dataType;                   /**< String description of the datatype that is written to file,
					       * obtained by calling arrayDataType() template function.*/
      MPI_Aint* displacements;                /**< Used in creation of an MPI_Struct in endMultiwrite.*/
      unsigned int endMultiwriteCounter;      /**< A counter used in endMultiwrite to synchronize threads.*/
      std::string fileName;                   /**< Name of the output file.*/
      bool fileOpen;                          /**< If true, a file has been successfully opened for writing.*/
      MPI_File fileptr;                       /**< MPI file pointer to the output file.*/
      bool initialized;                       /**< If true, VLSV Writer initialization is complete, does not tell if it was successful.*/
      int masterRank;                         /**< Rank of master process in communicator comm.*/
      bool multiwriteFinalized;               /**< If true, multiwrite array writing mode has finalized correctly. 
					       This variable is used to synchronize threads in endMultiwrite function..*/
      bool multiwriteInitialized;             /**< If true, multiwrite array writing mode has initialized correctly. 
					       This variable is used to synchronize threads in startMultiwrite function.*/
      
      std::vector<unsigned int> multiwriteOffsets; /**< Offset for each thread using VLSVWriter, used to load 
						    * data into an MPI struct in endMultiwrite.*/
      char* multiwriteOffsetPointer;          /**< Pointer to an array that is used to calculate offsets in an 
					       * MPI struct created in endMultiwrite.*/
      std::vector<std::list<Multi_IO_Unit> > multiwriteUnits; /**< Container for all multiwrite units for this process. 
							       * Each thread using VLSVWriter has its own list. This 
							       * allows vlsv::Writer::addMultiwriteUnit to be called without 
							       * thread synchronizations.*/   
      uint64_t myBytes;                       /**< Number of bytes this process is writing to the current array.*/
      int myrank;                             /**< Rank of this process in communicator comm.*/
      unsigned int N_multiwriteUnits;         /**< Total number of multiwrite units this process has. In multithreaded mode 
					       * this is equal to the sum of multiwrite units over all threads.*/
      int N_processes;                        /**< Number of processes in communicator comm.*/
      MPI_Offset offset;                      /**< MPI offset into output file for this process.*/
      MPI_Offset* offsets;                    /**< Array with N_processes elements. Used to scatter file offsets.*/
      MPI_Datatype* types;                    /**< Used in creation of an MPI_Struct in endMultiwrite.*/
      uint64_t vectorSize;                    /**< Number of elements in each data vector per array element,
					       * must have the same value on all participating processes.*/
      datatype::type vlsvType;                /**< Same as dataType but in an integer representation.*/
      muxml::MuXML* xmlWriter;                /**< Pointer to XML writer, used for writing a footer to the VLSV file.*/
   };

   template<typename T> inline
   bool Writer::addMultiwriteUnit(const T* array,const uint64_t& arrayElements) {
      // Check that startMultiwrite has initialized correctly:
      if (multiwriteInitialized == false) return false;
   
      // Cast away const-ness:
      T* arrayPtr = const_cast<T*>(array);
   
      // Each thread records their multiwrite units to per-thread storage,
      // so there is no need to synchronize access to vector multiwriteUnits:
      multiwriteUnits[0].push_back(Multi_IO_Unit(reinterpret_cast<char*>(arrayPtr),MPI_Type<T>(),arrayElements*vectorSize));
      return true;
   }

   /** Start an array writing process.
    * @param arraySize  Number of elements this MPI process will write to the output array. 
    * @param vectorSize Number of elements in each data vector, this value must have the 
    * same value on all participating MPI processes.
    * @return If true, multiwrite mode has initialized correctly and functions
    * vlsv::Writer::addMultiwriteUnit and VLSVWriter::endMultiwrite may be called.*/
   template<typename T> inline
   bool Writer::startMultiwrite(const uint64_t& arraySize,const uint64_t& vectorSize) {
      return startMultiwrite(getStringDatatype<T>(),arraySize,vectorSize,sizeof(T));
   }

   /** Write an array to the output file. This function is simply a wrapper to 
    * multiwrite functions, i.e. array is written to file with a single multiwrite unit.
    * @param tagName Name of the array, same as the XML tag name in output file.
    * @param attribs Other attributes for the output XML tag, given in [tag name,tag value] pairs.
    * @param arraySize Number of elements in array.
    * @param vectorSize Number of elements in vectors that comprise the array elements.
    * @param array Pointer to the output array.
    * @return If true, the array was successfully written to file.*/
   template<typename T> inline
   bool Writer::writeArray(const std::string& tagName,const std::map<std::string,std::string>& attribs,
			   const uint64_t& arraySize,const uint64_t& vectorSize,const T* array) {
      // Cast away const-ness of array pointer (required by MPI):
      T* arrayPtr = const_cast<T*>(array);
      return writeArray(tagName,attribs,getStringDatatype<T>(),arraySize,vectorSize,sizeof(T),reinterpret_cast<char*>(arrayPtr));
   }

   template<typename T> inline
   bool Writer::writeParameter(const std::string& parameterName,const T* const array) {
      std::map<std::string,std::string> attributes;
      attributes["name"] = parameterName;
   
      if (myrank == masterRank)
	return writeArray("PARAMETER",attributes,1,1,array);
      else
	return writeArray("PARAMETER",attributes,0,0,array);
   }

   template<typename T> inline
   bool Writer::writeWithReduction(const std::string& arrayName,const std::map<std::string,std::string>& attribs,
				   const uint64_t& arraySize,T* array,MPI_Op operation) {
      // Master process allocates a receive buffer for reduction:
      T* recvBuffer = NULL;
      if (myrank == masterRank) recvBuffer = new T[arraySize];

      // Reduce result to master:
      if (MPI_Reduce(array,recvBuffer,arraySize,MPI_Type<T>(),operation,masterRank,comm) != MPI_SUCCESS) {
	 delete [] recvBuffer; return false;
      }

      // Write result to file. Only master process has a non-zero array length, 
      // all other processes write a zero-length array:
      if (myrank == masterRank) {
	 writeArray(arrayName,attribs,1,arraySize,recvBuffer);
      } else {
	 writeArray(arrayName,attribs,0,0,recvBuffer);
      }
   
      delete [] recvBuffer; recvBuffer = NULL;
      return true;
   }

} // namespace vlsv
   
#endif
