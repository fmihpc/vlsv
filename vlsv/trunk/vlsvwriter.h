#ifndef VLSVWRITER2_H
#define VLSVWRITER2_H

#include <stdint.h>
#include <mpi.h>
#include <limits>
#ifdef THREADING
   #include <pthread.h>
#endif

#include "muxml.h"
#include "mpiconversion.h"
#include "vlsv_common.h"

namespace VLSV {
   struct WriteUnit {
      char* array;              /**< Pointer to data to be written.*/
      MPI_Datatype mpiType;     /**< MPI datatype of data that is written.*/
      uint64_t amount;          /**< How many elements are to be written.*/
      
      WriteUnit(char* array,const MPI_Datatype& mpiType,const uint64_t& amount);
   
    private:
   
      /** Private default constructor to prevent creation of empty multiwrite units.*/
      WriteUnit();
   };
}

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
class VLSVWriter {
 public:
   VLSVWriter();
   ~VLSVWriter();

   bool addMultiwriteUnit(char* array,const uint64_t& arrayElements,const int& threadID=0);
   bool close(const int& threadID=0);
   bool endMultiwrite(const std::string& tagName,const std::map<std::string,std::string>& attribs,const int& threadID=0);
   bool open(const std::string& fname,MPI_Comm comm,const int& masterProcessID,
	     const int& mpiThreadingLevel=MPI_THREAD_SINGLE,const int& N_threads=1,const int& threadID=0,const int& masterThreadID=0);
   bool startMultiwrite(const std::string& datatype,const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,const int& threadID);
   bool writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,const std::string& dataType,
		   const uint64_t& arraySize,const uint64_t& vectorSize,const uint64_t& dataSize,char* array,const int& threadID=0);
   
   // ***** TEMPLATE WRAPPER FUNCTIONS ***** //

   template<typename T> 
   bool addMultiwriteUnit(T* array,const uint64_t& arrayElements,const int& threadID=0);

   template<typename T>
   bool startMultiwrite(const uint64_t& arraySize,const uint64_t& vectorSize,const int& threadID=0);
   
   template<typename T> 
   bool writeArray(const std::string& arrayName,const std::map<std::string,std::string>& attribs,
		   const uint64_t& arraySize,const uint64_t& vectorSize,T* array,const int& threadID=0);

 private:

   uint64_t arraySize;                     /**< Number of array elements this process will write.*/
   int* blockLengths;                      /**< Used in creation of an MPI_Struct in endMultiwrite.*/
   uint64_t* bytesPerProcess;              /**< Array with N_processes elements. Used to gather myBytes.*/
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
   int masterThreadID;                     /**< ID of the master thread. In certain MPI threading levels
					    * only the master thread is allowed to call MPI routines.*/
   int mpiThreadingLevel;                  /**< Threading level supported by the underlying MPI library,
					    * obtained from MPI_Init_thread.*/
   bool multiwriteFinalized;               /**< If true, multiwrite array writing mode has finalized correctly. 
					    This variable is used to synchronize threads in endMultiwrite function..*/
   bool multiwriteInitialized;             /**< If true, multiwrite array writing mode has initialized correctly. 
					    This variable is used to synchronize threads in startMultiwrite function.*/

   std::vector<unsigned int> multiwriteOffsets; /**< Offset for each thread using VLSVWriter, used to load 
						 * data into an MPI struct in endMultiwrite.*/
   char* multiwriteOffsetPointer;          /**< Pointer to an array that is used to calculate offsets in an 
					    * MPI struct created in endMultiwrite.*/
   std::vector<std::list<VLSV::WriteUnit> > multiwriteUnits; /**< Container for all multiwrite units for this process. 
							      * Each thread using VLSVWriter has its own list. This 
							      * allows VLSVWriter::addMultiwriteUnit to be called without 
							      * thread synchronizations.*/   
   uint64_t myBytes;                       /**< Number of bytes this process is writing to the current array.*/
   int myrank;                             /**< Rank of this process in communicator comm.*/
   unsigned int N_multiwriteUnits;         /**< Total number of multiwrite units this process has. In multithreaded mode 
					    * this is equal to the sum of multiwrite units over all threads.*/
   int N_processes;                        /**< Number of processes in communicator comm.*/
   int N_threads;                          /**< Total number of threads using VLSV writer.*/
   MPI_Offset offset;                      /**< MPI offset into output file for this process.*/
   MPI_Offset* offsets;                    /**< Array with N_processes elements. Used to scatter file offsets.*/
   MPI_Datatype* types;                    /**< Used in creation of an MPI_Struct in endMultiwrite.*/
   uint64_t vectorSize;                    /**< Number of elements in each data vector per array element,
					    * must have the same value on all participating processes.*/
   VLSV::datatype vlsvType;                /**< Same as dataType but in an integer representation.*/
   MuXML* xmlWriter;                       /**< Pointer to XML writer, used for writing a footer to the VLSV file.*/
   
   // ***** VARIABLES ONLY USED IN MULTITHREADED MODE *****
   
   #ifdef THREADING
   pthread_barrier_t  barrier;             /**< Barrier, used to sync all threads using VLSVWriter.*/
   pthread_cond_t     closeCond;           /**< Conditional lock used in VLSVWriter::close to sync threads.*/
   pthread_mutex_t    closeLock;           /**< Mutex used in VLSVWriter::close to sync threads.
					    * Used together with VLSVWriter::closeCond.*/
   pthread_cond_t     endMultiwriteCond;   /**< Conditional lock used in VLSVWriter::endMultiwrite to sync threads.*/
   pthread_mutex_t    endMultiwriteLock;   /**< Mutex used in VLSVWriter::endMultiwrite to sync threads.
					    * Used together with VLSVWriter::endMultiwriteCond.*/
   pthread_cond_t     multiwriteStartCond; /**< Conditional lock used in VLSVWriter::multiwriteStart to sync threads.*/
   pthread_mutex_t    multiwriteStartLock; /**< Mutex used in VLSVWriter::multiwriteStart to sync threads. 
					    *Used together with VLSVWriter::multiwriteStartCond.*/
   pthread_cond_t     openCond;            /**< Conditional lock used in VLSVWriter::open to sync threads.*/
   pthread_mutex_t    openLock;            /**< Mutex used in VLSVWriter::open to sync threads.
					    * Used together with VLSVWriter::openCond.*/
   #endif
   
   /** Returns a string representation of an array that is to be written to file.
    * The correct C++ datatype can be deduced from the string value returned by this
    * function and from the byte size of the data in array. For example, datatype "float"
    * with byte size of 8 usually means that the values are doubles.
    * @return String representation of the datatype.
    */
   template<typename T> std::string arrayDataType();
   
   MPI_Datatype getMPIDatatype(VLSV::datatype dt,uint64_t dataSize);
   VLSV::datatype getVLSVDatatype(const std::string& s);
};

template<> inline std::string VLSVWriter::arrayDataType<int8_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int16_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int32_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<int64_t>() {return "int";}
template<> inline std::string VLSVWriter::arrayDataType<uint8_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint16_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint32_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<uint64_t>() {return "uint";}
template<> inline std::string VLSVWriter::arrayDataType<float>() {return "float";}
template<> inline std::string VLSVWriter::arrayDataType<double>() {return "float";}
template<> inline std::string VLSVWriter::arrayDataType<long double>() {return "float";}

template<typename T>
inline bool VLSVWriter::addMultiwriteUnit(T* array,const uint64_t& arrayElements,const int& threadID) {
   // Check that startMultiwrite has initialized correctly:
   if (multiwriteInitialized == false) return false;
   
   // Each thread records their multiwrite units to per-thread storage,
   // so there is no need to synchronize access to vector multiwriteUnits:
   multiwriteUnits[threadID].push_back(VLSV::WriteUnit(reinterpret_cast<char*>(array),MPI_Type<T>(),arrayElements*vectorSize));
   return true;
}

/** Start an array writing process.
 * @param arraySize  Number of elements this MPI process will write to the output array. 
 * In multithreaded mode all threads must call this function with the same value.
 * @param vectorSize Number of elements in each data vector, this value must have the 
 * same value on all participating MPI processes. In multithreaded mode all threads 
 * must call this function with the same value.
 * @param threadID Thread ID of the thread calling this function. If multithreaded mode
 * is not used, master thread ID should be used here. Defaults to value zero.
 * @return If true, multiwrite mode has initialized correctly and functions
 * VLSVWriter::addMultiwriteUnit and VLSVWriter::endMultiwrite may be called. In 
 * multithreaded mode every thread calling this function will return the same value.
 */
template<typename T>
inline bool VLSVWriter::startMultiwrite(const uint64_t& arraySize,const uint64_t& vectorSize,const int& threadID) {
   return startMultiwrite(arrayDataType<T>(),arraySize,vectorSize,sizeof(T),threadID);
}

/** Write an array to the output file. This function is simply a wrapper to 
 * multiwrite functions, i.e. array is written to file with a single multiwrite unit.
 * In multithreaded mode it is safe to call this function with all threads.
 * @param tagName Name of the array, same as the XML tag name in output file.
 * @param attribs Other attributes for the output XML tag, given in [tag name,tag value] pairs.
 * @param arraySize Number of elements in array.
 * @param vectorSize Number of elements in vectors that comprise the array elements.
 * @param array Pointer to the output array.
 * @param threadID Thread ID of the thread that calls this function. If multithreaded mode is 
 * not used, master thread ID should be used here. Defaults to value zero.
 * @return If true, the array was successfully written to file.
 */
template<typename T> 
inline bool VLSVWriter::writeArray(const std::string& tagName,const std::map<std::string,std::string>& attribs,
				   const uint64_t& arraySize,const uint64_t& vectorSize,T* array,const int& threadID) {
   return writeArray(tagName,attribs,arrayDataType<T>(),arraySize,vectorSize,sizeof(T),reinterpret_cast<char*>(array),threadID);
}

#endif
