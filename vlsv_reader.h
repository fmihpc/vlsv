/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013,2015 Finnish Meteorological Institute
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

#ifndef VLSV_READER_H
#define VLSV_READER_H

#include <stdint.h>
#include <list>
#include <set>
#include <fstream>

#include "muxml.h"
#include "vlsv_common.h"

namespace vlsv {

   class Reader {
    public:
      Reader();
      virtual ~Reader();
   
      virtual bool close();
      virtual bool getArrayAttributes(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribsIn,
                                      std::map<std::string,std::string>& attribsOut) const;
      virtual bool getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                                uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& byteSize);
      virtual const std::string getErrorString() const;
      virtual bool getFileName(std::string& openFile) const;
      virtual bool getUniqueAttributeValues(const std::string& tagName,const std::string& attribName,std::set<std::string>& output) const;
      virtual bool loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs);
      virtual bool open(const std::string& fname);
      virtual bool readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                             const uint64_t& begin,const uint64_t& amount,char* buffer);

      template<typename T>
      bool read(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                const uint64_t& begin,const uint64_t& amount,T*& buffer,bool allocateMemory=true);
      template<typename T>
      bool readParameter(const std::string& parameterName,T& value);
   
    protected:
      unsigned char endiannessFile;   /**< Endianness in VLSV file.*/
      unsigned char endiannessReader; /**< Endianness of computer which reads the data.*/
      error::type lastErrorCode;      /**< Code indicating last error that has occurred, if any.*/
      std::fstream filein;            /**< Input file stream.*/
      std::string fileName;           /**< Name of the input file.*/
      bool fileOpen;                  /**< If true, a file is currently open.*/
      bool swapIntEndianness;         /**< If true, endianness should be swapped on read data (not implemented yet).*/
      muxml::MuXML xmlReader;         /**< XML reader used to parse VLSV footer.*/
   
      /** Struct used to store information on the currently open array.*/
      struct ArrayOpen {
         std::streamoff offset;
         std::string tagName;
         std::string arrayName;
         datatype::type dataType;
         uint64_t arraySize;
         uint64_t vectorSize;
         uint64_t dataSize;
      } arrayOpen;
   };

   template<typename T> inline
   bool Reader::read(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                     const uint64_t& begin,const uint64_t& amount,T*& outBuffer,bool allocateMemory) {
      // Don't touch outBuffer if it has already been allocated:
      if (allocateMemory == true) outBuffer = NULL;
      if (amount == 0) return true;

      // Get array info:
      uint64_t arraySize;
      uint64_t vectorSize;
      datatype::type datatype;
      uint64_t dataSize;
      if (Reader::getArrayInfo(tagName,attribs,arraySize,vectorSize,datatype,dataSize) == false) {
         std::cerr << "vlsv::Reader failed to get array info" << std::endl;
         return false;
      }

      // Check that requested read is inside the array:
      if (begin > arraySize || (begin+amount) > arraySize) return false;

      // Read data into temporary buffer:
      char* buffer = new char[amount*vectorSize*dataSize];
      if (Reader::readArray(tagName,attribs,begin,amount,buffer) == false) {
         delete [] buffer; buffer = NULL;
         return false;
      }

      // Copy data from temporary buffer to output:
      if (allocateMemory == true) outBuffer = new T[amount*vectorSize];
      char* ptr = buffer;
      for (uint64_t i=0; i<amount; ++i) {
         for (uint64_t j=0; j<vectorSize; ++j) {
            convertValue<T>(outBuffer[i*vectorSize+j],ptr,datatype,dataSize,false);	 
            ptr += dataSize;
         }
      }

      delete [] buffer; buffer = NULL;
      return true;
   }

   template<typename T> inline
   bool Reader::readParameter(const std::string& parameterName,T& value) {
      std::list<std::pair<std::string,std::string> > attribs;
      const std::string name = "name";
      attribs.push_back(make_pair(name,parameterName));
      uint64_t arraySize;
      uint64_t vectorSize;
      datatype::type dataType;
      uint64_t dataSize;
      
      // Get array info containing parameter value:
      bool success = Reader::getArrayInfo("PARAMETER",attribs,arraySize,vectorSize,dataType,dataSize);
      if (success == false) {
         return success;
      }
      
      // Check that the array contains a single parameter value:
      if (arraySize != 1) return false;
      if (vectorSize != 1) return false;
      
      // Read parameter value into temporary buffer:
      char* buffer = new char[arraySize*vectorSize*dataSize];
      success = Reader::readArray("PARAMETER",attribs,0,1,buffer);

      // Convert read value into requested datatype:
      convertValue<T>(value,buffer,dataType,dataSize,false);
      delete [] buffer; buffer = NULL;
      return success;
   }

} // namespace vlsv
   
#endif
