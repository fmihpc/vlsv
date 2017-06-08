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

#include "portable_file_io.h"
#include "vlsv_reader.h"

using namespace std;

namespace vlsv {

   Reader::Reader() {
      endiannessReader = detectEndianness();
      fileOpen = false;
      swapIntEndianness = false;
   }

   Reader::~Reader() {
      filein.close();   
   }
   
   bool Reader::close() {
      filein.close();
      xmlReader.clear();
      fileOpen = false;
      return true;
   }

   /** Get attributes of the given XML tag.
    * @param tagName Name of the XML tag.
    * @param attribsIn Constraints that limit the search.
    * @param attribsOut Attributes of the XML tag, if one matched given constraints.
    * @return If true, an XML tag was found that mathes given constraints.*/
   bool Reader::getArrayAttributes(const string& tagName,const list<pair<string,string> >& attribsIn,map<string,string>& attribsOut) const {
      if (fileOpen == false) return false;
      muxml::XMLNode* node = xmlReader.find(tagName,attribsIn);
      if (node == NULL) return false;
      attribsOut = node->attributes;   
      return true;
   }

   /** Get metadata of given array.
    * @param tagName Name of the XML tag.
    * @param attribs Constraints that limit search.
    * @param arraySize Variable in which array size is written.
    * @param vectorSize Variable in which vector size of each array element is written.
    * @param dataType Variable in which the datatype stored to array is written.
    * @param dataSize Variable in which byte size of the datatype stored to array is written.
    * @return If true, an array was found that matched given search criteria and output variables 
    * contain meaningful values.*/
   bool Reader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
                             uint64_t& arraySize,uint64_t& vectorSize,datatype::type& dataType,uint64_t& dataSize) {
      if (fileOpen == false) return false;
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) return false;
   
      arraySize = atol(node->attributes["arraysize"].c_str());
      vectorSize = atol(node->attributes["vectorsize"].c_str());
      dataSize = atol(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "unknown") dataType = datatype::UNKNOWN;
      else if (node->attributes["datatype"] == "int") dataType = datatype::INT;
      else if (node->attributes["datatype"] == "uint") dataType = datatype::UINT;
      else if (node->attributes["datatype"] == "float") dataType = datatype::FLOAT;
      else {
         cerr << "vlsv::Reader ERROR: Unknown datatype '" << node->attributes["datatype"] << "' in tag!" << endl;
         return false;
      }
      return true;
   }

   const std::string Reader::getErrorString() const {
      return vlsv::getErrorString(lastErrorCode);
   }
   
   bool Reader::getFileName(std::string& openFile) const {
      if (fileOpen == false) {
         openFile = "";
      } else {
         openFile = fileName;
      }
      return fileOpen;
   }

   /** Get unique values of given XML tag attribute. This function can be used to query the names of 
    * all mesh variables, for example.
    * @param tagName Name of the XML tag whose attributes are included in the search.
    * @param attribName Name of the queried tag attribute.
    * @param output Set in which unique attribute values are written.
    * @return If true, output variable contain meaningful values.*/
   bool Reader::getUniqueAttributeValues(const string& tagName,const string& attribName,set<string>& output) const {
      if (fileOpen == false) return false;

      muxml::XMLNode* node = xmlReader.find("VLSV");
      for (multimap<string,muxml::XMLNode*>::const_iterator it=node->children.lower_bound(tagName); 
           it!=node->children.upper_bound(tagName); ++it) {
         map<string,string>::const_iterator tmp = it->second->attributes.find(attribName);
         if (tmp == it->second->attributes.end()) continue;
         output.insert(tmp->second);
      }
      return true;
   }

   bool Reader::loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
      if (fileOpen == false) return false;
   
      // Find tag corresponding to given array:
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) return false;

      // Copy array information from tag:
      arrayOpen.offset = atol(node->value.c_str());
      arrayOpen.tagName = tagName;
      arrayOpen.arraySize = atol(node->attributes["arraysize"].c_str());
      arrayOpen.vectorSize = atol(node->attributes["vectorsize"].c_str());
      arrayOpen.dataSize = atol(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "unknown") arrayOpen.dataType = datatype::UNKNOWN;
      else if (node->attributes["datatype"] == "int") arrayOpen.dataType = datatype::INT;
      else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = datatype::UINT;
      else if (node->attributes["datatype"] == "float") arrayOpen.dataType = datatype::FLOAT;
      else {
         cerr << "vlsv::Reader ERROR: Unknown datatype in tag!" << endl;
         return false;
      }   
      //if (arrayOpen.arraySize == 0) return false;
      if (arrayOpen.vectorSize == 0) return false;
      if (arrayOpen.dataSize == 0) return false;
   
      return true;
   }

   /** Open a VLSV file for reading. This function fails if a 
    * file is already open. 
    * @param fname File name.
    * @return If true, file was successfully opened.*/
   bool Reader::open(const std::string& fname) {
      bool success = true;
      lastErrorCode = error::NONE;
      if (fileOpen == true) {
         lastErrorCode = error::READ_FILE_ALREADY_OPEN;
         return false;
      }

      // If given filename includes path, chdir into that path:
      string fnameWithoutPath;
      string pathName;
      char cwd[1024];
      if (fileio::getcwd(cwd,sizeof(cwd)) != NULL) {
         const size_t position = fname.find_last_of("/");
         if (position == string::npos) {
            pathName = "";
            fnameWithoutPath = fname;
         } else {
            pathName = fname.substr(0,position);
            fnameWithoutPath = fname.substr(position+1);
            // Chdir to path containing input file and attempt to open the file, 
            // then chdir back to current working directory. 
            // Chdir returns zero value if it succeeds
            // Not done if the string is empty as chdir fials in that case.
            if (fileio::chdir(pathName.c_str()) != 0) {
               lastErrorCode = error::READ_CWD_FAIL;
               success = false;
            }
         }
      }

      filein.open(fnameWithoutPath.c_str(), fstream::in | fstream::binary);
      if (fileio::chdir(cwd) != 0) success = false;

      if (filein.good() == true) {
         fileName = fnameWithoutPath;
         fileOpen = true;
      } else {
         filein.close();
         success = false;
         lastErrorCode = error::READ_FILE_BAD;
      }

      if (success == false) return success;

      // Detect file endianness:
      char* ptr = reinterpret_cast<char*>(&endiannessFile);
      filein.read(ptr,1);
      if (filein.good() == false) {
         success = false;
         lastErrorCode = error::READ_FILE_ENDIANNESS;
         return success;
      }
      if (endiannessFile != endiannessReader) swapIntEndianness = true;

      // Read footer offset:
      uint64_t footerOffset;
      char buffer[16];
      filein.seekg(8);
      filein.read(buffer,8);
      if (filein.good() == false) {
         lastErrorCode = error::READ_FOOTER_OFFSET;
         success = false;
         return success;
      }
      footerOffset = convUInt64(buffer,swapIntEndianness);

      // Read footer XML tree:
      filein.seekg(footerOffset);
      if (filein.tellg() != footerOffset) {
         lastErrorCode = error::READ_NO_FOOTER;
         success = false;
      }
      
      if (success == true && xmlReader.read(filein) == false) {
         lastErrorCode = error::READ_FOOTER;
         success = false;
      }
      filein.clear();
      filein.seekg(16);

      return success;
   }

   /** Read given part of a given array from file.
    * @param tagName Name of the XML tag.
    * @param attribs List of attributes that uniquely determine the array.
    * @param begin Index of the first read array element.
    * @param amount How many array elements are read.
    * @param buffer Buffer in which data is copied.
    * @return If true, array was found and requested part was copied to buffer.*/
   bool Reader::readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			  const uint64_t& begin,const uint64_t& amount,char* buffer) {
      if (fileOpen == false) {
         cerr << "vlsv::Reader ERROR: readArray called but a file is not open!" << endl;
         return false;
      }
   
      // If zero-length read was requested, exit immediately:
      if (amount == 0) return true;
      
      // Find tag corresponding to given array:
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
         cerr << "vlsv::Reader ERROR: Failed to find tag='" << tagName << "' attribs:" << endl;
         for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
            cerr << '\t' << it->first << " = '" << it->second << "'" << endl;
         }
         return false;
      }
      
      // Copy array information from tag:
      arrayOpen.offset = atol(node->value.c_str());
      arrayOpen.tagName = tagName;
      arrayOpen.arraySize = atol(node->attributes["arraysize"].c_str());
      arrayOpen.vectorSize = atol(node->attributes["vectorsize"].c_str());
      arrayOpen.dataSize = atol(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "int") arrayOpen.dataType = datatype::INT;
      else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = datatype::UINT;
      else if (node->attributes["datatype"] == "float") arrayOpen.dataType = datatype::FLOAT;
      else {
         cerr << "vlsv::Reader ERROR: Unknown datatype in tag!" << endl;
         return false;
      }
   
      if (arrayOpen.arraySize == 0) return false;
      if (arrayOpen.vectorSize == 0) return false;
      if (arrayOpen.dataSize == 0) return false;
      
      // Sanity check on values:
      if (begin + amount > arrayOpen.arraySize) {
         cerr << "vlsv::Reader ERROR: Requested read exceeds array size. begin: " << begin;
         cerr << " amount: " << amount << " size: " << arrayOpen.arraySize << endl;
         return false;
      }

      // Read data from file:
      streamoff start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
      streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
      filein.clear();
      filein.seekg(start);
      filein.read(buffer,readBytes);
      
      // Check that we were able to read the requested amount of data:
      if (filein.gcount() != readBytes) {
         cerr << "vlsv::Reader ERROR: Failed to read requested amount of bytes!" << endl;      
         cerr << "tag name='" << tagName << "'" << endl;
         cerr << "attributes:" << endl;
         for (map<string,string>::const_iterator it=node->attributes.begin(); it!=node->attributes.end(); ++it) {
            cerr << '\t' << it->first << " = " << it->second << endl;
         }
         cerr << "array offset string '" << node->value.c_str() << "'" << endl;
         cerr << "start=" << start << " readBytes=" << readBytes << endl;
         cerr << "offset=" << arrayOpen.offset << " vectorsize=" << arrayOpen.vectorSize << " dataSize=" << arrayOpen.dataSize << endl;
         return false;
      }
      return true;
   }

} // namespace vlsv
