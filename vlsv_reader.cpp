/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2013, 2015 Finnish Meteorological Institute
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

   bool Reader::withError(const std::stringstream& ss) {
      errorString = ss.str();
      #ifndef NDEBUG
         cerr << errorString << endl;
      #endif
      return false;
   }

   /** Default constructor. Attempts to detect the endianness of the system.*/
   Reader::Reader() {
      endiannessReader = detectEndianness();
      fileOpen = false;
      swapIntEndianness = false;
   }

   /** Destructor. Closes the input file if one is open.*/
   Reader::~Reader() {
      filein.close();   
   }
   
   /** Close the input file if one is open.
    * @return Always returns true.*/
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
    * @return If true, an XML tag was found that matches given constraints.*/
   bool Reader::getArrayAttributes(const string& tagName,const list<pair<string,string> >& attribsIn,map<string,string>& attribsOut) {
      // Check that an input file is open:
      if (fileOpen == false) {
         stringstream ss;
         ss << "ERROR: Input file is not open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      muxml::XMLNode* node = xmlReader.find(tagName,attribsIn);
      if (node == NULL) {
         stringstream ss;
         ss << "ERROR: Could not find requested array '" << tagName << "' in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
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
      // Check that input file is open:
      if (fileOpen == false) {
         stringstream ss;
         ss << "ERROR: Input file is not open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
         stringstream ss;
         ss << "ERROR: Could not find requested array '" << tagName << "' in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }

      arraySize = atol(node->attributes["arraysize"].c_str());
      vectorSize = atol(node->attributes["vectorsize"].c_str());
      dataSize = atol(node->attributes["datasize"].c_str());
      if (node->attributes["datatype"] == "unknown") dataType = datatype::UNKNOWN;
      else if (node->attributes["datatype"] == "int") dataType = datatype::INT;
      else if (node->attributes["datatype"] == "uint") dataType = datatype::UINT;
      else if (node->attributes["datatype"] == "float") dataType = datatype::FLOAT;
      else {
         stringstream ss;
         ss << "ERROR: Unknown datatype '" << node->attributes["datatype"] << "' detected in ";
         ss << __FILE__ << ":" << __LINE__ << endl;
         return withError(ss);
      }
      return true;
   }

   /** Copy the name of the current input file to the given string variable.
    * If a file is not open, then an empty string "" is copied.
    * @param openFile Name of the input file is copied here.
    * @return If true, a file is currently open. Otherwise false is returned.*/
   bool Reader::getFileName(std::string& openFile) const {
      if (fileOpen == false) {
         openFile = "";
      } else {
         openFile = fileName;
      }
      return fileOpen;
   }

   /** Return a string describing the latest error, if any, that has occurred.
    * @return String describing the latest error.*/
   std::string Reader::getLastError() const {return errorString;}

   /** Get unique values of given XML tag attribute. This function can be used to query the names of 
    * all mesh variables, for example.
    * @param tagName Name of the XML tag whose attributes are included in the search.
    * @param attribName Name of the queried tag attribute.
    * @param output Set in which unique attribute values are written.
    * @return If true, output variable contain meaningful values.*/
   bool Reader::getUniqueAttributeValues(const string& tagName,const string& attribName,set<string>& output) {
      // Check that an input file is open:
      if (fileOpen == false) {
         stringstream ss;
         ss << "ERROR: Input file is not open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }

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
      // Check that a file is open:
      if (fileOpen == false) {
         stringstream ss;
         ss << "ERROR: file not open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }

      // Find tag corresponding to given array:
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
         stringstream ss;
         ss << "ERROR: Could not find requested array '" << tagName << "' in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }

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
         stringstream ss;
         ss << "ERROR: Unknown datatype in XML tag in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      if (arrayOpen.vectorSize == 0) {
         stringstream ss;
         ss << "ERROR: array has zero vector size in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      if (arrayOpen.dataSize == 0) {
         stringstream ss;
         ss << "ERROR: array has zero data size in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      return true;
   }

   /** Open a VLSV file for reading. This function fails if a 
    * file is already open, or if the file endianness or footer 
    * could not be read.
    * @param fname File name.
    * @return If true, file was successfully opened.*/
   bool Reader::open(const std::string& fname) {
      // Check that another file isn't already open:
      if (fileOpen == true) {
         stringstream ss;
         ss << "ERROR: File '" << fname << "' should be opened, but file '";
         ss << fileName << "' is currently open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      fileOpen = false;
   
      // If given filename includes path, chdir into that path:
      string fnameWithoutPath;
      string pathName;
      char cwd[1024];
      if (fileio::getcwd(cwd,sizeof(cwd)) != NULL) {
         #ifdef WINDOWS
            const size_t position = fname.find_last_of("\\");
         #else
            const size_t position = fname.find_last_of("/");
         #endif
         if (position == string::npos) {
            pathName = "";
            fnameWithoutPath = fname;
         } else {
            pathName = fname.substr(0,position);
            fnameWithoutPath = fname.substr(position+1);
         }
      }

      // Chdir to path containing input file and attempt to open the file, 
      // then chdir back to current working directory. 
      // Chdir returns zero value if it succeeds
      if (fileio::chdir(pathName.c_str()) != 0) {
		  stringstream ss;
		  ss << "ERROR: Could not chdir to input file dir '" << pathName << "' in " << __FILE__ << ":" << __LINE__;
          return withError(ss);
	  }
      filein.open(fnameWithoutPath.c_str(), fstream::in | fstream::binary);
      if (fileio::chdir(cwd) != 0) {
		  stringstream ss;
		  ss << "ERROR: Could not chdir to working dir '" << cwd << "' in " << __FILE__ << ":" << __LINE__;
          return withError(ss);
	  }

      if (filein.good() == true) {
         fileName = fnameWithoutPath;
      } else {
         filein.close();
		 stringstream ss;
		 ss << "ERROR: Input file not valid in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
   
      // Detect file endianness:
      char* ptr = reinterpret_cast<char*>(&endiannessFile);
      filein.read(ptr,1);
      if (filein.gcount() != 1) {
         stringstream ss;
         ss << "ERROR: Failed to read file endianness in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      if (endiannessFile != endiannessReader) swapIntEndianness = true;

      // Read footer offset:
      uint64_t footerOffset;
      char buffer[16];
      filein.seekg(8);
      filein.read(buffer,8);
      if (filein.gcount() != 8) {
         stringstream ss;
         ss << "ERROR: Failed to read footer offset in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      footerOffset = convUInt64(buffer,swapIntEndianness);
   
      // Read footer XML tree:
      filein.seekg(footerOffset);
      if (xmlReader.read(filein) == false) {
         stringstream ss;
         ss << "ERROR: Failed to read footer in " << __FILE__ << ":" << __LINE__ << endl;
         ss << "       MuXML says '" << xmlReader.getLastError() << "'";
         filein.close();
         return withError(ss);
      }
      filein.clear();
      filein.seekg(16);
      
      fileOpen = true;
      return true;
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
      // Check that a file is open:
      if (fileOpen == false) {
         stringstream ss;
         ss << "ERROR: readArray called but a file is not open in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
   
      // If zero-length read was requested, exit immediately:
      if (amount == 0) return true;
      
      // Find tag corresponding to given array:
      muxml::XMLNode* node = xmlReader.find(tagName,attribs);
      if (node == NULL) {
         stringstream ss;
         ss << "ERROR: Failed to find tag='" << tagName << "' attribs:" << endl;
         for (list<pair<string,string> >::const_iterator it=attribs.begin(); it!=attribs.end(); ++it) {
            ss << '\t' << it->first << " = '" << it->second << "'" << endl;
         }
         ss << "Error generated in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
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
         stringstream ss;
         ss << "ERROR: Unknown datatype in tag in " __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
   
      // Check that array contains sane values:
      if (arrayOpen.arraySize == 0) {
         stringstream ss;
         ss << "ERROR: Attempted to read from zero size array in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      if (arrayOpen.vectorSize == 0) {
         stringstream ss;
         ss << "ERROR: Attempted to read from array having zero vector size in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      if (arrayOpen.dataSize == 0) {
         stringstream ss;
         ss << "ERROR: Attempted to read from array having zero data size " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      
      // Sanity check on input values:
      if (begin + amount > arrayOpen.arraySize) {
         stringstream ss;
         ss << "ERROR: Starting element " << begin << " is larger than array size " << arrayOpen.arraySize;
         ss << " or amount to read " << begin+amount << " exceeds array size in ";
         ss << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }

      // Read data from file:
      streampos start      = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
      streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
      filein.clear();
      filein.seekg(start);
      if (filein.good() != true || filein.tellg() != start) {
         stringstream ss;
         ss << "ERROR: Failed to seek read position " << start << " in " << __FILE__ << ":" << __LINE__;
         return withError(ss);
      }
      filein.read(buffer,readBytes);

      // Check that we were able to read the requested amount of data:
      if (filein.gcount() != readBytes) {
         stringstream ss;
         ss << "ERROR: Failed to read requested " << readBytes << " B, only got " << filein.gcount() << " B" << endl;

         if (filein.eof() == true) ss << "fstream eof bit is set" << endl;
         if (filein.fail() == true) ss << "fstream fail and/or bad bit(s) are set" << endl;
         if (filein.bad() == true) ss << "fstream bad bit is set" << endl;
         ss << "current position in stream is " << (uint64_t)filein.tellg() << endl;

         filein.clear();
         filein.seekg(0);
         streampos fsize = filein.tellg();
         filein.seekg(0,ios::end);
         ss << "file size is " << filein.tellg()-fsize << " bytes " << endl;

         ss << "tag name='" << tagName << "'" << endl; 
         ss << "attributes:" << endl;
         for (map<string,string>::const_iterator it=node->attributes.begin(); it!=node->attributes.end(); ++it) {
            ss << '\t' << it->first << " = " << it->second << endl;
         }
         ss << "array offset string '" << node->value.c_str() << "'" << endl;
         ss << "begin =" << begin << endl;
         ss << "start =" << start << " readBytes=" << readBytes << endl;
         ss << "offset=" << arrayOpen.offset << " vectorsize=" << arrayOpen.vectorSize << " dataSize=" << arrayOpen.dataSize << endl;
         ss << "Error generated in " __FILE__ << ":" << __LINE__;

         return withError(ss);
      }
      return true;
   }

} // namespace vlsv