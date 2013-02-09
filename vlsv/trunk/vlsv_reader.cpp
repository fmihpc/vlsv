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
#include <unistd.h>

#include "vlsv_reader.h"

using namespace std;

VLSVReader::VLSVReader() {
   endiannessReader = detectEndianness();
   fileOpen = false;
   swapIntEndianness = false;
}

VLSVReader::~VLSVReader() {
   filein.close();   
}

bool VLSVReader::close() {
   filein.close();
   xmlReader.clear();
   return true;
}

/** Get attributes of the given XML tag.
 * @param tagName Name of the XML tag.
 * @param attribsIn Constraints that limit the search.
 * @param attribsOut Attributes of the XML tag, if one matched given constraints.
 * @return If true, an XML tag was found that mathes given constraints.*/
bool VLSVReader::getArrayAttributes(const string& tagName,const list<pair<string,string> >& attribsIn,map<string,string>& attribsOut) const {
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
bool VLSVReader::getArrayInfo(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			      uint64_t& arraySize,uint64_t& vectorSize,VLSV::datatype& dataType,uint64_t& dataSize) const {
   if (fileOpen == false) return false;
   muxml::XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) return false;
   
   arraySize = atol(node->attributes["arraysize"].c_str());
   vectorSize = atol(node->attributes["vectorsize"].c_str());
   dataSize = atol(node->attributes["datasize"].c_str());
   if (node->attributes["datatype"] == "unknown") dataType = VLSV::UNKNOWN;
   else if (node->attributes["datatype"] == "int") dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype '" << node->attributes["datatype"] << "' in tag!" << endl;
      return false;
   }
   return true;
}

bool VLSVReader::getFileName(std::string& openFile) const {
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
bool VLSVReader::getUniqueAttributeValues(const string& tagName,const string& attribName,set<string>& output) const {
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

bool VLSVReader::loadArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs) {
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
   if (node->attributes["datatype"] == "unknown") arrayOpen.dataType = VLSV::UNKNOWN;
   else if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }   
   if (arrayOpen.arraySize == 0) return false;
   if (arrayOpen.vectorSize == 0) return false;
   if (arrayOpen.dataSize == 0) return false;
   
   return true;
}

/** Open a VLSV file for reading. This function fails if a 
 * file is already open. 
 * @param fname File name.
 * @return If true, file was successfully opened.*/
bool VLSVReader::open(const std::string& fname) {
   bool success = true;
   if (fileOpen == true) {
      #ifndef NDEBUG
         cerr << "VLSVReader ERROR: File '" << fname << "' should be opened, but file '";
         cerr << fileName << "' is currently open." << endl;
      #endif
      return false;
   }
   
   // TEST
   string fnameWithoutPath;
   string pathName;
   char cwd[1024];
   if (getcwd(cwd,sizeof(cwd)) != NULL) {
      //cerr << "VLSVReader current directory is '" << cwd << "'" << endl;
      //cerr << "\t file '" << fname << "' should be opened" << endl;
      
      const size_t position = fname.find_last_of("/");
      if (position == string::npos) {
	 pathName = "";
	 fnameWithoutPath = fname;
      } else {
	 pathName = fname.substr(0,position);
	 fnameWithoutPath = fname.substr(position+1);
      }
      
      //cerr << "\t path '" << fname.substr(0,position) << "'" << endl;
      //cerr << "\t file '" << fname.substr(position+1) << "'" << endl;
   }
   // END TEST
   
   chdir(pathName.c_str());
   filein.open(fnameWithoutPath.c_str(), fstream::in);
   chdir(cwd);
   
   if (filein.good() == true) {
      fileName = fnameWithoutPath;
      fileOpen = true;
      // DEBUG
      //cerr << "\t\t SUCCESS" << endl;
      // END DEBUG
   } else {
      filein.close();
      success = false;
      // DEBUG
      //cerr << "\t\t FAILURE" << endl;
      // END DEBUG
   }
   if (success == false) {
      #ifndef NDEBUG
         cerr << "VLSVReader ERROR: File '" << fnameWithoutPath << "' could not be opened!" << endl;
      #endif
      return success;
   }
   
   // Detect file endianness:
   char* ptr = reinterpret_cast<char*>(&endiannessFile);
   filein.read(ptr,1);
   if (endiannessFile != endiannessReader) swapIntEndianness = true;

   // Read footer offset:
   uint64_t footerOffset;
   char buffer[16];
   filein.seekg(8);
   filein.read(buffer,8);
   footerOffset = convUInt64(buffer,swapIntEndianness);
   
   // Read footer XML tree:
   filein.seekg(footerOffset);
   xmlReader.read(filein);
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
bool VLSVReader::readArray(const std::string& tagName,const std::list<std::pair<std::string,std::string> >& attribs,
			   const uint64_t& begin,const uint64_t& amount,char* buffer) {
   if (fileOpen == false) {
      cerr << "VLSVReader ERROR: readArray called but a file is not open!" << endl;
      return false;
   }
   
   // If zero-length read was requested, exit immediately:
   if (amount == 0) return true;
   
   // Find tag corresponding to given array:
   muxml::XMLNode* node = xmlReader.find(tagName,attribs);
   if (node == NULL) {
      cerr << "VLSVReader ERROR: Failed to find tag='" << tagName << "' attribs:" << endl;
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
   if (node->attributes["datatype"] == "int") arrayOpen.dataType = VLSV::INT;
   else if (node->attributes["datatype"] == "uint") arrayOpen.dataType = VLSV::UINT;
   else if (node->attributes["datatype"] == "float") arrayOpen.dataType = VLSV::FLOAT;
   else {
      cerr << "VLSVReader ERROR: Unknown datatype in tag!" << endl;
      return false;
   }
   
   if (arrayOpen.arraySize == 0) return false;
   if (arrayOpen.vectorSize == 0) return false;
   if (arrayOpen.dataSize == 0) return false;
   
   // Sanity check on values:
   if (begin + amount > arrayOpen.arraySize) {
      cerr << "VLSVReader ERROR: Requested read exceeds array size. begin: " << begin << " amount: " << amount << " size: " << arrayOpen.arraySize << endl;
      return false;
   }

   streamoff start = arrayOpen.offset + begin*arrayOpen.vectorSize*arrayOpen.dataSize;
   streamsize readBytes = amount*arrayOpen.vectorSize*arrayOpen.dataSize;
   filein.clear();
   filein.seekg(start);
   filein.read(buffer,readBytes);
   if (filein.gcount() != readBytes) {
      cerr << "VLSVReader ERROR: Failed to read requested amount of bytes!" << endl;      
      cerr << "tag name='" << tagName << "'" << endl;
      cerr << "attributes:" << endl;
      for (map<string,string>::const_iterator it=node->attributes.begin(); it!=node->attributes.end(); ++it) {
	 cerr << '\t' << it->first << " = " << it->second << endl;
      }
      cerr << "array offset string '" << node->value.c_str() << "'" << endl;
      cerr << "start=" << start << " readBytes=" << readBytes << endl;
      cerr << "offset=" << arrayOpen.offset << " vectorsize=" << arrayOpen.vectorSize << " dataSize=" << arrayOpen.dataSize << endl;
      exit(1);
   }
   return true;
}
