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

#ifndef VLSV_COMMON_H
#define VLSV_COMMON_H

#include <cstdlib>
#include <iostream>
#include <stdint.h>

namespace vlsv {
   
   namespace datatype {
      const unsigned char ENDIANNESS_LITTLE = 0;
      const unsigned char ENDIANNESS_BIG    = 1;
      
      enum type {
	 UNKNOWN,                                            /**< Unknown or unsupported datatype.*/
	 INT,                                                /**< Signed integer datatype.*/
	 UINT,                                               /**< Unsinged integer datatype.*/
	 FLOAT                                               /**< Floating point datatype.*/
      };
   }
   
   namespace geometry {
      enum type {
	   UNKNOWN,                                          /**< Mesh has unknown or unsupported coordinate system.*/
	   CARTESIAN,                                        /**< Mesh uses Cartesian (x,y,z) coordinate system.*/
	   CYLINDRICAL,                                      /**< Mesh uses cylindrical (r,phi,z) coordinate system.*/
	   SPHERICAL                                         /**< Mesh uses spherical (r,theta,phi) coordinate system.*/
      };
      
      const std::string STRING_UNKNOWN = "unknown_geometry"; /**< Mesh has unknown or unsupported geometry.*/
      const std::string STRING_CARTESIAN = "cartesian";      /**< Cartesian mesh geometry.*/
      const std::string STRING_CYLINDRICAL = "cylindrical";  /**< Cylindrical mesh geometry.*/
      const std::string STRING_SPHERICAL = "spherical";      /**< Spherical mesh geometry.*/
   }
   
   namespace mesh {
      enum type {
	 UNKNOWN,
	 POINT,
	 QUAD,
	 QUAD_MULTI,
	 UCD_MULTI
      };
      
      const std::string STRING_UNKNOWN = "unknown";          /**< Unknown or unsupported mesh type.*/
      const std::string STRING_POINT = "point";              /**< Point mesh.*/
      const std::string STRING_QUAD = "quad";
      const std::string STRING_QUAD_MULTI = "multimesh";     /**< Multi-domain mesh with uniform cell size.*/
      const std::string STRING_UCD_MULTI = "multi_ucd";      /**< Unstructured non-uniform curvilinear multi-domain mesh.*/
   }
   
   template<typename T> T convertFloat(const char* const ptr);
   template<typename T> T convertInteger(const char* const ptr,const bool& swapEndianness=false);
   template<typename T> void convertValue(T& value,const char* const ptr,datatype::type dt,int dataSize,const bool& swapEndianness=false);
   
   /** Returns a string representation of an array that is to be written to file.
    * The correct C++ datatype can be deduced from the string value returned by this
    * function and from the byte size of the data in array. For example, datatype "float"
    * with byte size of 8 usually means that the values are doubles.
    * @return String representation of the datatype.*/
   template<typename T> std::string getStringDatatype();
   
   const std::string& getMeshGeometry(geometry::type geom);
   geometry::type getMeshGeometry(const std::string& s);
   datatype::type getVLSVDatatype(const std::string& s);
   
   // ********************************************* //
   // ***** DEFINITIONS OF TEMPLATE FUNCTIONS ***** //
   // ********************************************* //

   template<typename T> inline
   T convertFloat(const char* const ptr) {
      return *reinterpret_cast<const T*>(ptr);
   }
   
   template<typename T> inline 
   T convertInteger(const char* const ptr,const bool& swapEndianness) {
      if (swapEndianness == false) return *reinterpret_cast<const T*>(ptr);
      int index = 0;
      T tmp = 0;
      char* const ptrtmp = reinterpret_cast<char*>(&tmp);
      for (int i=sizeof(T)-1; i>=0; --i) {
	 ptrtmp[index] = ptr[i];
	 ++index;
      }
      return tmp;
   }
   
   /** Driver function for convertFloat and convertInteger. Value from given buffer 
    * is converted into datatype given with template parameter T. Endianness of 
    * integer datatypes is converted if necessary. Note that if vlsv::datatype is 
    * vlsv::UNKNOWN the contents of buffer are simply copied into output variable 'value'.
    * @param value Variable in which the value from buffer is copied.
    * @param buffer Byte array containing the desired value.
    * @param dt vlsv::datatype of the value in buffer.
    * @param dataSize Byte size of the value in buffer.
    * @param swapEndianness If true, endianness of integer datatypes is swapped before
    * the value is copied to output variable 'value'.*/
   template<typename T> inline
     void convertValue(T& value,const char* const buffer,datatype::type dt,int dataSize,const bool& swapEndianness=false) {
      char* valuePtr = NULL;
      // Switch according the native datatype of the value in buffer:
      switch (dt) {
       case datatype::UNKNOWN:
	 // Unknown datatype, just byte-copy buffer to 'value':
	 valuePtr = reinterpret_cast<char*>(&value);
	 for (int i=0; i<dataSize; ++i) valuePtr[i] = buffer[i];
	 break;
       case datatype::INT:
	 // Signed integer, switch according to byte size:
	 switch (dataSize) {
	  case sizeof(int8_t):
	    value = convertInteger<int8_t>(buffer,swapEndianness);
	    break;
	  case sizeof(int16_t):
	    value = convertInteger<int16_t>(buffer,swapEndianness);
	    break;
	  case sizeof(int32_t):
	    value = convertInteger<int32_t>(buffer,swapEndianness);
	    break;
	  case sizeof(int64_t):
	    value = convertInteger<int64_t>(buffer,swapEndianness);
	    break;
	  default:
	    std::cerr << "(VLSV) ERROR: Unsupported datatype in convertValue!" << std::endl; 
	    std::cerr << "\t Exiting at convertValue VLSV::INT." << std::endl;
	    exit(1);
	    break;
	 }
	 break;
       case datatype::UINT:
	 // Unsigned integer, switch according to byte size:
	 switch (dataSize) {
	  case sizeof(uint8_t):
	    value = convertInteger<uint8_t>(buffer,swapEndianness);
	    break;
	  case sizeof(uint16_t):
	    value = convertInteger<uint16_t>(buffer,swapEndianness);
	    break;
	  case sizeof(uint32_t):
	    value = convertInteger<uint32_t>(buffer,swapEndianness);
	    break;
	  case sizeof(uint64_t):
	    value = convertInteger<uint64_t>(buffer,swapEndianness);
	    break;
	  default:
	    std::cerr << "(VLSV) ERROR: Unsupported datatype in convertValue!" << std::endl; 
	    std::cerr << "\t Exiting at convertValue VLSV::UINT." << std::endl;
	    exit(1);
	    break;
	 }
	 break;
       case datatype::FLOAT:
	 // Floating point, switch according to byte size:
	 switch (dataSize) {
	  case sizeof(float):
	    value = convertFloat<float>(buffer);
	    break;
	  case sizeof(double):
	    value = convertFloat<double>(buffer);
	    break;
	  case sizeof(long double):
	    value = convertFloat<long double>(buffer);
	    break;
	  default:
	    std::cerr << "(VLSV) ERROR: Unsupported datatype in convertValue!" << std::endl;
	    std::cerr << "\t Exiting at convertValue VLSV::FLOAT." << std::endl;
	    exit(1);
	    break;
	 }
	 break;
       default:
	 // Error:
	 std::cerr << "(VLSV) ERROR: Unsupported datatype in convertValue!" << std::endl;
	 exit(1);
	 break;
      }
   }
   
   template<typename T> inline std::string getStringDatatype() {return "unknown";}
   template<> inline std::string getStringDatatype<bool>() {return "int";}
   template<> inline std::string getStringDatatype<char>() {return "uint";}
   template<> inline std::string getStringDatatype<int8_t>() {return "int";}
   template<> inline std::string getStringDatatype<int16_t>() {return "int";}
   template<> inline std::string getStringDatatype<int32_t>() {return "int";}
   template<> inline std::string getStringDatatype<int64_t>() {return "int";}
   template<> inline std::string getStringDatatype<uint8_t>() {return "uint";}
   template<> inline std::string getStringDatatype<uint16_t>() {return "uint";}
   template<> inline std::string getStringDatatype<uint32_t>() {return "uint";}
   template<> inline std::string getStringDatatype<uint64_t>() {return "uint";}
   template<> inline std::string getStringDatatype<float>() {return "float";}
   template<> inline std::string getStringDatatype<double>() {return "float";}
   template<> inline std::string getStringDatatype<long double>() {return "float";}

   unsigned char detectEndianness();

   int8_t convInt8(const char* const ptr,const bool& swapEndian=false);
   int16_t convInt16(const char* const ptr,const bool& swapEndian=false);
   int32_t convInt32(const char* const ptr,const bool& swapEndian=false);
   int64_t convInt64(const char* const ptr,const bool& swapEndian=false);
   uint8_t convUInt8(const char* const ptr,const bool& swapEndian=false);
   uint16_t convUInt16(const char* const ptr,const bool& swapEndian=false);
   uint32_t convUInt32(const char* const ptr,const bool& swapEndian=false);
   uint64_t convUInt64(const char* const ptr,const bool& swapEndian=false);
   float convReal4(const char* const ptr,const bool& swapEndian=false);
   double convReal8(const char* const ptr,const bool& swapEndian=false);

} // namespace vlsv
   
#endif
