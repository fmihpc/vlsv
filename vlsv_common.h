/** @file vlsv_common.h
 *  @par License
 *  This file is part of VLSV file format.
 *  @n@n
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  @n@n
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *  @n@n
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  @copyright Copyright 2011-2015 Finnish Meteorological Institute.
 *  @author Arto Sandroos
 */

#ifndef VLSV_COMMON_H
#define VLSV_COMMON_H

#include <cstdlib>
#include <iostream>
#include <stdint.h>

namespace vlsv {

   namespace celltype {
      enum type {
         UNKNOWN,
         VERTEX,
         LINE,
         TRIANGLE,
         QUAD,
         TETRA,
         PYRAMID,
         WEDGE,
         HEXAHEDRON,
         VOXEL
      };
   }
   
   /** Enumeration of all supported error codes.*/
   namespace error {
      enum type {
         NONE,                                           /**< No errors.*/
         UNKNOWN,                                        /**< Unknown or unsupported error.*/
         READ_CWD_FAIL,                                  /**< Reader failed to cwd to directory containing input file.*/
         READ_FILE_BAD,                                  /**< Reader failed to open input file.*/
         READ_FILE_ALREADY_OPEN,                         /**< Reader failed to open file because a file is alread open.*/
         READ_FILE_ENDIANNESS,                           /**< Reader failed to read file endianness.*/
         READ_NO_FOOTER,                                 /**< Input file has no footer.*/
         READ_FOOTER_OFFSET,                             /**< Reader failed to read footer offset.*/
         READ_FOOTER,                                    /**< Reader failed to read footer.*/
         SIZE
      };
   }
   
    /** Tells whether a datatype stored in a buffer or a file is a signed or unsigned integer, or a floating point number.
    * @brief Datatype description.*/
   namespace datatype {
      const uint8_t ENDIANNESS_LITTLE = 0;                   /**< Data in a buffer or a file has little endian encoding, only has an effect on integers.
                                                              * @brief Little endian encoding.*/
      const uint8_t ENDIANNESS_BIG    = 1;                   /**< Data in a buffer or a file has big endian encoding, only has an effect on integers.
                                                              * @brief Big endian encoding.*/

      enum type {
         UNKNOWN,                                            /**< @brief Unknown or unsupported datatype.*/
         INT,                                                /**< @brief Signed integer.*/
         UINT,                                               /**< @brief Unsinged integer.*/
         FLOAT                                               /**< @brief Floating point.*/
      };
   }
   
   namespace geometry {
      enum type {
	   UNKNOWN,                                          /**< Mesh has unknown or unsupported coordinate system.*/
	   CARTESIAN,                                        /**< Mesh uses Cartesian (x,y,z) coordinate system.*/
	   CYLINDRICAL,                                      /**< Mesh uses cylindrical (r,phi,z) coordinate system.*/
	   SPHERICAL,                                        /**< Mesh uses spherical (r,theta,phi) coordinate system.*/
	   UNSTRUCTURED                                      /**< Mesh is unstructured.*/
      };
      
      const std::string STRING_UNKNOWN = "unknown_geometry";  /**< Mesh has unknown or unsupported geometry.*/
      const std::string STRING_CARTESIAN = "cartesian";       /**< Cartesian mesh geometry.*/
      const std::string STRING_CYLINDRICAL = "cylindrical";   /**< Cylindrical mesh geometry.*/
      const std::string STRING_SPHERICAL = "spherical";       /**< Spherical mesh geometry.*/
      const std::string STRING_UNSTRUCTURED = "unstructured"; /**< Unstructured mesh geometry.*/
   }
   
   namespace mesh {
      enum type {
         UNKNOWN,
           POINT,
           QUAD,
           QUAD_MULTI,
           UCD_AMR,
           UCD_MULTI,
           UCD_GENERIC_MULTI
      };
      
      const std::string STRING_UNKNOWN = "unknown";          /**< Unknown or unsupported mesh type.*/
      const std::string STRING_POINT = "point";              /**< Point mesh.*/
      const std::string STRING_QUAD = "quad";
      const std::string STRING_QUAD_MULTI = "multimesh";     /**< Multi-domain mesh with uniform cell size.*/
      const std::string STRING_UCD_AMR = "amr_ucd";          /**< Unstructured refined mesh.*/
      const std::string STRING_UCD_MULTI = "multi_ucd";      /**< Unstructured non-uniform curvilinear multi-domain mesh.*/
      const std::string STRING_UCD_GENERIC_MULTI = "multi_ucd_generic"; /**< Generic unstructured multi-domain mesh.*/
   }
   
   namespace ucdgenericmulti {
      /** Bounding box contains information on the number of blocks per coordinate 
       * direction in the mesh, and the number of cells per coordinate direction in a block.
       * @brief Definition of elements in MESH_BBOX array.*/
      namespace bbox {
         /** @brief Description of MESH_BBOX entries.*/
         enum elements {
            X_BLOCKS,      /**< Number of blocks in x-direction.*/
            Y_BLOCKS,      /**< Number of blocks in y-direction.*/
            Z_BLOCKS,      /**< Number of blocks in z-direction.*/
            BLOCK_WIDTH_X, /**< Number of cells in block in x-direction.*/
            BLOCK_WIDTH_Y, /**< Number of cells in block in y-direction.*/
            BLOCK_WIDTH_Z, /**< Number of cells in block in z-direction.*/
            SIZE           /**< Size of MESH_BBOX array.*/
         };
      }
      
      /** Multi-domain mesh consists of one or more domains. The array MESH_DOMAIN_SIZES 
       * has an entry for each domain that tells how many real and ghost cells it contains.
       * @brief Definition of elements of an entry in MESH_DOMAIN_SIZES array.*/
      namespace domainsizes {
         /** @brief Definition of elements in MESH_DOMAIN_SIZES array.*/
         enum elements {
            TOTAL_BLOCKS, /**< Total number of blocks in domain.*/
            GHOST_BLOCKS, /**< Number of ghost blocks in domain.*/
            TOTAL_NODES,  /**< Total number of nodes in domain.*/
            GHOST_NODES,  /**< Number of ghost nodes in domain.*/
            SIZE          /**< Size of MESH_DOMAIN_SIZES entry for each domain.*/
         };
      }
      
      namespace offsets {
         /** Definition of elements in MESH_OFFSETS array.*/
         enum elements {
            ZONE_ENTRIES, /**< Number of zone connectivity entries in domain.*/
            NODE_ENTRIES, /**< Number of nodes in domain.*/
            SIZE          /**< Size of MESH_OFFSETS entry for each domain.*/
         };
      }
   }

   const std::string getErrorString(const vlsv::error::type& errorCode);
   
   template<typename T> T convertFloat(const char* const ptr);
   template<typename T> T convertInteger(const char* const ptr,const bool& swapEndianness=false);
   
   /** Returns a string representation of a basic datatype that 
    * is to be written to an array in output file.
    * The correct C++ basic datatype can be deduced from the string value returned by this
    * function and from the byte size of the data in array. For example, datatype "float"
    * with byte size of 8 usually means that the value is a "double".
    * @brief String representation of a basic datatype.
    * @return String representation of the datatype.*/
   template<typename T> std::string getStringDatatype();
   
   std::string getStringDatatype(const vlsv::datatype::type& dt);
   const std::string& getMeshGeometry(geometry::type geom);
   geometry::type getMeshGeometry(const std::string& s);
   datatype::type getVLSVDatatype(const std::string& s);
   
   // ********************************************* //
   // ***** DEFINITIONS OF TEMPLATE FUNCTIONS ***** //
   // ********************************************* //

   /** Given a pointer to a memory, convert the data to the floating point
    * datatype passed as the template parameter.
    * @brief Convert byte data to a floating point datatype.
    * @tparam T Basic C/C++ floating point datatype.
    * @param ptr Pointer to memory where the value starts.
    * @return Data in buffer converted to a floating point value.*/
   template<typename T> inline
   T convertFloat(const char* const ptr) {
      return *reinterpret_cast<const T*>(ptr);
   }
   
   /** Given a pointer to a memory, convert the data to the integer 
    * datatype passed as the template parameter.
    * @brief Convert byte data to a basic integer datatype.
    * @tparam T Basic C/C++ integer datatype.
    * @param ptr Pointer to memory where the value starts.
    * @param swapEndianness If true, the endianness is swapped before converting to basic datatype.
    * @return Data in buffer converted to an integer value.*/
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
    * vlsv::UNKNOWN the contents of buffer are simply byte-copied into output variable 'value'.
    * @brief Convert data in buffer to a basic datatype value.
    * @tparam Basic datatype that the buffer data is converted into.
    * @param value Variable in which the value from buffer is copied.
    * @param buffer Byte array containing the desired value.
    * @param dt vlsv::datatype of the value in buffer.
    * @param dataSize Byte size of the value in buffer.
    * @param swapEndianness If true, endianness of integer datatypes is swapped before
    * the value is copied to output variable 'value'.*/
   template<typename T> inline
   void convertValue(T& value,const char* const buffer,datatype::type dt,int dataSize,const bool& swapEndianness) {
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
#ifndef _WINDOWS
               case sizeof(long double):
#endif
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

   std::string printDataRate(const uint64_t& bytes,const double& t);
} // namespace vlsv
   
#endif
