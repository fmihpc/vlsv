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

#include <mesh_vtk.h>
#include <vtkType.h>

using namespace std;

namespace vlsvplugin {

   int getNumberOfVertices(uint32_t cellType) {
      switch (cellType) {
       case VTK_VERTEX:
         return 1;
         break;
       case VTK_LINE:
         return 2;
         break;
       case VTK_TRIANGLE:
         return 3;
         break;
       case VTK_QUAD:
         return 4;
         break;
       case VTK_TETRA:
         return 4;
         break;
       case VTK_VOXEL:
         return 8;
         break;
       case VTK_HEXAHEDRON:
         return 8;
         break;
       case VTK_WEDGE:
         return 6;
         break;
       case VTK_PYRAMID:
         return 5;
         break;
       default:
         cerr << "Unsupported VTK cell type, id# " << cellType << endl;
         exit(1);
         break;
      }
   }
   
   int getVtkCelltype(uint32_t cellType) {
      switch (cellType) {
       case vlsv::celltype::UNKNOWN:
         cerr << "UNKNOWN VTK cell type, exiting" << endl;
         exit(1);
         break;
       case vlsv::celltype::VERTEX:
         return VTK_VERTEX;
         break;
       case vlsv::celltype::LINE:
         return VTK_LINE;
         break;
       case vlsv::celltype::TRIANGLE:
         return VTK_TRIANGLE;
         break;
       case vlsv::celltype::QUAD:
         return VTK_QUAD;
         break;
       case vlsv::celltype::TETRA:
         return VTK_TETRA;
         break;
       case vlsv::celltype::PYRAMID:
         return VTK_PYRAMID;
         break;
       case vlsv::celltype::WEDGE:
         return VTK_WEDGE;
         break;
       case vlsv::celltype::HEXAHEDRON:
         return VTK_HEXAHEDRON;
         break;
       case vlsv::celltype::VOXEL:
         return VTK_VOXEL;
         break;
       default:
         cerr << "Unsupported VTK cell type in " << __FILE__ << ' ' << __LINE__ << endl;
         exit(1);
         break;
      }
   }
   
   int getVtkDatatype(const vlsv::datatype::type& datatype,const uint64_t& dataSize) {
      switch (datatype) {
       case vlsv::datatype::UNKNOWN:
         return VTK_DATATYPE_NOT_FOUND;
         break;
	 
       case vlsv::datatype::INT:
         switch (dataSize) {
          case sizeof(int8_t):
            return VTK_CHAR;
            break;
          case sizeof(int16_t):
            return VTK_SHORT;
            break;
          case sizeof(int32_t):
            return VTK_INT;
            break;
          case sizeof(int64_t):
            return VTK_LONG;
            break;
          default:
            return VTK_DATATYPE_NOT_FOUND;
            break;
         }
         break;
	 
       case vlsv::datatype::UINT:
         switch (dataSize) {
          case sizeof(uint8_t):
            return VTK_UNSIGNED_CHAR;
            break;
          case sizeof(uint16_t):
            return VTK_UNSIGNED_SHORT;
            break;
          case sizeof(uint32_t):
            return VTK_UNSIGNED_INT;
            break;
          case sizeof(uint64_t):
            return VTK_UNSIGNED_LONG;
            break;
          default:
            return VTK_DATATYPE_NOT_FOUND;
            break;
         }
         break;
         
       case vlsv::datatype::FLOAT:
         switch (dataSize) {
          case sizeof(float):
            return VTK_FLOAT;
            break;
          case sizeof(double):
            return VTK_DOUBLE;
            break;
#ifndef _WINDOWS
          case sizeof(long double):
            return VTK_DATATYPE_NOT_FOUND;
            break;
#endif
          default:
            return VTK_DATATYPE_NOT_FOUND;
            break;
         }
         break;
         
       default:
         return VTK_DATATYPE_NOT_FOUND;
         break;
      }
   }
   
} // namespace vlsvplugin
