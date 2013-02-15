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

   int getVtkDatatype(const vlsv::datatype::type& datatype,const uint64_t& dataSize) {
      switch (datatype) {
       case vlsv::datatype::UNKNOWN:
	 return VTK_DATATYPE_NOT_FOUND;
	 break;
	 
       case vlsv::datatype::INT:
	 switch (dataSize) {
	  case sizeof(char):
	    return VTK_CHAR;
	    break;
	  case sizeof(short):
	    return VTK_SHORT;
	    break;
	  case sizeof(int):
	    return VTK_INT;
	    break;
	  case sizeof(long int):
	    return VTK_LONG;
	    break;
	  default:
	    return VTK_DATATYPE_NOT_FOUND;
	    break;
	 }
	 break;
	 
       case vlsv::datatype::UINT:
	 switch (dataSize) {
	  case sizeof(char):
	    return VTK_UNSIGNED_CHAR;
	    break;
	  case sizeof(short):
	    return VTK_UNSIGNED_SHORT;
	    break;
	  case sizeof(int):
	    return VTK_UNSIGNED_INT;
	    break;
	  case sizeof(long int):
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
	  case sizeof(long double):
	    return VTK_DATATYPE_NOT_FOUND;
	    break;
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
