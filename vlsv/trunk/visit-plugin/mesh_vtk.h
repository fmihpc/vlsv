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

#ifndef MESH_VTK_H
#define MESH_VTK_H

#include <vlsv_common.h>
#include <vtkCellType.h>

namespace vlsvplugin {

   const int VTK_DATATYPE_NOT_FOUND = -999;

   int getNumberOfVertices(uint32_t cellType);
   int getVtkCelltype(uint32_t cellType);
   int getVtkDatatype(const vlsv::datatype::type& datatype,const uint64_t& dataSize);
   
} // namespace vlsvplugin
   
#endif
