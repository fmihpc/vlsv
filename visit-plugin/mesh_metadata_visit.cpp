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

#include <DebugStream.h>

#include "mesh_metadata_visit.h"

using namespace std;

namespace vlsvplugin {
   VisitMeshMetadata::VisitMeshMetadata(): MeshMetadata() {
      blockOrigin = 0;
   }
   
   VisitMeshMetadata::~VisitMeshMetadata() { }

   int VisitMeshMetadata::getBlockOrigin() const {return blockOrigin;}

   /** Get VisIt / VTK mesh type that corresponds to the VLSV mesh.
    * Type of VisIt / VTK mesh stored in the VLSV mesh.*/
   //avtMeshType VisitMeshMetadata::getAvtMeshType() const {return meshType;}
   
   /** Get string representation of VisIt / VTK mesh type that corresponds to the VLSV mesh.
    * @return String representation of VisIt / VTK mesh.*/
   //std::string VisitMeshMetadata::getAvtMeshTypeString() const {return meshTypeString;}
   
} // namespace vlsvplugin


