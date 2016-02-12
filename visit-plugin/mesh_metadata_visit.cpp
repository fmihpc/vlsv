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
   /** Default constructor.*/
   VisitMeshMetadata::VisitMeshMetadata(): MeshMetadata() {
      blockOrigin = 0;
   }
   
   /** Default virtual constructor.*/
   VisitMeshMetadata::~VisitMeshMetadata() { }

   /** Checks that the specified mesh in VLSV file is of the type the derived class reads.
    *  @param vlsv vlsv::Reader that has the input file open.
    *  @param attribs Attributes for the XML tag 'MESH'.
    *  @return If true, this metadata reader understands the mesh in attribs.*/
   bool VisitMeshMetadata::checkVlsvMeshType(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      const bool rvalue = MeshMetadata::checkVlsvMeshType(vlsvReader,attribs);
      if (rvalue == false) {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }
      return rvalue;
   }

   int VisitMeshMetadata::getBlockOrigin() const {return blockOrigin;}

   /** Read domain metadata for this mesh.
    *  This function calls readDomainMetadata in MeshMetadata (or in one of its subclasses)
    *  and writes additional debugging information to VisIt log files.
    *  @param vlsvReader vlsv::Reader that has the input file open.
    *  @return If true, domain metadata was successfully read.*/
   bool VisitMeshMetadata::readDomainMetadata(vlsv::Reader* vlsvReader) {
      debug2 << "VLSV\t\t VisitMeshMetadata::readDomains called" << endl;

      const bool rvalue = MeshMetadata::readDomainMetadata(vlsvReader);

      if (rvalue == true) {
         if (vlsvMeshType == vlsv::mesh::UCD_MULTI) {
            debug4 << "VLSV\t\t Mesh size in blocks: " << meshBoundingBox[0] << ' ' << meshBoundingBox[1] << ' ' << meshBoundingBox[2];
            debug4 << " block sizes: " << meshBoundingBox[3] << ' ' << meshBoundingBox[4] << ' ' << meshBoundingBox[5] << endl;
         } else if (vlsvMeshType == vlsv::mesh::QUAD_MULTI) {
            debug4 << "VLSV\t\t Mesh corner coordinates: " << meshCoordinates[0] << ' ' << meshCoordinates[1] << ' ' << meshCoordinates[2];
            debug4 << " cell sizes: " << meshCoordinates[3] << ' ' << meshCoordinates[4] << ' ' << meshCoordinates[5] << endl;
         }
      } else {
         debug3 << "VLSV\t\t " << MeshMetadata::getErrorString() << endl;
      }

      return rvalue;
   }

   /** Read names and types of the variables that belong to this mesh.
    *  This function calls readVariables function in MeshMetadata 
    *  (or in one of its subclasses) and writes additional debugging information to VisIt log files.
    *  @param vlsv vlsv::Reader that has the input file open.
    *  @param attribs Attributes for the XML tag 'MESH'.
    *  @return If true, variable metadata was read successfully.*/
   bool VisitMeshMetadata::readVariables(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      const bool rvalue = MeshMetadata::readVariables(vlsvReader,attribs);

      if (rvalue == true) {
         debug4 << "VLSV\t\t Found variables:" << endl;
         for (auto& it: variableMetadata) {
            debug4 << "\t\t\t '" << it.name << "' vectorsize = " << it.vectorSize << endl;
         }
      }

      return rvalue;
   }
   
} // namespace vlsvplugin


