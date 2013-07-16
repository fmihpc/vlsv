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

#ifndef MESH_METADATA_VISIT_POINT_H
#define MESH_METADATA_VISIT_POINT_H

#include <mesh_metadata_visit.h>

namespace vlsvplugin {
   class VisitPointMeshMetadata: public VisitMeshMetadata {
    public:
      VisitPointMeshMetadata();
      ~VisitPointMeshMetadata();

      const vlsv::geometry::type& getMeshGeometry() const;
      uint64_t getNumberOfGhostNodes(int domain) const;
      uint64_t getNumberOfGhostZones(int domain) const;
      uint64_t getNumberOfLocalNodes(int domain) const;
      uint64_t getNumberOfLocalZones(int domain) const;
      uint64_t getNumberOfTotalNodes(int domain) const;
      uint64_t getNumberOfTotalZones(int domain) const;
      
      bool read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs);
      
    protected:
      vlsv::geometry::type geometry;  /**< Mesh geometry (Cartesian, Cylindrical, etc.).*/
   };   
} // namespace vlsvplugin

#endif
