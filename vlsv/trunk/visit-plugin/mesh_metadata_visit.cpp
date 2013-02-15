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

#include "mesh_metadata.h"
#include "mesh_metadata_visit.h"

using namespace std;

namespace vlsvplugin {
   VisitMeshMetadata::VisitMeshMetadata(): MeshMetadata() { 
      blockOrigin = 0;
      N_domains = 1;
   }
   
   VisitMeshMetadata::~VisitMeshMetadata() { }
   
   int VisitMeshMetadata::getBlockOrigin() const {return blockOrigin;}
   
   avtMeshType VisitMeshMetadata::getMeshType() const {return meshType;}
   
   std::string VisitMeshMetadata::getMeshTypeString() const {return meshTypeString;}

   int VisitMeshMetadata::getNumberOfDomains() const {return N_domains;}
   
   int VisitMeshMetadata::getSpatialDimension() const {return spatialDimension;}
   
   int VisitMeshMetadata::getTopologicalDimension() const {return topologicalDimension;}

   bool VisitMeshMetadata::read(vlsv::Reader* vlsvReader,const std::map<std::string,std::string>& attribs) {
      return MeshMetadata::read(vlsvReader,attribs);
   }
   
} // namespace vlsvplugin


