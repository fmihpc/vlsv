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

#ifndef MESH_READER_VISIT_UCD_GENERIC_MULTI_H
#define MESH_READER_VISIT_UCD_GENERIC_MULTI_H

#include <unordered_map>
#include <mesh_reader.h>
#include <mesh_metadata_visit_ucd_generic_multi.h>

#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGrid.h>

namespace vlsvplugin {

   class VisitUCDGenericMultiMeshReader: public MeshReader {
    public:
      VisitUCDGenericMultiMeshReader();
      virtual ~VisitUCDGenericMultiMeshReader();
      
      virtual bool readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output);
      virtual bool readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output);
      
    protected:
      bool xPeriodic;
      bool yPeriodic;
      bool zPeriodic;
      bool nodeCoordinateArraysRead;

      virtual bool readCellVariable(vlsv::Reader* vlsvReader,VisitUCDGenericMultiMeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output);
      virtual bool readNodeVariable(vlsv::Reader* vlsvReader,VisitUCDGenericMultiMeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output);
      virtual bool readNodeCoordinateArrays(vlsv::Reader* vlsvReader,const std::string& meshName,uint64_t offset,uint64_t amount,vtkPoints* points);
   };

} // namespace plugin

#endif
