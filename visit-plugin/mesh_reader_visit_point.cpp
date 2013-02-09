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

#include <mesh_reader_visit_point.h>
#include <mesh_metadata_visit_point.h>
#include <mesh_vtk.h>

#include <list>
#include <typeinfo>

#include <vtkCellType.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <DebugStream.h>

using namespace std;

namespace vlsvplugin {

   VisitPointMeshReader::VisitPointMeshReader(): MeshReader() { }
   
   VisitPointMeshReader::~VisitPointMeshReader() { }
   
   bool VisitPointMeshReader::readMesh(VLSVReader* vlsv,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitPointMeshReader::read called" << endl;
      output = NULL;
      
      // Check that VLSVReader exists:
      if (vlsv == NULL) {
	 debug2 << "VLSV\t\t ERROR: VLSVReader is NULL" << endl;
	 return false;
      }
      
      // Check that metadata is not NULL:
      if (md == NULL) {
	 debug2 << "VLSV\t\t ERROR: MeshMetadata object is NULL" << endl;
	 return false;
      }
      
      // Check that given metadata is of correct type:
      const VisitPointMeshMetadata* const metadata = dynamic_cast<const VisitPointMeshMetadata*>(md);
      if (typeid(*md) != typeid(*metadata)) {
	 debug2 << "VLSV\t\t ERROR: Given mesh metadata object is not of type VisitPointMeshMedata" << endl;
	 return false;
      }
      /*
      // Read point mesh XML tag:      
      uint64_t arraySize,vectorSize,dataSize;
      VLSV::datatype datatype;
      
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",metadata->getName()));
      if (vlsv->getArrayInfo("MESH",attribs,arraySize,vectorSize,datatype,dataSize) == false) {
	 debug2 << "VLSV\t\t ERROR: VLSVReader failed to read array info." << endl;
	 return false;
      }
      */
      debug4 << "VLSV\t\t arraysize:  " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize: " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:   " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:   " << metadata->getDatatype() << endl;
      
      uint64_t readOffset = 0;
      uint64_t N_points = metadata->getArraySize();
      
      // Determine VTK datatype that corresponds to the one in VLSV file:
      vtkPoints* points = NULL;
      int vtkDatatype = vlsvplugin::getVtkDatatype(metadata->getDatatype(),metadata->getDataSize());
      if (vtkDatatype == VTK_DATATYPE_NOT_FOUND) {
	 debug4 << "VLSV\t\t ERROR: Could not find VTK datatype for VLSV datatype '" << metadata->getDatatype();
	 debug4 << "' datasize '" << metadata->getDataSize() << endl;
	 return false;
      }
      
      // Create new vtkPoints object and read data directly to its internal array:
      points = vtkPoints::New(vtkDatatype);
      points->SetNumberOfPoints(N_points);
      char* ptr = reinterpret_cast<char*>(points->GetVoidPointer(0));
      
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",metadata->getName()));
      if (vlsv->readArray("MESH",attribs,readOffset,N_points,ptr) == false) {
	 debug2 << "VLSV\t\t ERROR: VLSVReader failed to read array" << endl;
	 points->Delete();
	 return false;
      }
      
      // Create new unstructured grid that only contains single vertices (=points):
      vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
      ugrid->SetPoints(points);
      points->Delete();
      ugrid->Allocate(N_points);
      for (uint64_t i=0; i<N_points; ++i) {
	 vtkIdType onevertex = static_cast<vtkIdType>(i);
	 ugrid->InsertNextCell(VTK_VERTEX,1,&onevertex);
      }

      output = ugrid;
      return true;
   }

   bool VisitPointMeshReader::readVariable(VLSVReader* vlsv,MeshMetadata* md,const VariableMetadata& vmd,int domain,float*& output) {
      output = NULL;
      return false;
   }
   
} // namespace vlsvplugin
