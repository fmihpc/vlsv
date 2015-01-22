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
   
   bool VisitPointMeshReader::readMesh(vlsv::Reader* vlsvReader,MeshMetadata* md,int domain,void*& output) {
      debug2 << "VLSV\t VisitPointMeshReader::read called" << endl;
      output = NULL;
      
      // Check that VLSVReader exists:
      if (vlsvReader == NULL) {
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

      debug4 << "VLSV\t\t arraysize:  " << metadata->getArraySize() << endl;
      debug4 << "VLSV\t\t vectorsize: " << metadata->getVectorSize() << endl;
      debug4 << "VLSV\t\t datasize:   " << metadata->getDataSize() << endl;
      debug4 << "VLSV\t\t datatype:   " << metadata->getDatatype() << endl;
      
      uint64_t readOffset = 0;
      uint64_t N_points = metadata->getArraySize();
      
      // Determine VTK datatype that corresponds to the one in VLSV file:
      vtkPoints* points = NULL;
      int vtkDatatype = VTK_FLOAT;
      
      // Get mesh geometry:
      const vlsv::geometry::type& geometry = metadata->getMeshGeometry();
      
      // Create new vtkPoints object and read data directly to its internal array:
      points = vtkPoints::New(vtkDatatype);
      points->SetNumberOfPoints(N_points);
      float* ptr = reinterpret_cast<float*>(points->GetVoidPointer(0));
      
      list<pair<string,string> > attribs;
      attribs.push_back(make_pair("name",metadata->getName()));
      if (vlsvReader->read("MESH",attribs,readOffset,N_points,ptr,false) == false) {
         debug2 << "VLSV\t\t ERROR: VLSVReader failed to read array" << endl;
         points->Delete();
         return false;
      }

      // Transform to Cartesian coordinates if necessary:
      float x,y,z;
      float crds[3];
      const double* transform = metadata->getTransform();

      switch (geometry) {
       case vlsv::geometry::UNKNOWN:
         debug2 << "VLSV\t\t WARNING: Unknown coordinate system, assuming Cartesian" << endl;
         break;
       case vlsv::geometry::CARTESIAN:
         for (uint64_t p=0; p<N_points; ++p) {
            crds[0] = ptr[3*p+0];
            crds[1] = ptr[3*p+1];
            crds[2] = ptr[3*p+2];

            x = ptr[3*p+0];
            y = ptr[3*p+1];
            z = ptr[3*p+2];

            ptr[3*p+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            ptr[3*p+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            ptr[3*p+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
            debug5 << "VLSV\t\t Point #" << p << "\t coords " << x << '\t' << y << '\t' << z << endl;
         }
         break;
       case vlsv::geometry::CYLINDRICAL:
         for (uint64_t p=0; p<N_points; ++p) {
            x = ptr[3*p+0];
            y = ptr[3*p+1];
            z = ptr[3*p+2];
            
            crds[0] = x * cos(y);
            crds[1] = x * sin(y);
            crds[2] = z;
            
            ptr[3*p+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            ptr[3*p+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            ptr[3*p+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
         }
         break;
       case vlsv::geometry::SPHERICAL:
         for (uint64_t p=0; p<N_points; ++p) {
            x = ptr[3*p+0];
            y = ptr[3*p+1];
            z = ptr[3*p+2];
            
            crds[0] = x * sin(y) * cos(z);
            crds[1] = x * sin(y) * sin(z);
            crds[2] = x * cos(y);
            
            ptr[3*p+0] = transform[0]*crds[0] + transform[1]*crds[1] + transform[2 ]*crds[2] + transform[3 ];
            ptr[3*p+1] = transform[4]*crds[0] + transform[5]*crds[1] + transform[6 ]*crds[2] + transform[7 ];
            ptr[3*p+2] = transform[8]*crds[0] + transform[9]*crds[1] + transform[10]*crds[2] + transform[11];
            
            debug5 << "VLSV\t\t Point #" << p << "\t coords " << x << '\t' << y << '\t' << z << endl;
         }
         break;
       default:
         debug2 << "VLSV\t\t WARNING: Unknown coordinate system, assuming Cartesian" << endl;
         break;
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

   bool VisitPointMeshReader::readVariable(vlsv::Reader* vlsvReader,MeshMetadata* md,const VariableMetadata& vmd,int domain,void*& output) {
      output = NULL;
      return false;
   }
   
} // namespace vlsvplugin
