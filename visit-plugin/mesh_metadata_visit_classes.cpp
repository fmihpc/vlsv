/** This file is part of VLSV file format.
 * 
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

#include <mesh_metadata_visit_classes.h>

using namespace std;

namespace vlsvplugin {
   VisitPointMeshMetadata::VisitPointMeshMetadata(): PointMeshMetadata(),VisitMeshMetadata() { }

   VisitPointMeshMetadata::~VisitPointMeshMetadata() { }

   avtMeshType VisitPointMeshMetadata::getAvtMeshType() const {return AVT_POINT_MESH;}
   
   std::string VisitPointMeshMetadata::getAvtMeshTypeString() const {return "AVT_POINT_MESH";}
   
   const std::string& VisitPointMeshMetadata::getCorrectVlsvMeshType() const {return PointMeshMetadata::getCorrectVlsvMeshType();}

   

   VisitUCDMultiMeshMetadata::VisitUCDMultiMeshMetadata(): UCDMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitUCDMultiMeshMetadata::~VisitUCDMultiMeshMetadata() { }

   avtMeshType VisitUCDMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitUCDMultiMeshMetadata::getCorrectVlsvMeshType() const {return UCDMultiMeshMetadata::getCorrectVlsvMeshType();}
   
   

   VisitQuadMultiMeshMetadata::VisitQuadMultiMeshMetadata(): QuadMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitQuadMultiMeshMetadata::~VisitQuadMultiMeshMetadata() { }

   avtMeshType VisitQuadMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitQuadMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitQuadMultiMeshMetadata::getCorrectVlsvMeshType() const {return QuadMultiMeshMetadata::getCorrectVlsvMeshType();}

   

   VisitUCDAMRMetadata::VisitUCDAMRMetadata(): UCDAMRMetadata(),VisitMeshMetadata() { }

   VisitUCDAMRMetadata::~VisitUCDAMRMetadata() { }

   avtMeshType VisitUCDAMRMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDAMRMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}

   const std::string& VisitUCDAMRMetadata::getCorrectVlsvMeshType() const {return UCDAMRMetadata::getCorrectVlsvMeshType();}
   


   VisitUCDGenericMultiMeshMetadata::VisitUCDGenericMultiMeshMetadata(): UCDGenericMultiMeshMetadata(),VisitMeshMetadata() { }

   VisitUCDGenericMultiMeshMetadata::~VisitUCDGenericMultiMeshMetadata() { }

   avtMeshType VisitUCDGenericMultiMeshMetadata::getAvtMeshType() const {return AVT_UNSTRUCTURED_MESH;}

   std::string VisitUCDGenericMultiMeshMetadata::getAvtMeshTypeString() const {return "AVT_UNSTRUCTURED_MESH";}
   
   const std::string& VisitUCDGenericMultiMeshMetadata::getCorrectVlsvMeshType() const {
      return UCDGenericMultiMeshMetadata::getCorrectVlsvMeshType();
   }

   bool VisitUCDGenericMultiMeshMetadata::read(vlsv::Reader* vlsv,const std::map<std::string,std::string>& attribs) {
      return UCDGenericMultiMeshMetadata::read(vlsv,attribs);
   }
}