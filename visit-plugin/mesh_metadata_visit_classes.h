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

#pragma once

#ifndef MESH_METADATA_VISIT_CLASSES_H
#define MESH_METADATA_VISIT_CLASSES_H

//#include <mesh_metadata.h>
#include <mesh_metadata_visit.h>
#include <mesh_metadata_visit_point.h>
#include <mesh_metadata_visit_ucd_multi.h>
#include <mesh_metadata_visit_quad_multi.h>
#include <mesh_metadata_visit_ucd_amr.h>
#include <mesh_metadata_visit_ucd_generic_multi.h>

namespace vlsvplugin {

   class VisitPointMeshMetadata: public PointMeshMetadata,public VisitMeshMetadata {
   public:
      VisitPointMeshMetadata();
      virtual ~VisitPointMeshMetadata();

      virtual avtMeshType getAvtMeshType() const;
      virtual std::string getAvtMeshTypeString() const;

   protected:

   };

   class VisitUCDMultiMeshMetadata: public UCDMultiMeshMetadata,public VisitMeshMetadata {
   public:
      VisitUCDMultiMeshMetadata();
      virtual ~VisitUCDMultiMeshMetadata();

      virtual avtMeshType getAvtMeshType() const;
      virtual std::string getAvtMeshTypeString() const;

   protected:

   };

   class VisitQuadMultiMeshMetadata: public QuadMultiMeshMetadata,public VisitMeshMetadata {
   public:
      VisitQuadMultiMeshMetadata();
      virtual ~VisitQuadMultiMeshMetadata();

      virtual avtMeshType getAvtMeshType() const;
      virtual std::string getAvtMeshTypeString() const;

   protected:

   };

   class VisitUCDAMRMetadata: public UCDAMRMetadata,public VisitMeshMetadata {
   public:
      VisitUCDAMRMetadata();
      virtual ~VisitUCDAMRMetadata();

      virtual avtMeshType getAvtMeshType() const;
      virtual std::string getAvtMeshTypeString() const;

   protected:

   };

   class VisitUCDGenericMultiMeshMetadata: public UCDGenericMultiMeshMetadata,public VisitMeshMetadata {
   public:
      VisitUCDGenericMultiMeshMetadata();
      virtual ~VisitUCDGenericMultiMeshMetadata();

      virtual avtMeshType getAvtMeshType() const;
      virtual std::string getAvtMeshTypeString() const;

   protected:

   };

} // namespace vlsvplugin

#endif