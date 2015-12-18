/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2012, 2015 Finnish Meteorological Institute
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

#ifndef MUXML_H
#define MUXML_H

#include <ostream>
#include <map>
#include <list>
#include <utility>
#include <sstream>
#include <vector>

namespace muxml {

   /** Description of a node in XML tree. Root node is the only node 
    * with NULL parent.*/
   struct XMLNode {
      XMLNode* parent;                              /**< Pointer to parent.*/
      std::multimap<std::string,XMLNode*> children; /**< Names of, and pointers to, existing children.*/
      std::map<std::string,std::string> attributes; /**< Node's (XML tag's) attribute,value pairs.*/
      std::string value;                            /**< Value field of this XML tag.*/

      XMLNode(XMLNode* parent);
      ~XMLNode();
   };

   class MuXML {
    public:
      MuXML();
      ~MuXML();
      
      template<typename T> bool addAttribute(XMLNode* node,const std::string& attribName,const T& attribValue);
      template<typename T> XMLNode* addNode(XMLNode* parent,const std::string& nodeName,const T& nodeValue);
      template<typename T> bool changeValue(XMLNode* node,const T& value);
      void clear();
      XMLNode* find(const std::string& nodeName,const XMLNode* node = NULL) const;
      XMLNode* find(const std::string& nodeName,const std::list<std::pair<std::string,std::string> >& attribs,const XMLNode* node=NULL) const;
      
      std::string getAttributeValue(const XMLNode* node,const std::string& attribName);
      void getAttributes(const XMLNode* node,std::map<std::string,std::string>& attribs);
      std::string getLastError() const;
      std::string getNodeValue(const XMLNode* node);
      XMLNode* getRoot() const;

      void print(std::ostream& out,const int& level=0,const XMLNode* node=NULL) const;

      bool read(std::istream& in,XMLNode* parent=NULL,const int& level=0,const char& currentChar=' ');
   
    private:
      std::string errorString;   /**< String describing the last error that has occurred, if any.*/
      XMLNode* root;             /**< Pointer to root node.*/
   };

   // ***** START TEMPLATE FUNCTION DEFINITIONS ***** //

   /** Add a new attribute,value pair to an existing XML tag. The value, given in 
    * parameter attribValue, is converted into a std::string using std::stringstream.
    * If an attribute with name attribName already exists, then its value is
    * overwritten by the new value.
    * @param node XMLNode containing the XML tag.
    * @param attribName Name (unique) of the new attribute.
    * @param attribValue Value of the new attribute.
    * @return If true, attribute was added successfully.*/
   template<typename T> inline
   bool MuXML::addAttribute(XMLNode* node,const std::string& attribName,const T& attribValue) {
      if (node == NULL) return false;
      // Add attribute through stringstream:
      std::stringstream ss;
      ss << attribValue;
      (node->attributes)[attribName] = ss.str();
      return true;
   }

   /** Insert a new (nested) XML tag to an existing tag. New tag's value, 
    * given in parameter nodeValue, is converted into a std::string using std::stringstream.
    * @param parent XMLNode containing the existing tag.
    * @param nodeName Name (non-unique) of the inserted node.
    * @param nodeValue Value of the inserted node.
    * @return Upon success a pointer to an XMLNode containing the new tag is returned.
    * Otherwise a NULL pointer is returned.*/
   template<typename T> inline
   XMLNode* MuXML::addNode(XMLNode* parent,const std::string& nodeName,const T& nodeValue) {
      if (parent == NULL) return NULL;

      // Insert node:
      XMLNode* node = new XMLNode(parent);
      parent->children.insert(std::make_pair(nodeName,node));

      // Convert value using a stringstream:
      std::stringstream ss;
      ss << nodeValue;
      node->value = ss.str();   
      return node;
   }

   /** Set the value of an XML tag in the given node. Passed value 
    * is converted into a std::string using std::stringstream.
    * @param node XMLNode containing the tag.
    * @param value New value for the tag.
    * @return If true, new value was set successfully.*/
   template<typename T> inline
   bool MuXML::changeValue(XMLNode* node,const T& value) {
      if (node == NULL) return false;
      
      // Change value through stringstream:
      std::stringstream ss;
      ss << value;
      node->value = ss.str();

      return true;
   }

} // namespace muxml

#endif
