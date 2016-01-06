/** This file is part of VLSV file format.
 * 
 *  Copyright 2011, 2012, 2015 Finnish Meteorological Institute
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

#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "muxml.h"

using namespace std;
using namespace muxml;

/** Constructor for XMLNode. If this node is the root node, 
 * then a NULL pointer should be passed in variable parent.
 * @param parent Parent of this node.*/
XMLNode::XMLNode(XMLNode* parent): parent(parent) { }

/** Destructor for XMLNode. Recursively calls delete for all children.*/
XMLNode::~XMLNode() {
   // Delete all children:
   for (auto& it: children) {delete it.second; it.second=NULL;}
}

/** Constructor for MuXML. Creates the root node MuXML::root.*/
MuXML::MuXML() {
   root = new XMLNode(NULL);
}

/** Destructor for MuXML. Deletes the stored XML tree.*/
MuXML::~MuXML() {
   delete root;
   root = NULL;
}

/** Clear the stored XML tree, i.e., delete all existing nodes 
 * and create a new empty tree.*/
void MuXML::clear() {
   delete root;
   root = new XMLNode(NULL);
}

/** Find the XML tag with given name, starting the search from the given node.
 * @param nodeName Name of the searched tag. If a node has multiple children with the same 
 * name that match nodeName, then a pointer to the first matching children is returned. 
 * @param node Pointer to node where the search is started from.
 * @return Pointer to node containing the searched XML tag, if it was found.
 * Otherwise a NULL pointer is returned.*/
XMLNode* MuXML::find(const std::string& nodeName,const XMLNode* node) const {
   if (node == NULL) node = root;

   // Recursive depth-first find. First check if child's name matches the searched name. If it does, return 
   // pointer to child. Otherwise search through child's children to see if it has the searched node.
   for (auto it = node->children.begin(); it != node->children.end(); ++it) {
      if (it->first == nodeName) return it->second;
      XMLNode* tmp = find(nodeName,it->second);
      if (tmp != NULL) return tmp;
   }
   return NULL;
}

/** Find the XML tag with given name and whose attribute,value pairs match with those
 * given in variable attribs. The search is started from the XMLNode given in parameter node.
 * @param nodeName Name of the searched tag. If multiple children match the given nodeName and 
 * attribute,value pairs, then a pointer to the first matching child is returned.
 * @param attribs List containing attribute,value pairs that the searched tag must have.
 * @param node Pointer to node where the search is started from.
 * @return Pointer to a node that matches the search criteria. A NULL pointer is 
 * returned if the search was unsuccessful.*/
XMLNode* MuXML::find(const std::string& nodeName,const std::list<std::pair<std::string,std::string> >& attribs,const XMLNode* node) const {
   if (node == NULL) node = root;

   // Search each children
   for (auto it=node->children.begin(); it!=node->children.end(); ++it) {
      bool matchFound = true;
      
      if (it->first != nodeName) {
	     matchFound = false;
      }
      // Tag name matches, check that attributes match:
      if (matchFound == true) {
	     for (auto jt = attribs.begin(); jt!=attribs.end(); ++jt) {
            const auto tmp = it->second->attributes.find((*jt).first);
            if (tmp == it->second->attributes.end()) {matchFound = false; break;} // attribute name was not found
            if (tmp->second != (*jt).second) {matchFound = false; break;} // attribute value did not match
	     }
         if (matchFound == true) {
            return it->second;
         }
      }
      // Recursively check children's nodes:
      XMLNode* tmp = find(nodeName,attribs,it->second);
      if (tmp != NULL) return tmp;
   }
   return NULL;
}

/** Get the value of the given attribute.
 * @param node Node containing the XML tag.
 * @param attribName Name of the attribute.
 * @return Value of the attribute. If the attribute doesn't exist then an empty string "" is returned.*/
string MuXML::getAttributeValue(const XMLNode* node,const std::string& attribName) {
   if (node == NULL) {
      stringstream ss;
      ss << "ERROR: getAttributeValue called for NULL node in " << __FILE__ << ":" << __LINE__;
      errorString = ss.str();
      return string("");
   }
   const auto it=node->attributes.find(attribName);
   if (it == node->attributes.end()) return "";
   return it->second;
}

/** Get all attribute,value pairs of the given node (XML tag).
 * @param node Node containing the XML tag.
 * @param attribs Map where the attribute,value pairs are copied to.*/
void MuXML::getAttributes(const XMLNode* node,std::map<std::string,std::string>& attribs) {
   if (node == NULL) {
      stringstream ss;
      ss << "ERROR: getAttributes called for NULL node in " << __FILE__ << ":" << __LINE__;
      errorString = ss.str();
      attribs.clear();
      return;
   }
   attribs = node->attributes;
}

/** Get a string description of the latest error, if any, that has occurred.
 * @return Description of the latest error.*/
std::string MuXML::getLastError() const {return errorString;}

/** Get the value of the given node (XML tag).
 * @param node Node containing the XML tag.
 * @return Node's value field.*/
string MuXML::getNodeValue(const XMLNode* node) {
   if (node == NULL) {
      stringstream ss;
      ss << "ERROR: getNodeValue called for NULL node in " << __FILE__ << ":" << __LINE__;
      errorString = ss.str();
      return string("");
   }
   return node->value;
}

/** Get pointer to root node.
 * @return Pointer to root node.*/
XMLNode* MuXML::getRoot() const {return root;}

/** Recursively print the XML tree contents to given output stream.
 * @param out Output stream.
 * @param level How deep we are into the XML tree.
 * @param node Pointer to first printed node.*/
void MuXML::print(std::ostream& out,const int& level,const XMLNode* node) const {
   const int tab = 3;
   if (node == NULL) {
      node = root;
   }
   for (auto it=node->children.begin(); it!=node->children.end(); ++it) {
      // Indent
      for (auto i=0; i<level; ++i) out << ' ';
      
      // Write child's name and its attributes:
      out << '<' << it->first;
      for (auto jt=it->second->attributes.begin(); jt!=it->second->attributes.end(); ++jt) {
         out << ' ' << jt->first << "=\"" << jt->second << "\"";
      }

      // Write child's value:
      out << ">" << it->second->value;
      
      // Call print for the child:
      if (it->second->children.size() > 0) {
         out << endl;
         print(out,level+tab,it->second);
         for (int i=0; i<level; ++i) out << ' ';
      }
      
      // Write child's end tag:
      out << "</" << it->first << '>' << endl;
   }
}

/** Read in an XML tree from the given input stream. This function will recursively call itself.
 * @param in Input stream containing XML tags.
 * @param parent New XML tags are inserted into the XML tree as children of this node.
 * @param level The depth of parent in the XML tree.
 * @param currentChar Last character that was read from the input stream.
 * @return If true, XML tags were successfully read and inserted to the XML tree.
 * If an error occurred, false is returned and error description is written to 
 * MuXML::errorString.
 * @see MuXML::getLastError().*/
bool MuXML::read(std::istream& in,XMLNode* parent,const int& level,const char& currentChar) {
   in >> noskipws;
   if (parent == NULL) parent = root;
   if (in.good() == false) {
      stringstream ss;
      ss << "ERROR: Input stream is not good in " __FILE__ << ":" << __LINE__;
      errorString = ss.str();
      return false;
   }
   
   int index = 0;
   bool success = true;
   char c = currentChar;
   char buffer[1024];
   do {
      // Skip empty chars:
      while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
      
      // Start to read tag name and its attributes.
      // Tag is of the form '<tagname attrib1="value1" attrib2="value2" ...> tagvalue </tagname>'
      // (apostrophes not included). Note that tag can have any number of attribute="value" pairs.
      if (c == '<') { // tag name begins with '<'
	     in >> c; // forward from '<'
         // Skip (possible) empty spaces
	     while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 
         // If the next character is '/', then we are reading the end tag:
	     if (c == '/') {
	        while (in.eof() == false && c != '>') in >> c;
	        in >> c; // forward from '>'
	        if (in.good() == false) {
               stringstream ss;
               ss << "ERROR: input stream not good in " << __FILE__ << ":" << __LINE__;
               errorString = ss.str();
               success = false;
            }
	        return success;
	     }

         // Copy tag name to buffer:
	     index = 0;
	     while (c != ' ' && c != '>') {
	        buffer[index] = c;
	        ++index;
	        in >> c;
	        if (in.good() == false) {
               stringstream ss;
               ss << "ERROR: input stream not good in " << __FILE__ << ":" << __LINE__;
               errorString = ss.str();
               success = false;
               break;
            }
	     }
	     buffer[index] = '\0';
	     XMLNode* node = addNode(parent,buffer,"");

	     // Skip (possible) empty spaces between tag name and its attribute,value pairs:
	     while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;

         // Valid header cannot end until end tag has been read:
	     if (in.good() == false) {
            stringstream ss;
            ss << "ERROR: input stream not good in " << __FILE__ << ":" << __LINE__;
            errorString = ss.str();
            success = false;
            break;
         }

	     // If next char is '>' the tag ends. Otherwise read all attribute,value pairs:
	     if (c != '>') {
            // Attributes are given in form 'attribName="attribValue"' (parentheses included, apostrophes not included)
	        while (c != '>') {
               // Skip (possible) empty spaces
	           while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	       
               // Read attribute name to buffer:
	           index = 0;
	           while (c != '=' && in.good() == true) {
		          buffer[index] = c;
		          ++index;
		          in >> c;
	           }
	           buffer[index] = '\0';
	           string attribName = buffer;

               // Check that input stream did not end prematurely:
               if (in.good() == false) {
                  stringstream ss;
                  ss << "ERROR: Input stream ended before attribute '" << attribName << "' was read in ";
                  ss << __FILE__ << ":" << __LINE__;
                  errorString = ss.str();
                  return false;
               }

               // Next characters in stream are ="
	           in >> c; // Forward from '='
	           in >> c; // Forward from '"'
	           
               // Check that input stream did not end prematurely:
               if (in.good() == false) {
                  stringstream ss;
                  ss << "ERROR: Input stream ended before attribute '" << attribName << "' value was read in ";
                  ss << __FILE__ << ":" << __LINE__;
                  errorString = ss.str();
                  return false;
               }

               // Copy attribute value to buffer:
               index = 0;
	            while (c != '"' && in.good() == true) {
		            buffer[index] = c;
		            ++index;
		            in >> c;
	            }
               buffer[index] = '\0';
	            string attribValue = buffer;

               // Check that input stream did not end prematurely:
               if (in.good() == false) {
                  stringstream ss;
                  ss << "ERROR: Input stream ended before attribute '" << attribName << "' value was read in ";
                  ss << __FILE__ << ":" << __LINE__;
                  errorString = ss.str();
                  return false;
               }

	            in >> c; // Forward from '"' character that ends attribute's value

               // Check that input stream did not end prematurely, at least end tag 
               // should still be left in stream:
	            if (in.good() == false) {
                  stringstream ss;
                  ss << "ERROR: input stream not good in " << __FILE__ << ":" << __LINE__;
                  errorString = ss.str();
                  success = false; 
                  break;
               }

               // Insert the attribute,value pair to node:
	           addAttribute(node,attribName,attribValue);
	        }
            // When we reach this point the current character in stream is '>'
	        in >> c;
	     } else {
	        in >> c;
	     }

         // Check that input stream did not end prematurely:
         if (in.good() == false) {
            stringstream ss;
            ss << "ERROR: input stream ended before tag value read in " << __FILE__ << ":" << __LINE__;
            errorString = ss.str();
            return false;
         }

         // Skip (possible) empty spaces before tag value:
	     while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
	 
         // Check that input stream did not end prematurely:
         if (in.good() == false) {
            stringstream ss;
            ss << "ERROR: input stream ended before tag value read in " << __FILE__ << ":" << __LINE__;
            errorString = ss.str();
            return false;
         }

	     // Read tag's value:
	     index = 0;
	     while ((c != ' ' && c != '\t' && c != '\n' && c != '<') && in.good() == true) {
	        buffer[index] = c;
	        ++index;
	        in >> c;
	     }
         buffer[index] = '\0';
	     changeValue(node,buffer);

         // Check that input stream did not end before tag value was read:
	     if (in.good() == false) {
            stringstream ss;
            ss << "ERROR: input stream ended before tag value was read in " << __FILE__ << ":" << __LINE__;
            errorString = ss.str();
            success = false;
            break;
         }

         // Skip (possible) empty characters between tag value and end tag:
	     while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;

         // If the next character is '<' we have encountered a nested tag, read it recursively:
	     if (c == '<') {
	        if (read(in,node,level+1,c) == false) {
               success = false;
            }
	        in >> c;
	     }
	     while ((c == ' ' || c == '\t' || c == '\n') && in.good() == true) in >> c;
      }

      // Valid input stream cannot end if we are reading a nested tag:
      if (level > 0 && in.good() == false) {
         stringstream ss;
         ss << "ERROR: input stream not good at level " << level << " in " << __FILE__ << ":" << __LINE__;
         errorString = ss.str();
         return false;
      }
   } while (success == true);
   return true;
}
