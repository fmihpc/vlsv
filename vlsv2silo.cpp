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

#include <cstdlib>
#include <iostream>
#include <stdint.h>
#include <cmath>
#include <list>
#include <silo.h>
#include <sstream>
#include <dirent.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <iomanip>
#include <cstring>

#ifdef PROFILE
   #include <profile.hpp>
   int convertQuadMeshTotalID;
   //int quadMeshID;
   //int quadVariablesID;
   //int quadMeshReadID;
#endif

#include "vlsv_reader.h"

using namespace std;

// When making option lists SILO seems to only store a pointer to 
// variable giving the option's value. The pointers must be valid 
// when data is written to file. Thus, we need to define some global 
// variables here that store values put in option lists.

static int timestep = 0;
static float time_float = 0.0;
static double time_double = 0.0;

static string label_x = "x-coordinate"; /**< Default name for x-coordinate.*/
static string label_y = "y-coordinate"; /**< Default name for y-coordinate.*/
static string label_z = "z-coordinate"; /**< Default name for z-coordinate.*/
static string units_x = "m";            /**< Default unit for x-coordinate.*/
static string units_y = "m";            /**< Default unit for y-coordinate.*/
static string units_z = "m";            /**< Default unit for z-coordinate.*/

static float bbox[6];           // Simulation domain bounding box (x_min,y_min,z_min,dx,dy,dz) parameters
static int maxHash = 2097152-1; // =2^21, the (i,j,k) indices are packed into a 64-bit integer (3*21=63).

static DBfile* fileptr = NULL; // Pointer to file opened by SILO

static set<string> directories; /**< List of directories created to output SILO file. If one tries to 
				 * create a directory more than once, SILO will throw errors to the 
				 * console. Similarly attempting to cd into a non-existing directory 
				 * throws errors to console. In order to avoid these errors, currently 
				 * existing directories are stored in this set.*/

struct CurveData {
   string xlabel;
   string ylabel;
   string xunit;
   string yunit;
   map<double,double> data;
   
   CurveData();
};

CurveData::CurveData(): xlabel(""),ylabel(""),xunit(""),yunit("") { }

static map<string,CurveData> curves; /**< Container for all curve data stored to VLSV file(s).
				      * map<double,double> is for sorting data by ascending x values.*/

/** Convert an integer value written to a char array to int64_t. If the integer value 
 * in char array is uint64_t the conversion to int64_t may lead to unexpected results.
 * @param ptr Pointer to char array containing an integer value.
 * @param dataType Is the integer value signed or unsigned?
 * @param dataSize Byte size of the stored integer value.
 * @return Stored value converted to int64_t.*/
int64_t convInt(const char* ptr,const vlsv::datatype::type& dataType,const uint64_t& dataSize) {
   switch (dataType) {
    // Convert signed integer to int64_t:
    case (vlsv::datatype::INT):
      switch(dataSize) {
       case (sizeof(int8_t)):
	 return *reinterpret_cast<const int8_t*>(ptr);
	 break;
       case (sizeof(int16_t)):
	 return *reinterpret_cast<const int16_t*>(ptr);
	 break;
       case (sizeof(int32_t)):
	 return *reinterpret_cast<const int32_t*>(ptr);
	 break;
       case (sizeof(int64_t)):
	 return *reinterpret_cast<const int64_t*>(ptr);
	 break;
       default:
	 cerr << "ERROR Unsupported signed integer datatype passed to convInt!" << endl;
	 exit(1);
	 break;
      }
      break;
    // Convert unsigned integer to int64_t:  
    case (vlsv::datatype::UINT):
      switch(dataSize) {
       case (sizeof(uint8_t)):
	 return *reinterpret_cast<const uint8_t*>(ptr);
	 break;
       case (sizeof(uint16_t)):
	 return *reinterpret_cast<const uint16_t*>(ptr);
	 break;
       case (sizeof(uint32_t)):
	 return *reinterpret_cast<const uint32_t*>(ptr);
	 break;
       case (sizeof(uint64_t)):
	 return *reinterpret_cast<const uint64_t*>(ptr);
	 break;
       default:
	 cerr << "ERROR Unsupported unsigned integer datatype passed to convInt!" << endl;
	 exit(1);
	 break;
      }
      break;
    // Unknown 
    default:
      cerr << "ERROR Unsupported datatype passed to convInt!" << endl;
      exit(1);
      break;
   }
}

/** Convert an integer value written to a char array to uint64_t. If the integer value
 * in char array is signed, the conversion to uint64_t may lead to unexpected results.
 * @param ptr Pointer to char array containing an integer value.
 * @param dataType Is the integer value signed or unsigned?
 * @param dataSize Byte size of the stored integer value.
 * @return Stored value converted to uint64_t.*/
uint64_t convUInt(const char* ptr,const vlsv::datatype::type& dataType,const uint64_t& dataSize) {
   switch (dataType) {
    // Convert signed integer to uint64_t:
    case (vlsv::datatype::INT):
      switch (dataSize) {
       case sizeof(int8_t):
	 return *reinterpret_cast<const int8_t*>(ptr);
	 break;
       case sizeof(int16_t):
	 return *reinterpret_cast<const int16_t*>(ptr);
	 break;
       case sizeof(int32_t):
	 return *reinterpret_cast<const int32_t*>(ptr);
	 break;
       case sizeof(int64_t):
	 return *reinterpret_cast<const int64_t*>(ptr);
	 break;
       default:
	 cerr << "ERROR Unsupported signed integer datatype passed to convUInt!" << endl;
	 exit(1);
	 break;
      }
      break;
    // Convert unsigned integer to uint64_t:
    case (vlsv::datatype::UINT):
      switch (dataSize) {
       case sizeof(uint8_t):
	 return *reinterpret_cast<const uint8_t*>(ptr);
	 break;
       case sizeof(uint16_t):
	 return *reinterpret_cast<const uint16_t*>(ptr);
	 break;
       case sizeof(uint32_t):
	 return *reinterpret_cast<const uint32_t*>(ptr);
	 break;
       case sizeof(uint64_t):
	 return *reinterpret_cast<const uint64_t*>(ptr);
	 break;
       default:
	 cerr << "ERROR Unsupported unsigned integer datatype passed to convUInt!" << endl;
	 exit(1);
	 break;
      }
      break;
    // Unsupported datatype:
    default:
      cerr << "ERROR: Unsupported datatype passed to convUInt!" << endl; 
      exit(1);
      break;
   }
}

/** Get a SILO datatype corresponding to given VLSV datatype.
 * @param dataType Datatype. Either signed or unsigned integer, or floating point value.
 * @param dataSize Byte size of datatype.
 * @return Corresponding SILO datatype or -1 if conversion was not possible.*/
int SiloType(const vlsv::datatype::type& dataType,const uint64_t& dataSize) {
   switch (dataType) {
    case vlsv::datatype::UNKNOWN:
      return -1;
      break;
    case vlsv::datatype::INT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case vlsv::datatype::UINT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case vlsv::datatype::FLOAT:
      if (dataSize == 4) return DB_FLOAT;
      else if (dataSize == 8) return DB_DOUBLE;
      else return -1;
      break;
   }
   return -1;
}

const float EPS = 1.0e-6;

// ************************************************ //
// ***** FUNCTIONS RELATED TO MESH_QUAD_MULTI ***** //
// ************************************************ //

/** Struct for storing (i,j,k) indices of nodes in the mesh. This 
 * struct is, in practice, used to eliminate duplicate nodes from the 
 * mesh pieces written to SILO file.*/
struct NodeIndices {
   int32_t i;            /**< i-index of the node.*/
   int32_t j;            /**< j-index of the node.*/
   int32_t k;            /**< k-index of the node.*/

   /** Constructor.
    * @param i i-index of the node.
    * @param j j-index of the node.
    * @param k k-index of the node.*/
   NodeIndices(int32_t i,int32_t j,int32_t k): i(i),j(j),k(k) { }
};

/** Comparator object for struct NodeIndices, required for being
 * able to store NodeIndices to unordered_map.*/
struct NodesAreEqual {
   /** Compare given NodeIndices objects for equality. Objects 
    * are equal if their (i,j,k) indices are equal.
    * @param first First NodeIndices object to compare.
    * @param second Second NodeIndices object to compare.
    * @return If true, given objects are equal.*/
   bool operator()(const NodeIndices& first,const NodeIndices& second) const {
      if (first.i != second.i) return false;
      if (first.j != second.j) return false;
      if (first.k != second.k) return false;
      return true;
   }
};

/** Hash function implementation for object NodeIndices, required for 
 * being able to store NodeIndices to unordered_map.*/
struct NodeHash {
   /** Calculate hash value for given NodeIndices object.
    * @param Node NodeIndices object whose hash value is to be calculated.
    * @return Calculated hash value.*/
   uint64_t operator()(const NodeIndices& node) const {
      uint64_t result = 0;
      uint64_t tmp = node.i % maxHash;
      result = (result | tmp);
      
      tmp = node.j % maxHash;
      result = (result | (tmp << 21));
      
      tmp = node.k % maxHash;
      result = (result | (tmp << 42));
      return result;
   }
};

template<typename T>
void insertNodeCrds(T i,T j,T k,vector<float>& xcrds,vector<float>& ycrds,vector<float>& zcrds) {
   xcrds.push_back(bbox[0] + i*bbox[3]);
   ycrds.push_back(bbox[1] + j*bbox[4]);
   zcrds.push_back(bbox[2] + k*bbox[5]);
}

/** Iterate over all cells (local + ghosts) in the given multimesh piece, and eliminate 
 * duplicate nodes from the resulting node list. Additionally, unique (x,y,z) coordinate 
 * tuples are written to given vectors.
 * @param amount Number of cells (local + ghosts) in the multimesh piece.
 * @param vectorSize Size of data vector in VLSV file (should be three).
 * @param ptrMesh Pointer to array containing (i,j,k) index tuples that were read from VLSV file.
 * @param zoneList Vector in which the node list is written.
 * @param xcrds Vector in which x-coordinates of unique nodes are written.
 * @param xcrds Vector in which y-coordinates of unique nodes are written.
 * @param xcrds Vector in which z-coordinates of unique nodes are written.*/
template<typename T>
void eliminateDuplicateNodes(uint64_t amount,uint64_t vectorSize,T* ptrMesh,vector<int>& zoneList,
			     vector<float>& xcrds,vector<float>& ycrds,vector<float>& zcrds) {
   zoneList.clear();
   xcrds.clear();
   ycrds.clear();
   zcrds.clear();
   
   uint32_t counter = 0;
   unordered_map<NodeIndices,int,NodeHash,NodesAreEqual> nodeIndices;
   pair<unordered_map<NodeIndices,int,NodeHash,NodesAreEqual>::iterator,bool> result;
   
   // Iterate over all cells and attempt to insert the eight nodes associated with 
   // it to unordered_map nodeIndices. If the insertion succeeds, the node and its coordinates 
   // is also inserted to vectors zoneList,xcrds,ycrds,and zcrds. In practice this eliminates 
   // duplicate nodes.
   for (uint64_t cell=0; cell<amount; ++cell) {
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+0,ptrMesh[1]+1,ptrMesh[2]+0),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+0,ptrMesh[1]+1,ptrMesh[2]+0,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+0,ptrMesh[1]+0,ptrMesh[2]+0),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+0,ptrMesh[1]+0,ptrMesh[2]+0,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+1,ptrMesh[1]+0,ptrMesh[2]+0),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+1,ptrMesh[1]+0,ptrMesh[2]+0,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+1,ptrMesh[1]+1,ptrMesh[2]+0),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+1,ptrMesh[1]+1,ptrMesh[2]+0,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+0,ptrMesh[1]+1,ptrMesh[2]+1),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+0,ptrMesh[1]+1,ptrMesh[2]+1,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+0,ptrMesh[1]+0,ptrMesh[2]+1),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+0,ptrMesh[1]+0,ptrMesh[2]+1,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+1,ptrMesh[1]+0,ptrMesh[2]+1),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+1,ptrMesh[1]+0,ptrMesh[2]+1,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      result = nodeIndices.insert(make_pair(NodeIndices(ptrMesh[0]+1,ptrMesh[1]+1,ptrMesh[2]+1),counter));
      if (result.second == true)        {insertNodeCrds(ptrMesh[0]+1,ptrMesh[1]+1,ptrMesh[2]+1,xcrds,ycrds,zcrds); ++counter;}
      zoneList.push_back(result.first->second);
      ptrMesh += vectorSize;
   }
}

// ****************************************** //
// ***** FUNCTIONS RELATED TO MESH_QUAD ***** //
// ****************************************** //

template<typename REAL> struct NodeCrd {
   REAL x;
   REAL y;
   REAL z;
   NodeCrd(const REAL& x,const REAL& y,const REAL& z): x(x),y(y),z(z) { }
   NodeCrd(char* x0,char* y0,char* z0,
	   char* dx,char* dy,char* dz) {
      // Calculate node coordinates:
      x = *reinterpret_cast<REAL*>(x0) + *reinterpret_cast<REAL*>(dx);
      y = *reinterpret_cast<REAL*>(y0) + *reinterpret_cast<REAL*>(dy);
      z = *reinterpret_cast<REAL*>(z0) + *reinterpret_cast<REAL*>(dz);
      
      // Flush very small coordinate values to zero:
      if (fabs(x) < EPS) x = 0.0;
      if (fabs(y) < EPS) y = 0.0;
      if (fabs(z) < EPS) z = 0.0;
   }
  
   void xcopy(char* target,REAL scaling=1.0) const {
      const REAL tmp = x*scaling;
      const char* ptr = reinterpret_cast<const char*>(&tmp);
      for (unsigned int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
   }

   void ycopy(char* target,REAL scaling=1.0) const {
      const REAL tmp = y*scaling;
      const char* ptr = reinterpret_cast<const char*>(&tmp);
      for (unsigned int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
   }

   void zcopy(char* target,REAL scaling=1.0) const {
      const REAL tmp = z*scaling;
      const char* ptr = reinterpret_cast<const char*>(&tmp);
      for (unsigned int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
   }

   bool comp(const NodeCrd<REAL>& n) const {
      REAL EPS1,EPS2,EPS;
      EPS1 = 1.0e-6 * fabs(x);
      EPS2 = 1.0e-6 * fabs(n.x);
      if (x == 0.0) EPS1 = 1.0e-7;
      if (n.x == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(x - n.x) > EPS) return false;
      
      EPS1 = 1.0e-6 * fabs(y);
      EPS2 = 1.0e-6 * fabs(n.y);
      if (y == 0.0) EPS1 = 1.0e-7;
      if (n.y == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(y - n.y) > EPS) return false;
      
      EPS1 = 1.0e-6 * fabs(z);
      EPS2 = 1.0e-6 * fabs(n.z);
      if (z == 0.0) EPS1 = 1.0e-7;
      if (n.z == 0.0) EPS2 = 1.0e-7;
      EPS = max(EPS1,EPS2);
      if (fabs(z - n.z) > EPS) return false;
      return true;
   }
};

struct NodeComp {
    bool operator()(const NodeCrd<double>& a,const NodeCrd<double>& b) const {
      if (a.comp(b) == true) return false;
      double EPS = 0.5e-5 * (fabs(a.z) + fabs(b.z));      
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.y) + fabs(b.y));      
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.x) + fabs(b.x));      
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }
    
   bool operator()(const NodeCrd<float>& a,const NodeCrd<float>& b) const {
      if (a.comp(b) == true) return false;
      float EPS = 0.5e-5 * (fabs(a.z) + fabs(b.z));      
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.y) + fabs(b.y));      
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.x) + fabs(b.x));      
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }
   
   bool operator()(const NodeCrd<long double>& a,const NodeCrd<long double>& b) const {
      if (a.comp(b) == true) return false;
      float EPS = 0.5e-5 * (fabs(a.z) + fabs(b.z));
      if (a.z > b.z + EPS) return false;
      if (a.z < b.z - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.y) + fabs(b.y));
      if (a.y > b.y + EPS) return false;
      if (a.y < b.y - EPS) return true;
      
      EPS = 0.5e-5 * (fabs(a.x) + fabs(b.x));
      if (a.x > b.x + EPS) return false;
      if (a.x < b.x - EPS) return true;
      return false;
   }
};

// ****************************************** //
// ***** MISCELLANEOUS HELPER FUNCTIONS ***** //
// ****************************************** //

/** Create a SILO option list. By default simulation time and time step are inserted to 
 * the returned option list, as well as memory for extra options to be added later.
 * @param vlsvReader VLSV reader, used for reading input file.
 * @param N_extra How many additional options are needed.
 * @param insertTimes If true, simulation time and time step are inserted to the option list.
 * @return Pointer to created option list or a NULL pointer.*/
DBoptlist* getOptionList(vlsv::Reader& vlsvReader,const int& N_extra=0,const bool& insertTimes=true) {
   DBoptlist* optlist = NULL;
   if (insertTimes == false) {
      if (N_extra > 0) optlist = DBMakeOptlist(N_extra);
      return optlist;
   }
   
   list<pair<string,string> > attribs;
   uint64_t arraySize,vectorSize,dataSize;
   vlsv::datatype::type dataType;

   // Check if time value has been given in the file. If not, then return 
   // immediately with an empty list:
   int N_options = N_extra;
   set<string> parameterNames;
   if (vlsvReader.getUniqueAttributeValues("PARAMETER","name",parameterNames) == false) {
      if (N_options == 0) return NULL;
      return DBMakeOptlist(N_options);
   }
   if (parameterNames.find("time") != parameterNames.end()) ++N_options;
   if (parameterNames.find("timestep") != parameterNames.end()) ++N_options;
   
   if (N_options == 0) return NULL;
   optlist = DBMakeOptlist(N_options);
   
   // Read time value from file:
   if (parameterNames.find("time") != parameterNames.end()) {
      attribs.push_back(make_pair("name","time"));
      vlsvReader.getArrayInfo("PARAMETER",attribs,arraySize,vectorSize,dataType,dataSize);      
      char* valbuffer = new char[arraySize*vectorSize*dataSize];
      vlsvReader.readArray("PARAMETER",attribs,0,arraySize,valbuffer);
   
      // Write time to SILO file either as float or double, 
      // depending on the format given in the file:
      if (dataSize == sizeof(float)) {
	 time_float = *reinterpret_cast<float*>(valbuffer);
	 DBAddOption(optlist,DBOPT_TIME,&time_float);
      } else if (dataSize == sizeof(double)) {
	 time_double = *reinterpret_cast<double*>(valbuffer);
	 DBAddOption(optlist,DBOPT_DTIME,&time_double);
      } else if (dataSize == sizeof(long double)) {
	 time_double = *reinterpret_cast<long double*>(valbuffer);
	 DBAddOption(optlist,DBOPT_DTIME,&time_double);
      } else {
	 cerr << "TIME given in unknown datatype, aborting" << endl; exit(1);
      }
      delete [] valbuffer; valbuffer = NULL;
   }
   
   if (parameterNames.find("timestep") != parameterNames.end()) {
      attribs.clear();
      attribs.push_back(make_pair("name","timestep"));
      vlsvReader.getArrayInfo("PARAMETER",attribs,arraySize,vectorSize,dataType,dataSize);
      char* valbuffer = new char[arraySize*vectorSize*dataSize];
      vlsvReader.readArray("PARAMETER",attribs,0,arraySize,valbuffer);
      
      if (dataType == vlsv::datatype::INT) timestep = convInt(valbuffer,dataType,dataSize);
      else if (dataType == vlsv::datatype::UINT) timestep = convUInt(valbuffer,dataType,dataSize);
      else {
	 cerr << "timestep given in unsupported datatype!" << endl;
      }
      
      DBAddOption(optlist,DBOPT_CYCLE,&timestep);
      delete [] valbuffer; valbuffer = NULL;
   }   
   return optlist;
}

/** Parse coordinate names and units from given XML tag attribute list, and 
 * add them to SILO option list. Note that the option list must have space for 
 * six options.
 * @param optlist Option list where coordinate names and units are added.
 * @parma attribsOut Map containing XML tag attributes.*/
void parseCoordinateNames(DBoptlist* optlist,map<string,string>& attribsOut) {
   // Default values for coordinate names & units:
   label_x = "x-coordinate";
   label_y = "y-coordinate";
   label_z = "z-coordinate";
   units_x = "m";
   units_y = "m";
   units_z = "m";
   
   // Coordinate names & units in XML tag overwrite the default ones:
   if (attribsOut.find("xlabel") != attribsOut.end()) label_x = attribsOut["xlabel"];
   if (attribsOut.find("ylabel") != attribsOut.end()) label_y = attribsOut["ylabel"];
   if (attribsOut.find("zlabel") != attribsOut.end()) label_z = attribsOut["zlabel"];
   if (attribsOut.find("xunit") != attribsOut.end()) units_x = attribsOut["xunit"];
   if (attribsOut.find("yunit") != attribsOut.end()) units_y = attribsOut["yunit"];
   if (attribsOut.find("zunit") != attribsOut.end()) units_z = attribsOut["zunit"];

   // Add coordinate names & units to option list:
   DBAddOption(optlist,DBOPT_XLABEL,const_cast<char*>(label_x.c_str()));
   DBAddOption(optlist,DBOPT_YLABEL,const_cast<char*>(label_y.c_str()));
   DBAddOption(optlist,DBOPT_ZLABEL,const_cast<char*>(label_z.c_str()));
   DBAddOption(optlist,DBOPT_XUNITS,const_cast<char*>(units_x.c_str()));
   DBAddOption(optlist,DBOPT_YUNITS,const_cast<char*>(units_y.c_str()));
   DBAddOption(optlist,DBOPT_ZUNITS,const_cast<char*>(units_z.c_str()));
}

/** Create a directory to output SILO file for the given variable.
 * @param vlsvReader VLSV reader used to read input file.
 * @param meshName Name of the mesh, same as the name of resulting multimesh.
 * @param varName Name of the variable.
 * @param outputDir Variable in which the name of the output directory is copied.*/
void createVariableDirectory(vlsv::Reader& vlsvReader,const string& meshName,const string& varName,string& outputDir) {
   // By default variables are written to root directory:
   outputDir = "/";
   
   // Search for XML attribute "dir" for the given variable. If it is defined, 
   // create a subdirectory with a same name as the attribute value. Otherwise exit:
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("mesh",meshName));
   attributes.push_back(make_pair("name",varName));
   map<string,string> arrayAttributes;
   
   if (vlsvReader.getArrayAttributes("VARIABLE",attributes,arrayAttributes) == false) return;
   map<string,string>::const_iterator it = arrayAttributes.find("dir");
   if (it == arrayAttributes.end()) return;
   
   const string& dir = it->second;
   set<string>::const_iterator dirExists = directories.find(dir);
   if (dirExists == directories.end()) {
      string tmp;
      string s = dir;
      string current = "";
      size_t pos = s.find('/');
      while (pos != string::npos) {
	 tmp = s.substr(0,pos);
	 if (pos != 0) current = current + '/' + tmp;
	 
	 // This exception will make this work correctly if the
	 // first character in directory name is '/':
	 if (tmp.size() == 0) {
	    s = s.substr(pos+1);
	    pos = s.find('/');
	    continue;
	 }
	 
	 // If current directory does not exist, insert it:
	 if (directories.find(current) == directories.end()) {
	    directories.insert(current);
	    DBMkDir(fileptr,const_cast<char*>(current.c_str()));
	 }
	 
	 s = s.substr(pos+1);
	 pos = s.find('/');
      }
      tmp = s.substr(0,pos);
      if (tmp.size() > 0) current = current + '/' + tmp;
      
      // If directory name was given without a trailing '/', e.g. "/dir1/dir2",
      // this if inserts the final directory:
      if (tmp.size() > 0) {
	 if (directories.find(current) == directories.end()) {
	    directories.insert(current);
	    DBMkDir(fileptr,const_cast<char*>(current.c_str()));
	 }
      }
      outputDir = current;
   } else {
      // Variable directory already exists, just copy its name:
      outputDir = *dirExists;
   }
}

/** Write a multimesh variable from VLSV file to output SILO file.
 * @param vlsvReader VLSV reader used to read the input file.
 * @param meshName Name of the mesh, not mesh piece. Multimesh consists of one or more mesh pieces 
 * each having a unique name.
 * @param varName Name of the variable.
 * @param meshDir Directory in the SILO file where the variable data is written.
 * @param meshNumber Which mesh piece variable data is to be converted.
 * @param variableOffsets Array containing offsets to variable data in VLSV file for each mesh piece.
 * @param N_cells Total number of cells in given mesh piece (local + ghosts).
 * @param N_ghosts Number of ghost cells in given mesh piece.
 * @param ghostLocalIDs For each ghost cell, index to variable data in mesh piece containing the data.
 * @param ghostDomains For each ghost cell, number of mesh piece containing the data.
 * @param fullName Variable in which the full name of the converted multimesh variable is copied (does not include filename).
 * @return If true, variable was converted successfully.*/
bool convertMultimeshVariable(vlsv::Reader& vlsvReader,const string& meshName,const string& varName,const string& meshDir,
			      int meshNumber,unsigned int* variableOffsets,uint64_t N_cells,uint64_t N_ghosts,
			      unsigned int* ghostLocalIDs,unsigned int* ghostDomains,string& fullName) {
   bool success = true;
   fullName = "";
   
   // Writing a unstructured grid variable is a rather straightforward process. The
   // only complication here is that some of the variables are actually vectors, i.e.
   // vectorSize > 1 (vectorSize == 1 for scalars). Format in which vectors are stored in VLSV
   // differ from format in which they are written to SILO files.
   vlsv::datatype::type dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("mesh",meshName));
   attributes.push_back(make_pair("name",varName));
   if (vlsvReader.getArrayInfo("VARIABLE",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      // This will exit if this variable does not belong to given mesh
      return false;
   }   
   // Check that variable is a scalar or 3D vector:
   if (vectorSize != 1 && vectorSize != 3) return false;

   // Read variable data. We do not actually need to care if
   // the data is given as floats, doubles, or long doubles. We
   // just need to tell SILO what kind of data the buffer contains.
   // Here we read all data local to this domain (N_cells-N_ghosts cells):
   char* buffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("VARIABLE",attributes,0,arraySize,buffer) == false) {
      success = false;
      delete [] buffer; buffer = NULL;
      return success;
   }
   
   // Vector variables need to be copied to temporary arrays before
   // writing to SILO file:
   char** components = new char*[vectorSize];
   for (uint64_t i=0; i<vectorSize; ++i) components[i] = new char[N_cells*dataSize];
   
   // Copy vector variables from local cells:
   for (uint64_t block=0; block<N_cells-N_ghosts; ++block) {
      const uint64_t compOffset  = block*dataSize;
      const uint64_t localOffset = (variableOffsets[meshNumber] + block)*vectorSize*dataSize;
      for (uint64_t j=0; j<vectorSize; ++j) {
	 const uint64_t bufferOffset = localOffset + j*dataSize;
	 for (uint64_t i=0; i<dataSize; ++i) {
	    components[j][compOffset+i] = buffer[bufferOffset + i];
	 }
      }
   }
   
   // Copy vector variables from ghost cells:
   const uint64_t ghostOffset = N_cells-N_ghosts;
   for (uint64_t block=0; block<N_ghosts; ++block) {
      // Offset to components[j] where ghost cell value is stored:
      const uint64_t compOffset = (ghostOffset+block)*dataSize;
      // Offset to neighbour domain cell where ghost cell value is read:
      const uint64_t nbrOffset  = (variableOffsets[ghostDomains[block]] + ghostLocalIDs[block])*vectorSize*dataSize;
      // Copy data:
      for (uint64_t j=0; j<vectorSize; ++j) {
	 const uint64_t bufferOffset = nbrOffset + j*dataSize;
	 for (uint64_t i=0; i<dataSize; ++i) {
	    components[j][compOffset+i] = buffer[bufferOffset + i];
	 }
      }
   }
   
   delete [] buffer; buffer = NULL;
   
   // SILO requires one variable name per (vector) component, but we only have one.
   // That is, for electric field SILO would like to get "Ex","Ey", and "Ez", but we
   // only have "E" in VLSV file. Use the VLSV variable name for all components:
   vector<string> varNames(vectorSize);
   vector<char*> varNamePtrs(vectorSize);
   for (uint64_t i=0; i<vectorSize; ++i) {
      stringstream ss;
      ss << varName << (i+1);
      varNames[i] = ss.str();
      varNamePtrs[i] = const_cast<char*>(varNames[i].c_str());
   }
   
   // Create a subdirectory to SILO file if array tag contains attribute "dir":
   map<string,string> arrayAttributes;
   const string rootDir = "/" + meshDir;
   fullName = rootDir + '/';
   if (vlsvReader.getArrayAttributes("VARIABLE",attributes,arrayAttributes) == true) {
      map<string,string>::const_iterator it = arrayAttributes.find("dir");
      if (it != arrayAttributes.end()) {
	 const string& dir = it->second;
	 
	 // If dir is not yet in SILO file, create it. The directories
	 // need to be created one level at a time:
	 if (directories.find(dir) == directories.end()) {
	    string tmp;
	    string s = dir;
	    string current = '/' + meshDir;
	    
	    size_t pos = s.find('/');
	    while (pos != string::npos) {
	       tmp = s.substr(0,pos);
	       if (pos != 0) {
		  current = current + '/' + tmp;
	       }
	       
	       // This exception will make this work correctly if the
	       // first character in directory name is '/':
   	       if (tmp.size() == 0) {
		  s = s.substr(pos+1);
		  pos = s.find('/');
		  continue;
	       }
	       
	       // If current directory does not exist, insert it:
	       if (directories.find(current) == directories.end()) {
		  directories.insert(current);
		  if (DBMkDir(fileptr,const_cast<char*>(current.c_str())) < 0) {
		     cerr << "ERROR failed to create dir '" << current << "'" << endl;
		  }
	       }
	       
	       s = s.substr(pos+1);
	       pos = s.find('/');
	    }
	    tmp = s.substr(0,pos);
	    if (tmp.size() > 0) current = current + '/' + tmp;
	    
	    // If directory name was given without a trailing '/', e.g. "/dir1/dir2",
	    // this if-statement inserts the final directory:
	    if (tmp.size() > 0) {
	       if (directories.find(current) == directories.end()) {
		  directories.insert(current);
		  if (DBMkDir(fileptr,const_cast<char*>(current.c_str())) < 0) {
		     cerr << "ERROR failed to create dir '" << current << "'" << endl;
		  }
	       }
	    }
	 }
	 
	 // cd into dir:
	 const string dirName = rootDir + dir;
	 if (DBSetDir(fileptr,const_cast<char*>(dirName.c_str())) < 0) {
	    cerr << "ERROR failed to cd into dir '" << dirName << "'" << endl;
	 }
	 fullName = dirName + '/';
      }
   }   
   fullName += varName;

   // Make an option list. If physical units of the quantity were given
   // in tag attributes, write units to SILO:
   unsigned int N_options = 0;
   if (arrayAttributes.find("unit") != arrayAttributes.end()) ++N_options;
   DBoptlist* optlist = getOptionList(vlsvReader,N_options);
            
   if (arrayAttributes.find("unit") != arrayAttributes.end()) {
      DBAddOption(optlist,DBOPT_UNITS,const_cast<char*>(arrayAttributes["unit"].c_str()));
   }
   
   // Write the unstructured mesh variable to SILO:
   const string fullMeshName = '/' + meshDir + '/' + meshDir;
   if (DBPutUcdvar(fileptr,varName.c_str(),fullMeshName.c_str(),vectorSize,&(varNamePtrs[0]),components,N_cells,NULL,0,
		   SiloType(dataType,dataSize),DB_ZONECENT,optlist) < 0) success = false;
   if (optlist != NULL) DBFreeOptlist(optlist);
   
   // Deallocate memory:
   for (uint64_t i=0; i<vectorSize; ++i) {delete [] components[i]; components[i] = NULL;}
   delete [] components; components = NULL;

   // cd back into root dir:
   if (DBSetDir(fileptr,const_cast<char*>(rootDir.c_str())) < 0) {
      cerr << "Failed to return to root dir in SILO!" << endl;
      success = false;
   }
   
   return success;
}

bool convertMeshVariable(vlsv::Reader& vlsvReader,const string& meshName,const string& varName) {   
   bool success = true;

   // Writing a unstructured grid variable is a rather straightforward process. The 
   // only compilation here is that some of the variables are actually vectors, i.e. 
   // vectorSize > 1 (vectorSize == 1 for scalars). Format in which vectors are stored in VLSV 
   // differ from format in which they are written to SILO files.
   vlsv::datatype::type dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("mesh",meshName));
   attributes.push_back(make_pair("name",varName));
   if (vlsvReader.getArrayInfo("VARIABLE",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      return false;
   }

   // Read variable data. We do not actually need to care if 
   // the data is given as floats, doubles, or long doubles. We 
   // just need to tell SILO what kind of data the buffer contains:
   char* buffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("VARIABLE",attributes,0,arraySize,buffer) == false) success = false;
   if (success == false) {
      delete [] buffer;
      return success;
   }
      
   // Vector variables need to be copied to temporary arrays before 
   // writing to SILO file:
   char** components = new char*[vectorSize];
   for (uint64_t i=0; i<vectorSize; ++i) {
      components[i] = new char[arraySize*dataSize];
      for (uint64_t j=0; j<arraySize; ++j) for (uint64_t k=0; k<dataSize; ++k) 
	components[i][j*dataSize+k] = buffer[j*vectorSize*dataSize + i*dataSize + k];
   }

   // SILO requires one variable name per (vector) component, but we only have one. 
   // That is, for electric field SILO would like to get "Ex","Ey", and "Ez", but we 
   // only have "E" in VLSV file. Use the VLSV variable name for all components.
   vector<string> varNames(vectorSize);
   vector<char*> varNamePtrs(vectorSize);
   for (uint64_t i=0; i<vectorSize; ++i) {
      stringstream ss;
      ss << varName << (i+1);
      varNames[i] = ss.str();
      varNamePtrs[i] = const_cast<char*>(varNames[i].c_str());
   }

   // Create a subdirectory to SILO file if array tag contains attribute "dir":
   map<string,string> arrayAttributes;
   const string rootDir = "/";
   if (vlsvReader.getArrayAttributes("VARIABLE",attributes,arrayAttributes) == true) {
      map<string,string>::const_iterator it = arrayAttributes.find("dir");
      if (it != arrayAttributes.end()) {
	 const string& dir = it->second;
	 
	 // If dir is not yet in SILO file, create it. The directories 
	 // need to be created one level at a time:
	 if (directories.find(dir) == directories.end()) {
	    string tmp;
	    string s = dir;
	    string current = "";
	    size_t pos = s.find('/');
	    while (pos != string::npos) {
	       tmp = s.substr(0,pos);
	       if (pos != 0) current = current + '/' + tmp;
	       
	       // This exception will make this work correctly if the 
	       // first character in directory name is '/':
	       if (tmp.size() == 0) {
		  s = s.substr(pos+1);
		  pos = s.find('/');
		  continue;
	       }

	       // If current directory does not exist, insert it:
	       if (directories.find(current) == directories.end()) {
		  directories.insert(current);
		  DBMkDir(fileptr,const_cast<char*>(current.c_str()));
	       }
	       
	       s = s.substr(pos+1);
	       pos = s.find('/');
	    }
	    tmp = s.substr(0,pos);
	    if (tmp.size() > 0) current = current + '/' + tmp;
	    
	    // If directory name was given without a trailing '/', e.g. "/dir1/dir2",
	    // this if inserts the final directory:
	    if (tmp.size() > 0) {
	       if (directories.find(current) == directories.end()) {
		  directories.insert(current);
		  DBMkDir(fileptr,const_cast<char*>(current.c_str()));
	       }
	    }
	 }
	 
	 // cd into dir:
	 DBSetDir(fileptr,const_cast<char*>(dir.c_str()));
      }
   }
   
   // Make an option list. If physical units of the quantity were given 
   // in tag attributes, write units to SILO:
   unsigned int N_options = 0;
   if (arrayAttributes.find("unit") != arrayAttributes.end()) ++N_options;
   DBoptlist* optlist = getOptionList(vlsvReader,N_options);
   
   if (arrayAttributes.find("unit") != arrayAttributes.end()) {
      DBAddOption(optlist,DBOPT_UNITS,const_cast<char*>(arrayAttributes["unit"].c_str()));
   }
   
   // Mesh name must have full path name so that VisIt will find them:
   const string mesh = "/" + meshName;

   // Write the unstructured mesh variable to SILO:
   if (DBPutUcdvar(fileptr,varName.c_str(),mesh.c_str(),vectorSize,&(varNamePtrs[0]),components,arraySize,NULL,0,
		   SiloType(dataType,dataSize),DB_ZONECENT,optlist) < 0) success = false;
   if (optlist != NULL) DBFreeOptlist(optlist);
   for (uint64_t i=0; i<vectorSize; ++i) {delete [] components[i]; components[i] = NULL;}
   delete [] components; components = NULL;
   delete [] buffer; buffer = NULL;

   // cd back into root dir:
   if (DBSetDir(fileptr,const_cast<char*>(rootDir.c_str())) < 0) {
      cerr << "Failed to return to root dir in SILO!" << endl;
      success = false;
   }
   return success;
}

bool convertPointMesh(vlsv::Reader& vlsvReader,const string& meshName) {
   bool success = true;
   //cerr << "Converting point mesh '" << meshName << "'" << endl;
   
   // Fetch mesh coordinate array info and do sanity check on values:
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("name",meshName));
   uint64_t arraySize,vectorSize,dataSize;
   vlsv::datatype::type dataType;
   if (vlsvReader.getArrayInfo("MESH",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Array MESH info could not be obtained!" << endl;
      return false;
   }
   if (dataType != vlsv::datatype::FLOAT) {
      cerr << "Mesh coordinates are not floating point values!" << endl;
      return false; // Coordinates must be floating point values
   }
   if (dataSize != sizeof(float) && (dataSize != sizeof(double) && dataSize != sizeof(long double))) {
      cerr << "Mesh coordinates have unsupported floating point byte size!" << endl;
      cerr << "\t byte size obtained: " << dataSize << endl;
      return false;
   }
   if (vectorSize < 1 || vectorSize > 3) {
      cerr << "Mesh dimensionality must be between 1 and 3!" << endl;
      return false;
   }
   
   if (arraySize == 0 || vectorSize == 0 || dataSize == 0) return true;

   // Read array XML attributes:
   map<string,string> attribsOut;
   if (vlsvReader.getArrayAttributes("MESH",attributes,attribsOut) == false) {
      cerr << "\t\tERROR: Failed to obtain XML tags of MESH array!" << endl;
   }
   // If XML tag has attribute 'meshinfo', read mesh information from that mesh instead of this point mesh:
   if (attribsOut.find("meshinfo") != attribsOut.end()) {
      list<pair<string,string> > tmp;
      const string meshInfoName = attribsOut["meshinfo"];
      tmp.push_back(make_pair("name",meshInfoName));
      if (vlsvReader.getArrayAttributes("MESH",tmp,attribsOut) == false) {
	 cerr << "\t\tERROR: Could not fetch mesh info from '" << meshInfoName << "' for point mesh '" << meshName << "'" << endl;
	 attribsOut.clear();
      }
   }
   
   const uint64_t N_points = arraySize;
   const uint64_t N_dims   = vectorSize;
   
   // Read all points from file:
   char* inbuffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("MESH",attributes,0,arraySize,inbuffer) == false) {
      cerr << "Failed to read point mesh '" << meshName << "' coordinates!" << endl;
      delete[] inbuffer; inbuffer = NULL;
      return false;
   }
   
   // Create x,y,z coordinate arrays for output (three pointers per possible floating point type).
   // Only the pointers with the same datatype as in VLSV file have arrays allocated to them.
   char* coordinateArrays[3];   
   float* xout4 = NULL; float* yout4 = NULL; float* zout4 = NULL;
   double* xout8 = NULL; double* yout8 = NULL; double* zout8 = NULL;
   long double* xout12 = NULL; long double* yout12 = NULL; long double* zout12 = NULL;
   switch (dataSize) {
    case (sizeof(float)):
      xout4 = new float[N_points];
      yout4 = new float[N_points];
      zout4 = new float[N_points];
      break;
    case (sizeof(double)):
      xout8 = new double[N_points];
      yout8 = new double[N_points];
      zout8 = new double[N_points];
      break;
    case (sizeof(long double)):
      xout12 = new long double[N_points];
      yout12 = new long double[N_points];
      zout12 = new long double[N_points];
      break;
   }

   // Apply scale factors to coordinate points:
   float xscale4 = 1.0; float yscale4 = 1.0; float zscale4 = 1.0;
   double xscale8 = 1.0; double yscale8 = 1.0; double zscale8 = 1.0;
   long double xscale12 = 1.0; long double yscale12 = 1.0; long double zscale12 = 1.0;
   
   switch (dataSize) {
    case (sizeof(float)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale4 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale4 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale4 = atof(attribsOut["zscaling"].c_str());
      
      for (uint64_t point=0; point<N_points; ++point) {
	 xout4[point] = reinterpret_cast<float*>(inbuffer)[point*N_dims+0] * xscale4;
	 yout4[point] = reinterpret_cast<float*>(inbuffer)[point*N_dims+1] * yscale4;
	 zout4[point] = reinterpret_cast<float*>(inbuffer)[point*N_dims+2] * zscale4;	 
      }
      
      coordinateArrays[0] = reinterpret_cast<char*>(xout4);
      coordinateArrays[1] = reinterpret_cast<char*>(yout4);
      coordinateArrays[2] = reinterpret_cast<char*>(zout4);      
      break;
    case (sizeof(double)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale8 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale8 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale8 = atof(attribsOut["zscaling"].c_str());

      for (uint64_t point=0; point<N_points; ++point) {
	 xout8[point] = reinterpret_cast<double*>(inbuffer)[point*N_dims+0] * xscale8;
	 yout8[point] = reinterpret_cast<double*>(inbuffer)[point*N_dims+1] * xscale8;
	 zout8[point] = reinterpret_cast<double*>(inbuffer)[point*N_dims+2] * xscale8;
      }
      
      coordinateArrays[0] = reinterpret_cast<char*>(xout8);
      coordinateArrays[1] = reinterpret_cast<char*>(yout8);
      coordinateArrays[2] = reinterpret_cast<char*>(zout8);
      break;
    case (sizeof(long double)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale12 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale12 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale12 = atof(attribsOut["zscaling"].c_str());
      
      for (uint64_t point=0; point<N_points; ++point) {
	 xout12[point] = reinterpret_cast<long double*>(inbuffer)[point*N_dims+0] * xscale12;
	 yout12[point] = reinterpret_cast<long double*>(inbuffer)[point*N_dims+1] * yscale12;
	 zout12[point] = reinterpret_cast<long double*>(inbuffer)[point*N_dims+2] * zscale12;
      }
      
      coordinateArrays[0] = reinterpret_cast<char*>(xout12);
      coordinateArrays[1] = reinterpret_cast<char*>(yout12);
      coordinateArrays[2] = reinterpret_cast<char*>(zout12);
      break;
   }
   delete [] inbuffer; inbuffer = NULL;

   // Make an option list and insert simulation time and timestep, if they are available:
   DBoptlist* optlist = getOptionList(vlsvReader,6);

   // Get coordinate names & units from VLSV file if available, otherwise use default values:
   parseCoordinateNames(optlist,attribsOut);

   if (DBPutPointmesh(fileptr,meshName.c_str(),N_dims,coordinateArrays,N_points,SiloType(dataType,dataSize),optlist) < 0) {
      cerr << "Failed to write the point mesh to file!" << endl;
      success = false;
   }
   if (optlist != NULL) DBFreeOptlist(optlist);

   // Deallocate memory and exit:
   delete [] xout4; delete [] yout4; delete [] zout4;
   delete [] xout8; delete [] yout8; delete [] zout8;
   delete [] xout12; delete [] yout12; delete [] zout12;
   return success;
}

/** Convert given multimesh piece from input VLSV file to SILO format.
 * @param vlsvReader VLSV reader used to read the input file.
 * @param meshName Name of the multimesh.
 * @param siloMeshName Name of given mesh piece in output SILO file.
 * @param start Offset to VLSV file array 'MESH' that contains (i,j,k) indices of cells for the given mesh piece.
 * @param amount Number of (i,j,k) tuples to read, i.e. total number of cells (local+ghosts) in the given mesh piece.
 * @param N_ghosts Number of ghost cells in given mesh piece.
 * @return If true, given mesh piece was successfully written to output file.*/
bool convertQuadMesh2(vlsv::Reader& vlsvReader,const string& meshName,const string& siloMeshName,
		      const uint64_t& start,const uint64_t& amount,const uint64_t& N_ghosts) {
   bool success = true;
   
   vlsv::datatype::type dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("mesh",meshName));
   
   // Get info on array containing mesh bounding box parameters:
   if (vlsvReader.getArrayInfo("MESH_BBOX",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Could not obtain info on array 'MESH_BBOX'!" << endl;
      return false;
   }
   // Read mesh bounding box data:
   char* meshBboxIn = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("MESH_BBOX",attributes,0,arraySize,meshBboxIn) == false) {
      cerr << "Failed to read array 'MESH_BBOX' !" << endl;
      delete [] meshBboxIn; meshBboxIn = NULL;
      return false;
   }
   if (arraySize != 6) {
      cerr << "Array 'MESH_BBOX' has incorrect size!" << endl;
      delete [] meshBboxIn; meshBboxIn = NULL;
      return false;
   }
   if (dataType != vlsv::datatype::FLOAT) {
      cerr << "Incorrect datatype in array 'MESH_BBOX' !" << endl;
      delete [] meshBboxIn; meshBboxIn = NULL;
      return false;
   }
   
   // Convert read bounding box parameters to floats:
   float* meshbbox4 = NULL;
   double* meshbbox8 = NULL;
   long double* meshbbox12 = NULL;
   switch (dataSize) {
    case (sizeof(float)):
      meshbbox4 = reinterpret_cast<float*>(meshBboxIn);
      for (int i=0; i<6; ++i) bbox[i] = meshbbox4[i];
      break;
    case (sizeof(double)):
      meshbbox8 = reinterpret_cast<double*>(meshBboxIn);
      for (int i=0; i<6; ++i) bbox[i] = meshbbox8[i];
      break;
    case (sizeof(long double)):
      meshbbox12 = reinterpret_cast<long double*>(meshBboxIn);
      for (int i=0; i<6; ++i) bbox[i] = meshbbox12[i];
      break;
    default:
      cerr << "Unknown floating point datatype in array 'MESH_BBOX' !" << endl;
      delete [] meshBboxIn; meshBboxIn = NULL;
      return false;
      break;
   }
   delete [] meshBboxIn; meshBboxIn = NULL;
   
   // Get info on array containing cell (i,j,k) indices:
   attributes.clear();
   attributes.push_back(make_pair("name",meshName));
   attributes.push_back(make_pair("type",vlsv::mesh::STRING_QUAD_MULTI));
   if (vlsvReader.getArrayInfo("MESH",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Array MESH info could not be obtained!" << endl;
      return false;
   }
   if (dataType != vlsv::datatype::UINT && dataType != vlsv::datatype::INT) {
      cerr << "Array 'MESH' does not contain integer indices!" << endl;
      return false;
   }
   if (dataSize != sizeof(uint16_t) && dataSize != sizeof(uint32_t) && dataSize != sizeof(uint64_t)) {
      cerr << "Array 'MESH' has unsupported (un)signed integer datatype!" << endl;
      cerr << "Byte size is " << dataSize << endl;
      cerr << "Supported sizes are: " << sizeof(uint16_t) << '\t' << sizeof(uint32_t) << '\t' << sizeof(uint64_t) << endl;
      return false;
   }
   
   // Get array 'MESH' XML attributes:
   map<string,string> attribsOut;
   if (vlsvReader.getArrayAttributes("MESH",attributes,attribsOut) == false) {
      cerr << "\t\tERROR: Failed to obtain XML tags of MESH array!" << endl;
   }
   
   // Create pointers to input data that are used to convert 
   // supported unsigned integer datatypes to internally used datatype:
   char* indicesIn = new char[amount*vectorSize*dataSize];
   int16_t* ptrMeshi2 = reinterpret_cast<int16_t*>(indicesIn);
   int32_t* ptrMeshi4 = reinterpret_cast<int32_t*>(indicesIn);
   int64_t* ptrMeshi8 = reinterpret_cast<int64_t*>(indicesIn);
   uint16_t* ptrMeshui2 = reinterpret_cast<uint16_t*>(indicesIn);
   uint32_t* ptrMeshui4 = reinterpret_cast<uint32_t*>(indicesIn);
   uint64_t* ptrMeshui8 = reinterpret_cast<uint64_t*>(indicesIn);
   
   // Read (i,j,k) indices:
   if (vlsvReader.readArray("MESH",attributes,start,amount,indicesIn) == false) {
      cerr << "Failed to read array 'MESH' !" << endl;
      delete [] indicesIn; indicesIn = NULL;
      return false;
   }

   // Create arrays for storing x,y,z coordinates of cells:
   size_t outputArraySize = static_cast<size_t>(ceil(amount * 1.2));
   vector<float> xcrds;
   vector<float> ycrds;
   vector<float> zcrds;
   vector<int> zoneList;
   xcrds.reserve(outputArraySize);
   ycrds.reserve(outputArraySize);
   zcrds.reserve(outputArraySize);
   zoneList.reserve(arraySize*8);
   
   switch (dataType) {
    // Eliminate duplicate nodes when (i,j,k) indices were given as signed integers:
    case (vlsv::datatype::INT):
      switch (dataSize) {
       case (sizeof(int16_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshi2,zoneList,xcrds,ycrds,zcrds);
	 break;
       case (sizeof(int32_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshi4,zoneList,xcrds,ycrds,zcrds);
	 break;
       case (sizeof(int64_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshi8,zoneList,xcrds,ycrds,zcrds);
	 break;
       default:
	 cerr << "ERROR: Unsupported signed integer datatype in convertQuadMesh2!" << endl; exit(1);
	 break;
      }
      break;
    // Eliminate duplicate nodes when (i,j,k) indices were given as unsigned integers:
    case (vlsv::datatype::UINT):
      switch (dataSize) {
       case (sizeof(uint16_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshui2,zoneList,xcrds,ycrds,zcrds);
	 break;
       case (sizeof(uint32_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshui4,zoneList,xcrds,ycrds,zcrds);
	 break;
       case (sizeof(uint64_t)):
	 eliminateDuplicateNodes(amount,vectorSize,ptrMeshui8,zoneList,xcrds,ycrds,zcrds);
	 break;
       default:
	 cerr << "ERROR: Unsupported unsigned integer datatype in convertQuadMesh2!" << endl; exit(1);
	 break;
      }
      break;
    // Indices were given with an unsupported datatype:
    default:
      cerr << "ERROR: Unsupported datatype in convertQuadMesh2!" << endl; exit(1);
      break;
   }
   delete [] indicesIn; indicesIn = NULL;
   
   // Write unstructured mesh to SILO file:
   const int N_dims  = 3;                      // Number of dimensions
   const int N_zones = amount;
   const int N_nodes = xcrds.size();           // Total number of nodes (>N_zones)
   int shapeTypes[] = {DB_ZONETYPE_HEX};       // Hexahedrons only
   int shapeSizes[] = {8};                     // Each hexahedron has 8 nodes
   int shapeCnt[] = {N_zones};                 // Only 1 shape type (hexahedron)
   const int N_shapes = 1;                     //  -- "" --

   void* coords[3];                            // Pointers to coordinate arrays
   coords[0] = &(xcrds[0]);
   coords[1] = &(ycrds[0]);
   coords[2] = &(zcrds[0]);

   // Write zone list into silo file:
   const string zoneListName = siloMeshName + "Zones";
   if (DBPutZonelist2(fileptr,zoneListName.c_str(),N_zones,N_dims,&(zoneList[0]),8*N_zones,0,0,N_ghosts,shapeTypes,shapeSizes,shapeCnt,N_shapes,NULL) < 0) success = false;
   
   // Make an option list (do not insert simulation time or timestep):
   DBoptlist* optlist = getOptionList(vlsvReader,7,false);
   
   // Get coordinate names & units from VLSV file if available, otherwise use default values:
   parseCoordinateNames(optlist,attribsOut);
   int connectivity = 1;
   DBAddOption(optlist,DBOPT_TIME,&connectivity);
   
   // Write UCD grid to SILO file:
   if (DBPutUcdmesh(fileptr,siloMeshName.c_str(),N_dims,NULL,coords,N_nodes,N_zones,zoneListName.c_str(),NULL,SiloType(vlsv::datatype::FLOAT,sizeof(float)),optlist) < 0) success = false;
   if (optlist != NULL) DBFreeOptlist(optlist);
   return success;
}

/** Convert multimesh from VLSV file to SILO format.
 * @param vlsvReader VLSV reader that is used to read the input file.
 * @param meshName Name of the multimesh.
 * @param fileName Name of input file.
 * @param outputFileName Name of output file.
 * @return If true, multimesh and associated variables were written successfully to output file.*/
bool convertMultimesh(vlsv::Reader& vlsvReader,const string& meshName,const string& fileName,const string& outputFileName) {
   bool success = true;
   
   // Read array 'MESH' attributes first to check if the mesh is stored in another VLSV file:
   vlsv::datatype::type meshDataType;
   uint64_t meshArraySize,meshVectorSize,meshDataSize;
   list<pair<string,string> > attributes;
   map<string,string> meshAttributes;
   attributes.push_back(make_pair("name",meshName));
   if (vlsvReader.getArrayAttributes("MESH",attributes,meshAttributes) == false) {
      cerr << "ERROR: Array MESH attributes could not be obtained!" << endl;
      return false;
   }
   
   // Read mesh using a separate vlsv::Reader, this makes it easier to handle cases 
   // where the mesh only exists in some VLSV files. If array MESH_ZONES contains an attribute 
   // 'file', the mesh is read from that file. Otherwise mesh is read from variable fileName:
   vlsv::Reader meshReader;
   bool writeMeshOut = true;
   if (meshAttributes.find("file") != meshAttributes.end()) {
      // Mesh already exists in an older SILO file, do not write it again:
      if (meshReader.open(meshAttributes["file"]) == false) {
	 cerr << "ERROR: Could not open file '" << meshAttributes["file"] << "' !" << endl;
	 return false;
      }
      writeMeshOut = false;
      // Replace meshAttributes["file"] suffix ".vlsv" with ".silo":
      string fileout = meshAttributes["file"];
      size_t pos = fileout.rfind(".vlsv");
      if (pos != string::npos) fileout.replace(pos,5,".silo");
      meshAttributes["file"] = fileout;
   } else {
      // Mesh needs to be re-written:
      if (meshReader.open(fileName) == false) {
	 cerr << "ERROR: Could not open file '" << fileName << "' !" << endl;
	 return false;
      }
      meshAttributes["file"] = outputFileName;
   }
   
   // Read the number of meshes, usually equal to the number of 
   // MPI processes in the simulation:
   attributes.clear();
   attributes.push_back(make_pair("mesh",meshName));
   if (meshReader.getArrayInfo("MESH_ZONES",attributes,meshArraySize,meshVectorSize,meshDataType,meshDataSize) == false) {
      cerr << "ERROR: Array MESH_ZONES info could not be obtained!" << endl;
      return false;
   }
   
   char* N_meshes = new char[meshArraySize*meshVectorSize*meshDataSize];
   if (meshReader.readArray("MESH_ZONES",attributes,0,meshArraySize,N_meshes) == false) {
      cerr << "ERROR: Could not read array 'MESH_ZONES' !" << endl;
      delete [] N_meshes; N_meshes = NULL;
      return false;
   }
   
   // Allocate array that contains the names of mesh pieces. This 
   // needs to be char** array so that it can be passed to SILO writer:
   char** meshNames = new char* [meshArraySize];
   for (uint64_t m=0; m<meshArraySize; ++m) meshNames[m] = NULL;

   // Array that contains, for each variable, the full paths to 
   // multimesh variable arrays, which will be written to same 
   // directories as the corresponding mesh pieces:
   map<string,vector<string> > multivars;

   // Get total number of cells, and number of ghosts cells, in each mesh piece
   // and calculate offsets into variable arrays based on those numbers:
   uint64_t* N_cells  = new uint64_t[meshArraySize];
   uint64_t* N_ghosts = new uint64_t[meshArraySize];
   unsigned int* variableOffsets = new unsigned int[meshArraySize+1];
   unsigned int* ghostOffsets = new unsigned int[meshArraySize+1];
   variableOffsets[0] = 0;
   ghostOffsets[0] = 0;
   for (uint64_t m=0; m<meshArraySize; ++m) {
      const uint64_t byteSize = meshVectorSize*meshDataSize;
      switch (meshDataType) {
       case (vlsv::datatype::INT):
	 N_cells[m] = convInt(N_meshes+m*byteSize,meshDataType,meshDataSize);
	 N_ghosts[m] = convInt(N_meshes+m*byteSize+meshDataSize,meshDataType,meshDataSize);
	 break;
       case (vlsv::datatype::UINT):
	 N_cells[m] = convUInt(N_meshes+m*byteSize,meshDataType,meshDataSize);
	 N_ghosts[m] = convUInt(N_meshes+m*byteSize+meshDataSize,meshDataType,meshDataSize);
	 break;
       default:
	 cerr << "Unsupported datatype in array 'MESH_ZONES'" << endl;
	 delete [] N_meshes; N_meshes = NULL;
	 delete [] N_cells; N_cells = NULL;
	 delete [] N_ghosts; N_ghosts = NULL;
	 delete [] variableOffsets; variableOffsets = NULL;
	 delete [] ghostOffsets; ghostOffsets = NULL;
	 return false;
	 break;
      }
      variableOffsets[m+1] = variableOffsets[m] + (N_cells[m]-N_ghosts[m]);
      ghostOffsets[m+1] = ghostOffsets[m] + N_ghosts[m];
   }

   // Convert mesh piece m, and its variables, to SILO format:
   uint64_t cellOffset = 0;
   
   // Number of next mesh piece written to SILO file. This is different from 
   // loop counter m below if one or more mesh pieces have zero cells.
   uint64_t currentMeshNumber = 0;
   for (uint64_t m=0; m<meshArraySize; ++m) {
      // Skip empty multimesh pieces:
      if (N_cells[m] == 0) continue;
      
      // Read ghost cells local ids and domain array info for this mesh piece:
      vlsv::datatype::type domainDatatype,localidDatatype;
      uint64_t domainArraySize,domainVectorSize,domainDataSize;
      uint64_t localidArraySize,localidVectorSize,localidDataSize;
      if (meshReader.getArrayInfo("MESH_GHOST_DOMAINS",attributes,domainArraySize,domainVectorSize,domainDatatype,domainDataSize) == false) success = false;
      if (meshReader.getArrayInfo("MESH_GHOST_LOCALIDS",attributes,localidArraySize,localidVectorSize,localidDatatype,localidDataSize) == false) success = false;
      if (success == false) cerr << "ERROR" << endl;
      if (domainVectorSize != 1 || (domainDatatype != vlsv::datatype::INT && domainDatatype != vlsv::datatype::UINT)) success = false;
      if (localidVectorSize != 1 || (localidDatatype != vlsv::datatype::INT && localidDatatype != vlsv::datatype::UINT)) success = false;
      if (success == false) {
	 cerr << "ERROR occurred while obtaining arrays 'MESH_GHOST_DOMAINS' and/or 'MESH_GHOST_LOCALIDS' !" << endl;
	 exit(1);
      }
      
      // Read ghost localid & domain arrays, and convert them to uints:
      char* bfr = new char[N_ghosts[m]*localidDataSize];
      unsigned int* ghostLocalIDs = NULL;
      if (meshReader.readArray("MESH_GHOST_LOCALIDS",attributes,ghostOffsets[m],N_ghosts[m],bfr) == false) success = false;
      if (success == true) {
	 if (localidDatatype != vlsv::datatype::UINT || localidDataSize != sizeof(unsigned int)) {
	    // File has wrong datatype, convert to uints:
	    ghostLocalIDs = new unsigned int[N_ghosts[m]];
	    for (uint64_t i=0; i<N_ghosts[m]; ++i) {
	       ghostLocalIDs[i] = convUInt(bfr+i*localidDataSize,localidDatatype,localidDataSize);
	    }
	 } else {
	    // File has correct datatype, just swap pointers:
	    ghostLocalIDs = reinterpret_cast<unsigned int*>(bfr);
	    bfr = NULL;
	 }
      }
      delete [] bfr; bfr = NULL;

      bfr = new char[N_ghosts[m]*domainDataSize];
      unsigned int* ghostDomains = NULL;
      if (meshReader.readArray("MESH_GHOST_DOMAINS",attributes,ghostOffsets[m],N_ghosts[m],bfr) == false) success = false;
      if (success == true) {
	 if (domainDatatype != vlsv::datatype::UINT || domainDataSize != sizeof(unsigned int)) {
	    // File has wrong datatype, convert to uints:
	    ghostDomains = new unsigned int[N_ghosts[m]];
	    for (uint64_t i=0; i<N_ghosts[m]; ++i) {
	       ghostDomains[i] = convUInt(bfr+i*domainDataSize,domainDatatype,domainDataSize);
	    }
	 } else {
	    // File has correct datatype, just swap pointers:
	    ghostDomains = reinterpret_cast<unsigned int*>(bfr);
	    bfr = NULL;
	 }
      }
      delete [] bfr; bfr = NULL;

      // Check that everything is ok:
      if (success == false) {
	 cerr << "ERROR has occurred while reading arrays 'MESH_GHOST_DOMAINS' and/or 'MESH_GHOST_LOCALIDS' !" << endl;
	 delete [] ghostLocalIDs;
	 delete [] ghostDomains;
	 delete [] meshNames;
	 delete [] N_meshes;
	 exit(1);
      }
      
      // Create a new directory for mesh m:
      stringstream ss;
      ss.fill('0');
      ss << "mesh" << setw(8) << currentMeshNumber;
      string meshDir;
      ss >> meshDir;
      
      if (DBMkDir(fileptr,meshDir.c_str()) < 0) {
	 cerr << "ERROR occurred while creating a directory to SILO file!" << endl;
	 success = false; break;
      }

      // Create full SILO pathname (file+dir) for mesh piece m and store it to meshNames:
      char* mn = new char[fileName.size() + 1 + 2*meshDir.size() + 5];
      mn[0] = '\0';
      strcat(mn,meshAttributes["file"].c_str());
      strcat(mn,":");
      strcat(mn,"/");
      strcat(mn,meshDir.c_str());
      strcat(mn,"/");
      strcat(mn,meshDir.c_str());
      meshNames[currentMeshNumber] = mn;
      
      // Change current SILO directory to meshDir:
      if (DBSetDir(fileptr,meshDir.c_str()) < 0) {
	 cerr << "ERROR occurred while trying to chdir to '" << meshDir << "'" << endl;
	 success = false; break;
      }
      
      // Write mesh piece m:
      if (writeMeshOut == true) {
	 if (convertQuadMesh2(meshReader,meshName,meshDir,cellOffset,N_cells[m],N_ghosts[m]) == false) {
	    cerr << "ERROR occurred while converting multimesh #" << m << endl;
	    success = false; break;
	 }
      }

      // Write all variables in mesh piece:
      set<string> variableNames;
      if (vlsvReader.getUniqueAttributeValues("VARIABLE","name",variableNames) == false) {
	 success = false; break;
      }  
      string fullVariableName;
      for (set<string>::const_iterator it=variableNames.begin(); it!=variableNames.end(); ++it) {
	 if (convertMultimeshVariable(vlsvReader,meshName,*it,meshDir,m,variableOffsets,
				      N_cells[m],N_ghosts[m],ghostLocalIDs,ghostDomains,fullVariableName) == true) {
	    multivars[*it].push_back(fullVariableName);
	 }
      }
      if (success == false) break;

      // Change current SILO directory back to root:
      if (DBSetDir(fileptr,"..") < 0) {
	 cerr << "ERROR occurred while trying to chdir to root" << endl;
	 success = false; break;
      }
      cellOffset += N_cells[m];
      
      // Clear memory:
      directories.clear();      
      delete [] ghostLocalIDs; ghostLocalIDs = NULL;
      delete [] ghostDomains; ghostDomains = NULL;
      
      // Increase mesh number:
      ++currentMeshNumber;
   }   
   const uint64_t N_validMeshPieces = currentMeshNumber;
   
   // Write multimesh to SILO file:
   int* meshTypes = NULL;
   if (success == true) {
      meshTypes = new int[N_validMeshPieces];
      for (uint64_t m=0; m<N_validMeshPieces; ++m) meshTypes[m] = DB_UCDMESH;
      DBoptlist* optlist = getOptionList(vlsvReader,1,true);
      int connectivity = 1;
      DBAddOption(optlist,DBOPT_TV_CONNECTIVITY,&connectivity);
      if (DBPutMultimesh(fileptr,meshName.c_str(),N_validMeshPieces,meshNames,meshTypes,optlist) < 0) {
	 cerr << "ERROR occurred while writing multimesh!" << endl;
	 success = false;
      }
      delete [] meshTypes; meshTypes = NULL;
      if (optlist != NULL) DBFreeOptlist(optlist);
   }

   // Create requested SILO directories for multimesh variables. Directory names 
   // are given in XML attributes called 'dir'. Directory names default to root '/':
   directories.clear();
   set<string> variableNames;
   map<string,string> variableDirectories;
   if (vlsvReader.getUniqueAttributeValues("VARIABLE","name",variableNames) == false) success = false;
   if (success == true) for (set<string>::const_iterator it=variableNames.begin(); it!=variableNames.end(); ++it) {
      createVariableDirectory(vlsvReader,meshName,*it,variableDirectories[*it]);
   }
   
   // Write variables as SILO multivars:
   for (map<string,vector<string> >::iterator it=multivars.begin(); it!=multivars.end(); ++it) {
      // Allocate array for multivar paths, and copy path names to it:
      char** varNames = new char* [it->second.size()];
      int* varTypes = new int[it->second.size()];
      for (size_t i=0; i<it->second.size(); ++i) {
	 varNames[i] = new char[(it->second)[i].size()+1];
	 varNames[i][0] = '\0';
	 strcat(varNames[i],(it->second)[i].c_str());
	 varTypes[i] = DB_UCDVAR;
      }
      
      // cd to variable directory:
      map<string,string>::const_iterator dir = variableDirectories.find(it->first);
      if (dir == variableDirectories.end()) {
	 cerr << "ERROR could not find directory for variable '" << it->first << "'" << endl;
      } else {
	 if (DBSetDir(fileptr,dir->second.c_str()) < 0) {
	    cerr << "ERROR occurred while attempting to cd into SILO dir '" << dir->second << "'" << endl;
	 }
      }
      
      // Write multivar to SILO:
      DBoptlist* optList = DBMakeOptlist(1);
      DBAddOption(optList,DBOPT_MMESH_NAME,const_cast<char*>(meshName.c_str()));
      
      if (DBPutMultivar(fileptr,it->first.c_str(),it->second.size(),varNames,varTypes,optList) < 0) {
	 cerr << "ERROR failed to write multimesh variable '" << it->first << "'" << endl;
	 success = false;
      }
      
      // cd back to root:
      if (DBSetDir(fileptr,"/") < 0) {
	 cerr << "ERROR occurred while attempting to cd back into SILO root dir after writing multivariable" << endl;
      }
      
      // Deallocate memory:
      if (optList != NULL) DBFreeOptlist(optList);      
      for (size_t i=0; i<it->second.size(); ++i) {delete varNames[i]; varNames[i] = NULL;}
      delete [] varNames; varNames = NULL;
      delete [] varTypes; varTypes = NULL;
   }
   
   // Deallocate memory and exit:
   meshReader.close();
   for (uint64_t m=0; m<meshArraySize; ++m) {delete meshNames[m]; meshNames[m] = NULL;}
   delete [] meshNames; meshNames = NULL;
   delete [] N_meshes; N_meshes = NULL;
   delete [] variableOffsets; variableOffsets = NULL;
   delete [] ghostOffsets; ghostOffsets = NULL;
   delete [] N_cells; N_cells = NULL;
   delete [] N_ghosts; N_ghosts = NULL;
   return success;
}

bool convertQuadMesh(vlsv::Reader& vlsvReader,const string& meshName) {
   bool success = true;
   
   // First task is to push all unique node coordinates into a map.
   // This is not too difficult for unrefined grids, since each spatial cell stores
   // its bottom lower left corner coordinate and size. For refined grid the situation
   // is more complex, as there are more unique nodes than the lower left corners:
   map<NodeCrd<float>,uint64_t,NodeComp> nodes4;
   map<NodeCrd<double>,uint64_t,NodeComp> nodes8;
   map<NodeCrd<long double>,uint64_t,NodeComp> nodes12;
   
   vlsv::datatype::type dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("name",meshName));
   attributes.push_back(make_pair("type",vlsv::mesh::STRING_QUAD));
   if (vlsvReader.getArrayInfo("MESH",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Array MESH info could not be obtained!" << endl;
      return false;
   }   
   if (dataType != vlsv::datatype::FLOAT) {
      cerr << "Mesh coordinates are not floating point values!" << endl;
      return false; // Coordinates must be floating point values
   }
   if (dataSize != sizeof(float) && (dataSize != sizeof(double) && dataSize != sizeof(long double))) {
      cerr << "Mesh coordinates have unsupported floating point byte size!" << endl;
      cerr << "\t byte size obtained: " << dataSize << endl;
      return false;
   }
   
   #ifdef PROFILE
      profile::start(convertQuadMeshTotalID);
   #endif
   
   // Get array XML attributes:
   map<string,string> attribsOut;
   if (vlsvReader.getArrayAttributes("MESH",attributes,attribsOut) == false) {
      cerr << "\t\tERROR: Failed to obtain XML tags of MESH array!" << endl;
   }
   // If XML tag has attribute 'meshinfo', read mesh information from that mesh instead of this quad mesh:
   if (attribsOut.find("meshinfo") != attribsOut.end()) {
      list<pair<string,string> > tmp;
      const string meshInfoName = attribsOut["meshinfo"];
      tmp.push_back(make_pair("name",meshInfoName));
      if (vlsvReader.getArrayAttributes("MESH",tmp,attribsOut) == false) {
	 cerr << "\t\tERROR: Could not fetch mesh info from '" << meshInfoName << "' for point quad '" << meshName << "'" << endl;
	 attribsOut.clear();
      }
   }
   
   // Read the coordinate array one node (of a spatial cell) at a time
   // and create a map which only contains each existing node once.
   // Pointer ptr points correctly to all floating point types:
   char* ptr = new char[vectorSize*dataSize];
   const int ds = dataSize;
   
   // Create a pointer to value zero that is of the same datatype 
   // as the floating point value in file:
   float       ZERO4 = 0.0;
   double      ZERO8 = 0.0;
   long double ZERO12 = 0.0;
   char* zeroPtr = NULL;
   switch (dataSize) {
    case (sizeof(float)):
      zeroPtr = reinterpret_cast<char*>(&ZERO4);
      break;
    case (sizeof(double)):
      zeroPtr = reinterpret_cast<char*>(&ZERO8);
      break;
    case (sizeof(long double)):
      zeroPtr = reinterpret_cast<char*>(&ZERO12);
      break;
   }
   
   // Insert all eight nodes of each cell into map nodes.
   // NOTE: map is a unique associative container - given a suitable comparator, map 
   // will filter out duplicate nodes:
   #ifdef PROFILE
      profile::start("data read+insert");
   #endif
   switch (dataSize) {
    case (sizeof(float)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds),0));
	 nodes4.insert(make_pair(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds),0));
      }
      break;
    case (sizeof(double)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds),0));
	 nodes8.insert(make_pair(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds),0));
      }
      break;
    case (sizeof(long double)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds),0));
	 nodes12.insert(make_pair(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds),0));
      }
      break;
   }
   #ifdef PROFILE
      profile::stop("data read+insert");
   #endif
   
   // Copy unique node x,y,z coordinates into separate arrays, 
   // which will be passed to SILO writer:
   int N_nodes = 0;
   if (dataSize == sizeof(float)) N_nodes = nodes4.size();
   if (dataSize == sizeof(double)) N_nodes = nodes8.size();
   if (dataSize == sizeof(long double)) N_nodes = nodes12.size();

   // Scale factors for coordinates:
   float xscale4 = 1.0; float yscale4 = 1.0; float zscale4 = 1.0;
   double xscale8 = 1.0; double yscale8 = 1.0; double zscale8 = 1.0;
   long double xscale12 = 1.0; long double yscale12 = 1.0; long double zscale12 = 1.0;
   
   uint64_t counter = 0;
   char* xcrds = new char[N_nodes*dataSize];
   char* ycrds = new char[N_nodes*dataSize];
   char* zcrds = new char[N_nodes*dataSize];
   switch (dataSize) {
    case (sizeof(float)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale4 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale4 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale4 = atof(attribsOut["zscaling"].c_str());
      for (map<NodeCrd<float>,uint64_t>::iterator it=nodes4.begin(); it!=nodes4.end(); ++it) {
	 it->second = counter;
	 it->first.xcopy(xcrds + counter*dataSize,xscale4);
	 it->first.ycopy(ycrds + counter*dataSize,yscale4);
	 it->first.zcopy(zcrds + counter*dataSize,zscale4);
	 ++counter;
      }
      break;
    case (sizeof(double)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale8 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale8 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale8 = atof(attribsOut["zscaling"].c_str());
      for (map<NodeCrd<double>,uint64_t>::iterator it=nodes8.begin(); it!=nodes8.end(); ++it) {
	 it->second = counter;
	 it->first.xcopy(xcrds + counter*dataSize,xscale8);
	 it->first.ycopy(ycrds + counter*dataSize,yscale8);
	 it->first.zcopy(zcrds + counter*dataSize,zscale8);
	 ++counter;
      }
      break;
    case (sizeof(long double)):
      if (attribsOut.find("xscaling") != attribsOut.end()) xscale12 = atof(attribsOut["xscaling"].c_str());
      if (attribsOut.find("yscaling") != attribsOut.end()) yscale12 = atof(attribsOut["yscaling"].c_str());
      if (attribsOut.find("zscaling") != attribsOut.end()) zscale12 = atof(attribsOut["zscaling"].c_str());
      for (map<NodeCrd<long double>,uint64_t>::iterator it=nodes12.begin(); it!=nodes12.end(); ++it) {
	 it->second = counter;
	 it->first.xcopy(xcrds + counter*dataSize,xscale12);
	 it->first.ycopy(ycrds + counter*dataSize,yscale12);
	 it->first.zcopy(zcrds + counter*dataSize,zscale12);
	 ++counter;
      }
      break;
   }

   // Read through the coordinate array again and create a node list. Each 3D spatial cell is 
   // associated with 8 nodes, and most of these nodes are shared with neighbouring cells. In 
   // order to get VisIt display the data correctly, the duplicate nodes should not be used. 
   // Here we create a list of indices into xcrds,ycrds,zcrds arrays, with eight entries per cell:
   int* nodeList = new int[8*arraySize];
   map<NodeCrd<float>,uint64_t,NodeComp>::const_iterator it4;
   map<NodeCrd<double>,uint64_t,NodeComp>::const_iterator it8;
   map<NodeCrd<long double>,uint64_t,NodeComp>::const_iterator it12;
   
   #ifdef PROFILE
      profile::start("data read+find");
   #endif
   switch (dataSize) {
    case (sizeof(float)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr)); nodeList[i*8+0] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr)); nodeList[i*8+1] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr)); nodeList[i*8+2] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr)); nodeList[i*8+3] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds)); nodeList[i*8+4] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds)); nodeList[i*8+5] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds)); nodeList[i*8+6] = it4->second;
	 it4 = nodes4.find(NodeCrd<float>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds)); nodeList[i*8+7] = it4->second;
      }
      break;
    case (sizeof(double)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr)); nodeList[i*8+0] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr)); nodeList[i*8+1] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr)); nodeList[i*8+2] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr)); nodeList[i*8+3] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds)); nodeList[i*8+4] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds)); nodeList[i*8+5] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds)); nodeList[i*8+6] = it8->second;
	 it8 = nodes8.find(NodeCrd<double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds)); nodeList[i*8+7] = it8->second;
      }
      break;
    case (sizeof(long double)):
      for (uint64_t i=0; i<arraySize; ++i) {
	 if (vlsvReader.readArray("MESH",attributes,i,1,ptr) == false) {success = false;}
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr, zeroPtr)); nodeList[i*8+0] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr, zeroPtr)); nodeList[i*8+1] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds, zeroPtr)); nodeList[i*8+2] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds, zeroPtr)); nodeList[i*8+3] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr, zeroPtr,ptr+5*ds)); nodeList[i*8+4] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds, zeroPtr,ptr+5*ds)); nodeList[i*8+5] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds,ptr+3*ds,ptr+4*ds,ptr+5*ds)); nodeList[i*8+6] = it12->second;
	 it12 = nodes12.find(NodeCrd<long double>(ptr+0*ds,ptr+1*ds,ptr+2*ds, zeroPtr,ptr+4*ds,ptr+5*ds)); nodeList[i*8+7] = it12->second;
      }
      break;
   }
   #ifdef PROFILE
      profile::stop("data read+find");
      profile::start("data write");
   #endif
   
   // Write the unstructured mesh to SILO file:
   const int N_dims  = 3;                      // Number of dimensions
   const int N_zones = arraySize;              // Total number of zones (=spatial cells)
   int shapeTypes[] = {DB_ZONETYPE_HEX};       // Hexahedrons only
   int shapeSizes[] = {8};                     // Each hexahedron has 8 nodes
   int shapeCnt[] = {N_zones};                 // Only 1 shape type (hexahedron)
   const int N_shapes = 1;                     //  -- "" --
   
   void* coords[3];                            // Pointers to coordinate arrays
   coords[0] = xcrds;
   coords[1] = ycrds;
   coords[2] = zcrds;
   
   // Write zone list into silo file:
   const string zoneListName = meshName + "Zones";
   if (DBPutZonelist2(fileptr,zoneListName.c_str(),N_zones,N_dims,nodeList,8*arraySize,0,0,0,shapeTypes,shapeSizes,shapeCnt,N_shapes,NULL) < 0) success = false;

   // Make an option list and insert simulation time and timestep, if they are available:
   DBoptlist* optlist = getOptionList(vlsvReader,6);
   
   // Get coordinate names & units from VLSV file if available, otherwise use default values:
   parseCoordinateNames(optlist,attribsOut);
   
   // Write grid into silo file:
   if (DBPutUcdmesh(fileptr,meshName.c_str(),N_dims,NULL,coords,N_nodes,N_zones,zoneListName.c_str(),NULL,SiloType(dataType,dataSize),optlist) < 0) success = false;
   if (optlist != NULL) DBFreeOptlist(optlist);
   
   // Deallocate memory:
   nodes4.clear();
   nodes8.clear();
   nodes12.clear();
   delete [] ptr; ptr = NULL;
   delete [] nodeList; nodeList = NULL;
   delete [] xcrds; xcrds = NULL;
   delete [] ycrds; ycrds = NULL;
   delete [] zcrds; zcrds = NULL;
   #ifdef PROFILE
      profile::stop("data write");
      profile::start("variables total");
   #endif
   
   // Write all variables of this mesh into silo file:
   set<string> variableNames;
   if (vlsvReader.getUniqueAttributeValues("VARIABLE","name",variableNames) == false) {
      #ifdef PROFILE
         profile::stop("variables total");
         profile::stop(convertQuadMeshTotalID);
      #endif
      return false;
   }
   for (set<string>::const_iterator it=variableNames.begin(); it!=variableNames.end(); ++it) {
      if (convertMeshVariable(vlsvReader,meshName,*it) == false) success = false;
   }

   #ifdef PROFILE
      profile::stop("variables total");
      profile::stop(convertQuadMeshTotalID);
   #endif
   return success;
}

/** Append curve (x,y) data pair to map curves. VLSV file array that contains 
 * data for y points must be of type vlsv::datatype::FLOAT. Furthermore, array and vector sizes 
 * must equal unity, i.e. only single x-value and a single y-value. Finally, (x,y) values 
 * are written as doubles to the SILO file irrespective of their actual floating point types.
 * @param vlsvReader VLSV file reader that is used to read data.
 * @param varName Name of the array containing the curve data in VLSV file.
 * @param isTimeSeries If true, curve is a time series, i.e. time is used as x-value.
 * @param If true, curve data was successfully extracted from the file.*/
bool appendCurveValue(vlsv::Reader& vlsvReader,const string& varName,bool isTimeSeries) {
   bool success = true;

   vlsv::datatype::type xDataType;
   uint64_t xArraySize,xVectorSize,xDataSize;
   list<pair<string,string> > attribs;
   map<string,string> attribsOut;
   
   // Buffer for reading data:
   char buffer[sizeof(long double)];
   
   // Get x value from VLSV file:
   attribs.push_back(make_pair("name","time"));
   if (vlsvReader.getArrayInfo("PARAMETER",attribs,xArraySize,xVectorSize,xDataType,xDataSize) == false) return false;
   if (xDataType != vlsv::datatype::FLOAT) return false;
   if (xArraySize != 1 || xVectorSize != 1) return false;
   if (vlsvReader.readArray("PARAMETER",attribs,0,xArraySize,buffer) == false) return false;

   // Convert x value to double:
   double x = NAN;
   switch (xDataSize) {
    case (sizeof(float)):
      x = *reinterpret_cast<float*>(buffer);
      break;
    case (sizeof(double)):
      x = *reinterpret_cast<double*>(buffer);
      break;
    case (sizeof(long double)):
      x = *reinterpret_cast<long double*>(buffer);
      break;
   }
   
   // Get info of array containing the y value:
   attribs.clear();
   attribs.push_back(make_pair("name",varName));
   if (vlsvReader.getArrayInfo("TIMESERIES",attribs,xArraySize,xVectorSize,xDataType,xDataSize) == false) return false;
   if (xDataType != vlsv::datatype::FLOAT) return false;
   if (xArraySize != 1 || xVectorSize != 1) return false;
   
   // Get y-value array attributes:
   if (vlsvReader.getArrayAttributes("TIMESERIES",attribs,attribsOut) == false) return false;
   
   // Get y value from VLSV file:
   if (vlsvReader.readArray("TIMESERIES",attribs,0,xArraySize,buffer) == false) return false;

   // Convert y value to double:
   double y = NAN;
   switch (xDataSize) {
    case (sizeof(float)):
      y = *reinterpret_cast<float*>(buffer);
      break;
    case (sizeof(double)):
      y = *reinterpret_cast<double*>(buffer);
      break;
    case (sizeof(long double)):
      y = *reinterpret_cast<long double*>(buffer);
      break;
   }

   // Store (x,y) value pair to map curves. Note that the VLSV files are not read in 
   // chronological order. Thus, the (x,y) pairs are also read in unspecified order 
   // and need to be sorted by ascending x-value. The sort is done here:
   //curves[varName][x] = y;
   if (attribsOut.find("xlabel") != attribsOut.end()) curves[varName].xlabel = attribsOut["xlabel"];
   if (attribsOut.find("ylabel") != attribsOut.end()) curves[varName].ylabel = attribsOut["ylabel"];
   if (attribsOut.find("xunit") != attribsOut.end()) curves[varName].xunit = attribsOut["xunit"];
   if (attribsOut.find("yunit") != attribsOut.end()) curves[varName].yunit = attribsOut["yunit"];
   curves[varName].data[x] = y;
   
   return success;
}

bool convertSILO(const string& fname) {
   bool success = true;
   
   // Open VLSV file for reading:
   vlsv::Reader vlsvReader;
   if (vlsvReader.open(fname) == false) {
      cerr << "Failed to open '" << fname << "'" << endl;
      return false;
   }
   
   // Open SILO file for writing:
   string fileout = fname;
   size_t pos = fileout.rfind(".vlsv");
   if (pos != string::npos) fileout.replace(pos,5,".silo");
   
   fileptr = DBCreate(fileout.c_str(),DB_CLOBBER,DB_LOCAL,"VLSV data file",DB_PDB);
   if (fileptr == NULL) return false;
   
   // Get the names of all meshes in vlsv file, and write into silo file:
   set<string> meshNames;
   if (vlsvReader.getUniqueAttributeValues("MESH","name",meshNames) == false) {
      DBClose(fileptr);
      return false;
   }
   
   // Convert each mesh into SILO format:
   for (set<string>::const_iterator it=meshNames.begin(); it!=meshNames.end(); ++it) {
      map<string,string> attribsOut;
      list<pair<string,string> > attribsIn;
      attribsIn.push_back(make_pair("name",*it));
      vlsvReader.getArrayAttributes("MESH",attribsIn,attribsOut);
      if (attribsOut.find("type") == attribsOut.end()) {
	 cout << "\t Skipping mesh '" << *it << "' because it has unspecified type" << endl;
	 continue;
      }
      
      if (attribsOut["type"] == vlsv::mesh::STRING_QUAD) convertQuadMesh(vlsvReader,*it);
      else if (attribsOut["type"] == vlsv::mesh::STRING_POINT) convertPointMesh(vlsvReader,*it);
      else if (attribsOut["type"] == vlsv::mesh::STRING_QUAD_MULTI) convertMultimesh(vlsvReader,*it,fname,fileout);
      else {
	 cout << "\t Skipping mesh '" << *it << "' because it has unknown type" << endl;
      }
   }
   vlsvReader.close();
   DBClose(fileptr); fileptr = NULL;
   return success;
}

/** Convert curve data from VLSV files into SILO format. Curve data is here understood to 
 * consist of (x,y) value pairs, one pair per VLSV file.
 * @param fname Name of a VLSV file to read.
 * @param lastFile If true, all VLSV files have been read and curve data should be written to SILO file.
 * @return If true, curve data was converted successfully.*/
bool convertCurveSILO(const string& fname,bool lastFile) {
   bool success = true;
   
   // Open VLSV file for reading:
   vlsv::Reader vlsvReader;
   if (lastFile == false) {
      if (vlsvReader.open(fname) == false) {
	 cerr << "Failed to open '" << fname << "'" << endl;
	 return false;
      }

      // Get all unique attribute 'name' values from arrays that have name TIMESERIES:
      set<string> parameterNames;
      if (vlsvReader.getUniqueAttributeValues("TIMESERIES","name",parameterNames) == false) {
	 vlsvReader.close();
	 return false;
      }

      // Append (x,y) pair from each TIMESERIES array to curves:
      for (set<string>::const_iterator it=parameterNames.begin(); it!=parameterNames.end(); ++it) {
	 appendCurveValue(vlsvReader,*it,true);
      }
   
      vlsvReader.close();
   }

   // Write curve data to SILO file:
   if (lastFile == true) {
      cerr << "Writing file curves.silo" << endl;
      string outName = "curves.silo";
      fileptr = DBCreate(outName.c_str(),DB_CLOBBER,DB_LOCAL,"VLSV curve data",DB_PDB);
      if (fileptr == NULL) return false;

      for (map<string,CurveData>::const_iterator it=curves.begin(); it!=curves.end(); ++it) {
	 vector<double> x(it->second.data.size());
	 vector<double> y(it->second.data.size());
	 size_t index = 0;
	 for (map<double,double>::const_iterator jt=it->second.data.begin(); jt!=it->second.data.end(); ++jt) {
	    x[index] = jt->first;
	    y[index] = jt->second;
	    ++index;
	 }
	 
	 DBoptlist* optlist = DBMakeOptlist(4);
	 DBAddOption(optlist,DBOPT_XLABEL,const_cast<char*>(it->second.xlabel.c_str()));
	 DBAddOption(optlist,DBOPT_YLABEL,const_cast<char*>(it->second.ylabel.c_str()));
	 DBAddOption(optlist,DBOPT_XUNITS,const_cast<char*>(it->second.xunit.c_str()));
	 DBAddOption(optlist,DBOPT_YUNITS,const_cast<char*>(it->second.yunit.c_str()));
	 
	 if (DBPutCurve(fileptr,it->first.c_str(),&(x[0]),&(y[0]),DB_DOUBLE,x.size(),optlist) < 0) {
	    cerr << "Failed to write curve '" << it->first << "'" << endl; success = false;
	 }	 
	 DBFreeOptlist(optlist);
      }
      DBClose(fileptr); fileptr = NULL;
   }
   return success;
}

int main(int argn,char* args[]) {
   
   if (argn < 2) {
      cout << endl;
      cout << "USAGE: ./vlsv2silo <input file mask(s)>" << endl;
      cout << "Each VLSV in the current directory is compared against the given file mask(s)," << endl;
      cout << "and if match is found, that file is converted into SILO format." << endl;
      cout << endl;
      return 1;
   }
   
   #ifdef PROFILE
      MPI_Init(&argn,&args);
      convertQuadMeshTotalID = profile::initializeTimer("Quad mesh total");
   #endif
   
   // Convert file masks into strings:
   vector<string> masks;
   for (int i=1; i<argn; ++i) masks.push_back(args[i]);

   // Compare directory contents against each mask:
   set<string> inputFiles;
   const string directory = ".";
   const string suffix = ".vlsv";
   for (size_t mask=0; mask<masks.size(); ++mask) {
      cout << "Comparing mask '" << masks[mask] << "'" << endl;
      unsigned int filesConverted = 0;
      DIR* dir = opendir(directory.c_str());
      if (dir == NULL) continue;
      
      struct dirent* entry = readdir(dir);
      while (entry != NULL) {
	 const string entryName = entry->d_name;
	 // Compare entry name against given mask and file suffix ".vlsv":
	 if (entryName.find(masks[mask]) == string::npos || entryName.find(suffix) == string::npos) {
	    entry = readdir(dir);
	    continue;
	 }
	 inputFiles.insert(entryName);
	 entry = readdir(dir);
      }
      closedir(dir);

      // Convert each file in set inputFiles. The files should be in 
      // chronological order:
      for (set<string>::iterator it=inputFiles.begin(); it!=inputFiles.end(); ++it) {
	 cout << "\t converting '" << *it << "'" << endl;
	 directories.clear();
	 convertSILO(*it);
	 ++filesConverted;
      }
      
      // Convert curve(s) stored to VLSV files, these include time series:
      for (set<string>::iterator it=inputFiles.begin(); it!=inputFiles.end(); ++it) {
	 convertCurveSILO(*it,false);
      }
      convertCurveSILO(string(""),true);
      
      if (filesConverted == 0) cout << "\t no matches found" << endl;
   }
   
   #ifdef PROFILE
      profile::print(MPI_COMM_SELF);
      MPI_Finalize();
   #endif
   return 0;
}
