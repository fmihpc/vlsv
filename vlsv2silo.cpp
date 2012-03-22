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

#ifdef PROFILE
   #include <profile.hpp>
   int convertQuadMeshTotalID;
   //int quadMeshID;
   //int quadVariablesID;
   //int quadMeshReadID;
#endif

#include "vlsvreader.h"

using namespace std;

// When making option lists SILO seems to only store a pointer to 
// variable giving the option's value. The pointers must be valid 
// when data is written to file. Thus, we need to define some global 
// variables here that store values put in option lists.

static int timestep = 0;
static float time_float = 0.0;
static double time_double = 0.0;

static string label_x = "x-coordinate";
static string label_y = "y-coordinate";
static string label_z = "z-coordinate";
static string units_x = "m";
static string units_y = "m";
static string units_z = "m";

static DBfile* fileptr = NULL; // Pointer to file opened by SILO

static set<string> directories; /**< List of directories created to output SILO file. If one tries to 
				 * create a directory more than once, SILO will throw errors to the 
				 * console. Similarly attempting to cd into a non-existing directory 
				 * throws errors to console. In order to avoid these errors, currently 
				 * existing directories are stored in this set.*/

static map<string,map<double,double> > curves; /**< Container for all curve data stored to VLSV file(s).
						* map<double,double> is for sorting data by ascending x values.*/

int64_t convInt(const char* ptr,const VLSV::datatype& dataType,const uint64_t& dataSize) {
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
      exit(1);
      break;
   }
}

uint64_t convUInt(const char* ptr,const VLSV::datatype& dataType,const uint64_t& dataSize) {
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
      exit(1);
      break;
   }
}

int SiloType(const VLSV::datatype& dataType,const uint64_t& dataSize) {
   switch (dataType) {
    case VLSV::INT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case VLSV::UINT:
      if (dataSize == 2) return DB_SHORT;
      else if (dataSize == 4) return DB_INT;
      else if (dataSize == 8) return DB_LONG;
      else return -1;
      break;
    case VLSV::FLOAT:
      if (dataSize == 4) return DB_FLOAT;
      else if (dataSize == 8) return DB_DOUBLE;
      else return -1;
      break;
   }
   return -1;
}

const float EPS = 1.0e-6;

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
      for (int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
   }

   char ycopy(char* target,REAL scaling=1.0) const {
      const REAL tmp = y*scaling;
      const char* ptr = reinterpret_cast<const char*>(&tmp);
      for (int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
   }

   char zcopy(char* target,REAL scaling=1.0) const {
      const REAL tmp = z*scaling;
      const char* ptr = reinterpret_cast<const char*>(&tmp);
      for (int i=0; i<sizeof(REAL); ++i) target[i] = ptr[i];
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

/** SILO format does not seem to like have both DBOPT_CYCLE (=timestep) and 
 * DBOPT_TIME/DBOPT_DTIME (time in float/double) defined in the option list, 
 * so we just check if time has been given in the file and add its value to 
 * the option list.
 */
DBoptlist* getOptionList(VLSVReader& vlsvReader,const int& N_extra=0) {
   DBoptlist* optlist = NULL;
   list<pair<string,string> > attribs;
   uint64_t arraySize,vectorSize,dataSize;
   VLSV::datatype dataType;

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
      
      if (dataType == VLSV::INT) timestep = convInt(valbuffer,dataType,dataSize);
      else if (dataType == VLSV::UINT) timestep = convUInt(valbuffer,dataType,dataSize);
      else {
	 cerr << "timestep given in unsupported datatype!" << endl;
      }
      
      DBAddOption(optlist,DBOPT_CYCLE,&timestep);
      delete [] valbuffer; valbuffer = NULL;
   }   
   return optlist;
}

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

bool convertMeshVariable(VLSVReader& vlsvReader,const string& meshName,const string& varName) {   
//   #ifdef PROFILE
//      profile::start(quadVariablesID);
//   #endif
   bool success = true;

   // Writing a unstructured grid variable is a rather straightforward process. The 
   // only compilation here is that some of the variables are actually vectors, i.e. 
   // vectorSize > 1 (vectorSize == 1 for scalars). Format in which vectors are stored in VLSV 
   // differ from format in which they are written to SILO files.
   VLSV::datatype dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("mesh",meshName));
   attributes.push_back(make_pair("name",varName));
   if (vlsvReader.getArrayInfo("VARIABLE",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
//      #ifdef PROFILE
//         profile::stop(quadVariablesID);
//      #endif
      return false;
   }

   // Read variable data. We do not actually need to care if 
   // the data is given as floats, doubles, or long doubles. We 
   // just need to tell SILO what kind of data the buffer contains:
   char* buffer = new char[arraySize*vectorSize*dataSize];
   if (vlsvReader.readArray("VARIABLE",attributes,0,arraySize,buffer) == false) success = false;
   if (success == false) {
      delete [] buffer;
//      #ifdef PROFILE
//         profile::stop(quadVariablesID);
//      #endif
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
//   #ifdef PROFILE
//      profile::stop(quadVariablesID);
//   #endif
   return success;
}

bool convertPointMesh(VLSVReader& vlsvReader,const string& meshName) {
   bool success = true;
   //cerr << "Converting point mesh '" << meshName << "'" << endl;
   
   // Fetch mesh coordinate array info and do sanity check on values:
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("name",meshName));
   uint64_t arraySize,vectorSize,dataSize;
   VLSV::datatype dataType;
   if (vlsvReader.getArrayInfo("MESH",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Array MESH info could not be obtained!" << endl;
      return false;
   }
   if (dataType != VLSV::FLOAT) {
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

bool convertQuadMesh(VLSVReader& vlsvReader,const string& meshName) {
   bool success = true;
   const float EPS = 1.0e-7;
   
   // First task is to push all unique node coordinates into a map.
   // This is not too difficult for unrefined grids, since each spatial cell stores
   // its bottom lower left corner coordinate and size. For refined grid the situation
   // is more complex, as there are more unique nodes than the lower left corners:
   map<NodeCrd<float>,uint64_t,NodeComp> nodes4;
   map<NodeCrd<double>,uint64_t,NodeComp> nodes8;
   map<NodeCrd<long double>,uint64_t,NodeComp> nodes12;
   
   VLSV::datatype dataType;
   uint64_t arraySize,vectorSize,dataSize;
   list<pair<string,string> > attributes;
   attributes.push_back(make_pair("name",meshName));
   attributes.push_back(make_pair("type",VLSV::MESH_QUAD));
   if (vlsvReader.getArrayInfo("MESH",attributes,arraySize,vectorSize,dataType,dataSize) == false) {
      cerr << "Array MESH info could not be obtained!" << endl;
      return false;
   }   
   if (dataType != VLSV::FLOAT) {
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
   char* zeroPtr;
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
 * data for y points must be of type VLSV::FLOAT. Furthermore, array and vector sizes 
 * must equal unity, i.e. only single x-value and a single y-value. Finally, (x,y) values 
 * are written as doubles to the SILO file irrespective of their actual floating point types.
 * @param vlsvReader VLSV file reader that is used to read data.
 * @param varName Name of the array containing the curve data in VLSV file.
 * @param isTimeSeries If true, curve is a time series, i.e. time is used as x-value.
 * @param If true, curve data was successfully extracted from the file.*/
bool appendCurveValue(VLSVReader& vlsvReader,const string& varName,bool isTimeSeries) {
   bool success = true;

   VLSV::datatype xDataType;
   uint64_t xArraySize,xVectorSize,xDataSize;
   list<pair<string,string> > attribs;
   
   // Buffer for reading data:
   char buffer[sizeof(long double)];
   
   // Get x value from VLSV file:
   attribs.push_back(make_pair("name","time"));
   if (vlsvReader.getArrayInfo("PARAMETER",attribs,xArraySize,xVectorSize,xDataType,xDataSize) == false) return false;
   if (xDataType != VLSV::FLOAT) return false;
   if (xArraySize != 1 || xVectorSize != 1) return false;
   if (vlsvReader.readArray("PARAMETER",attribs,0,xArraySize,buffer) == false) return false;

   // Convert x value to double:
   double x;
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
   
   // Get y value from VLSV file:
   attribs.clear();
   attribs.push_back(make_pair("name",varName));
   if (vlsvReader.getArrayInfo("TIMESERIES",attribs,xArraySize,xVectorSize,xDataType,xDataSize) == false) return false;
   if (xDataType != VLSV::FLOAT) return false;
   if (xArraySize != 1 || xVectorSize != 1) return false;
   if (vlsvReader.readArray("TIMESERIES",attribs,0,xArraySize,buffer) == false) return false;

   // Convert y value to double:
   double y;
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
   curves[varName][x] = y;
   return success;
}

bool convertSILO(const string& fname) {
   bool success = true;
   
   // Open VLSV file for reading:
   VLSVReader vlsvReader;
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
      
      if (attribsOut["type"] == VLSV::MESH_QUAD) convertQuadMesh(vlsvReader,*it);
      else if (attribsOut["type"] == VLSV::MESH_POINT) convertPointMesh(vlsvReader,*it);
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
   VLSVReader vlsvReader;
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

      for (map<string,map<double,double> >::const_iterator it=curves.begin(); it!=curves.end(); ++it) {
	 vector<double> x(it->second.size());
	 vector<double> y(it->second.size());
	 size_t index = 0;
	 for (map<double,double>::const_iterator jt=it->second.begin(); jt!=it->second.end(); ++jt) {
	    x[index] = jt->first;
	    y[index] = jt->second;
	    ++index;
	 }
	 
	 if (DBPutCurve(fileptr,it->first.c_str(),&(x[0]),&(y[0]),DB_DOUBLE,x.size(),NULL) < 0) {
	    cerr << "Failed to write curve '" << it->first << "'" << endl; success = false;
	 }
      }
      DBClose(fileptr); fileptr = NULL;
   }
   return success;
}

int main(int argn,char* args[]) {
   
   if (argn < 2) {
      cout << endl;
      cout << "USAGE: ./vlsv2vtk <input file mask(s)>" << endl;
      cout << "Each VLSV in the current directory is compared against the given file mask(s)," << endl;
      cout << "and if match is found, that file is converted into SILO format." << endl;
      cout << endl;
      return 1;
   }
   
   #ifdef PROFILE
      MPI_Init(&argn,&args);
      convertQuadMeshTotalID = profile::initializeTimer("Quad mesh total");
      //quadMeshID             = profile::initializeTimer("mesh convert");
      //quadVariablesID        = profile::initializeTimer("variables");
      //quadMeshReadID         = profile::initializeTimer("data read");
   #endif
   
   // Convert file masks into strings:
   vector<string> masks;
   for (int i=1; i<argn; ++i) masks.push_back(args[i]);

   // Compare directory contents against each mask:
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
	 cout << "\t converting '" << entryName << "'" << endl;
	 directories.clear();
	 convertSILO(entryName);
	 ++filesConverted;
	 entry = readdir(dir);
      }
      closedir(dir);

      // Convert curve(s) stored to VLSV files, these include time series:
      dir = opendir(directory.c_str());
      if (dir == NULL) continue;
      
      entry = readdir(dir);
      while (entry != NULL) {
	 const string entryName = entry->d_name;
	 // Compare entry name against given mask and file suffix ".vlsv":
	 if (entryName.find(masks[mask]) == string::npos || entryName.find(suffix) == string::npos) {
	    entry = readdir(dir);
	    continue;
	 }
	 directories.clear();
	 convertCurveSILO(entryName,false);
	 entry = readdir(dir);	 
      }
      closedir(dir);
      convertCurveSILO(string(""),true);
      
      if (filesConverted == 0) cout << "\t no matches found" << endl;
   }
   
   #ifdef PROFILE
      profile::print(MPI_COMM_SELF);
      MPI_Finalize();
   #endif
   return 0;
}
