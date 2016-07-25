/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2016 Finnish Meteorological Institute
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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include <mpi.h>

#include "vlsv_writer.h"

using namespace std;

/** Struct for storing data read from a matrix market input file.*/
struct Data {
   uint64_t N_rows;                /**< Number of rows in matrix.*/
   uint64_t N_cols;                /**< Number of columns in matrix.*/
   uint64_t N_nz;                  /**< Number of non-zeroes in (sparse) matrix.*/
   uint64_t N_values;              /**< Number of variables.*/

   vector<uint64_t> rows;          /**< Row (y) node coordinates.*/
   vector<uint64_t> cols;          /**< Column (x) node coordinates.*/
   vector<uint64_t> globalIDs;     /**< Global IDs of matrix elements.*/
   vector<vector<float> > values;  /**< Variable data.*/

   Data(): N_rows(0),N_cols(0),N_nz(0),N_values(0) { }

};

/** Split a string into one or mode substrings based on the given delimiter.
 * @param s String to be split.
 * @param delim Delimiter, e.g., ' ' (whitespace).
 * @return Vector of found substrings.*/
vector<string> splitString(const string& s,const char& delim) {
   vector<string> elems;
   istringstream ss(s);
   string item;
   while (getline(ss,item,delim)) {
      elems.push_back(item);
   }
   return elems;
}

/** Read the given file and store sparse matrix data to matrixData.
 * @param fname Name of the input file that contains matrix data in Matrix Market format.
 * @param matrixData Struct where read data is stored.
 * @return If true, input file was read successfully and values in matrixData are valid.*/
bool readMatrixMarket(const string& fname,Data& matrixData) {
   bool success = true;
   
   // Open input file
   ifstream in(fname.c_str(), ios_base::in);
   if (in.good() == false) {
      cerr << "Error opening input file" << endl;
      return false;
   }

   // First line should contain '%%MatrixMarket' 'matrix' <format> <field> <symmetry>
   string line;
   getline(in,line);
     {
        istringstream ss(line);
        string title,object,format,field,symmetry;
        if (!(ss >> title >> object >> format >> field >> symmetry)) {
           cerr << "Input file not recognized as a Matrix Market file" << endl;
           return false;
        }
        
        if (title != "%%MatrixMarket") {
           cerr << "Input file not recognized as a Matrix Market file" << endl;
           return false;
        }
        
        if (object != "matrix") {
           cerr << "Input file does not contain a matrix" << endl;
           return false;
        }
        
        if (object == "array") {
           cerr << "'array' format for matrices is not currently supported" << endl;
           return false;
        }
        
        if (field == "complex") {
           cerr << "Complex data is not currently supported" << endl;
           return false;
        }
        
        if (symmetry != "general") {
           cerr << "Matrix has unsupported symmetry '" << symmetry << "'" << endl;
           return false;
        }
     }
   line.clear();

   // Skip comment section. Upon exit line should contain matrix size.
   while (in.good() == true) {
      getline(in,line);
      if (line[0] != '%') break;
      line.clear();
   }
   
   // Read matrix rows, columns, and amount of non-zeroes
   if (in.good() == false) {
      cerr << "Unexpected end of file" << endl;
      return false;
   }
     {
        istringstream ss(line);
        if (!(ss >> matrixData.N_rows >> matrixData.N_cols >> matrixData.N_nz)) {
           cerr << "Error reading matrix size" << endl;
           return false;
        }
     }

   // Create x,y coordinate arrays:
   matrixData.rows.resize(matrixData.N_rows+1);
   for (uint64_t r=0; r<matrixData.N_rows; ++r) {
      matrixData.rows[r] = matrixData.N_rows - r;
   }
   matrixData.cols.resize(matrixData.N_cols+1);
   for (uint64_t c=0; c<matrixData.N_cols; ++c) {
      matrixData.cols[c] = c;
   }
   matrixData.globalIDs.resize(matrixData.N_nz);

   // Read and store (row,column,value) tuples. Assume that each line has 
   // the same number of data values (used to allocate array values).
   line.clear();
   getline(in,line);
   if (in.good() == false) {
      cerr << "Unexpected end of file" << endl;
      return false;
   }

   uint64_t counter = 0;
   while (in.good() == true) {
      const vector<string> elems = splitString(line,' ');
      
      // Allocate memory for all data values
      if (matrixData.values.size() == 0 && elems.size() > 2) {         
         matrixData.N_values = elems.size() - 2;
         matrixData.values.resize(matrixData.N_values);
         for (uint64_t i=0; i<matrixData.N_values; ++i) {
            matrixData.values[i].resize(matrixData.N_nz);            
         }
      }

      // Lambda for calculating cell's global index
      auto calcIndex = [&](const uint64_t& row,const uint64_t& col) {
         return row*matrixData.N_cols + col;
      };

      istringstream ss(line);
      uint64_t row,col;
      if (!(ss >> row >> col)) {
         cerr << "Error occurred while reading matrix rows and columns" << endl;
         return false;
      }
      matrixData.globalIDs[counter] = calcIndex(row-1,col-1);

      for (uint64_t i=0; i<matrixData.N_values; ++i) {
         float value;
         if (!(ss >> value)) {
            cerr << "Error occurred while reading matrix data" << endl;
            return false;
         }
         matrixData.values[i][counter] = value;
      }
      
      line.clear();
      getline(in,line);
      ++counter;
   }

   cout << "Successfully read " << counter << " / " << matrixData.N_nz << " lines" << endl;
   cout << "\t Found " << matrixData.values.size() << " variables" << endl;
   return success;
}

/** Write sparse matrix data from matrixData to output VLSV file. The output file 
 * name will be the same as the input file name, except that the suffix is replaced 
 * by ".vlsv".
 * @param fname Input file name.
 * @param matrixData Struct containing the matrix data.
 * @return If true, a VLSV file containing matrix data was successfully written.*/
bool writeVlsv(const string& fname,Data& matrixData) {
   bool success = true;

   const string meshName = "matrix";
   
   string fname_out = fname;
   
   // Strip the possible path from input file name
   if (fname.find_last_of('/') <= fname.size())
     fname_out = fname.substr(fname.find_last_of('/')+1,string::npos);

   // Replace the externsion with '.vlsv'
   fname_out = fname_out.substr(0,fname_out.find_last_of('.')) + ".vlsv";
   
   vlsv::Writer vlsv;
   if (vlsv.open(fname_out,MPI_COMM_SELF,0) == false) {
      cerr << "Error opening output file" << endl;
      return false;
   }
   
   // Write Global IDs
   map<string,string> xmlAttribs;
   xmlAttribs["name"] = meshName;
   xmlAttribs["type"] = vlsv::mesh::STRING_UCD_MULTI;
   xmlAttribs["spatial_dimension"] = "2";

   if (vlsv.writeArray("MESH",xmlAttribs,matrixData.N_nz,1,matrixData.globalIDs.data()) == false) success = false;

   // Write mesh bounding box
   uint64_t bbox[6];
   bbox[0] = matrixData.N_cols;
   bbox[1] = matrixData.N_rows;
   bbox[2] = 1;
   bbox[3] = 1;
   bbox[4] = 1;
   bbox[5] = 1;
   xmlAttribs.clear();
   xmlAttribs["mesh"] = meshName;

   if (vlsv.writeArray("MESH_BBOX",xmlAttribs,6,1,bbox) == false) success = false;
   
   cerr << "rows,cols " << matrixData.N_cols << ' ' << matrixData.N_rows << endl;
   
   // Write node coordinates
   if (vlsv.writeArray("MESH_NODE_CRDS_X",xmlAttribs,matrixData.cols.size(),1,matrixData.cols.data()) == false) success = false;
   if (vlsv.writeArray("MESH_NODE_CRDS_Y",xmlAttribs,matrixData.rows.size(),1,matrixData.rows.data()) == false) success = false;
   
   uint64_t z_crds[2];
   z_crds[0] = 0; z_crds[1] = 1;
   if (vlsv.writeArray("MESH_NODE_CRDS_Z",xmlAttribs,2,1,z_crds) == false) success = false;
   
   // Write number of cells in each multimesh domain
   uint64_t N_zones[2];
   N_zones[0] = matrixData.N_nz;
   N_zones[1] = 0;   
   if (vlsv.writeArray("MESH_DOMAIN_SIZES",xmlAttribs,1,2,N_zones) == false) success = false;
   if (vlsv.writeArray("MESH_GHOST_DOMAINS",xmlAttribs,0,2,N_zones) == false) success = false;
   if (vlsv.writeArray("MESH_GHOST_LOCALIDS",xmlAttribs,0,2,N_zones) == false) success = false;
   
   // Write variables
   for (size_t i=0; i<matrixData.values.size(); ++i) {
      stringstream varName;
      varName << "Var" << i;
      xmlAttribs["name"] = varName.str();
      if (vlsv.writeArray("VARIABLE",xmlAttribs,matrixData.values[i].size(),1,matrixData.values[i].data()) == false) success = false;
   }
   
   vlsv.close();
   return success;
}


int main(int argn,char* args[]) {
   
   if (argn != 2) {
      cerr << endl;
      cerr << "USAGE:" << endl;
      cerr << "./conv_mtx_vlsv <mtx input file>" << endl;
      cerr << endl;
      return 1;
   }
   
   Data matrixData;
   if (readMatrixMarket(args[1],matrixData) == false) {
      cerr << "Failed to read input file" << endl;
      return 1;
   }
   
   MPI_Init(&argn,&args);
   if (writeVlsv(args[1],matrixData) == false) {
      cerr << "Failed to write VLSV file" << endl;
      return 1;
   }
   MPI_Finalize();
   
   return 0;
}


