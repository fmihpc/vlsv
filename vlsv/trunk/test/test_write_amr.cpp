#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <map>
#include <time.h>
#include <cmath>

#include "../vlsv_writer.h"
#include "../vlsv_amr.h"

using namespace std;

uint64_t calculateGlobalID(uint64_t i,uint64_t j,uint64_t k,uint64_t Nx,uint64_t Ny,uint64_t Nz) {
   return k*Ny*Nx + j*Nx + i;
}

void refineBlock(const size_t& index,vector<uint64_t>& globalIDs) {
   uint32_t i,j,k,refLevel;
   vlsv::calculateCellIndices(globalIDs[index],refLevel,i,j,k);
   
   i *= 2;
   j *= 2;
   k *= 2;
   const uint64_t newGlobalID = vlsv::calculateGlobalID(refLevel+1,i,j,k);

   globalIDs[index] = newGlobalID;
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i+1,j  ,k  ));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i  ,j+1,k  ));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i+1,j+1,k  ));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i  ,j  ,k+1));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i+1,j  ,k+1));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i  ,j+1,k+1));
   globalIDs.push_back(vlsv::calculateGlobalID(refLevel+1,i+1,j+1,k+1));
}

bool writeMesh(vlsv::Writer& vlsv,const std::string& geometryString,const int& dims) {
   bool success = true;
   srand(time(NULL));
   
   string meshName;
   if (dims == 2 && geometryString == vlsv::geometry::STRING_CARTESIAN) meshName = "cart_amr_2d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_CARTESIAN) meshName = "cart_amr_3d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_CYLINDRICAL) meshName = "cyl_amr_3d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_SPHERICAL) meshName = "sphe_amr_3d";
   else {
      cerr << "unknown geometry and/or mesh dimensionality" << endl; exit(1);
   }
   
   map<string,string> attributes;
   attributes["mesh"] = meshName;
   
   // Number of unrefined blocks in each coordinate direction:
   uint32_t Nx0,Ny0,Nz0;
   if (geometryString == vlsv::geometry::STRING_CARTESIAN) {
      Nx0 = 5;
      Ny0 = 5;
      Nz0 = 1;
   } else if (geometryString == vlsv::geometry::STRING_CYLINDRICAL) {
      Nx0 = 4;
      Ny0 = 8;
      Nz0 = 1;
   } else if (geometryString == vlsv::geometry::STRING_SPHERICAL) {
      Nx0 = 4;
      Ny0 = 2;
      Nz0 = 8;
   } else {
      cerr << "unknown geometry" << endl; exit(1);
   }
   
   uint32_t Nx_cells = 1;
   uint32_t Ny_cells = 1;
   uint32_t Nz_cells = 1;
   uint32_t N_blocks = Nx0*Ny0*Nz0;
   uint32_t N_ghosts = 0;
   uint32_t maxRefinementLevel = 8;
   
   vlsv::initMesh(Nx0,Ny0,Nz0,maxRefinementLevel);
   
   // Write mesh bounding box:
   uint32_t bbox[6];
   bbox[0] = Nx0;
   bbox[1] = Ny0;
   bbox[2] = Nz0;
   bbox[3] = Nx_cells;
   bbox[4] = Ny_cells;
   bbox[5] = Nz_cells;
   if (vlsv.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;

   // Write unrefined node coordinates:
   float* x_crds = new float[Nx0*Nx_cells+1];
   float* y_crds = new float[Ny0*Ny_cells+1];
   float* z_crds = new float[Nz0*Nz_cells+1];
   
   if (geometryString == vlsv::geometry::STRING_CARTESIAN) {
      for (int i=0; i<Nx0*Nx_cells+1; ++i) x_crds[i] = 1.0*i;
      for (int i=0; i<Ny0*Ny_cells+1; ++i) y_crds[i] = 1.0*i;
      for (int i=0; i<Nz0*Nz_cells+1; ++i) z_crds[i] = 1.0*i;
   } else if (geometryString == vlsv::geometry::STRING_CYLINDRICAL) {
      for (int i=0; i<Nx0*Nx_cells+1; ++i) x_crds[i] = 1.0*i;
      const double d_phi = 2*M_PI/(Ny0*Ny_cells);
      for (int i=0; i<Ny0*Ny_cells+1; ++i) y_crds[i] = i*d_phi;
      for (int i=0; i<Nz0*Nz_cells+1; ++i) z_crds[i] = 1.0*i;
   } else if (geometryString == vlsv::geometry::STRING_SPHERICAL) {
      for (int i=0; i<Nx0*Nx_cells+1; ++i) x_crds[i] = 1.0*i;
      const double d_theta = 30.0/(Ny0*Ny_cells) * M_PI/180.0;
      const double d_phi   = 2.0*M_PI/(Nz0*Nz_cells);
      for (int i=0; i<Ny0*Ny_cells+1; ++i) y_crds[i] = 75.0*M_PI/180.0 + i*d_theta;
      for (int i=0; i<Nz0*Nz_cells+1; ++i) z_crds[i] = i*d_phi;
   } else {
      cerr << "Unknown geometry" << endl; exit(1);
   }

   if (vlsv.writeArray("MESH_NODE_CRDS_X",attributes,Nx0*Nx_cells+1,1,x_crds) == false) success = false;
   if (vlsv.writeArray("MESH_NODE_CRDS_Y",attributes,Ny0*Ny_cells+1,1,y_crds) == false) success = false;
   if (vlsv.writeArray("MESH_NODE_CRDS_Z",attributes,Nz0*Nz_cells+1,1,z_crds) == false) success = false;
   
   delete [] x_crds;
   delete [] y_crds;
   delete [] z_crds;

   vector<uint64_t> globalIDs;
   for (int k=0; k<Nz0; ++k) for (int j=0; j<Ny0; ++j) for (int i=0; i<Nx0; ++i) {
      uint64_t globalID = calculateGlobalID(i,j,k,Nx0,Ny0,Nz0);
      globalIDs.push_back(globalID);
   }

   int N_refinements = 5;
   for (int i=0; i<N_refinements; ++i) {
      size_t refinedIndex = globalIDs.size() * 1.0*rand()/RAND_MAX;
      refineBlock(refinedIndex,globalIDs);
   }
   
   stringstream ss;
   ss << maxRefinementLevel;
   attributes["name"] = meshName;
   attributes["type"] = vlsv::mesh::STRING_UCD_AMR;
   attributes["max_refinement_level"] = ss.str();
   attributes["geometry"] = geometryString;
   if (geometryString == vlsv::geometry::STRING_CYLINDRICAL) attributes["yperiodic"] = "yes";
   if (geometryString == vlsv::geometry::STRING_SPHERICAL) attributes["zperiodic"] = "yes";
   if (vlsv.writeArray("MESH",attributes,globalIDs.size(),1,&(globalIDs[0])) == false) success = false;

   // Write number of total blocks, and the number of ghost blocks in domain:
   uint64_t domainSize[2];
   domainSize[0] = globalIDs.size();
   domainSize[1] = N_ghosts;
     
   attributes.clear();
   attributes["mesh"] = meshName;
   if (vlsv.writeArray("MESH_DOMAIN_SIZES",attributes,1,2,domainSize) == false) success = false;

   // Write ghost zones (not applicable here):
   uint64_t dummy;
   if (vlsv.writeArray("MESH_GHOST_LOCALIDS",attributes,N_ghosts,1,&dummy) == false) success = false;
   if (vlsv.writeArray("MESH_GHOST_DOMAINS",attributes,N_ghosts,1,&dummy) == false) success = false;

   vector<uint64_t> cellIDs;
   for (size_t block=0; block<globalIDs.size(); ++block) {
      for (uint32_t i=0; i<Nx_cells*Ny_cells*Nz_cells; ++i) {
	 cellIDs.push_back(globalIDs[block]);
      }
   }

   // Write global IDs as a variable:
   attributes.clear();
   if (dims == 2 && geometryString == vlsv::geometry::STRING_CARTESIAN) attributes["name"] = "cart_globalID_2d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_CARTESIAN) attributes["name"] = "cart_globalID_3d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_CYLINDRICAL) attributes["name"] = "cyl_globalID_3d";
   else if (dims == 3 && geometryString == vlsv::geometry::STRING_SPHERICAL) attributes["name"] = "sphe_globalID_3d";

   attributes["mesh"] = meshName;
   if (vlsv.writeArray("VARIABLE",attributes,cellIDs.size(),1,&(cellIDs[0])) == false) success = false;

   return success;   
}

int main(int argn,char* args[]) {
   bool success = true;

   // Init MPI:
   MPI_Init(&argn,&args);
   
   vlsv::Writer vlsv;
   if (vlsv.open("amr.vlsv",MPI_COMM_WORLD,0) == false) {
      success = false;
      MPI_Finalize();	
      return 1;
   }

   if (writeMesh(vlsv,vlsv::geometry::STRING_CARTESIAN,3) == false) success = false;
   if (writeMesh(vlsv,vlsv::geometry::STRING_CYLINDRICAL,3) == false) success = false;
   if (writeMesh(vlsv,vlsv::geometry::STRING_SPHERICAL,3) == false) success = false;
   
   if (vlsv.close() == false) {
	success = false;
   }
   
   MPI_Finalize();
   if (success == false) return 1;
   return 0;
}




