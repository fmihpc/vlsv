/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2014 Finnish Meteorological Institute
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
#include <cmath>
#include <map>
#include <algorithm>

#include <vlsv_writer.h>

#include "amr_mesh.h"

using namespace std;

AmrMesh::AmrMesh(const uint32_t& Nx0,const uint32_t& Ny0,const uint32_t& Nz0,const uint32_t& xCells,
		 const uint32_t& yCells,const uint32_t& zCells,const uint8_t& maxRefLevel) {
   bbox[0] = Nx0;
   bbox[1] = Ny0;
   bbox[2] = Nz0;
   N_blocks0       = bbox[0]*bbox[1]*bbox[2];

   bbox[3]    = xCells;
   bbox[4]    = yCells;
   bbox[5]    = zCells;
   refLevelMaxAllowed = maxRefLevel;
   
   initialized = false;
}

AmrMesh::~AmrMesh() { }

/** Get a const iterator pointing to the first existing block in mesh.
 * @return Iterator pointing to the first existing block.*/
const std::unordered_map<uint64_t,uint8_t>::iterator AmrMesh::begin() {
   return globalIDs.begin();
}

/** Attempt to coarsen the given block. Coarsen will not succeed if it 
 * would create a block with more than one refinement level difference 
 * between it and its neighbors.
 * @param globalID Global ID of the block.
 * @return If true, block was coarsened.*/
bool AmrMesh::coarsen(const uint64_t& globalID) {
   if (globalIDs.find(globalID) == globalIDs.end()) return false;
   
   uint32_t refLevel,i,j,k;
   getIndices(globalID,refLevel,i,j,k);
   if (refLevel == 0) return false;

   vector<uint64_t> nbrs;
   getSiblingNeighbors(globalID,nbrs);

   // Check if block can be coarsened:
   for (size_t n=0; n<nbrs.size(); ++n) {
      vector<uint64_t> children;
      getChildren(nbrs[n],children);
      for (size_t c=0; c<children.size(); ++c) {
	 if (globalIDs.find(children[c]) != globalIDs.end()) return false;
      }
   }

   // Remove the block and its siblings:
   vector<uint64_t> siblings;
   getSiblings(globalID,siblings);
   for (size_t s=0; s<siblings.size(); ++s) {
      globalIDs.erase(siblings[s]);
   }
   
   // Insert parent:
   globalIDs.insert(make_pair(getParent(globalID),refLevel-1));
   return true;
}

/** Get a const iterator pointing past the last existing block.
 * @param Iterator pointing past the last existing block.*/
const std::unordered_map<uint64_t,uint8_t>::iterator AmrMesh::end() {
   return globalIDs.end();
}

/** Finalize the class. Deallocates all memory.
 * @return If true, class finalized correctly.*/
bool AmrMesh::finalize() {return true;}

/** Get global IDs of block's children. Note that the children may 
 * or may not exists -- this function simply calculates the global IDs.
 * @param globalID Global ID of the block.
 * @param children Vector where global IDs of children are inserted.*/
void AmrMesh::getChildren(const uint64_t& globalID,std::vector<uint64_t>& children) {
   children.clear();
   
   uint32_t refLevel,i,j,k;
   getIndices(globalID,refLevel,i,j,k);
   if (refLevel+1 > refLevelMaxAllowed) return;
   
   i *= 2;
   j *= 2;
   k *= 2;
   
   children.push_back(getGlobalID(refLevel+1,i  ,j  ,k  ));
   children.push_back(getGlobalID(refLevel+1,i+1,j  ,k  ));
   children.push_back(getGlobalID(refLevel+1,i  ,j+1,k  ));
   children.push_back(getGlobalID(refLevel+1,i+1,j+1,k  ));
   children.push_back(getGlobalID(refLevel+1,i  ,j  ,k+1));
   children.push_back(getGlobalID(refLevel+1,i+1,j  ,k+1));
   children.push_back(getGlobalID(refLevel+1,i  ,j+1,k+1));
   children.push_back(getGlobalID(refLevel+1,i+1,j+1,k+1));
}

/** Get the global ID of a block with the given indices and refinement level. 
 * The refinement level must be equal or greater than zero, and less than or equal 
 * to the maximum refinement level.
 * @param refLevel Block's refinement level.
 * @param i i-index of the block at given refinement level.
 * @param j j-index of the block at given refinement level.
 * @param k k-index of the block at given refinement level.
 * @return Global ID of the block.*/
uint64_t AmrMesh::getGlobalID(const uint32_t& refLevel,const uint32_t& i,const uint32_t& j,const uint32_t& k) {
   const uint32_t multiplier = pow(2,refLevel);
   return offsets[refLevel] + k*bbox[1]*bbox[0]*multiplier*multiplier + j*bbox[0]*multiplier + i;
}

/** Get i,j,k indices of the block with given global ID, and its refinement level.
 * @param globalID Global ID of the block.
 * @param refLevel Block's refinement level is written here.
 * @param i Block's i-index is written here.
 * @param j Block's j-index is written here.
 * @param k Block's k-index is written here.*/
void AmrMesh::getIndices(const uint64_t& globalID,uint32_t& refLevel,uint32_t& i,uint32_t& j,uint32_t& k) {
   refLevel   = upper_bound(offsets.begin(),offsets.end(),globalID)-offsets.begin()-1;
   const uint64_t cellOffset = offsets[refLevel];

   const uint64_t multiplier = pow(2,refLevel);
   const uint64_t Nx = bbox[0] * multiplier;
   const uint64_t Ny = bbox[1] * multiplier;

   uint64_t index = globalID - cellOffset;
   k = index / (Ny*Nx);
   index -= k*Ny*Nx;
   j = index / Nx;
   i = index - j*Nx;
}

/** Get global IDs of block's neighbors. The neighbor IDs are calculated at the 
 * same refinement level as the block. Thus, some of the returned neighbors may 
 * not actually exist.
 * @param globalID Global ID of the block.
 * @param neighboIDs Vector where neighbors global IDs are written to.*/
void AmrMesh::getNeighbors(const uint64_t& globalID,std::vector<uint64_t>& neighborIDs) {
   neighborIDs.clear();

   uint32_t i,j,k,refLevel;
   getIndices(globalID,refLevel,i,j,k);

   const uint32_t Nx_max = bbox[0] * pow(2,refLevel);
   const uint32_t Ny_max = bbox[1] * pow(2,refLevel);
   const uint32_t Nz_max = bbox[2] * pow(2,refLevel);
   for (int k_off=-1; k_off<2; ++k_off) {
      if (k+k_off >= Nz_max) continue;
      for (int j_off=-1; j_off<2; ++j_off) {
	 if (j+j_off >= Ny_max) continue;
	 for (int i_off=-1; i_off<2; ++i_off) {
	    if (i+i_off >= Nx_max) continue;
	    if (i_off == 0 && (j_off == 0 && k_off == 0)) continue;
	    neighborIDs.push_back(getGlobalID(refLevel,i+i_off,j+j_off,k+k_off));
	 }
      }
   }
}

/** Get global ID of block's parent. If the block is at refinement level 0, 
 * i.e., at the base grid level, block's global ID is returned instead.
 * @param globalID Global ID of the block.
 * @return Global ID of block's parent.*/
uint64_t AmrMesh::getParent(const uint64_t& globalID) {
   uint32_t refLevel,i,j,k;
   getIndices(globalID,refLevel,i,j,k);

   if (refLevel == 0) return globalID;
   
   i /= 2;
   j /= 2;
   k /= 2;
   return getGlobalID(refLevel-1,i,j,k);
}

/** Get global IDs of all neighbors of this block and its siblings 
 * at the same refinement level as the block.
 * If the block is not at the boundary of the simulation domain, the 
 * returned vector should contain 56 neighbor IDs.
 * @param globalID Global ID of the block.
 * @param nbrs Vector where the neighbor IDs are written to.*/
void AmrMesh::getSiblingNeighbors(const uint64_t& globalID,std::vector<uint64_t>& nbrs) {
   nbrs.clear();

   uint32_t refLevel,i,j,k;
   getIndices(globalID,refLevel,i,j,k);

   i -= (i % 2);
   j -= (j % 2);
   k -= (k % 2);
   
   for (int k_off=-1; k_off<3; ++k_off) {
     for (int j_off=-1; j_off<3; ++j_off) {
       for (int i_off=-1; i_off<3; ++i_off) {
	  int cntr=0;
	  if (i_off == 0 || i_off == 1) ++cntr;
	  if (j_off == 0 || j_off == 1) ++cntr;
	  if (k_off == 0 || k_off == 1) ++cntr;
	  if (cntr == 3) continue;
	  nbrs.push_back(getGlobalID(refLevel,i+i_off,j+j_off,k+k_off));
       }
     }
   }
}

/** Get global IDs of block's siblings at the same refinement level as the block. 
 * Note that some of the siblings may not exists, for example, 
 * if the sibling has been refined. Returned vector also contains the ID of this block.
 * @param globalID Global ID of the block.
 * @param siblings Vector where sibling's global IDs are written to.*/
void AmrMesh::getSiblings(const uint64_t& globalID,std::vector<uint64_t>& siblings) {
   siblings.clear();

   uint32_t refLevel,i,j,k;
   getIndices(globalID,refLevel,i,j,k);

   i -= (i % 2);
   j -= (j % 2);
   k -= (k % 2);
   
   siblings.push_back(getGlobalID(refLevel,i  ,j  ,k  ));
   siblings.push_back(getGlobalID(refLevel,i+1,j  ,k  ));
   siblings.push_back(getGlobalID(refLevel,i  ,j+1,k  ));
   siblings.push_back(getGlobalID(refLevel,i+1,j+1,k  ));
   siblings.push_back(getGlobalID(refLevel,i  ,j  ,k+1));
   siblings.push_back(getGlobalID(refLevel,i+1,j  ,k+1));
   siblings.push_back(getGlobalID(refLevel,i  ,j+1,k+1));
   siblings.push_back(getGlobalID(refLevel,i+1,j+1,k+1));
}

/** Initialize the mesh.
 * @param xmin
 * @param xmax
 * @param ymin
 * @param ymax
 * @param zmin
 * @param zmax
 * @return If true, mesh was successfully initialized and is ready for use.*/
bool AmrMesh::initialize(const double& xmin,const double& xmax,const double& ymin,const double& ymax,
			 const double& zmin,const double& zmax,const uint8_t refLevel) {   
   if (initialized == true) return true;
   initialized = false;
   if (refLevel > refLevelMaxAllowed) return false;

   // Calculate block global ID offsets for each refinement level:
   offsets.resize(refLevelMaxAllowed+1);
   offsets[0] = 0;
   for (size_t i=1; i<refLevelMaxAllowed+1; ++i) {
      offsets[i] = offsets[i-1] + N_blocks0 * pow(8,i-1);
   }

   // Insert all blocks at given refinement level to mesh:
   const uint32_t factor = pow(2,refLevel);
   
   size_t counter = 0;
   for (uint32_t k=0; k<bbox[2]*factor; ++k) 
     for (uint32_t j=0; j<bbox[1]*factor; ++j)
       for (uint32_t i=0; i<bbox[0]*factor; ++i) {
	  globalIDs.insert(make_pair(getGlobalID(refLevel,i,j,k),refLevel));
       }
   
   meshLimits[0] = xmin;
   meshLimits[1] = xmax;
   meshLimits[2] = ymin;
   meshLimits[3] = ymax;
   meshLimits[4] = zmin;
   meshLimits[5] = zmax;
   
   initialized = true;
   return initialized;
}

/** Refine the block. This function will additionally refine block's neighbors, 
 * if it is necessary to maintain maximum difference of one refinement level 
 * between neighboring blocks.
 * @param globalID Global ID of the block.
 * @return If true, block was refined.*/
bool AmrMesh::refine(const uint64_t& globalID) {
   if (globalIDs.find(globalID) == globalIDs.end()) return false;

   vector<uint64_t> nbrs;
   getNeighbors(globalID,nbrs);

   uint32_t i,j,k,refLevel;
   getIndices(globalID,refLevel,i,j,k);
   if (refLevel == refLevelMaxAllowed) {
      return false;
   }

   i *= 2;
   j *= 2;
   k *= 2;
   globalIDs.erase(globalID);
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i  ,j  ,k  ),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i+1,j  ,k  ),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i  ,j+1,k  ),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i+1,j+1,k  ),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i  ,j  ,k+1),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i+1,j  ,k+1),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i  ,j+1,k+1),refLevel+1));
   globalIDs.insert(make_pair(getGlobalID(refLevel+1,i+1,j+1,k+1),refLevel+1));

   // Enforce that neighbors have at maximum one refinement level difference:
   for (size_t n=0; n<nbrs.size(); ++n) {
      // Get parent of the neighbor:
      const uint64_t parentID = getParent(nbrs[n]);
      
      // If the parent exists, it is at two refinement levels higher than 
      // the block that was refined above. Refine parent of neighbor:
      if (parentID == nbrs[n]) continue;
      unordered_map<uint64_t,uint8_t>::iterator parent = globalIDs.find(parentID);
      if (parent != globalIDs.end()) {
	 refine(parentID);
      }
   }

   return true;
}

/** Get the number of blocks in the mesh.
 * @return Number of blocks in the mesh.*/
size_t AmrMesh::size() const {
   return globalIDs.size();
}

/** Write the mesh to given file.
 * @param fileName Name of the output file.
 * @return If true, mesh was successfully written.*/
bool AmrMesh::write(const std::string fileName) {
   if (initialized == false) return false;
   bool success = true;

   const string meshName = "amr_mesh";
   
   vlsv::Writer vlsv;
   if (vlsv.open(fileName,MPI_COMM_WORLD,0) == false) {
      success = false;
      return false;
   }

   // Write block global IDs:
   map<string,string> attributes;
   attributes["name"] = meshName;
   attributes["type"] = vlsv::mesh::STRING_UCD_AMR;
   stringstream ss;
   ss << (uint32_t)refLevelMaxAllowed;
   attributes["max_refinement_level"] = ss.str();
   attributes["geometry"] = vlsv::geometry::STRING_CARTESIAN;
     {
	vector<uint64_t> buffer;
	for (unordered_map<uint64_t,uint8_t>::iterator it=globalIDs.begin(); it!=globalIDs.end(); ++it) {
	   buffer.push_back(it->first);
	}
	if (vlsv.writeArray("MESH",attributes,globalIDs.size(),1,&(buffer[0])) == false) success = false;
     }

   // Write mesh bounding box:
   attributes.clear();
   attributes["mesh"] = meshName;
   if (vlsv.writeArray("MESH_BBOX",attributes,6,1,bbox) == false) success = false;
   
   // Write domain sizes:
   uint64_t domainSize[2];
   domainSize[0] = globalIDs.size();
   domainSize[1] = 0;
   if (vlsv.writeArray("MESH_DOMAIN_SIZES",attributes,1,2,domainSize) == false) success = false;

   // Write ghost zone data (not applicable here):
   uint64_t dummy;
   if (vlsv.writeArray("MESH_GHOST_LOCALIDS",attributes,domainSize[1],1,&dummy) == false) success = false;
   if (vlsv.writeArray("MESH_GHOST_DOMAINS",attributes,domainSize[1],1,&dummy) == false) success = false;

   // Write node coordinates:
   vector<float> coords;
   coords.resize(bbox[0]*bbox[3]+1);
   float dx = (meshLimits[1]-meshLimits[0])/(bbox[0]*bbox[3]);
   for (uint32_t i=0; i<bbox[0]*bbox[3]+1; ++i) {
      coords[i] = meshLimits[0] + i*dx;
   }
   if (vlsv.writeArray("MESH_NODE_CRDS_X",attributes,coords.size(),1,&(coords[0])) == false) success = false;
   
   coords.resize(bbox[1]*bbox[4]+1);
   dx = (meshLimits[3]-meshLimits[2])/(bbox[1]*bbox[4]);
   for (uint32_t i=0; i<bbox[1]*bbox[4]+1; ++i) {
      coords[i] = meshLimits[2] + i*dx;
   }
   if (vlsv.writeArray("MESH_NODE_CRDS_Y",attributes,coords.size(),1,&(coords[0])) == false) success = false;
   
   coords.resize(bbox[2]*bbox[5]+1);
   dx = (meshLimits[5]-meshLimits[4])/(bbox[2]*bbox[5]);
   for (uint32_t i=0; i<bbox[2]*bbox[5]+1; ++i) {
      coords[i] = meshLimits[4] + i*dx;
   }
   if (vlsv.writeArray("MESH_NODE_CRDS_Z",attributes,coords.size(),1,&(coords[0])) == false) success = false;
   
   return success;
}

