/** This file is part of VLSV file format.
 * 
 *  Copyright 2011-2016 Finnish Meteorological Institute
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

#ifndef VLSV_COMMON_MPI_H
#define VLSV_COMMON_MPI_H

#include <mpi.h>
#include "vlsv_common.h"

namespace vlsv {
   uint64_t getMaxBytesPerRead();
   uint64_t getMaxBytesPerWrite();

   bool broadcast(const std::string& input,std::string& output,MPI_Comm comm,const int& masterRank);
   bool checkSuccess(const bool& myStatus,MPI_Comm comm);
   MPI_Datatype getMPIDatatype(datatype::type dt,uint64_t dataSize);
}

#endif
