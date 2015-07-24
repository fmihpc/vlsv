/** @file mpiconversion.h
 *  @brief Functions for converting basic C/C++ datatypes to MPI primitive datatypes (MPI_INT etc).
 * 
 *  These functions are meant to be 
 *  used when passing passing datatypes to MPI commands. For example, instead of
 *  @verbatim MPI_Send(...,MPI_INT,...) @endverbatim you would write 
 *  @verbatim MPI_Send(...,MPI_Type<int>(),...). @endverbatim
 * 
 *  This makes it painless to typedef datatypes in codes. For 
 *  example, you might usually want to typedef float as 'Real', but use doubles 
 *  in some rare cases where you need more accuracy. In this example you would 
 *  put MPI_Type<Real>() everywhere where MPI_Datatype is needed.
 * 
 *  @par License
 *  This file is part of VLSV file format.
 *  @n@n
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  @n@n
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *  @n@n
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *  
 * @copyright Copyright 2011-2013,2015 Finnish Meteorological Institute.
 * @author Arto Sandroos
 * */

#ifndef MPI_CONVERSION_H
#define MPI_CONVERSION_H
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

/** Return the MPI primitive datatype corresponding to the C++ datatype 
 * given as template parameter. The default function will print an error 
 * message to stderr and exit.
 * @tparam T C++ basic datatype that is to be converted into MPI datatype.
 * @return MPI primitive datatype.*/
template<typename T> inline MPI_Datatype MPI_Type() {
   std::cerr << "(mpiconversion.h): NULL datatype returned, byte size of original is " << sizeof(T) << std::endl;
   exit(1);
   return 0;
}

/** Overloaded template function that converts "char" to MPI_CHAR.
 * @tparam char C++ "char" datatype.
 * @return MPI_CHAR.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<char>() {return MPI_CHAR;}

// *****  SIGNED INTEGER TYPES ***** //

/** Overloaded template function that converts "signed char" to MPI_CHAR.
 * @tparam char C++ "signed char" datatype.
 * @return MPI_CHAR.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<signed char>() {return MPI_CHAR;}

/** Overloaded template function that converts "short" to MPI_SHORT.
 * @tparam char C++ "short" datatype.
 * @return MPI_SHORT.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<short>() {return MPI_SHORT;}

/** Overloaded template function that converts "int" to MPI_INT.
 * @tparam char C++ "int" datatype.
 * @return MPI_INT.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<int>() {return MPI_INT;}

/** Overloaded template function that converts "long int" to MPI_CHAR.
 * @tparam char C++ "long int" datatype.
 * @return MPI_LONG.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<long int>() {return MPI_LONG;}

/** Overloaded template function that converts "long long int" to MPI_LONG_LONG.
 * @tparam char C++ "long long int" datatype.
 * @return MPI_LONG_LONG.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<long long int>() {return MPI_LONG_LONG;}

// ***** UNSIGNED INTEGER TYPES ***** //

/** Overloaded template function that converts "unsigned char" to MPI_UNSIGNED_CHAR.
 * @tparam char C++ "unsigned char" datatype.
 * @return MPI_UNSIGNED_CHAR.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<unsigned char>() {return MPI_UNSIGNED_CHAR;}

/** Overloaded template function that converts "unsigned short" to MPI_UNSIGNED_SHORT.
 * @tparam char C++ "unsigned short" datatype.
 * @return MPI_UNSIGNED_SHORT.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<unsigned short>() {return MPI_UNSIGNED_SHORT;}

/** Overloaded template function that converts "unsigned int" to MPI_UNSIGNED.
 * @tparam char C++ "unsigned int" datatype.
 * @return MPI_UNSIGNED.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<unsigned int>() {return MPI_UNSIGNED;}

/** Overloaded template function that converts "unsigned long int" to MPI_UNSIGNED_LONG.
 * @tparam char C++ "unsigned long int" datatype.
 * @return MPI_UNSIGNED_LONG.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<unsigned long int>() {return MPI_UNSIGNED_LONG;}

/** Overloaded template function that converts "unsigned long long int" to MPI_UNSIGNED_LONG_LONG.
 * @tparam char C++ "unsigned long long int" datatype.
 * @return MPI_UNSIGNED_LONG_LONG.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<unsigned long long int>() {return MPI_UNSIGNED_LONG_LONG;}

// ***** FLOATING POINT TYPES ***** //

/** Overloaded template function that converts "float" to MPI_FLOAT.
 * @tparam char C++ "float" datatype.
 * @return MPI_FLOAT.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<float>() {return MPI_FLOAT;}

/** Overloaded template function that converts "double" to MPI_DOUBLE.
 * @tparam char C++ "double" datatype.
 * @return MPI_DOUBLE.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<double>() {return MPI_DOUBLE;}

/** Overloaded template function that converts "long double" to MPI_LONG_DOUBLE.
 * @tparam char C++ "long double" datatype.
 * @return MPI_LONG_DOUBLE.
 * @see MPI_Type(). */
template<> inline MPI_Datatype MPI_Type<long double>() {return MPI_LONG_DOUBLE;}

#endif
