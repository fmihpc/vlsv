#ifndef MPI_CONVERSION_H
#define MPI_CONVERSION_H
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <stdint.h>

// Overloaded templates which return the corresponding data type
// for C++ native data types. For example, if float has been 
// typedef'd as Real, then MPI_Type<Real>() should return MPI_FLOAT.
// If you later on change the typedef to double, MPI_Type<Real>() 
// still works.
template<typename T> inline MPI_Datatype MPI_Type() {
   std::cerr << "(mpiconversion.h): NULL datatype returned, byte size of original is " << sizeof(T) << std::endl;
   exit(1);
   return 0;
}

template<> inline MPI_Datatype MPI_Type<char>() {return MPI_CHAR;}

// Signed integer types
template<> inline MPI_Datatype MPI_Type<int8_t>() {return MPI_CHAR;}
template<> inline MPI_Datatype MPI_Type<int16_t>() {return MPI_SHORT;}
template<> inline MPI_Datatype MPI_Type<int32_t>() {return MPI_INT;}
template<> inline MPI_Datatype MPI_Type<int64_t>() {return MPI_LONG;}
template<> inline MPI_Datatype MPI_Type<long long int>() {return MPI_LONG_LONG;}

// Unsigned integer types
template<> inline MPI_Datatype MPI_Type<uint8_t>() {return MPI_UNSIGNED_CHAR;}
template<> inline MPI_Datatype MPI_Type<uint16_t>() {return MPI_UNSIGNED_SHORT;}
template<> inline MPI_Datatype MPI_Type<uint32_t>() {return MPI_UNSIGNED;}
template<> inline MPI_Datatype MPI_Type<uint64_t>() {return MPI_UNSIGNED_LONG;}
template<> inline MPI_Datatype MPI_Type<unsigned long long int>() {return MPI_UNSIGNED_LONG_LONG;}

// Floating point types
template<> inline MPI_Datatype MPI_Type<float>() {return MPI_FLOAT;}
template<> inline MPI_Datatype MPI_Type<double>() {return MPI_DOUBLE;}
template<> inline MPI_Datatype MPI_Type<long double>() {return MPI_LONG_DOUBLE;}

#endif
