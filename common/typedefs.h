//=============================================================================
// Basic types.
//
// Copyright 2001 and onwards Kirill Shabunov
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//=============================================================================

#ifndef TYPEDEFS_H

#define TYPEDEFS_H

#include <stdint.h>

typedef unsigned char byte;

#include <limits.h>

#if (CHAR_BIT!=8)
#error Type char is assumed to be 8 bit.
#endif

#if (USHRT_MAX!=0xFFFF)
#error Type short is assumed to be 16 bit. Redefine int16 type in common/typedefs.h
#endif

typedef signed char int8;
typedef signed short int16;
typedef unsigned char uint8;
typedef unsigned short uint16;

#if (UINT_MAX==0xFFFFFFFFL)
typedef signed int int32;
typedef unsigned int uint32;
#elif (ULONG_MAX==0xFFFFFFFFL)
typedef signed long int32;
typedef unsigned long uint32;
#else
#error Unable to define unsigned 32 bit. Redefine int32 type in typedefs.h
#endif

#endif // #ifndef TYPEDEFS_H