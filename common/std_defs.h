//=============================================================================
// Standard defines.
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

#ifndef STD_DEFS_H

#define STD_DEFS_H

//-----------------------------------------------------------------------------
// Configuration switches.

//-----------------------------------------------------------------------------
// Includes.

//-----------------------------------------------------------------------------
// Defines

// Memory
#define CHK_FREE(p) {if ((p) != NULL) { free(p); (p) = NULL; }}

// General math
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Float / double
#define FLOAT_EPS 1e-10
#define FLOAT_EPS2ONE (1.0 - FLOAT_EPS)
#define SGNF(x) (((x) >= 0.0) ? 1.0 : -1.0)
#define ABSF(x) (((x) >= 0.0) ? (x) : -(x))
#define FLOAT_EQ(x, y) (ABSF((x) - (y)) < FLOAT_EPS * ABSF(x))
#define FLOAT_EQE(x, y, e) (ABSF((x) - (y)) < ((e) * ABSF(x)))


//-----------------------------------------------------------------------------
// Typedefs.

//-----------------------------------------------------------------------------
// Global data.

#ifndef STD_DEFS_C


#endif // #ifndef STD_DEFS_C

//-----------------------------------------------------------------------------
// Prototypes.

#endif // #ifndef STD_DEFS_H