//=============================================================================
// Common include for reliability recalculations.
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

#ifndef FORMATS_H

#define FORMATS_H

//-----------------------------------------------------------------------------
// Configuration switches.

//-----------------------------------------------------------------------------
// Includes.

#include <math.h>

//#include "format_gm1.h"
//#include "format_gm2.h"
//#include "format_gm3.h"

//#include "format_eps1.h" // eps, exact.
//#include "format_eps2.h"
//#include "format_eps3.h"
//#include "format_eps4.h" // eps., approx.
//#include "format_eps4h.h" // eps., approx., hard input.
#include "format_eps5.h" // eps, exact, clipping.
//#include "format_eps5h.h" // eps, exact, hard input, clipping.
//#include "format_eps6.h"
//#include "format_eps6a.h"
//#include "format_eps6b.h"
//#include "format_eps7.h"
//#include "format_eps8.h"

//#include "format_rho1.h" // LLR, exact.
//#include "format_rho2.h"
//#include "format_rho3.h"
//#include "format_rho3h.h"
//#include "format_rho4.h"
//#include "format_rho5.h" // LLR, min/sum.
//#include "format_rho6.h"
//#include "format_rho6a.h"
//#include "format_rho7.h"
//#include "format_rho7a.h"
//#include "format_rho8.h" // rho, exact, clipping.

//#include "format_q2.h"

#ifdef FORMAT_GAMMA

#include "gm_util.h"

#endif // FORMAT_GAMMA


//-----------------------------------------------------------------------------
// Defines

#ifdef FORMAT_EPS
// Convert from y to the epsilon format.

#ifdef HARD_INPUT
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) e[i] = (y[i] > 0.0) ? -0.6 : 0.6; \
}
#endif // if HARD_INPUT

#ifndef HARD_INPUT

#ifdef FORMAT_USE_CLIPPING
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) { \
      e[i] = exp(y[i]); \
      e[i] = (1.0 - e[i]) / (1.0 + e[i]); \
      FORMAT_CLIP(e[i]); \
   } \
}
#endif // if FORMAT_USE_CLIPPING

#ifndef FORMAT_USE_CLIPPING
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) { \
      e[i] = exp(y[i]); \
      e[i] = (1.0 - e[i]) / (1.0 + e[i]); \
   } \
}
#endif // if not FORMAT_USE_CLIPPING

#endif // if not HARD_INPUT

#endif // if FORMAT_EPS

#ifdef FORMAT_GAMMA
// Convert y to the gamma format.
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) { \
      e[i] = exp(y[i]); \
      e[i] = (1.0 - e[i]) / (1.0 + e[i]); \
      e[i] = ((e[i] > 0.0) ? 1.0 : -1.0) - e[i]; \
   } \
}
#endif // if FORMAT_GAMMA

#ifdef FORMAT_RHO
// We use likelihoods of 0.
#ifdef HARD_INPUT
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) e[i] = (y[i] > 0.0) ? -1.0 : 1.0; \
}
#else // if HARD_INPUT
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) e[i] = -y[i]; \
}
#endif // else if HARD_INPUT
#endif // if FORMAT_RHO

#ifdef FORMAT_RHO4
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   for (i = 0; i < n; i++) e[i] = quantize(-y[i]); \
}
#endif // if FORMAT_RHO4

#ifdef FORMAT_Q2
#define VY_TO_FORMAT(y, e, n) { \
   int i; \
   double d1; \
   for (i = 0; i < n; i++) { \
      d1 = exp(y[i]); \
      d1 = (1.0 - d1) / (1.0 + d1); \
      e[i] = quantize(d1); \
   } \
}
#endif // if FORMAT_Q2

//-----------------------------------------------------------------------------
// Typedefs.

//-----------------------------------------------------------------------------
// Global data.

#ifndef FORMATS_C


#endif // #ifndef FORMATS_C

//-----------------------------------------------------------------------------
// Prototypes.

#endif // #ifndef FORMATS_H