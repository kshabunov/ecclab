//=============================================================================
// Defines for reliability recalculations.
// Exact for epsilon + clipping: |e| < 1 - CLIP_EPS
// ev = e1 * e2.
// eu = (e1 + e2) / (1 + e1 * e2).
// Metric: ln(2 * p).
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

#ifndef FORMAT_EPS5_H

#define FORMAT_EPS5_H

#define FORMAT_EPS
#define FORMAT_USE_CLIPPING

#define CLIP_EPS 1e-8
#define CLIP_T (1.0 - CLIP_EPS)

#define FORMAT_CLIP(e) { \
   if ((e) > CLIP_T) (e) = CLIP_T; \
   if ((e) < -CLIP_T) (e) = -CLIP_T; \
}

#define STRONGEST_EST0 (1.0)

#define YLDEC0 1.0
#define YLDEC1 (-1.0)

#define EPS2LOGP0(e) log(1.0 + (e))
#define EPS2LOGP1(e) log(1.0 - (e))

// True if estimate votes for 0.
#define CHECK_EST0(y) ((y) > 0.0)

// True if e1 is a stronger estimate of 0 than est2.
#define STRONGER_EST0(e1, e2) ((e1) > (e2))

// Invert the estimate.
#define INV_EST(y) (-(y))

// ln(Pr{y = 0 | y}).
#define EST_TO_LNP0(y) EPS2LOGP0(y)

// ln(Pr{y = 1 | y}).
#define EST_TO_LNP1(y) EPS2LOGP1(y)

// ln(Pr{y = 0 | y}) for y voting for 0.
#define EST0_TO_LNP0(y) EPS2LOGP0(y)
// ln(Pr{y = 0 | y}) for y voting for 1.
#define EST1_TO_LNP0(y) EPS2LOGP0(y)

// ln(Pr{y = 1 | y}) for y voting for 0.
#define EST0_TO_LNP1(y) EPS2LOGP1(y)
// ln(Pr{y = 1 | y}) for y voting for 1.
#define EST1_TO_LNP1(y) EPS2LOGP1(y)

// EST0_TO_LNP1(y) - EST0_TO_LNP0(y).
#define EST0_ADD_LNP1_SUB_LNP0(y) log((1.0 - (y)) / (1.0 + (y)))

// Convert estimate to log likelihood ratio.
#define EST2LLR(e) log((1.0 + (e)) / (1.0 - (e)))

#define XOR_EST(y1, y2) ((y1) * (y2))

#define ADD_EST(y1, y2, y3) { \
   y3 = ((y1) + (y2)) / (1.0 + (y1) * (y2)); \
   FORMAT_CLIP(y3); \
}

#define XOR_YLDEC(y1, y2) ((y1) * (y2))

// if (y == YLDEC0) return e; else return INV_EST(e);
#define EST_XOR_YLDEC(e, y) ((e) * (y))

//-----------------------------------------------------------------------------
// List items typedefs.

typedef int xlitem;
typedef double ylitem;
typedef double slitem;

#endif // #ifndef FORMAT_EPS5_H