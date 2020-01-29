//=============================================================================
// Interface for dtrm_glp_inner.c
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

#ifndef DTRM_GLP_INNER_H

#define DTRM_GLP_INNER_H

//-----------------------------------------------------------------------------
// Defines

#define FLSIZ_MULT 4


//-----------------------------------------------------------------------------
// Includes.

#include "../common/typedefs.h"
#include "../formats/formats.h"


//-----------------------------------------------------------------------------
// Typedefs.

typedef struct xlist_item {
  xlitem x;
  struct xlist_item* p;
} xlist_item;

typedef struct {
  int c_r; // RM r parameter.
  int c_m; // RM m parameter.
  int c_n; // Code length.
  int c_k; // Code dimension.
  int peak_lsiz; // Peak list size.
  uint32 *node_table;
  int node_table_len;
  int p_num; // # of permutations.
  int *pxarr; // permutations for inf. seq. (after decoding).
  int *pyarr; // perm. for ch. output vector (before decoding).
  int ret_s_sum; // If 1 - return SUM s_i for list.
  uint8 *mem_buf;
} decoder_type;


//-----------------------------------------------------------------------------
// Prototypes.


// Decoding procedure.
int
rm_dec(
  decoder_type *dd, // Decoder instance data.
  double *y_in, // Decoder input.
  int *x_dec, // Decoded information sequence.
  double *s_dec // Metric of x_dec (set NULL, if not needed).
);

#endif // #ifndef DTRM_GLP_INNER_H