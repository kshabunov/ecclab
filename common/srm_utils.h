//=============================================================================
// Common utilities for RM-like codes.
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

#ifndef SRM_UTILS_H

#define SRM_UTILS_H

//-----------------------------------------------------------------------------
// Includes.

#include "typedefs.h"

//-----------------------------------------------------------------------------
// Prototypes.

// Calculate dimension of RM(r, m) code.
// \sum_0^r {m \choose i}.
int calc_rm_k(int m, int r);

// Encode (m, m) code (for BSC).
// (Or decode, since the inverse matrix is the same.)
// Modified generator matrix.
int mrm_enc_mm_bsc(
   int m,
   int *x, // Information vector.
   int *y // Output codeword.
);

// Encode RM code (for BSC).
// Modified generator matrix.
int mrm_enc_bsc(
   int m, // RM m parameter.
   int r, // RM r parameter.
   int *x, // Information vector.
   int *y // Output codeword.
);

// Encode RM code with BPSK (0 --> -1, 1 --> +1).
// Modified generator matrix.
int mrm_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   int *x, // Information vector.
   double *y // Output codeword.
);

// Calculate dimension of SubRM(r, m) tr/P(node_table) code.
int calc_srm_k(int m, int r, uint32 *node_table);

// Encode SubRM(r, m) tr0/P(node_table) code (for BSC).
// Modified generator matrix.
int smrm_enc_bsc(
   int m,
   int r,
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
);

// Encode SubRM tr0/P(node_table) code with BPSK (0 --> -1, 1 --> +1).
// Modified generator matrix.
int smrm_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   double *y // Output codeword.
);

// Calculate dimension of SubRM(r, m) pr0/P(node_table) code.
int calc_par0_srm_k(int m, int r, uint32 *node_table);

// Encode SubRM(r, m) par0/P(node_table) code (for BSC).
int smrm_par0_enc_bsc(
   int m,
   int r,
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
);

// Encode SubRM par0/P(node_table) code with BPSK (0 --> -1, 1 --> +1).
int smrm_par0_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   double *y // Output codeword.
);

#endif // #ifndef SRM_UTILS_H