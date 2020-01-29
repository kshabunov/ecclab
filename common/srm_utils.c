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

#define SRM_UTILS_C

//-----------------------------------------------------------------------------
// Configuration switches.

//-----------------------------------------------------------------------------
// Includes.

#include <stdlib.h>
#include "srm_utils.h"
#include "../interfaces/ui_utils.h"

//-----------------------------------------------------------------------------
// Internal defines.

//-----------------------------------------------------------------------------
// Internal typedefs.

//-----------------------------------------------------------------------------
// Global data.

// For calc_srm_k(), smrm_enc_bsc_p().
int enc_utils_ntp;

//-----------------------------------------------------------------------------
// Functions.

// Calculate dimension of RM(r, m) code.
// \sum_0^r {m \choose i}.
int calc_rm_k(int m, int r) {
   if (r == 0) return 1;
   if (m == r) return 1 << m;
   return calc_rm_k(m - 1, r - 1) + calc_rm_k(m - 1, r);
}

// Encode (m, m) code (for BSC).
// (Or decode, since the inverse matrix is the same.)
// Modified generator matrix.
int mrm_enc_mm_bsc(
   int m,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2;
   int i;

   if (m == 0) {
      y[0] = x[0];
      return 1;
   }

   n2 = n >> 1;

   mrm_enc_mm_bsc(m - 1, x, y);
   mrm_enc_mm_bsc(m - 1, x + n2, y + n2);
   
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];

   return n;
}

// Encode RM code (for BSC).
// Modified generator matrix.
int mrm_enc_bsc(
   int m, // RM m parameter.
   int r, // RM r parameter.
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2, kv, ku;
   int i;

   if (r == 0) {
      for (i = 0; i < n; i++) y[i] = x[0];
      return 1;
   }

   n2 = n >> 1;

   if (r == m) {
      kv = mrm_enc_bsc(m - 1, m - 1, x, y);
      ku = mrm_enc_bsc(m - 1, m - 1, x + kv, y + n2);
   }
   else {
      kv = mrm_enc_bsc(m - 1, r - 1, x, y);
      ku = mrm_enc_bsc(m - 1, r, x + kv, y + n2);
   }
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];
   return kv + ku;
}

// Encode RM code with BPSK (0 --> -1, 1 --> +1).
// Modified generator matrix.
int mrm_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   int *x, // Information vector.
   double *y // Output codeword.
) 
{
   int n = 1 << m;
   int *c;
   int i;

   c = (int *)malloc(n * sizeof(int));
   if (c == NULL) {
      err_msg("mrm_enc_bpsk(): short of memory.");
      return 1;
   }
   mrm_enc_bsc(m, r, x, c);
   for (i = 0; i < n; i++) y[i] = (double)((c[i] << 1) - 1);
   free(c);

   return 0;
}

// Inner recursive procedure for calc_srm_k().
int calc_srm_k_p(int m, int r, uint32 *node_table) {
   if (r == 0) {
      if (node_table[enc_utils_ntp++] == 0) return 0;
      else return 1;
   }
   if (m == r) {
      if (node_table[enc_utils_ntp++] == 0) return 0;
      else return 1 << m;
   }
   return calc_srm_k_p(m - 1, r - 1, node_table) + calc_srm_k_p(m - 1, r, node_table);
}

// Calculate dimension of SubRM(r, m) tr0/P(node_table) code.
int calc_srm_k(int m, int r, uint32 *node_table) {
   enc_utils_ntp = 0;
   return calc_srm_k_p(m, r, node_table);
}

//
// Utils for tr0/P(node_table) codes, modified gen. matrix.
//

// Inner recursive procedure for smrm_enc_bsc().
int smrm_enc_bsc_p(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2, kv, ku;
   int i;

   if (r == 0) {
      if (node_table[enc_utils_ntp++] == 0) {
         for (i = 0; i < n; i++) y[i] = 0;
         return 0;
      }
      else {
         for (i = 0; i < n; i++) y[i] = x[0];
         return 1;
      }
   }
   if (r == m) {
      if (node_table[enc_utils_ntp++] == 0) {
         for (i = 0; i < n; i++) y[i] = 0;
         return 0;
      }
      else return mrm_enc_mm_bsc(m, x, y);
   }

   n2 = n >> 1;

   kv = smrm_enc_bsc_p(m - 1, r - 1, node_table, x, y);
   ku = smrm_enc_bsc_p(m - 1, r, node_table, x + kv, y + n2);
   
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];

   return kv + ku;
}

// Encode SubRM(r, m) tr0/P(node_table) code (for BSC).
// Modified generator matrix.
int smrm_enc_bsc(
   int m,
   int r,
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   enc_utils_ntp = 0;
   return smrm_enc_bsc_p(m, r, node_table, x, y);
}

// Encode SubRM tr0/P(node_table) code with BPSK (0 --> -1, 1 --> +1).
// Modified generator matrix.
int smrm_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   double *y // Output codeword.
) 
{
   int n = 1 << m;
   int *c;
   int i;

   c = (int *)malloc(n * sizeof(int));
   if (c == NULL) {
      err_msg("smrm_enc_bpsk(): short of memory.");
      return 1;
   }
   smrm_enc_bsc(m, r, node_table, x, c);
   for (i = 0; i < n; i++) y[i] = (double)((c[i] << 1) - 1);
   free(c);

   return 0;
}

//
// Utils for par0/P(node_table) codes
//

// Inner recursive procedure for calc_par0_srm_k().
int calc_par0_srm_k_p(int m, int r, uint32 *node_table) {
   if (r == 0) {
      if (node_table[enc_utils_ntp++] == 0) return 0;
      else return 1;
   }
   return calc_par0_srm_k_p(m - 1, r - 1, node_table)
      + calc_par0_srm_k_p(m - 1, (r == m) ? r - 1 : r, node_table);
}

// Calculate dimension of SubRM(r, m) pr0/P(node_table) code.
int calc_par0_srm_k(int m, int r, uint32 *node_table) {
   enc_utils_ntp = 0;
   return calc_par0_srm_k_p(m, r, node_table);
}

// Inner recursive procedure for smrm_par0_enc_bsc().
int smrm_par0_enc_bsc_p(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   int n = 1 << m;
   int n2, kv, ku;
   int i;

   if (r == 0) {
      if (node_table[enc_utils_ntp++] == 0) {
         for (i = 0; i < n; i++) y[i] = 0;
         return 0;
      }
      else {
         for (i = 0; i < n; i++) y[i] = x[0];
         return 1;
      }
   }
   n2 = n >> 1;

   kv = smrm_par0_enc_bsc_p(m - 1, r - 1, node_table, x, y);
   ku = smrm_par0_enc_bsc_p(m - 1, (r == m) ? r - 1 : r, node_table, x + kv, y + n2);
   
   for (i = 0; i < n2; i++) y[i] ^= y[i + n2];

   return kv + ku;
}

// Encode SubRM(r, m) par0/P(node_table) code (for BSC).
int smrm_par0_enc_bsc(
   int m,
   int r,
   uint32 *node_table,
   int *x, // Information vector.
   int *y // Output codeword.
) 
{
   enc_utils_ntp = 0;
   return smrm_par0_enc_bsc_p(m, r, node_table, x, y);
}

// Encode SubRM par0/P(node_table) code with BPSK (0 --> -1, 1 --> +1).
int smrm_par0_enc_bpsk(
   int m, // RM m parameter.
   int r, // RM r parameter.
   uint32 *node_table,
   int *x, // Information vector.
   double *y // Output codeword.
) 
{
   int n = 1 << m;
   int *c;
   int i;

   c = (int *)malloc(n * sizeof(int));
   if (c == NULL) {
      err_msg("smrm_par0_enc_bpsk(): short of memory.");
      return 1;
   }
   smrm_par0_enc_bsc(m, r, node_table, x, c);
   for (i = 0; i < n; i++) y[i] = (double)((c[i] << 1) - 1);
   free(c);

   return 0;
}
