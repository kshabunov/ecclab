//=============================================================================
// ML decoder for RM(1, m) codes.
//
// Copyright 1999 and onwards Kirill Shabunov
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

//-----------------------------------------------------------------------------
// Includes.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../common/typedefs.h"
#include "../common/std_defs.h"
#include "../common/count_ops.h"
#include "../formats/formats.h"
#include "rm1_ml.h"

//-----------------------------------------------------------------------------
// Functions.

int rm1_dec_sh(
   int c_m,
   double y_in[],
   int x_dec[]
)
{
   int c_n = 1 << c_m;
   int m, n, n2, i, j;
   ylitem *ylist, *yp1, *yp2, *yp3, *yp4, y1;
   slitem *slist, s1, s2;
   uint8 *mem_buf;

   // Allocate memory.
   mem_buf = (uint8 *)malloc(c_n * (c_m + 1) * sizeof(ylitem) + c_n * 2 * sizeof(slitem));
   if (mem_buf == NULL) return 1;
   slist = (slitem *)mem_buf;
   ylist = (ylitem *)(mem_buf + c_n * 2 * sizeof(slitem));

   // Calc. initial vector.

#ifdef FORMAT_EPS
   // Convert y_in to the epsilon format.
   for (i = 0; i < c_n; i++) {
      ylist[i] = exp(y_in[i]);
      ylist[i] = (1.0 - ylist[i]) / (1.0 + ylist[i]);
   }
#endif // FORMAT_EPS

#ifdef FORMAT_GAMMA
   // Convert y_in to the gamma format.
   for (i = 0; i < c_n; i++) {
      ylist[i] = exp(y_in[i]);
      ylist[i] = (1.0 - ylist[i]) / (1.0 + ylist[i]);
      ylist[i] = ((ylist[i] > 0.0) ? 1.0 : -1.0) - ylist[i];
   }
#endif // FORMAT_GAMMA

#ifdef FORMAT_RHO
   // We use likelihoods of 0.
   for (i = 0; i < c_n; i++) ylist[i] = -y_in[i];
#ifdef HARD_INPUT
   for (i = 0; i < c_n; i++) ylist[i] = (ylist[i] >= 0.0) ? 1.0 : -1.0;
#endif // HARD_INPUT
#endif // FORMAT_RHO

#ifdef FORMAT_Q
   for (i = 0; i < c_n; i++) ylist[i] = quantize(-y_in[i]);
#endif // FORMAT_Q

   for (i = 0; i < c_n; i++) slist[i] = 0.0;

   // All end nodes except the last.
   n = c_n;
   yp1 = ylist;
   for (m = c_m; m >= 1; m--) {
      n2 = n / 2;
      for (i = 0; i < c_n; i += n) {
         yp2 = yp1 + n2;
         yp3 = yp1 + c_n;
         yp4 = yp2 + c_n;
         s1 = 0.0;
         s2 = 0.0;
         for (j = 0; j < n2; j++) {
            y1 = XOR_EST(yp1[j], yp2[j]);
            if (CHECK_EST0(y1)) {
               s1 += EST0_TO_LNP0(y1);
               s2 += EST0_TO_LNP1(y1);
            }
            else {
               s1 += EST1_TO_LNP0(y1);
               s2 += EST1_TO_LNP1(y1);
            }
            ADD_EST(yp1[j], yp2[j], yp3[j]);
            ADD_EST(INV_EST(yp1[j]), yp2[j], yp4[j]);
         }
         slist[i + n2] = slist[i] + s2; // Likelihood corresp. to 1.
         slist[i] += s1; // Likelihood corresp. to 0.
         yp1 += n;
      } // for i.
      n /= 2;
   }

   // Last end node (0, 0).
   s2 = slist[0] + EST_TO_LNP0(yp1[0]);
   j = 0;
   for (i = 0; i < c_n; i++) {
      s1 = slist[i] + EST_TO_LNP0(yp1[i]);
      if (s1 > s2) { s2 = s1; j = i << 1; }
      s1 = slist[i] + EST_TO_LNP1(yp1[i]);
      if (s1 > s2) { s2 = s1; j = (i << 1) + 1; }
   }

   // "Unpack" j to get information bits.
   for (i = c_m; i >= 0; i--) {
      x_dec[i] = j & 1;
      j >>= 1;
   }

   // Free memory.
   free(mem_buf);

   return 0;
}

int rm1_dec_lst(
   int c_m,
   double y_in[],
   int lsiz, // size of the list (currently <= 4096).
   int x_dec[]
)
{
   int c_n = 1 << c_m;
   int m, n, n2, i, j, ja[4096];
   ylitem *ylist, *yp1, *yp2, *yp3, *yp4, y1;
   slitem *slist, s1, s2, sa[4096];
   uint8 *mem_buf;

   // Allocate memory.
   mem_buf = (uint8 *)malloc(c_n * (c_m + 1) * sizeof(ylitem) + c_n * 2 * sizeof(slitem));
   if (mem_buf == NULL) return 1;
   slist = (slitem *)mem_buf;
   ylist = (ylitem *)(mem_buf + c_n * 2 * sizeof(slitem));

   // Calc. initial vector.

#ifdef FORMAT_EPS
   // Convert y_in to the epsilon format.
   for (i = 0; i < c_n; i++) {
      ylist[i] = exp(y_in[i]);
      ylist[i] = (1.0 - ylist[i]) / (1.0 + ylist[i]);
   }
#endif // FORMAT_EPS

#ifdef FORMAT_GAMMA
   // Convert y_in to the gamma format.
   for (i = 0; i < c_n; i++) {
      ylist[i] = exp(y_in[i]);
      ylist[i] = (1.0 - ylist[i]) / (1.0 + ylist[i]);
      ylist[i] = ((ylist[i] > 0.0) ? 1.0 : -1.0) - ylist[i];
   }
#endif // FORMAT_GAMMA

#ifdef FORMAT_RHO
   // We use likelihoods of 0.
   for (i = 0; i < c_n; i++) ylist[i] = -y_in[i];
#ifdef HARD_INPUT
   for (i = 0; i < c_n; i++) ylist[i] = (ylist[i] >= 0.0) ? 1.0 : -1.0;
#endif // HARD_INPUT
#endif // FORMAT_RHO

#ifdef FORMAT_Q
   for (i = 0; i < c_n; i++) ylist[i] = quantize(-y_in[i]);
#endif // FORMAT_Q

   for (i = 0; i < c_n; i++) slist[i] = 0.0;

   // All end nodes except the last.
   n = c_n;
   yp1 = ylist;
   for (m = c_m; m >= 1; m--) {
      n2 = n / 2;
      for (i = 0; i < c_n; i += n) {
         yp2 = yp1 + n2;
         yp3 = yp1 + c_n;
         yp4 = yp2 + c_n;
         s1 = 0.0;
         s2 = 0.0;
         for (j = 0; j < n2; j++) {
            y1 = XOR_EST(yp1[j], yp2[j]);
            if (CHECK_EST0(y1)) {
               s1 += EST0_TO_LNP0(y1);
               s2 += EST0_TO_LNP1(y1);
            }
            else {
               s1 += EST1_TO_LNP0(y1);
               s2 += EST1_TO_LNP1(y1);
            }
            ADD_EST(yp1[j], yp2[j], yp3[j]);
            ADD_EST(INV_EST(yp1[j]), yp2[j], yp4[j]);
         }
         slist[i + n2] = slist[i] + s2; // Likelihood corresp. to 1.
         slist[i] += s1; // Likelihood corresp. to 0.
         yp1 += n;
      } // for i.
      n /= 2;
   }

   // Last end node (0, 0).
   for (j = 0; j < lsiz; j++) {
      sa[j] = -1e300;
      ja[j] = 1000;
      for (i = 0; i < c_n; i++) {
         s1 = slist[i] + EST_TO_LNP0(yp1[i]);
         if (s1 > sa[j]) {
            sa[j] = s1; ja[j] = i << 1; 
            m = i;
         }
         s1 = slist[i] + EST_TO_LNP1(yp1[i]);
         if (s1 > sa[j]) {
            sa[j] = s1; ja[j] = (i << 1) + 1; 
            m = i;
         }
      }
      slist[m]= -1e300;
   }

   // "Unpack" ja[] to get information bits.
   if (lsiz > c_n) lsiz = c_n;
   if (lsiz > 4096) lsiz = 4096;
   for (j = 0; j < lsiz; j++) {
      for (i = (j + 1) * (c_m + 1) - 1; i >= j * (c_m + 1); i--) {
         x_dec[i] = ja[j] & 1;
         ja[j] >>= 1;
      }
   }

   // Free memory.
   free(mem_buf);

   return 0;
}

// Same as rm1_dec_sh() except the input is considered
// to be already in the internal format.
int rm1_dec_sh_finp(
   int c_m,
   double y_in[],
   int x_dec[]
)
{
   int c_n = 1 << c_m;
   int m, n, n2, i, j;
   ylitem *ylist, *yp1, *yp2, *yp3, *yp4, y1;
   slitem *slist, s1, s2;
   uint8 *mem_buf;

   // Allocate memory.
   mem_buf = (uint8 *)malloc(c_n * (c_m + 1) * sizeof(ylitem) + c_n * 2 * sizeof(slitem));
   if (mem_buf == NULL) return 1;
   slist = (slitem *)mem_buf;
   ylist = (ylitem *)(mem_buf + c_n * 2 * sizeof(slitem));

   // Init the list.
   for (i = 0; i < c_n; i++) ylist[i] = y_in[i];
   for (i = 0; i < c_n; i++) slist[i] = 0.0;

   // All end nodes except the last.
   n = c_n;
   yp1 = ylist;
   for (m = c_m; m >= 1; m--) {
      n2 = n / 2;
      for (i = 0; i < c_n; i += n) {
         yp2 = yp1 + n2;
         yp3 = yp1 + c_n;
         yp4 = yp2 + c_n;
         s1 = 0.0;
         s2 = 0.0;
         for (j = 0; j < n2; j++) {
            y1 = XOR_EST(yp1[j], yp2[j]);
            ADD_OPS(0, 1);
            if (CHECK_EST0(y1)) {
               s1 += EST0_TO_LNP0(y1);
               s2 += EST0_TO_LNP1(y1);
               ADD_OPS(0, 2 * 2);
               ADD_OPS(1, 1 * 2);
            }
            else {
               s1 += EST1_TO_LNP0(y1);
               s2 += EST1_TO_LNP1(y1);
               ADD_OPS(0, 2 * 2);
               ADD_OPS(1, 1 * 2);
            }
            // ADD_OPS(2, 1);
            ADD_EST(yp1[j], yp2[j], yp3[j]);
            ADD_EST(INV_EST(yp1[j]), yp2[j], yp4[j]);
            ADD_OPS(0, 3 * 2);
            ADD_OPS(1, 1 * 2);
            ADD_OPS(3, 1);
         }
         slist[i + n2] = slist[i] + s2; // Likelihood corresp. to 1.
         slist[i] += s1; // Likelihood corresp. to 0.
         ADD_OPS(0, 1 * 2);
         yp1 += n;
      } // for i.
      n /= 2;
   }

   // Last end node (0, 0).
   s2 = slist[0] + EST_TO_LNP0(yp1[0]);
   j = 0;
   for (i = 0; i < c_n; i++) {
      s1 = slist[i] + EST_TO_LNP0(yp1[i]);
      ADD_OPS(0, 2);
      ADD_OPS(1, 1);
      if (s1 > s2) { s2 = s1; j = i << 1; }
      ADD_OPS(2, 1);
      s1 = slist[i] + EST_TO_LNP1(yp1[i]);
      ADD_OPS(0, 2);
      ADD_OPS(1, 1);
      if (s1 > s2) { s2 = s1; j = (i << 1) + 1; }
      ADD_OPS(2, 1);
   }

   // if (j != 0) printf("%d, %d, %g\n", j, c_m, slist[0]);
   // "Unpack" j to get information bits.
   for (i = c_m; i >= 0; i--) {
      x_dec[i] = j & 1;
      j >>= 1;
   }

   // Free memory.
   free(mem_buf);

   return 0;
}

// Max dot product.
int rm1_dec_sh_dotp(
   int c_m,
   double y_in[],
   int x_dec[]
)
{
   int c_n = 1 << c_m;
   int m, n, n2, i, j;
   ylitem *ylist, *yp1, *yp2, *yp3, *yp4;
   double s1, s2;
   uint8 *mem_buf;

   // Allocate memory.
   mem_buf = (uint8 *)malloc(c_n * (c_m + 1) * sizeof(ylitem));
   if (mem_buf == NULL) return 1;
   ylist = (ylitem *)(mem_buf);

   // Init the list.
   for (i = 0; i < c_n; i++) ylist[i] = y_in[i];

   // All end nodes except the last.
   n = c_n;
   yp1 = ylist;
   for (m = c_m; m >= 1; m--) {
      n2 = n / 2;
      for (i = 0; i < c_n; i += n) {
         yp2 = yp1 + n2;
         yp3 = yp1 + c_n;
         yp4 = yp2 + c_n;
         for (j = 0; j < n2; j++) {
            ADD_EST(yp1[j], yp2[j], yp3[j]);
            ADD_EST(INV_EST(yp1[j]), yp2[j], yp4[j]);
            ADD_OPS(0, 2);
            ADD_OPS(3, 1);
         }
         yp1 += n;
      } // for i.
      n /= 2;
   }

   // Last end node (0, 0).
   s2 = yp1[0];
   j = 0;
   for (i = 0; i < c_n; i++) {
      s1 = yp1[i];
      if (s1 > s2) { s2 = s1; j = i << 1; }
      ADD_OPS(2, 1);
      s1 = -yp1[i];
      ADD_OPS(3, 1);
      if (s1 > s2) { s2 = s1; j = (i << 1) + 1; }
      ADD_OPS(2, 1);
   }

   // if (j != 0) printf("%d, %d, %g\n", j, c_m, slist[0]);
   // "Unpack" j to get information bits.
   for (i = c_m; i >= 0; i--) {
      x_dec[i] = j & 1;
      j >>= 1;
   }

   // Free memory.
   free(mem_buf);

   return 0;
}
