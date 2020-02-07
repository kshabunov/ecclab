//=============================================================================
// dtrm0
// Recursive decoding on triangle down to (0, m) and (m, m).
// Needs channel sigma (define global DEC_NEEDS_SIGMA).
// Parameters: RM_m, RM_r
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

#include <stdlib.h>
#include <string.h>
#include "../common/typedefs.h"
#include "../interfaces/codec.h"
#include "../common/spf_par.h"
#include "../interfaces/ui_utils.h"
#include "../common/srm_utils.h"
#include "../formats/formats.h"

#ifndef FORMAT_EPS
#error Works only in eps format.
#endif

//-----------------------------------------------------------------------------
// Internal defines.

// Max # of SNR values.
#define VARL_NUM_MAX        50

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
   int c_n; // Code length.
   int c_k; // Code dimension.
   int c_m; // RM m.
   int c_r; // RM r.
   int csnrn;
   double sg22; // 2 / sg^2 for calculation of epsilons.
} cdc_inst_type;


//-----------------------------------------------------------------------------
// Global data.

//-----------------------------------------------------------------------------
// Internal functions.

int
tri0_dec(
   double *y_in,
   int r,
   int m,
   int n,
   double *y_dec,
   int *x_dec,
   int *x_dec_len,
   double *aux_buf
)
{
   int i, n2, cur_x_dec_len;
   double p1, p2;
   double *y1, *y2, *y_dec1, *y_dec2;

   // RM(0, m) case.
   if (r == 0) {

      (*x_dec_len) = 1; // one information bit.

#ifdef FORMAT_EPS4_H
      p1 = 0.0;
      for (i = 0; i < n; i++) {
         p1 += y_in[i];
      }
      if (p1 < 0.0) {
#else // FORMAT_EPS4_H
      p1 = 1.0;
      p2 = 1.0;
      for (i = 0; i < n; i++) {
         p1 *= 1.0 - y_in[i];
         p2 *= 1.0 + y_in[i];
      }
      if (p1 > p2) {
#endif // else FORMAT_EPS4_H
         x_dec[0] = 1;
         for (i = 0; i < n; i++) y_dec[i] = -1.0;
      }
      else {
         x_dec[0] = 0;
         for (i = 0; i < n; i++) y_dec[i] = 1.0;
      }

      return 0;

   } // if (r == 0).

   // RM(m, m) case.
   if (r == m) {

      int tmp[1024];

      (*x_dec_len) = n; // n information bits.

      for (i = 0; i < n; i++) {
         if (y_in[i] > 0) {
            tmp[i] = 0;
            y_dec[i] = 1.0;
         }
         else {
            tmp[i] = 1;
            y_dec[i] = -1.0;
         }
      }
      mrm_enc_mm_bsc(m, tmp, x_dec);

      return 0;

   } // if (r == m).

   // Further reduction.

   n2 = n / 2;
   y1 = y_in;
   y2 = y_in + n2;
   y_dec1 = y_dec;
   y_dec2 = y_dec + n2;

   // Decode u1.
   // [x_dec1, y_dec1] = dtrm00_she(r - 1, m - 1, y1 .* y2);
   for (i = 0; i < n2; i++) aux_buf[i] = y1[i] * y2[i];
   tri0_dec(
      aux_buf, r - 1, m - 1, n2,
      y_dec1, x_dec, &cur_x_dec_len, aux_buf + n2
   );
   (*x_dec_len) = cur_x_dec_len;

   // Decode u2.
   // [x_dec2, y_dec2] = dtrm00_she(r, m - 1, add_est_eps(y1 .* y_dec1, y2, n2));
   for (i = 0; i < n2; i++) y1[i] *= y_dec1[i];
   VADD_EST(y1, y2, aux_buf, n2);
   // ADD_EST_EPS(y1, y2, aux_buf, n2);
   tri0_dec(
      aux_buf, r, m - 1, n2,
      y_dec2, x_dec + cur_x_dec_len, &cur_x_dec_len, aux_buf + n2
   );
   (*x_dec_len) += cur_x_dec_len;

   // y_dec = [y_dec1 .* y_dec2, y_dec2];
   for (i = 0; i < n2; i++) y_dec1[i] *= y_dec2[i];

   return 0;
}

//-----------------------------------------------------------------------------
// codec interface functions.

int
cdc_init(
   char param_str[],
   void **cdc
)
{
   cdc_inst_type *dd = NULL;
   char *str1 = NULL, *token;
   int c_n = 0, c_k = 0, c_m = 0, c_r = 0;

   // Allocate decoder instance.
   dd = (cdc_inst_type *)malloc(sizeof(cdc_inst_type));
   if (dd == NULL) {
      show_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   memset(dd, 0, sizeof(cdc_inst_type));

   // Allocate temporary string for parsing.
   str1 = (char *)malloc(SP_STR_MAX);
   if (str1 == NULL) {
      err_msg("cdc_init: Short of memory!");
      goto ret_err;
   }
   str1[0] = 0;

   // ----- Parse parameters.
   strcpy(str1, param_str);

   token = strtok(str1, tk_seps_prepared);
   while (token != NULL) {
      TRYGET_INT_TOKEN(token, "RM_m", c_m);
      TRYGET_INT_TOKEN(token, "RM_r", c_r);
      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   free(str1); str1 = NULL;

   dd->c_m = c_m;
   dd->c_r = c_r;

   // Check if all necessary parameters are specified.
   if ((c_m == 0) || (c_r == 0)) {
      err_msg("cdc_init: RM_m or RM_r are not specified.");
      goto ret_err;
   }

   // Complete settings.
   dd->c_n = c_n = 1 << c_m;
   dd->c_k = c_k = calc_rm_k(c_m, c_r);

#ifdef FORMAT_Q
   if (dec_tables_init()) {
      err_msg("cdc_init: error in dec_tables_init().");
      goto ret_err;
   }
#endif // FORMAT_Q

   (*cdc) = dd;

   return 0;

ret_err:

   if (str1 != NULL) free(str1);
   if (dd != NULL) {
      free(dd);
      (*cdc) = NULL;
   }
   return 1;
}

void
cdc_close(
   void *cdc
)
{
   cdc_inst_type *dd;

   if (cdc == NULL) return;
   dd = (cdc_inst_type *)cdc;
   free(dd);
#ifdef FORMAT_Q
   dec_tables_free();
#endif // FORMAT_Q
}

int cdc_get_n(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->c_n;
}

int cdc_get_k(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return 0;
   dd = (cdc_inst_type *)cdc;
   return dd->c_k;
}

void
cdc_set_sg(
   void *cdc,
   double noise_sg
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   dd->sg22 = 2.0 / (noise_sg * noise_sg);
}

int
enc_bpsk(
   void *cdc,
   int x[],
   double y[]
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   return mrm_enc_bpsk(dd->c_m, dd->c_r, x, y);
}

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int x_dec[]
)
{
   cdc_inst_type *dd;
   double *dec_input, *y_dec, *aux_buf;
   // int max_peak_lsiz;
   void *mem_buf;
   int i;

   dd = (cdc_inst_type *)cdc;

   mem_buf = (double *)malloc(dd->c_n * 3 * sizeof(double));
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }

   dec_input = (double *)mem_buf;
   y_dec = dec_input + dd->c_n;
   aux_buf = y_dec + dd->c_n;
   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;
   VY_TO_FORMAT(dec_input, dec_input, dd->c_n);

   tri0_dec(dec_input, dd->c_r, dd->c_m, dd->c_n, y_dec, x_dec, &i, aux_buf);

   free(mem_buf);

   return 0;
}
