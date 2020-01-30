//=============================================================================
// Recursive decoding on triangle, global list with permutations.
// Needs channel sigma (define global DEC_NEEDS_SIGMA).
// May use up to mCr permutations.
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

//-----------------------------------------------------------------------------
// Includes.

#include <stdlib.h>
#include <string.h>
#include "../common/typedefs.h"
#include "../common/std_defs.h"
#include "../interfaces/codec.h"
#include "../common/spf_par.h"
#include "../interfaces/ui_utils.h"
#include "../common/srm_utils.h"
#include "../dtrm_glp/dtrm_glp_inner.h"

//-----------------------------------------------------------------------------
// Internal defines.

// Max # of SNR values.
#define VARL_NUM_MAX        50

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
   decoder_type dc;
   int c_n; // Code length.
   int c_k; // Code dimension.
   int c_m; // RM m.
   int c_r; // RM r.
   int csnrn;
   double sg22; // 2 / sg^2 for calculation of epsilons.
   int use_var_list;
   int lsiz[VARL_NUM_MAX];
   double dist_t[VARL_NUM_MAX];
   int dist_t_n;
} cdc_inst_type;


//-----------------------------------------------------------------------------
// Global data.

//-----------------------------------------------------------------------------
// Internal functions.

int get_membuf_size(decoder_type *dd) {
   int flsiz = MAX(dd->peak_lsiz, dd->p_num) * FLSIZ_MULT;
   int mem_buf_size = 0;

   mem_buf_size += dd->c_n * sizeof(ylitem); // y_in
   mem_buf_size += flsiz * sizeof(int); // frind
   mem_buf_size += flsiz * sizeof(int); // lorder
   mem_buf_size += (dd->c_k * flsiz + 1) * sizeof(xlist_item); // xlist
   mem_buf_size += flsiz * sizeof(xlist_item *); // lind2xl
   mem_buf_size += dd->c_m * flsiz * sizeof(ylitem *); // ylist
   mem_buf_size += dd->c_m * sizeof(ylitem *) + (dd->c_n - 1) * MAX(dd->peak_lsiz, dd->p_num) * 4 * sizeof(ylitem); // yitems
   mem_buf_size += flsiz * sizeof(slitem); // slist
   mem_buf_size += flsiz * sizeof(int); // plist.
   mem_buf_size += dd->c_n * sizeof(xlitem); // xtmp.
   mem_buf_size += flsiz * sizeof(int); // parent.

   return mem_buf_size;
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
   int i, j, i1, i2, pn;
   int *p1, *px, *py;
   int p_info[10000], inf_coeff[2048];
   int p_num = 0;

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
      if (strcmp(token, "border_node_mask") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.node_table_len = strlen(token);
         dd->dc.node_table = (uint32 *)malloc(dd->dc.node_table_len * 4);
         if (dd->dc.node_table == NULL) {
            err_msg("cdc_init: Short of memory.");
            goto ret_err;
         }
         for (i = 0; i < dd->dc.node_table_len; i++)
            dd->dc.node_table[i] = (token[i] == '0') ? 0 : 1;
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "border_node_lsize") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.peak_lsiz = atoi(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "permutations") == 0) {
         if (c_m == 0) {
            err_msg("cdc_init(): parameters file: RM_m must be specified before permutations.");
            goto ret_err;
         }
         token = strtok(NULL, tk_seps_prepared);
         if (token[0] == GROUP_CHAR_BEGIN) {
            token = strtok(NULL, tk_seps_prepared);
            i1 = 0;
            p_num = 0;
            while (token[0] != GROUP_CHAR_END) {
               for (i = 0; i < c_m; i++) {
                  p_info[i1++] = atoi(token);
                  token = strtok(NULL, tk_seps_prepared);
                  if (token == NULL) {
                     err_msg("cdc_init(): parameters file: permutations.");
                     goto ret_err;
                  }
               }
               p_num++;
            }
         }
         else {
            err_msg("cdc_init(): parameters file: permutations.");
            goto ret_err;
         }
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      TRYGET_ONOFF_TOKEN(token, "use_variable_list", dd->use_var_list);
      TRYGET_GRINTDOUBLE_TOKEN(token, "distance_threshold",
         dd->lsiz, dd->dist_t, dd->dist_t_n, VARL_NUM_MAX,
         "cdc_init: distance_threshold.");

      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   free(str1); str1 = NULL;

   dd->c_m = dd->dc.c_m = c_m;
   dd->c_r = dd->dc.c_r = c_r;

   // Check if all necessary parameters are specified.
   if ((c_m == 0) || (c_r == 0)) {
      err_msg("cdc_init: RM_m or RM_r are not specified.");
      goto ret_err;
   }
   if (dd->dc.node_table_len == 0) {
      err_msg("cdc_init: No border nodes map.");
      goto ret_err;
   }
   if (dd->dc.peak_lsiz == 0) {
      err_msg("cdc_init: No list size.");
      goto ret_err;
   }
   if ((dd->use_var_list != 0) && (dd->dist_t_n == 0)) {
      err_msg("cdc_init: No distance_threshold.");
      goto ret_err;
   }

   // Complete settings.
   dd->c_n = dd->dc.c_n = c_n = 1 << c_m;
   for (i = 0; i < dd->dc.node_table_len; i++)
      dd->dc.node_table[i] = (dd->dc.node_table[i] == 0) ? 0 : dd->dc.peak_lsiz;
   dd->c_k = dd->dc.c_k = c_k = calc_srm_k(c_m, c_r, dd->dc.node_table);

   // Generate permutations.
   if (p_num == 0) {
      // permutations parameter was not given => only one permutation
      // (0 1 ... m-1) (no permutation).
      for (i = 0; i < c_m; i++) p_info[i] = i;
      p_num = 1;
   }
   dd->dc.p_num = p_num;
   if (p_num > 0) {
      dd->dc.pxarr = (int *)malloc(p_num * c_k * sizeof(int));
      dd->dc.pyarr = (int *)malloc(p_num * c_n * sizeof(int));
      if ((dd->dc.pxarr == NULL) || (dd->dc.pxarr == NULL)) {
         err_msg("cdc_init: Short of memory.");
         goto ret_err;
      }
      px = dd->dc.pxarr;
      py = dd->dc.pyarr;
      p1 = p_info;
   }
   // Calculate inf_coeff[]. It contains all integers with weight <= c_r
   // in reverse lexicographic order.
   i1 = 0;
   for (i = c_n - 1; i >= 0; i--) {
      // i2 <-- weight(i);
      i2 = i & 1;
      for (j = 1; j < c_m; j++) i2 += (i >> j) & 1;
      if (i2 <= c_r) inf_coeff[i1++] = i;
   }
   // Generate permutation arrays.
   for (pn = 0; pn < p_num; pn++) {
      // Next permutation for y.
      for (i = 0; i < c_n; i++) {
         i1 = 0;
         for (j = 0; j < c_m; j++) i1 |= ((i >> j) & 1) << p1[j];
         py[i] = i1;
      }
      // for (i = 0; i < c_n; i++) msg_printf("%d ", py[i]); // DDD
      // msg_printf("\n"); // DDD
      py += c_n;
      // Next permutation for x.
      for (i = 0; i < c_k; i++) {
         i1 = 0;
         i2 = inf_coeff[i];
         for (j = 0; j < c_m; j++) i1 |= ((i2 >> j) & 1) << p1[j];
         for (j = 0; j < c_k; j++) if (inf_coeff[j] == i1) break;
         px[i] = j;
      }
      // for (i = 0; i < c_k; i++) msg_printf("%d ", px[i]); // DDD
      // msg_printf("\n"); // DDD
      px += c_k;
      p1 += c_m;
   }

#ifdef FORMAT_Q
   if (dec_tables_init()) {
      err_msg("cdc_init: error in dec_tables_init().");
      goto ret_err;
   }
#endif // FORMAT_Q

   // Allocate decoder memory buffer
   dd->dc.mem_buf = (uint8 *)malloc(get_membuf_size(&(dd->dc)));
   if (dd->dc.mem_buf == NULL) {
      err_msg("cdc_init: Short of memory.");
      goto ret_err;
   }

   (*cdc) = dd;

   return 0;

ret_err:

   if (str1 != NULL) free(str1);
   if (dd != NULL) {
      if (dd->dc.mem_buf != NULL) free(dd->dc.mem_buf);
      if (dd->dc.pxarr != NULL) free(dd->dc.pxarr);
      if (dd->dc.pyarr != NULL) free(dd->dc.pyarr);
      if (dd->dc.node_table != NULL) free(dd->dc.node_table);
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
   if (dd->dc.mem_buf != NULL) free(dd->dc.mem_buf);
   if (dd->dc.node_table != NULL) free(dd->dc.node_table);
   if (dd->dc.pxarr != NULL) free(dd->dc.pxarr);
   if (dd->dc.pyarr != NULL) free(dd->dc.pyarr);
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

void
cdc_set_csnrn(
   void *cdc,
   int csnrn
)
{
   cdc_inst_type *dd;

   dd = (cdc_inst_type *)cdc;
   dd->csnrn = csnrn;
   if ((dd->use_var_list != 0) && (dd->dist_t_n <= csnrn)) {
      err_msg("cdc_set_csnrn: dist_t_n <= csnrn.");
   }
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
   return smrm_enc_bpsk(dd->c_m, dd->c_r, dd->dc.node_table, x, y);
}

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int x_dec[]
)
{
   cdc_inst_type *dd;
   double *dec_input, *y_dec;
   // int max_peak_lsiz;
   void *mem_buf;
   int i, j;
   double d1;

   dd = (cdc_inst_type *)cdc;

   mem_buf = (double *)malloc(dd->c_n * 2 * sizeof(double));
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return -1;
   }

   dec_input = (double *)mem_buf;
   // dec_input <-- 2 * c_out / sg^2.
   for (i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;

   if (dd->use_var_list == 0)
      rm_dec(&(dd->dc), dec_input, x_dec, NULL);
   else {
      y_dec = ((double *)mem_buf) + dd->c_n;
      // max_peak_lsiz = dd->dc.peak_lsiz;
      for (i = 0; i < dd->dist_t_n; i++) {
         dd->dc.peak_lsiz = dd->lsiz[i];
         for (j = 0; j < dd->dc.node_table_len; j++)
            dd->dc.node_table[j] = (dd->dc.node_table[j] == 0) ? 0 : dd->lsiz[i];
         rm_dec(&(dd->dc), dec_input, x_dec, NULL);
         enc_bpsk(cdc, x_dec, y_dec);
         d1 = 0.0;
         // d1 <-- dist(dec_input, y_dec).
         for (j = 0; j < dd->c_n; j++)
            d1 += (c_out[j] - y_dec[j]) * (c_out[j] - y_dec[j]);
         d1 *= dd->sg22 / 2.0;
         if (d1 < dd->dist_t[i]) break;
      }
      // msg_printf("i = %d, d1 = %g.\n", i, d1);
      // msg_flush();
      // dd->dc.peak_lsiz = max_peak_lsiz;
   }

   free(mem_buf);

   return 0;
}
