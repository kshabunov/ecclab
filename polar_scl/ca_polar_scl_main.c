//=============================================================================
// CRC-aided (CA) Polar SCL codec.
// Implements codec.h
// Needs channel sigma (define global DEC_NEEDS_SIGMA).
//
// Copyright 2019 and onwards Kirill Shabunov.
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
#include "../common/crc.h"
#include "../common/srm_utils.h"
#include "polar_scl_inner.h"

//-----------------------------------------------------------------------------
// Internal defines.

// #define DBG

// Max # of SNR values.
#define VARL_NUM_MAX        50

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
   decoder_type dc;
   int c_n; // Code length.
   int c_k; // Code dimension.
   int eff_k; // Effective code dimension (without CRC).
   int c_m; // Polar m.
   int *crc;
   int crc_len;
   int csnrn;
   double sg22; // 2 / sg^2 for calculation of epsilons.
   int lsiz[VARL_NUM_MAX];
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

int get_k_from_node_table(
  int c_n,
  uint32 *node_table
) {
  int i, k = 0;
  for (i = 0; i < c_n; i++) {
    k += node_table[i];
  }
  return k;
}

int get_best_in_list(
  int *x,
  int x_len,
  int lsiz,
  int *crc,
  int crc_len
) {
  if (crc_len <= 0) {
    return 0;
  }
  else {
    int sum_len = x_len + crc_len;
    int i;
    int *tmp = (int *)malloc(x_len * sizeof(int));
    if (tmp == NULL) {
      return 0;
    }
    for (i = 0; i < lsiz; i++) {
      if (is_valid_crc(x + i * sum_len, x_len, crc, crc_len, tmp)) {
        free(tmp);
        return i;
      }
    }
    free(tmp);
    return -1; // no candidates with valid CRC found
  }
}

#ifdef DBG
void print_bin_vect(char rem[], int *x, int x_len) {
  int i;
  printf("%s: ", rem);
  for (i = 0; i < x_len; i++) {
    printf("%d", x[i]);
  }
  printf("\n");
}
#endif

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
   int c_n = 0, c_k = 0, c_m = 0;
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

      TRYGET_INT_TOKEN(token, "c_m", c_m);
      if (strcmp(token, "info_bits_mask") == 0) {
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
      if (strcmp(token, "ca_polar_crc") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->crc_len = strlen(token);
         dd->crc = (int *)malloc(dd->crc_len * sizeof(int));
         if (dd->crc == NULL) {
            err_msg("cdc_init: Short of memory.");
            goto ret_err;
         }
         for (i = 0; i < dd->crc_len; i++)
            dd->crc[i] = (token[i] == '0') ? 0 : 1;
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "list_size") == 0) {
         token = strtok(NULL, tk_seps_prepared);
         dd->dc.peak_lsiz = atoi(token);
         token = strtok(NULL, tk_seps_prepared);
         continue;
      }
      if (strcmp(token, "permutations") == 0) {
         // TODO: Permutations that were used for RM codes are not usable with Polar, but let's leave it here for now.
         if (c_m == 0) {
            err_msg("cdc_init(): parameters file: c_m must be specified before permutations.");
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

      SPF_SKIP_UNKNOWN_PARAMETER(token);
   }

   free(str1); str1 = NULL;

   dd->c_m = dd->dc.c_m = c_m;

   // Check if all necessary parameters are specified.
   if (c_m == 0) {
      err_msg("cdc_init: c_m is not specified.");
      goto ret_err;
   }
   if (dd->dc.node_table_len == 0) {
      err_msg("cdc_init: No info_bits_mask.");
      goto ret_err;
   }
   if (dd->dc.peak_lsiz == 0) {
      err_msg("cdc_init: No list size.");
      goto ret_err;
   }

   // Complete settings.
   dd->c_n = dd->dc.c_n = c_n = 1 << c_m;
   dd->c_k = dd->dc.c_k = c_k = get_k_from_node_table(c_n, dd->dc.node_table);
   dd->eff_k = c_k - dd->crc_len;
   dd->dc.ret_list = (dd->crc_len > 0) ? 1 : 0;
   for (i = 0; i < dd->dc.node_table_len; i++)
      dd->dc.node_table[i] = (dd->dc.node_table[i] == 0) ? 0 : dd->dc.peak_lsiz;

   // Generate permutations.
   // TODO: Permutations.
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
      if ((dd->dc.pxarr == NULL) || (dd->dc.pyarr == NULL)) {
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
      if (i2 <= c_m) inf_coeff[i1++] = i;
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

   return RC_OK;

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
   return RC_ERROR;
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
   if (dd->crc != NULL) free(dd->crc);
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

   if (cdc == NULL) return RC_ERROR;
   dd = (cdc_inst_type *)cdc;
   return dd->c_n;
}

int cdc_get_k(void *cdc) 
{
   cdc_inst_type *dd;

   if (cdc == NULL) return RC_ERROR;
   dd = (cdc_inst_type *)cdc;
   return dd->eff_k;
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
   int *x_crc;

   dd = (cdc_inst_type *)cdc;

   if (dd->crc_len > 0) {
      x_crc = (int *)malloc(dd->c_k * sizeof(int));
      if (x_crc == NULL) {
         err_msg("enc_bpsk: Short of memory.");
         return RC_ERROR;
      }
      calc_crc(x, dd->eff_k, dd->crc, dd->crc_len, x_crc + dd->crc_len);
      memcpy(x_crc, x, dd->eff_k * sizeof(int));
      smrm_par0_enc_bpsk(dd->c_m, dd->c_m, dd->dc.node_table, x_crc, y);
      free(x_crc);
      return dd->eff_k;
   }
   else {
      return smrm_par0_enc_bpsk(dd->c_m, dd->c_m, dd->dc.node_table, x, y);
   }
}

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int x_dec[]
)
{
   cdc_inst_type *dd;
   double *dec_input;
   int *x_candidates;
   uint8 *mem_buf, *mem_buf_ptr;
   int mem_buf_size;
   int best_ind;
   int is_erasure = 0;

   dd = (cdc_inst_type *)cdc;

   mem_buf_size = dd->c_n * sizeof(double); // dec_input
   if (dd->dc.ret_list) {
      mem_buf_size += dd->dc.peak_lsiz * dd->c_k * sizeof(int); // x_candidates
   }
   mem_buf = (uint8 *)malloc(mem_buf_size);
   if (mem_buf == NULL) {
      err_msg("dec_bpsk: Short of memory.");
      return RC_ERROR;
   }
   mem_buf_ptr = mem_buf;

   dec_input = (double *)mem_buf_ptr;
   mem_buf_ptr += dd->c_n * sizeof(double);

   if (dd->dc.ret_list) {
      x_candidates = (int *)mem_buf_ptr;
      mem_buf_ptr += dd->dc.peak_lsiz * dd->c_k * sizeof(int);
   }

   // dec_input <-- 2 * c_out / sg^2.
   for (int i = 0; i < dd->c_n; i++) dec_input[i] = c_out[i] * dd->sg22;

   if (dd->dc.ret_list) {
      polar_dec(&(dd->dc), dec_input, x_candidates, NULL);
      best_ind = get_best_in_list(x_candidates, dd->eff_k, dd->dc.peak_lsiz, dd->crc, dd->crc_len);
#ifdef DBG
      printf("best_ind: %d\n", best_ind);
#endif
      if (best_ind < 0) {
        is_erasure = 1;
        best_ind = 0;
      }
      memcpy(x_dec, x_candidates + best_ind * dd->c_k, dd->eff_k * sizeof(int));
   }
   else {
      polar_dec(&(dd->dc), dec_input, x_dec, NULL);
   }

   free(mem_buf);

   return is_erasure ? RC_DEC_ERASURE : RC_OK;
}
