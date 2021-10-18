//=============================================================================
// Load simulation results file data to an Octave/Matlab structure variable.
//
// Copyright 2021 and onwards Kirill Shabunov.
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
#include <string.h>
#include "../interfaces/ui_utils.h"
#include "../common/std_defs.h"
#include "../common/spf_par.h"
#include "mex.h"

//-----------------------------------------------------------------------------
// Internal defines.

// Max # of SNR values.
#define SNR_NUM_MAX        50

#define MEX_TRYGET_GRDOUBLE_TOKEN(tk, tk_samp, x, n, nmax, erstr) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      int i = 0; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      if ((tk)[0] == GROUP_CHAR_BEGIN) { \
         (tk) = strtok(NULL, tk_seps_prepared); \
         while (((tk)[0] != GROUP_CHAR_END) && (i < (nmax))) { \
            (x)[i++] = atof(tk); \
            (tk) = strtok(NULL, tk_seps_prepared); \
         } \
      } \
      else { \
         err_msg(erstr); \
      } \
      (n) = i; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

//-----------------------------------------------------------------------------
// Functions.

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  FILE *fp;
  char *str1, *token; // For old res file parsing.
  int c_n = 0, c_k = 0;
  int snr_num = 0;
  double snr_db[SNR_NUM_MAX];
  int trn[SNR_NUM_MAX];
  int ml_trn[SNR_NUM_MAX];
  int en_bit[SNR_NUM_MAX];
  int en_bl[SNR_NUM_MAX];
  int er_n[SNR_NUM_MAX];
  int enml_bl[SNR_NUM_MAX];
  int trn_n = 0;
  int ml_trn_n = 0;
  int en_bit_n = 0;
  int en_bl_n = 0;
  int er_n_n = 0;
  int enml_bl_n = 0;
  int extr_mlb = 0;
  int n;
  const char *keys[] = {"n", "k", "EbNo", "BER", "WER", "ERR", "MLER"};
  char err_reading_msg[] = "Error reading SRF file.";
  mxArray *mxArr, *res;
  double *d;

  // Check arguments.
  if (nrhs < 1) {
    msg_printf("Usage: %s('path/to/file.srf')\n", mexFunctionName());
    err_msg("No input file name.\n");
    return;
  }

  // Get data from the SRF file.
  char *infName = mxArrayToString(prhs[0]);
  if (spf_read_preparse(infName, &str1) != 0) {
    err_msg(err_reading_msg);
    return;
  }
  token = strtok(str1, tk_seps_prepared);
  while (token != NULL) {
    TRYGET_INT_TOKEN(token, code_n_token, c_n)
    TRYGET_INT_TOKEN(token, code_k_token, c_k)
    MEX_TRYGET_GRDOUBLE_TOKEN(token, SNR_token, snr_db, snr_num, SNR_NUM_MAX, err_reading_msg)
    if (tryread_group_int_param(&token, tr_num_token, trn, sizeof(int), &trn_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    if (tryread_group_int_param(&token, en_bit_token, en_bit, sizeof(int), &en_bit_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    if (tryread_group_int_param(&token, en_bl_token, en_bl, sizeof(int), &en_bl_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    if (tryread_group_int_param(&token, er_n_token, er_n, sizeof(int), &er_n_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    if (tryread_group_int_param(&token, ml_tr_num_token, ml_trn, sizeof(int), &ml_trn_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    if (tryread_group_int_param(&token, enml_bl_token, enml_bl, sizeof(int), &enml_bl_n, SNR_NUM_MAX, err_reading_msg) == RC_OK) continue;
    SPF_SKIP_UNKNOWN_PARAMETER(token)
  }
  free(str1);

  // Check validity of the data.
  if ((snr_num != trn_n) || (en_bit_n != trn_n) || (en_bl_n != trn_n) || (er_n_n > 0 && er_n_n != trn_n)
      || (enml_bl_n != trn_n) || (ml_trn_n != trn_n)
      || (c_n == 0) || (c_k == 0)) {
    err_msg("Invalid SRF file.");
    return;
  }

  res = mxCreateStructMatrix(1, 1, sizeof(keys) / sizeof(char *), keys);

  mxSetFieldByNumber(res, 0, 0, mxCreateDoubleScalar(c_n));
  mxSetFieldByNumber(res, 0, 1, mxCreateDoubleScalar(c_k));

  mxArr = mxCreateDoubleMatrix(1, snr_num, mxREAL);
  d = mxGetPr(mxArr);
  for (n = 0; n < snr_num; n++) d[n] = snr_db[n];
  mxSetFieldByNumber(res, 0, 2, mxArr);

  mxArr = mxCreateDoubleMatrix(1, snr_num, mxREAL);
  d = mxGetPr(mxArr);
  for (n = 0; n < snr_num; n++) {
    d[n] = (trn[n] > 0) ? (double)(en_bit[n]) / (trn[n] * c_k) : 0.0;
  }
  mxSetFieldByNumber(res, 0, 3, mxArr);

  mxArr = mxCreateDoubleMatrix(1, snr_num, mxREAL);
  d = mxGetPr(mxArr);
  for (n = 0; n < snr_num; n++) {
    d[n] = (trn[n] > 0) ? (double)(en_bl[n]) / trn[n] : 0.0;
  }
  mxSetFieldByNumber(res, 0, 4, mxArr);

  mxArr = mxCreateDoubleMatrix(1, snr_num, mxREAL);
  d = mxGetPr(mxArr);
  for (n = 0; n < snr_num; n++) {
    d[n] = (trn[n] > 0 && er_n_n > 0) ? (double)(er_n[n]) / trn[n] : 0.0;
  }
  mxSetFieldByNumber(res, 0, 5, mxArr);

  mxArr = mxCreateDoubleMatrix(1, snr_num, mxREAL);
  d = mxGetPr(mxArr);
  for (n = 0; n < snr_num; n++) {
    d[n] = (ml_trn[n] > 0) ? (double)(enml_bl[n]) / ml_trn[n] : 0.0;
  }
  mxSetFieldByNumber(res, 0, 6, mxArr);

  if (nlhs > 0) {
    plhs[0] = res;
  }
}