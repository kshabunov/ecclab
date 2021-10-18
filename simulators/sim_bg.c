//=============================================================================
// Simulation in the binary Gaussian channel
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

#define SIM_BG_C

//-----------------------------------------------------------------------------
// Configuration switches.

// Use the stdlib rand() function as the uniform random numbers generator.
// #define USE_STDLIB_RND

//-----------------------------------------------------------------------------
// Includes.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/std_defs.h"
#include "../common/spf_par.h"
#include "../interfaces/simul.h"
#include "../interfaces/codec.h"
#include "../interfaces/ui_utils.h"

//-----------------------------------------------------------------------------
// Internal defines.

// Max # of SNR values.
#define SNR_NUM_MAX        50

#define MAX_TRIALS_PER_SNR 1e10

// Minimum number of observed block errors to estimate required number of trials.
#define MIN_EN_FOR_TRN_ESTIMATION 20

// Max length of file names.
#define FN_LEN_MAX         1000

// Defines for pseudorandom numbers generator.
#ifdef USE_STDLIB_RND
#define RND ((double)rand() / RAND_MAX)
#else // USE_STDLIB_RND
// #define rnd_A     16807UL
// #define rnd_M     2147483647UL
#define rnd_A     16807.0
#define rnd_M     2147483647.0 // (2^31 - 1) - a prime number.
#define rnd_D     4.656612875e-10
// #define RND ((double)(rnd_seed = (rnd_A * rnd_seed) % rnd_M) * rnd_D)
#define RND ((rnd_seed = fmod(rnd_A * rnd_seed, rnd_M)) * rnd_D)
#endif // else USE_STDLIB_RND

//-----------------------------------------------------------------------------
// Internal typedefs.

typedef struct {
  int trn;
  int en_bit;
  int en_bl;
  int er_n;
  int trn_ml;
  int enml_bl;
} sim_point;

typedef struct {
  void *dc_inst; // Codec instance.
  char rf_name[FN_LEN_MAX]; // Results file name.
  int ret_int; // Return interval.
  int code_n;
  int code_k;
  int use_rndcw; // If 1 use random codeword.
  double fixedR; // If set, simulate as if the code has this rate (error rate vs 1/R * Es/N0).
  int snr_num;
  int csnrn;
  double snr_db[SNR_NUM_MAX];
  int min_trn;
  int min_en_bl;
  int trn_req[SNR_NUM_MAX];
  int do_ml; // If 1 evaluate ML LB.
  int do_ml_hd; // If 1 evaluate ML LB for hard dec. decoder.
  sim_point pt[SNR_NUM_MAX];
  sim_point pt_saved[SNR_NUM_MAX];
} sim_bg_inst;

//-----------------------------------------------------------------------------
// Global data.

// For pseudorandom numbers generator.
// unsigned long rnd_seed = 1;
double rnd_seed = 1.0;

// Simulation parameters file keywords.
char res_file_token[] = "res_file";
char snr_val_trn_token[] = "SNR_val_trn";
char ml_lb_token[] = "ml_lb";
char random_codeword_token[] = "random_codeword";
char fixedR_token[] = "fixed_R";

//-----------------------------------------------------------------------------
// Functions.

static double db2val(double x) {
  return exp(log(10.0) * x / 10.0);
}

static void set_rnd_seed(long s) {
#ifdef USE_STDLIB_RND
  srand((unsigned)s);
#else
  rnd_seed = (double) s;
#endif
}

/* Pauses for a specified number of milliseconds. */
static void sleep(clock_t wait) {
  clock_t goal;
  goal = wait + clock();
  while (goal > clock());
}

// Copy and add Gaussian noise using Marsaglia polar method.
static void copy_add_noise(
  double const ys[],
  int n, // Length of y[] (must be even!).
  double sg, // Standard deviation (sigma).
  double yd[]
) {
  double v1, v2, r;
  int i = 0;

  while (i < n) {
    do {
      v1 = 2.0 * RND - 1.0;
      v2 = 2.0 * RND - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);
    r = sqrt((-2.0 * log(r)) / r);
    yd[i] = ys[i] + v1 * r * sg;
    i++;
    yd[i] = ys[i] + v2 * r * sg;
    i++;
  }
}

void init_sim_point(sim_point *p) {
  memset(p, 0, sizeof(sim_point));
}

void copy_sim_point(sim_point *dst, sim_point *src) {
  memcpy(dst, src, sizeof(sim_point));
}

void add_sim_points(sim_point *dst, sim_point *src1, sim_point *src2) {
  dst->trn = src1->trn + src2->trn;
  dst->en_bit = src1->en_bit + src2->en_bit;
  dst->en_bl = src1->en_bl + src2->en_bl;
  dst->er_n = src1->er_n + src2->er_n;
  dst->trn_ml = src1->trn_ml + src2->trn_ml;
  dst->enml_bl = src1->enml_bl + src2->enml_bl;
}

void sub_sim_points(sim_point *dst, sim_point *src1, sim_point *src2) {
  dst->trn = src1->trn - src2->trn;
  dst->en_bit = src1->en_bit - src2->en_bit;
  dst->en_bl = src1->en_bl - src2->en_bl;
  dst->er_n = src1->er_n - src2->er_n;
  dst->trn_ml = src1->trn_ml - src2->trn_ml;
  dst->enml_bl = src1->enml_bl - src2->enml_bl;
}

// Allocate and init a simulation instance.
int sim_init(
  sim_init_params *sp, // Simulation parameters.
  void **inst // Simulation instance.
) {
  char *sp_str; // Simulation parameters string.
  sim_bg_inst *sim; // Simulator instance.
  // Temporary variables.
  char *token, *str1;

  // Read and preparse simulation parameters.
  if (spf_read_preparse(sp->spf_name, &sp_str) != RC_OK) return RC_ERROR;

  // Allocate simulator instance.
  sim = (sim_bg_inst *)malloc(sizeof(sim_bg_inst));
  if (sim == NULL) {
    err_msg("sim_bg init error: short of memory!");
    return RC_ERROR;
  }
  memset(sim, 0, sizeof(sim_bg_inst));

  // Allocate temporary string for parsing.
  str1 = (char *)malloc(SP_STR_MAX);
  if (str1 == NULL) {
    err_msg("sim_bg init error: short of memory!");
    return RC_ERROR;
  }
  str1[0] = 0;

  // Parse simulation parameters.
  strcpy(str1, sp_str);
  token = strtok(str1, tk_seps_prepared);
  while (token != NULL) {

    if (strcmp(token, res_file_token) == 0) {
      token = strtok(NULL, tk_seps_prepared);
      if (strlen(token) >= FN_LEN_MAX) {
        err_msg("sim_bg init error: output file name is too long.");
        return RC_ERROR;
      }
      strcpy(sim->rf_name, token);
      token = strtok(NULL, tk_seps_prepared);
      continue;
    }
    if (sim->snr_num < 1) {
      if (strcmp(token, snr_val_trn_token) == 0) {
        int i = 0;
        token = strtok(NULL, tk_seps_prepared);
        if (token[0] == GROUP_CHAR_BEGIN) {
          token = strtok(NULL, tk_seps_prepared);
          while ((token[0] != GROUP_CHAR_END) && (i < SNR_NUM_MAX)) {
            sim->snr_db[i] = atof(token);
            token = strtok(NULL, tk_seps_prepared);
            sim->trn_req[i++] = atoi(token);
            token = strtok(NULL, tk_seps_prepared);
          }
        }
        else {
          err_msg("sim_bg init error: cannot read SNR_val_trn.");
          return RC_ERROR;
        }
        sim->snr_num = i;
        token = strtok(NULL, tk_seps_prepared);
        continue;
      }
    }
    if (sim->snr_num < 1) {
      TRYGET_GRDOUBLE_TOKEN(token, "EbNo_values", sim->snr_db, sim->snr_num, SNR_NUM_MAX, "sim_bg init error: too many Eb/No values!");
    }
    TRYGET_INT_TOKEN(token, "min_trials_per_snr", sim->min_trn);
    TRYGET_INT_TOKEN(token, "min_errors_per_snr", sim->min_en_bl);
    TRYGET_ONOFF_TOKEN(token, ml_lb_token, sim->do_ml);
    TRYGET_ONOFF_TOKEN(token, random_codeword_token, sim->use_rndcw);
    TRYGET_ONOFF_TOKEN(token, "ml_lb_hard", sim->do_ml_hd);
    TRYGET_FLOAT_TOKEN(token, fixedR_token, sim->fixedR);

    SPF_SKIP_UNKNOWN_PARAMETER(token);
  }

  // Init codec.
  strcpy(str1, sp_str);
  if (cdc_init(str1, &(sim->dc_inst))) {
    err_msg("sim_bg init error: cannot init codec.");
    return RC_ERROR;
  }

  sim->code_n = cdc_get_n(sim->dc_inst);
  sim->code_k = cdc_get_k(sim->dc_inst);

  // Check if all necessary parameters are specified.
  if (sim->snr_num < 1) {
    err_msg("sim_bg init error: no SNR values are specified.");
    return RC_ERROR;
  }
  if (sim->min_trn <= 0) {
    sim->min_trn = 1;
  }
  if (sim->min_en_bl <= 0) {
    sim->min_en_bl = 1;
  }
  if ((sim->code_n == 0) || (sim->code_k == 0)) {
    err_msg("sim_bg init error: code_n or code_k are not specified.");
    return RC_ERROR;
  }

  // Assign other settings and init counters in the simulation instance.
  for (int i = 0; i < sim->snr_num; i++) {
    if (sim->trn_req[i] < 1) {
      sim->trn_req[i] = sim->min_trn;
    }
  }
  sim->csnrn = 0;
  sim->ret_int = sp->ret_int;
  if (sp->dont_randomize == 0) set_rnd_seed((int)time(NULL));
  if (sim->do_ml_hd) sim->do_ml = 1;

  free(sp_str);
  free(str1);

  (*inst) = sim;

  return RC_OK;
}

// Delete simulation instance.
void sim_close(
  void *inst // Simulation instance.
) {
  sim_bg_inst *sim;

  sim = (sim_bg_inst *)inst;
  cdc_close(sim->dc_inst);
  free(sim);
}

static void update_trn_req(
  sim_bg_inst *sim,
  int csnrn
) {
  if (sim->pt[csnrn].en_bl < sim->min_en_bl) {
    double estimate;
    if (sim->pt[csnrn].en_bl > MIN_EN_FOR_TRN_ESTIMATION) {
      estimate = (double)sim->pt[csnrn].trn * sim->min_en_bl / sim->pt[csnrn].en_bl;
      sim->trn_req[csnrn] = (int)MAX(sim->min_trn, MIN(estimate + 0.5, MAX_TRIALS_PER_SNR));
      return;
    }
    if (sim->pt[csnrn].trn > sim->trn_req[csnrn] / 2) {
      estimate = sim->pt[csnrn].trn * 1.5;
      sim->trn_req[csnrn] = (int)MIN(estimate + 0.5, MAX_TRIALS_PER_SNR);
    }
  }
}

// Main simulation cycle routine.
// Return: 0 - completed, >0 - not completed, <0 - error.
int sim_run(
  void *inst // Simulation instance.
) {
  sim_bg_inst *sim;
  int c_n; // Code length.
  int c_k; // Code dimension.
  double c_R; // Code rate.
  int csnrn; // Current SNR value number.
  // int trn; // (Current) trials number.
  double noise_sg; // Current value of sigma.
  int *x; // Information vector.
  int *x_dec; // Decoded information vector.
  double *c_in; // Channel input.
  double *c_out; // Channel output.
  time_t start_time;
  int en = 0;
  int i;

  sim = (sim_bg_inst *)inst;

  c_n = sim->code_n;
  c_k = sim->code_k;
  c_R = sim->fixedR ? sim->fixedR : ((double)c_k / (double)c_n);
  csnrn = sim->csnrn; // Current SNR value number.

  time(&start_time); // Remember start time.

  // Allocate buffers.
  x = (int *)malloc(c_n * sizeof(int));
  x_dec = (int *)malloc(c_n * sizeof(int));
  c_in = (double *)malloc(c_n * sizeof(double));
  c_out = (double *)malloc(c_n * sizeof(double));
  if ((x == NULL) || (x_dec == NULL) || (c_in == NULL) || (c_out == NULL)) {
    err_msg("sim_bg run error: short of memory.");
    return RC_ERROR;
  }

  // Main simulation loop.
  while (csnrn < sim->snr_num) {

    noise_sg = 1 / sqrt(2 * c_R * db2val(sim->snr_db[csnrn]));

#ifdef DEC_NEEDS_CSNRN
    // Set current SNR # in the codec.
    cdc_set_csnrn(sim->dc_inst, csnrn);
#endif // DEC_NEEDS_CSNRN

#ifdef DEC_NEEDS_SIGMA
    // Set sg in the codec.
    cdc_set_sg(sim->dc_inst, noise_sg);
#endif // DEC_NEEDS_SIGMA

    while (sim->pt[csnrn].trn < sim->trn_req[csnrn]) {

      sim->pt[csnrn].trn++;

      // Generate random or zero code word.
      if (sim->use_rndcw) {
        for (i = 0; i < c_k; i++) x[i] = (RND > 0.5) ? 1 : 0;
        enc_bpsk(sim->dc_inst, x, c_in);
      }
      else {
        for (i = 0; i < c_k; i++) x[i] = 0;
        for (i = 0; i < c_n; i++) c_in[i] = -1.0;
      }

      // Generate channel output for the simulated codeword.
      // (Channel simulation.)
      // c_out <-- c_in + noise.
      copy_add_noise(c_in, c_n, noise_sg, c_out);

      // Decode.
      int rc = dec_bpsk(sim->dc_inst, c_out, x_dec);
      if (rc == RC_ERROR) {
        err_msg("sim_bg run error: error while decoding.");
        return RC_ERROR;
      }

      // Count the number of incorrect information bits.
      en = 0;
      for (i = 0; i < c_k; i++) if (x_dec[i] != x[i]) en++;
      // Adjust error counters.
      if (en) {
        sim->pt[csnrn].en_bit += en;
        sim->pt[csnrn].en_bl++;
      }

      // Count erasures.
      if (rc == RC_DEC_ERASURE) {
        sim->pt[csnrn].er_n++;
      }

      // Check if the ML decoding would fail too.
      if (sim->do_ml && rc != RC_DEC_ERASURE) {
        sim->pt[csnrn].trn_ml++;
        double d1 = 0.0, d2 = 0.0;
        if (sim->do_ml_hd)
          for (i = 0; i < c_n; i++) c_out[i] = (c_out[i] > 0.0) ? 1.0 : -1.0;
        // d1 <-- dist(c_in, c_out).
        for (i = 0; i < c_n; i++) d1 += (c_in[i] - c_out[i]) * (c_in[i] - c_out[i]);
        // d2 <-- dist[encode(x_dec), c_out].
        enc_bpsk(sim->dc_inst, x_dec, c_in);
        for (i = 0; i < c_n; i++) d2 += (c_in[i] - c_out[i]) * (c_in[i] - c_out[i]);
        if (d1 > d2) sim->pt[csnrn].enml_bl++;
      }

      update_trn_req(sim, csnrn);

      // Check elapsed time.
      if (time(NULL) - start_time > sim->ret_int) {
        free(c_out);
        free(c_in);
        free(x_dec);
        free(x);
        return RC_SIMUL_NOT_COMPLETED;
      }
    }

    sim->csnrn = ++csnrn;
  }

  free(c_out);
  free(c_in);
  free(x_dec);
  free(x);

  return RC_OK;
}

// Return short simulation status string.
void
sim_state_str(
  void *inst, // Simulation instance.
  char ds[] // Destination string.
) {
  sim_bg_inst *sim;
  int n;

  sim = (sim_bg_inst *)inst;
  n = sim->csnrn;
  if (sim->do_ml) {
    sprintf(ds,
      "SNR: %3.2f (%d / %d), trn: %d / %d, en: %d, ML en: %d.",
      sim->snr_db[n], n + 1, sim->snr_num,
      sim->pt[n].trn, sim->trn_req[n],
      sim->pt[n].en_bl, sim->pt[n].enml_bl
    );
  }
  else {
    sprintf(ds,
      "SNR: %3.2f (%d / %d), trn: %d / %d, en: %d.",
      sim->snr_db[n], n + 1, sim->snr_num,
      sim->pt[n].trn, sim->trn_req[n],
      sim->pt[n].en_bl
    );
  }
}

// Return current simulation results string.
void
sim_cur_res_str(
  void *inst, // Simulation instance.
  char ds[] // Destination string.
) {
  sim_bg_inst *sim;
  int n, snrs_to_show;
  char str1[1024];

  sim = (sim_bg_inst *)inst;
  snrs_to_show = (sim->csnrn < sim->snr_num) ? sim->csnrn + 1 : sim->snr_num;
  if (sim->do_ml) {
    sprintf(ds, "SNR\t ep_bit\t\t ep_bl\t\t epml_bl");
    for (n = 0; n < snrs_to_show; n++) {
      sprintf(str1,
        "\n%3.2f\t %.3e\t %.3e\t %.3e",
        sim->snr_db[n],
        (double)(sim->pt[n].en_bit) / (sim->pt[n].trn * sim->code_k),
        (double)(sim->pt[n].en_bl) / sim->pt[n].trn,
        sim->pt[n].trn_ml ? (double)(sim->pt[n].enml_bl) / sim->pt[n].trn_ml : 0
      );
      strcat(ds, str1);
    }
  }
  else {
    sprintf(ds, "SNR\t ep_bit\t\t ep_bl");
    for (n = 0; n < snrs_to_show; n++) {
      sprintf(str1,
        "\n%3.2f\t %.3e\t %.3e",
        sim->snr_db[n],
        (double)(sim->pt[n].en_bit) / (sim->pt[n].trn * sim->code_k),
        (double)(sim->pt[n].en_bl) / sim->pt[n].trn
      );
      strcat(ds, str1);
    }
  }
}

// Save simulation results.
int
sim_save_res(
  void *inst // Simulation instance.
) {
  sim_bg_inst *sim;
  FILE *fp;
  char bsyfn[FN_LEN_MAX]; // Busy flag file name.
  char *str1, *token; // For old res file parsing.
  int snr_num = 0;
  double snr_db[SNR_NUM_MAX];
  sim_point pt_file[SNR_NUM_MAX];
  int trn_n = 0;
  int ml_trn_n = 0;
  int en_bit_n = 0;
  int en_bl_n = 0;
  int er_n_n = 0;
  int enml_bl_n = 0;
  char err_reading_str[] = "sim_bg saving error: error reading old res file.";
  int i, j, n, snrs_to_save;

  sim = (sim_bg_inst *)inst;

  // Check if the file is busy.
  strcpy(bsyfn, sim->rf_name);
  strcat(bsyfn, ".bsy");
  n = 0;
  while ((fp = fopen(bsyfn, "r")) != NULL) {
    fclose(fp);
    if (n++ > 30) {
      err_msg("sim_bg saving error: output file is busy too long.");
      return RC_ERROR;
    }
    sleep(1000);
  }

  // Create the busy flag.
  if ((fp = fopen(bsyfn, "w")) == NULL) {
    err_msg("sim_bg saving error: cannot create busy flag file.");
    return RC_ERROR;
  }
  fprintf(fp, "Busy");
  fclose(fp);

  // Try to read and preparse old res file.
  if (spf_tryread_preparse(sim->rf_name, &str1) == 0) {

    // Get old data.
    token = strtok(str1, tk_seps_prepared);
    while (token != NULL) {
      int rc;
      // TRYGET_INT_TOKEN(token, code_n_token, code_n);
      // TRYGET_INT_TOKEN(token, code_k_token, code_k);
      TRYGET_GRDOUBLE_TOKEN(token, SNR_token, snr_db, snr_num, SNR_NUM_MAX, err_reading_str);
      rc = tryread_group_int_param(
        &token, tr_num_token, &pt_file[0].trn, sizeof(sim_point), &trn_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      rc = tryread_group_int_param(
        &token, en_bit_token, &pt_file[0].en_bit, sizeof(sim_point), &en_bit_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      rc = tryread_group_int_param(
        &token, en_bl_token, &pt_file[0].en_bl, sizeof(sim_point), &en_bl_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      rc = tryread_group_int_param(
        &token, er_n_token, &pt_file[0].er_n, sizeof(sim_point), &er_n_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      rc = tryread_group_int_param(
        &token, ml_tr_num_token, &pt_file[0].trn_ml, sizeof(sim_point), &ml_trn_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      rc = tryread_group_int_param(
        &token, enml_bl_token, &pt_file[0].enml_bl, sizeof(sim_point), &enml_bl_n, SNR_NUM_MAX, err_reading_str
      );
      if (rc == RC_OK) continue;
      if (rc == RC_ERROR) return RC_ERROR;
      SPF_SKIP_UNKNOWN_PARAMETER(token);
    }
    free(str1);

    // Check validity of the data.
    if ((snr_num != trn_n) || (en_bit_n != trn_n) || (en_bl_n != trn_n) || (er_n_n != trn_n)
        || (enml_bl_n != trn_n) || (ml_trn_n != trn_n)) {
      err_msg("sim_bg saving error: invalid old res file.");
      return RC_ERROR;
    }

    // Update the data.
    for (i = 0; (i <= sim->csnrn) && (i < sim->snr_num); i++) {
      if (sim->pt[i].trn == sim->pt_saved[i].trn) continue;
      n = 0;
      while ((snr_db[n] < sim->snr_db[i]) && (n < snr_num)) n++;
      if ((snr_db[n] > sim->snr_db[i]) || (n >= snr_num)) {
        // Make room for the new at n if required.
        for (j = snr_num; j > n; j--) {
          snr_db[j] = snr_db[j - 1];
          copy_sim_point(&pt_file[j], &pt_file[j - 1]);
        }
        snr_num++;
        // Add new at n.
        snr_db[n] = sim->snr_db[i];
        init_sim_point(&pt_file[n]);
      }
      // Update at n.
      sim_point tmppt;
      sub_sim_points(&tmppt, &sim->pt[i], &sim->pt_saved[i]);
      add_sim_points(&pt_file[n], &pt_file[n], &tmppt);
    }
    snrs_to_save = snr_num;
  }
  else {

    // Assume there was no old res file.
    snrs_to_save = (sim->csnrn < sim->snr_num) ? sim->csnrn + 1 : sim->snr_num;
    for (n = 0; n < snrs_to_save; n++) {
      snr_db[n] = sim->snr_db[n];
      copy_sim_point(&pt_file[n], &sim->pt[n]);
    }
  }
  for (i = 0; (i <= sim->csnrn) && (i < sim->snr_num); i++) {
    copy_sim_point(&sim->pt_saved[i], &sim->pt[i]);
  }

  // Save updated data.

  fp = fopen(sim->rf_name, "wt");
  if (fp == NULL) {
    err_msg("sim_bg saving error: cannot open file to save results.");
    return RC_ERROR;
  }

  fprintf(fp, "code_n %d\ncode_k %d\n\n", sim->code_n, sim->code_k);

  fprintf(fp, "SNR { ");
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%g ", snr_db[n]);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", tr_num_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].trn);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", en_bit_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].en_bit);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", en_bl_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].en_bl);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", er_n_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].er_n);
  fprintf(fp, "}\n");

  fprintf(fp, "WER { ");
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%.3e ", (double)(pt_file[n].en_bl) / pt_file[n].trn);
  fprintf(fp, "}\n");

  fprintf(fp, "ERR { ");
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%.3e ", (double)(pt_file[n].er_n) / pt_file[n].trn);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", ml_tr_num_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].trn_ml);
  fprintf(fp, "}\n");

  fprintf(fp, "%s { ", enml_bl_token);
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%d ", pt_file[n].enml_bl);
  fprintf(fp, "}\n");

  fprintf(fp, "MLER { ");
  for (n = 0; n < snrs_to_save; n++) fprintf(fp, "%.3e ", pt_file[n].trn_ml ? (double)pt_file[n].enml_bl / pt_file[n].trn_ml : 0.0);
  fprintf(fp, "}\n");

  fprintf(fp, "\n");
  fprintf(fp, "%% SNR   BER         WER         ERR         ML LB");
  for (n = 0; n < snrs_to_save; n++) {
    fprintf(fp,
      "\n%% %3.2f  %.3e  %.3e  %.3e  %.3e",
      snr_db[n],
      (double)pt_file[n].en_bit / ((double)pt_file[n].trn * sim->code_k),
      (double)pt_file[n].en_bl / pt_file[n].trn,
      (double)pt_file[n].er_n / pt_file[n].trn,
      pt_file[n].trn_ml ? (double)pt_file[n].enml_bl / pt_file[n].trn_ml : 0.0
    );
  }

  // Delete busy flag file.
  remove(bsyfn);

  fclose(fp);
  return RC_OK;
}

// Control simulation parameters "on the fly".
void sim_control(
  void *inst, // Simulation instance.
  int ctrl_code
) {
  sim_bg_inst *sim;
  int csnrn; // Current SNR value number.

  sim = (sim_bg_inst *)inst;
  csnrn = sim->csnrn;

  switch (ctrl_code) {
    case SIM_CTRL_CUR_MORE :
      sim->trn_req[csnrn] += sim->trn_req[csnrn] / 2;
      break;
    case SIM_CTRL_CUR_LESS :
      sim->trn_req[csnrn] -= sim->trn_req[csnrn] / 3;
      if (sim->trn_req[csnrn] < sim->pt[csnrn].trn) sim->trn_req[csnrn] = sim->pt[csnrn].trn;
      break;
    case SIM_CTRL_CUR_NEXT :
      sim->trn_req[csnrn] = sim->pt[csnrn].trn;
  }
}