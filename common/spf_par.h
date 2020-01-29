//=============================================================================
// Simulation parameters file parsing utilities interface.
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

#ifndef SPF_PAR_H

#define SPF_PAR_H

//-----------------------------------------------------------------------------
// Includes.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


//-----------------------------------------------------------------------------
// Defines

#define SP_STR_MAX         140000 // Max length of params file (with includes).
#define LF_CHAR            10    // LF character.
#define COMMENT_CHAR       '%'   // Comment starts from this character.
#define GROUP_CHAR_BEGIN   '{'   // Group starts from this character.
#define GROUP_CHAR_END     '}'   // Group starts from this character.
// This character is used to separate tokens in processed params string.
#define TOKEN_SEP_CHAR     '\1'

#define TRYGET_INT_TOKEN(tk, tk_samp, x) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      (tk) = strtok(NULL, tk_seps_prepared); \
      (x) = atoi(tk); \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define TRYGET_FLOAT_TOKEN(tk, tk_samp, x) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      (tk) = strtok(NULL, tk_seps_prepared); \
      (x) = atof(tk); \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define TRYGET_ONOFF_TOKEN(tk, tk_samp, x) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      (tk) = strtok(NULL, tk_seps_prepared); \
      (x) = (!strcmp(tk, switch_on_token)) ? 1 : 0; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define TRYGET_GRINT_TOKEN(tk, tk_samp, x, n, nmax, erstr) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      int i = 0; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      if ((tk)[0] == GROUP_CHAR_BEGIN) { \
         (tk) = strtok(NULL, tk_seps_prepared); \
         while (((tk)[0] != GROUP_CHAR_END) && (i < (nmax))) { \
            (x)[i++] = atoi(tk); \
            (tk) = strtok(NULL, tk_seps_prepared); \
         } \
      } \
      else { \
         err_msg(erstr); \
         return 1; \
      } \
      (n) = i; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define TRYGET_GRDOUBLE_TOKEN(tk, tk_samp, x, n, nmax, erstr) { \
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
         return 1; \
      } \
      (n) = i; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define TRYGET_GRINTDOUBLE_TOKEN(tk, tk_samp, x, y, n, nmax, erstr) { \
   if (strcmp((tk), (tk_samp)) == 0) { \
      int i = 0; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      if ((tk)[0] == GROUP_CHAR_BEGIN) { \
         (tk) = strtok(NULL, tk_seps_prepared); \
         while (((tk)[0] != GROUP_CHAR_END) && (i < (nmax))) { \
            (x)[i] = atoi(tk); \
            (tk) = strtok(NULL, tk_seps_prepared); \
            (y)[i++] = atof(tk); \
            (tk) = strtok(NULL, tk_seps_prepared); \
         } \
      } \
      else { \
         err_msg(erstr); \
         return 1; \
      } \
      (n) = i; \
      (tk) = strtok(NULL, tk_seps_prepared); \
      continue; \
   } \
}

#define SPF_SKIP_UNKNOWN_PARAMETER(tk) { \
   (tk) = strtok(NULL, tk_seps_prepared); \
   if ((tk)[0] == GROUP_CHAR_BEGIN) { \
      (tk) = strtok(NULL, tk_seps_prepared); \
      while (((tk)[strlen(tk) - 1] != GROUP_CHAR_END) && ((tk) != NULL)) \
         (tk) = strtok(NULL, tk_seps_prepared); \
   } \
   (tk) = strtok(NULL, tk_seps_prepared); \
}


//-----------------------------------------------------------------------------
// Typedefs.

//-----------------------------------------------------------------------------
// Global data.

#ifndef SPF_PAR_C

extern char token_seps[];
extern char tk_seps_prepared[];

extern char include_token[];
extern char switch_on_token[];

extern char code_n_token[];
extern char code_k_token[];
extern char code_G_binary_token[];
extern char code_name_token[];

// Simulation results file keywords.
extern char SNR_token[];
extern char tr_num_token[];
extern char en_bit_token[];
extern char en_bl_token[];
extern char ml_tr_num_token[];
extern char enml_bl_token[];

#endif // #ifndef SPF_PAR_C

//-----------------------------------------------------------------------------
// Prototypes.

// Copy parameters string from sstr to dstr without comments.
int
strcpy_strip_comments(
   char sstr[], // Source string.
   int sstr_len, // Source string length.
   int write_term_flag, // if 1 - write 0 at the end.
   char dstr[] // destination string.
);

// Open, read and preparse simulation parameters file to sp_str.
// Return: 0 - OK, !0 - error.
int
spf_read_preparse(
   char spf_name[], // Simulation parameters file name.
   char **out_str // Simulation parameters string.
);

// Try to open, read and preparse simulation parameters file to sp_str.
// Return: 0 - OK, !0 - error.
int
spf_tryread_preparse(
   char spf_name[], // Simulation parameters file name.
   char **out_str // Simulation parameters string.
);

#endif // #ifndef SPF_PAR_H