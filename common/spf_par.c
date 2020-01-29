//=============================================================================
// Simulation parameters file parsing utilities.
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

#define SPF_PAR_C

//-----------------------------------------------------------------------------
// Includes.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spf_par.h"
#include "../interfaces/ui_utils.h"

//-----------------------------------------------------------------------------
// Global data.

char token_seps[] = " ,\t\n\r";
char tk_seps_prepared[] = {TOKEN_SEP_CHAR, 0};

char include_token[] = "include";
char switch_on_token[] = "on";

char code_n_token[] = "code_n";
char code_k_token[] = "code_k";
char code_G_binary_token[] = "code_G_binary";
char code_name_token[] = "code_name";

// Simulation results file keywords.
char SNR_token[] = "SNR";
char tr_num_token[] = "tr_num";
char en_bit_token[] = "en_bit";
char en_bl_token[] = "en_bl";
char ml_tr_num_token[] = "ml_tr_num";
char enml_bl_token[] = "enml_bl";

//char _token[] = "";

//-----------------------------------------------------------------------------
// Functions.

// Copy parameters string from sstr to dstr.
// While copying remove comments.
// Return: resulting string length.
int
strcpy_strip_comments(
   char sstr[], // Source string.
   int sstr_len, // Source string length.
   int write_term_flag, // if 1 - write 0 at the end.
   char dstr[] // destination string.
)
{
   int dstr_len = 0;
   int i2 = 0; // Starting not in comment mode.
   int i;

   if (sstr_len <0) sstr_len = strlen(sstr);

   for (i = 0; i < sstr_len; i++) {
      if (i2 == 0) {
         if (sstr[i] != COMMENT_CHAR) dstr[dstr_len++] = sstr[i];
         else {
            i2 = 1;
         }
      }
      else {
         if (sstr[i] == LF_CHAR) {
            dstr[dstr_len++] = sstr[i];
            i2 = 0;
         }
      }
   }
   if (write_term_flag == 1) dstr[dstr_len] = 0;
   return dstr_len;
}

// Open, read and preparse simulation parameters file to sp_str.
// Return: 0 - OK, !0 - error.
int
spf_read_preparse_ext(
   char spf_name[], // Simulation parameters file name.
   char **out_str, // Simulation parameters string.
   int try_flag // If 1 don't give error message when can't open file.
)
{
   FILE *sp_file; // Simulation parameters file.
   char *sp_str, *str1, *str2, *token;
   int sp_str_len, str2_len, i1;

   // Read simulation parameters file to str1.
   if ((sp_file = fopen(spf_name, "rb")) == NULL) {
      if (try_flag != 1) err_msg("spf_read_preparse(): can't open file.");
      return 1;
   }
   str1 = (char *)malloc(SP_STR_MAX * 3);
   if (str1 == NULL) {
      err_msg("spf_read_preparse(): short of memory.");
      fclose(sp_file);
      return 1;
   }
   str2 = str1 + SP_STR_MAX;
   if ((i1 = fread(str1, 1, SP_STR_MAX, sp_file)) == 0) {
      err_msg("spf_read_preparse(): can't read file.");
      free(str1);
      fclose(sp_file);
      return 1;
   }
   fclose(sp_file);
   if (i1 >= SP_STR_MAX) {
      err_msg("spf_read_preparse(): file is too huge.");
      free(str1);
      return 1;
   }
   str1[i1] = 0;

   // Allocate future sp_str.
   sp_str = (char *)malloc(SP_STR_MAX);
   if (sp_str == NULL) {
      err_msg("spf_read_preparse(): short of memory!");
      free(str1);
      return 1;
   }
   sp_str_len = 0;
   sp_str[0] = 0;

   // Remove comments from str1 and copy the result to str2.
   str2_len = strcpy_strip_comments(str1, i1, 1, str2);

   // ----- Preparse parameters.

   token = strtok(str2, token_seps);
   while (token != NULL) {

      // Includes.
      if (strcmp(token, include_token) == 0) {
         token = strtok(NULL, token_seps);
         if ((sp_file = fopen(token, "rb")) == NULL) {
            err_msg("spf_read_preparse(): can't open include file.");
            free(sp_str);
            free(str1);
            return 1;
         }
         if ((i1 = fread(str1, 1, SP_STR_MAX, sp_file)) == 0) {
            err_msg("spf_read_preparse(): can't read include file.");
            fclose(sp_file);
            free(sp_str);
            free(str1);
            return 1;
         }
         fclose(sp_file);
         str1[i1] = 0;
         i1 = strcpy_strip_comments(str1, i1, 1, str1);
         if (i1 + str2_len >= SP_STR_MAX) {
            err_msg("spf_read_preparse(): file is too huge.");
            free(sp_str);
            free(str1);
            return 1;
         }
         // Insert the file instead of the include command.
         memmove(token + strlen(token) + 1 + i1, token + strlen(token) + 1, str2_len + 1);
         memcpy(token + strlen(token) + 1, str1, i1);
         str2_len += i1;
         //str2_len += strcpy_strip_comments(str1, i1, str2 + str2_len);
         token = strtok(NULL, token_seps);
         continue;
      }
      
      // Group begins.
      if (token[0] == GROUP_CHAR_BEGIN) {
         sp_str[sp_str_len++] = GROUP_CHAR_BEGIN;
         sp_str[sp_str_len++] = TOKEN_SEP_CHAR;
         sp_str[sp_str_len] = 0;
         if (strlen(token) > 1) {
            strcat(sp_str, token + 1);
            sp_str_len += strlen(token) - 1;
            sp_str[sp_str_len++] = TOKEN_SEP_CHAR;
            sp_str[sp_str_len] = 0;
         }
         token = strtok(NULL, token_seps);
         continue;
      }

      // Group ends.
      if (token[strlen(token) - 1] == GROUP_CHAR_END) {
         if (strlen(token) > 1) {
            token[strlen(token) - 1] = 0;
            strcat(sp_str, token);
            sp_str_len += strlen(token);
            sp_str[sp_str_len++] = TOKEN_SEP_CHAR;
            sp_str[sp_str_len] = 0;
         }
         sp_str[sp_str_len++] = GROUP_CHAR_END;
         sp_str[sp_str_len++] = TOKEN_SEP_CHAR;
         sp_str[sp_str_len] = 0;
         token = strtok(NULL, token_seps);
         continue;
      }

      // All other tokens.
      strcat(sp_str, token);
      sp_str_len += strlen(token);
      sp_str[sp_str_len++] = TOKEN_SEP_CHAR;
      sp_str[sp_str_len] = 0;
      
      token = strtok(NULL, token_seps);
   }

   free(str1);
   (*out_str) = sp_str;

   return 0;
}

// Open, read and preparse simulation parameters file to sp_str.
// Return: 0 - OK, !0 - error.
int
spf_read_preparse(
   char spf_name[], // Simulation parameters file name.
   char **out_str // Simulation parameters string.
) {
   return spf_read_preparse_ext(spf_name, out_str, 0);
}

// Try to open, read and preparse simulation parameters file to sp_str.
// Return: 0 - OK, !0 - error.
int
spf_tryread_preparse(
   char spf_name[], // Simulation parameters file name.
   char **out_str // Simulation parameters string.
) {
   return spf_read_preparse_ext(spf_name, out_str, 1);
}

