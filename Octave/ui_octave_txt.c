//=============================================================================
// Text mode user interface utilities implementation for Octave.
//
// Copyright 2021 and onwards Kirill Shabunov
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

#define UI_TXT_C

//-----------------------------------------------------------------------------
// Includes.

#include <stdarg.h>
#include <stdio.h>
#include "mex.h"
#include "../interfaces/ui_utils.h"

//-----------------------------------------------------------------------------
// Functions.

void show_msg(char msg[]) {
  mexPrintf("%s\n", msg);
}

void err_msg(char msg[]) {
  mexErrMsgTxt(msg);
}

void msg_printf(
    char *msgf, // Message format string.
    ...
)
{
  char str_buf[MSG_MAX_STR_LEN];
  va_list args;

  va_start(args, msgf);
  vsprintf(str_buf, msgf, args);
  va_end(args);

  mexPrintf("%s", str_buf);
}

void msg_flush(void) {
  return;
}
