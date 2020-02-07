//=============================================================================
// Header file for rm1_dec.c
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

#ifndef RM1_DEC_H

#define RM1_DEC_H

//-----------------------------------------------------------------------------
// Prototypes.

int rm1_dec_sh(
   int c_m,
   double y_in[], // in y format.
   int x_dec[]
);

int rm1_dec_lst(
   int c_m,
   double y_in[], // in y format.
   int lsiz, // size of the list (currently <= 4096).
   int x_dec[]
);

// Same as rm1_dec_sh() except the input is considered to be already
// converted to the internal format.
int rm1_dec_sh_finp(
   int c_m,
   double y_in[],
   int x_dec[]
);

// Max dot product.
int rm1_dec_sh_dotp(
   int c_m,
   double y_in[],
   int x_dec[]
);

#endif // #ifndef RM1_DEC_H