//=============================================================================
// Codec interface.
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

#ifndef CODEC_H

#define CODEC_H

//-----------------------------------------------------------------------------
// Defines

// #define DEC_NEEDS_SIGMA
// #define DEC_NEEDS_CSNRN

//-----------------------------------------------------------------------------
// Prototypes.

int
cdc_init(
   char param_str[],
   void **cdc
);

int cdc_get_n(void *cdc);
int cdc_get_k(void *cdc);

#ifdef DEC_NEEDS_CSNRN
void
cdc_set_csnrn(
   void *cdc,
   int csnrn
);
#endif // DEC_NEEDS_CSNRN

#ifdef DEC_NEEDS_SIGMA
void
cdc_set_sg(
   void *cdc,
   double noise_sg
);
#endif // DEC_NEEDS_SIGMA

int
enc_bpsk(
   void *cdc,
   int x[],
   double y[]
);

int
dec_bpsk(
   void *cdc,
   double c_out[],
   int xd[]
);

void
cdc_close(
   void *cdc
);

#endif // #ifndef CODEC_H