//=============================================================================
// Simulation interface.
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

#ifndef SIMUL_H

#define SIMUL_H

//-----------------------------------------------------------------------------
// Configuration switches.

//-----------------------------------------------------------------------------
// Includes.

//-----------------------------------------------------------------------------
// Defines

// For sim_control().
#define SIM_CTRL_CUR_MORE     1
#define SIM_CTRL_CUR_LESS     2
#define SIM_CTRL_CUR_NEXT     3

//-----------------------------------------------------------------------------
// Typedefs.

typedef struct {
   char *spf_name; // Simulation parameters file name.
   int ret_int; // Maximal interval to return from sim_run() in sec.
   int dont_randomize; // If !0 -- Don't randomize random numbers generator.
} sim_init_params;

//-----------------------------------------------------------------------------
// Global data.

#ifndef SIMUL_C


#endif // #ifndef SIMUL_C

//-----------------------------------------------------------------------------
// Prototypes.

// Allocate and init a simulation instance.
// Return: 0 - OK, !0 - error.
int sim_init(
   sim_init_params *sp, // Simulation parameters.
   void **inst // Simulation instance.
);

// Delete simulation instance.
void sim_close(
   void *inst // Simulation instance.
);

// Main simulation cycle routine.
// Return: 0 - completed, >0 - not completed, <0 - error.
int sim_run(
   void *inst // Simulation instance.
);

// Return short simulation status string.
void sim_state_str(
   void *inst, // Simulation instance.
   char ds[] // Destination string.
);

// Return current simulation results string.
void sim_cur_res_str(
   void *inst, // Simulation instance.
   char ds[] // Destination string.
);

// Save simulation results.
int sim_save_res(
   void *inst // Simulation instance.
);

// Control simulation parameters "on the fly".
void sim_control(
   void *inst, // Simulation instance.
   int ctrl_code
);

#endif // #ifndef SIMUL_H