//=============================================================================
// Main function to run simulation in text mode.
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

#define MAIN_TXT_C

//-----------------------------------------------------------------------------
// Includes.

#ifdef WIN32
#include <conio.h>
#include <windows.h>
#endif // if WIN32
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../common/std_defs.h"
#include "../interfaces/ui_utils.h"
#include "../interfaces/simul.h"

//-----------------------------------------------------------------------------
// Internal defines.

#define DEFAULT_RET_INT    10
#define DEFAULT_SAVE_INT   600

//-----------------------------------------------------------------------------
// Functions.

int main(int argc, char **argv) {

   sim_init_params sp;
   void *sim_inst; // Simulation instance.
   char msg_str[10000]; // String buffer for messages.
   time_t start_time;
   int save_int = DEFAULT_SAVE_INT;
   int i, rc;
#ifdef WIN32
   char ch1;
   HANDLE hMain;
#endif // WIN32

   // Check arguments.
   if (argc < 2) {
      msg_printf("Usage : %s <simulation_parameters_file> [options]\n", argv[0]);
      msg_printf("   -nr - don't randomize.\n");
      msg_printf("   -si <int> - saving interval in sec. (default %d).\n", DEFAULT_SAVE_INT);
      msg_printf("   -ri <int> - return interval in sec. (default %d).\n", DEFAULT_RET_INT);
      return RC_ERROR;
   }

   memset(&sp, 0, sizeof(sim_init_params));
   sp.spf_name = argv[1];
   sp.ret_int = DEFAULT_RET_INT;

   for (i = 2; i < argc; i++) {
      if (strcmp(argv[i], "-nr") == 0) sp.dont_randomize = 1;
      if (strcmp(argv[i], "-si") == 0) save_int = atoi(argv[++i]);
      if (strcmp(argv[i], "-ri") == 0) sp.ret_int = atoi(argv[++i]);
   }

#ifdef WIN32
   // Set low priority so when we run a long simulation
   // it will less affect system performance.
   hMain = GetCurrentProcess();
   SetPriorityClass(hMain, IDLE_PRIORITY_CLASS);
#endif // if WIN32

   // Init simulation instance.
   if (sim_init(&sp, &sim_inst) != RC_OK) {
      err_msg("Can't initialize simulation instance.");
      return RC_ERROR;
   }

   show_msg("Simulation start.");

   time(&start_time);

   // Call simulation procedure.
   while ((rc = sim_run(sim_inst)) == RC_SIMUL_NOT_COMPLETED) {
#ifdef WIN32
      if (_kbhit()) {
         ch1 = _getch();
         switch (ch1) {
         case 's' : sim_save_res(sim_inst); break;
         case ' ' :
            msg_printf("\n");
            sim_cur_res_str(sim_inst, msg_str);
            show_msg(msg_str);
            break;
         case '>' : sim_control(sim_inst, SIM_CTRL_CUR_MORE); break;
         case '<' : sim_control(sim_inst, SIM_CTRL_CUR_LESS); break;
         case 'N' : sim_control(sim_inst, SIM_CTRL_CUR_NEXT); break;
         }
      }
#endif // WIN32
      sim_state_str(sim_inst, msg_str);
      if (strlen(msg_str) < 75) {
         for (i = strlen(msg_str); i < 75; i++) msg_str[i] = ' ';
         msg_str[75] = 0;
      }
      msg_printf("\r%s", msg_str);
      msg_flush();
      // Check if it's time to save intermediate results.
      if (time(NULL) - start_time > save_int) {
         sim_save_res(sim_inst);
         time(&start_time);
      }
   }
   msg_printf("\n");

   show_msg("Simulation is over.");
   sim_cur_res_str(sim_inst, msg_str);
   show_msg(msg_str);

   sim_save_res(sim_inst);

   sim_close(sim_inst);

   return RC_OK;
}