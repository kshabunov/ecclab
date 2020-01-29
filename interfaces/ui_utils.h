//=============================================================================
// User interface utilities interface.
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

#ifndef UI_UTILS_H

#define UI_UTILS_H

//-----------------------------------------------------------------------------
// Defines

#define MSG_MAX_STR_LEN    2048

//-----------------------------------------------------------------------------
// Prototypes.

void show_msg(char msg[]);

void err_msg(char msg[]);

void msg_printf(
    char *msgf, // Message format string.
    ...
);

void msg_flush(void);

#endif // #ifndef UI_UTILS_H