//=============================================================================
// Operations counting stuff.
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

#ifndef COUNT_OPS_H
#define COUNT_OPS_H 

#ifdef COUNT_OPS
// 0 -- float, float.
// 1 -- float, const.
// 2 -- float compare.
// 3 -- float inverse.
// 4 -- float * (+/-1).
// 5 -- binary, binary.
#define ADD_OPS(i, n) gl_ops[i] += (n);
#else
#define ADD_OPS(i, n)
#endif

#ifdef COUNT_OPS
extern int gl_ops[];
#endif


#endif // #ifndef COUNT_OPS_H