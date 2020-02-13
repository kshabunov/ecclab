//=============================================================================
// Polar SCL.
// Can return the list or the best candidate.
// Internal format can be changed.
// Derived from dtrm_glp
//
// Copyright 2019 and onwards Kirill Shabunov.
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


//-----------------------------------------------------------------------------
// Configuration switches.

//#define DBG

//-----------------------------------------------------------------------------
// Includes.

#ifdef DBG
#include <stdio.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../common/typedefs.h"
#include "../common//std_defs.h"
#include "polar_scl_inner.h"


//-----------------------------------------------------------------------------
// Internal defines.

#define NEAR_ZERO (1e-300)

// List access defines.
#define GETX(i) (lind2xl[i]->x)
#define PUSHX(v, i) { next_xle_ptr->x = v; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }
#define PUSHZEROX(i) { next_xle_ptr->x = 0; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }
#define PUSHONEX(i) { next_xle_ptr->x = 1; next_xle_ptr->p = lind2xl[i]; lind2xl[i] = next_xle_ptr++; }

#ifdef DBG2
#define YLIST_POP(i) (printf("pop i: %d, yitems[i]: %p, yfrindp[i]: %d\n", i, yitems[i], yfrindp[i]), yitems[i] + (1 << (i)) * (yfrindp[i]++));
#else // DBG2
#define YLIST_POP(i) (yitems[i] + (1 << (i)) * (yfrindp[i]++))
#endif // DBG2
#define YLIST_RESET(i) { yfrindp[i] = 0; }

#define YLISTP(i, j) ylist[(i) * dd->c_m + (j)]
#define YLISTPP(i, j) (ylist + (i) * dd->c_m + (j))


//-----------------------------------------------------------------------------
// Internal typedefs.

//-----------------------------------------------------------------------------
// Global data.

uint32 *node_table; // The table of list sizes at the border nodes.
int node_counter;

slitem *slist;
xlist_item *xlist;
int *plist; // Permutation index list.
int *parent; // index of the parent list element.
xlist_item **lind2xl; // list index to xlist

ylitem **yitems; // peak_lsiz * 2 * c_n
ylitem **ylist; // current buffer - YLISTP(list_index, m)
int yfrindp[32]; // yfrindp[i] points to the next available row for ylist in yitems[i].

int *lorder; // list ordering
int cur_lsiz; // Current list size.

int *frind; // Stack containing indexes of free list cells.
int frindp; // frind pointer. Points to the next available element.

xlist_item *next_xle_ptr; // pointer to the next available xlist element

xlitem *xtmp;


//-----------------------------------------------------------------------------
// Functions.

void vgetx(int listInd, xlitem *res_x, int n) {
  int i;
  xlist_item *xle = lind2xl[listInd];
  for (i = 0; i < n; i++) {
    res_x[n - i - 1] = xle->x;
    xle = xle->p;
  }
}

int branch(int i0, xlist_item *xle0) {
  int i1 = frind[frindp++];
  parent[i1] = i0;
  lind2xl[i1] = xle0;
  lorder[cur_lsiz++] = i1;
  slist[i1] = slist[i0];
  plist[i1] = plist[i0];
  return i1;
}

#ifdef DBG
int vgetx_all(int listInd, xlitem *res_x) {
  int i, n;
  xlitem x;
  xlist_item *xle = lind2xl[listInd];

  i = 0;
  while (xle->p != NULL) {
    res_x[i++] = xle->x;
    xle = xle->p;
  }
  n = i;
  for (i = 0; i < n / 2; i++) {
    x = res_x[i];
    res_x[i] = res_x[n - i - 1];
    res_x[n - i - 1] = x;
  }
  return n;
}
void print_list(char rem[]) {
  int i, j, n;
  xlitem tmp[4096];
  if (rem != NULL) {
    printf("%s: ", rem);
  }
  for (i = 0; i < cur_lsiz; i++) {
    printf("%3.1f:", slist[lorder[i]]);
    n = vgetx_all(lorder[i], tmp);
    for (j = 0; j < n; j++) {
      printf("%1d", tmp[j]);
    }
    printf(" ");
  }
  printf("\n");
#ifdef DBG2
  printf("Y buffers: ");
  for (i = 0; i < 4; i++) {
    printf("%d:%p:%d, ", i, yitems[i], yfrindp[i]);
  }
  printf("\n");
#endif // DBG2
}
#endif // DBG

#define SWAP(a, b) { tmp = v[a]; v[a] = v[b]; v[b] = tmp; }

void qpartition(int *v, int len, int k) {
  int i, st, tmp;

  for (st = i = 0; i < len - 1; i++) {
    if (slist[v[i]] < slist[v[len - 1]]) continue;
    SWAP(i, st);
    st++;
  }

  SWAP(len - 1, st);

  if (k == st) return;
  if (st > k) qpartition(v, st, k);
  else qpartition(v + st, len - st, k - st);
}

// For qsort().
static int lst_comp(const void *i, const void *j)
{
   if (slist[*(int *)i] < slist[*(int *)j]) return 1;
   return -1;
}

void polar0_branch(
   decoder_type *dd, // Decoder instance data.
   int peak_lsiz
) {
  int cur_lsiz_old = cur_lsiz;
  int cur_ind0, cur_ind1;
  ylitem y;
  xlitem x;
  ylitem *yp;
  int i;

#ifdef DBG
  printf("rm00 branch.\n");
  print_list("a");
#endif

  for (i = 0; i < cur_lsiz_old; i++) {

    cur_ind0 = lorder[i];
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0]);

    yp = YLISTP(cur_ind0, 0);
    if (CHECK_EST0(yp[0])) {
      y = yp[0];
      x = 0;
    }
    else {
      y = INV_EST(yp[0]);
      x = 1;
    }

    PUSHX(x, cur_ind0);
    slist[cur_ind0] += EST0_TO_LNP0(y);

    PUSHX(1 - x, cur_ind1);
    slist[cur_ind1] += EST0_TO_LNP1(y);
  }

#ifdef DBG
  print_list("b");
#endif

  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) {
    qpartition(lorder, cur_lsiz, peak_lsiz);
    while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
  }

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, 0), YLISTPP(cur_ind0, 0), dd->c_m * sizeof(ylitem *));
    }
    x = GETX(cur_ind1);
    yp = YLISTP(cur_ind1, 0) = YLIST_POP(0);
    yp[0] = x ? YLDEC1 : YLDEC0;
  }

#ifdef DBG
  print_list("c");
#endif

}

void polar0_skip(
  decoder_type *dd // Decoder instance data.
) {
  int cur_ind0;
  ylitem *yp;
  int i;

#ifdef DBG
  printf("polar0 skip.\n");
  print_list("a");
#endif

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    yp = YLISTP(cur_ind0, 0);
    slist[cur_ind0] += EST_TO_LNP0(yp[0]);
    yp = YLISTP(cur_ind0, 0) = YLIST_POP(0);
    yp[0] = YLDEC0;
  }

#ifdef DBG
  print_list("b");
#endif
}

void polar_dec_inner(
  decoder_type *dd, // Decoder instance data.
  int m
) {

  int n = 1 << m;
  int n2 = n / 2;
  int cur_ind0;
  ylitem *yp1, *yp2, *vp, *up;
  ylitem y1;
  int i, j;

  if (m == 0) {
    if (node_table[node_counter] == 0) {
      polar0_skip(dd);
    }
    else {
      polar0_branch(dd, node_table[node_counter]);
    }
    node_counter++;
    return;
  }

#ifdef DBG2
  printf("innr a: m=%d, n=%3d.\n", m, n);
  print_list("innr a");
#endif

  // Calculate y_v = y_1 xor y_2.
  for (i = 0; i < cur_lsiz; i++) {
     cur_ind0 = lorder[i];
     yp1 = YLISTP(cur_ind0, m);
     yp2 = yp1 + n2;
     vp = YLISTP(cur_ind0, m - 1) = YLIST_POP(m - 1);
     VXOR_EST(yp1, yp2, vp, n2);
  }

  polar_dec_inner(dd, m - 1);

#ifdef DBG2
  printf("innr b: m=%d, n=%3d.\n", m, n);
  print_list("innr b");
#endif

  // Calculate y_u = y_1 xor v + y_2.
  // Also save ref to v in y.
  for (i = 0; i < cur_lsiz; i++) {
     cur_ind0 = lorder[i];
     yp1 = YLISTP(cur_ind0, m);
     yp2 = yp1 + n2;
     vp = YLISTP(cur_ind0, m) = YLISTP(cur_ind0, m - 1); // y <-- v
     up = YLISTP(cur_ind0, m - 1) = YLIST_POP(m - 1);
     for (j = 0; j < n2; j++) {
        y1 = EST_XOR_YLDEC(yp1[j], vp[j]); // y1 <-- y_1[j] xor v[j].
        ADD_EST(y1, yp2[j], up[j]); // y_u <-- y1 + y_2[j];
     }
  }

#ifdef DBG2
  printf("innr c: m=%d, n=%3d.\n", m, n);
#endif

  polar_dec_inner(dd, m - 1);

#ifdef DBG2
  printf("innr d: m=%d, n=%3d.\n", m, n);
#endif

  // y_dec <-- (u xor v | u).
  for (i = 0; i < cur_lsiz; i++) {
     cur_ind0 = lorder[i];
     vp = YLISTP(cur_ind0, m);
     yp1 = YLISTP(cur_ind0, m) = YLIST_POP(m);
     yp2 = yp1 + n2;
     up = YLISTP(cur_ind0, m - 1);
     for (j = 0; j < n2; j++) {
        yp1[j] = XOR_YLDEC(vp[j], up[j]);
        yp2[j] = up[j];
     }
  }

  YLIST_RESET(m - 1);
}

// Decoding procedure.
int
polar_dec(
   decoder_type *dd, // Decoder instance data.
   double *y_input, // Decoder input.
   int *x_dec, // Decoded information sequence.
   double *s_dec // Metric of x_dec (set NULL, if not needed).
)
{
  int i, j, i1;
  int c_n = dd->c_n;
  int c_k = dd->c_k;
  int c_m = dd->c_m;
  int peak_lsiz; // Peak list size.
  int flsiz; // Size of allocated list.
  ylitem *y_in;
  int n2 = dd->c_n / 2; // Half of the current length (n).
  int cur_ind0;
  ylitem *vp, *up;
  int *p1; // For permutations.
  slitem s1;
  ylitem y1;
  int *xp;
  uint8 *mem_buf_ptr;
  int mem_buf_size;

  // ----- Shortcuts assignment.

  node_table = dd->node_table;
  peak_lsiz = dd->peak_lsiz;
  flsiz = MAX(dd->peak_lsiz, dd->p_num) * FLSIZ_MULT;

  // ----- Memory allocation.

  mem_buf_ptr = dd->mem_buf;

  y_in = (ylitem *)mem_buf_ptr;
  mem_buf_ptr += c_n * sizeof(ylitem);

  frind = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  lorder = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  xlist = (xlist_item *)mem_buf_ptr;
  mem_buf_ptr += (c_k * flsiz + 1) * sizeof(xlist_item);

  lind2xl = (xlist_item **)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(xlist_item *);

  ylist = (ylitem **)mem_buf_ptr;
  //mem_buf_ptr += dd->c_m * MAX(dd->peak_lsiz, dd->p_num) * sizeof(ylitem *);
  mem_buf_ptr += dd->c_m * flsiz * sizeof(ylitem *);

  yitems = (ylitem **)mem_buf_ptr;
  mem_buf_ptr += dd->c_m * sizeof(ylitem *);
#ifdef DBG2
  printf("\nStart buffers: %ld\n", mem_buf_ptr - dd->mem_buf);
#endif
  for (i = 0; i < dd->c_m; i++) {
    yitems[i] = (ylitem *)mem_buf_ptr;
    mem_buf_ptr += (1 << i) * MAX(dd->peak_lsiz, dd->p_num) * 4 * sizeof(ylitem);
#ifdef DBG2
    printf("Buffer %d: yitems[i]=%p, %ld\n", i, yitems[i], mem_buf_ptr - dd->mem_buf);
#endif
  }

  slist = (slitem *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(slitem);

  plist = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  parent = (int *)mem_buf_ptr;
  mem_buf_ptr += flsiz * sizeof(int);

  xtmp = (xlitem *)mem_buf_ptr;
  mem_buf_ptr += c_n * sizeof(xlitem);

#ifdef DBG2
  printf("\nAssigned buffer size: %ld\n", mem_buf_ptr - dd->mem_buf);
#endif

  // ----- Initial assignments.

  VY_TO_FORMAT(y_input, y_in, c_n);

  for (i = 0; i < flsiz; i++) frind[i] = i;
  for (i = 0; i < 32; i++) YLIST_RESET(i);

  node_counter = 0;

  //---------- Root node.

  // Calculate y_v = y_1 xor y_2.
  // Branch with permutations at (c_m, c_m).
  // List size is going to be dd->p_num.
  // All path metrics are 0 so far.
  for (i = 0; i < dd->p_num; i++) {
    slist[i] = 0.0;
    lorder[i] = i;
    lind2xl[i] = xlist;
    p1 = dd->pyarr + c_n * i;
    vp = YLISTP(i, c_m - 1) = YLIST_POP(c_m - 1);
    // TODO: Permutations are not usable with general Polar, but let's leave it as it is for now.
    for (j = 0; j < n2; j++) vp[j] = XOR_EST(y_in[p1[j]], y_in[p1[j + n2]]);
    plist[i] = i;
  }
  cur_lsiz = frindp = dd->p_num;
  next_xle_ptr = xlist + 1;
  xlist[0].x = 0;
  xlist[0].p = NULL;

#ifdef DBG2
  printf("root 1: m=%d, n=%3d.\n", c_m, c_n);
  print_list("root a");
#endif

  // Decode v.
  polar_dec_inner(dd, c_m - 1);

#ifdef DBG2
  print_list("root b");
#endif

  // Calculate y_u = y_1 xor v + y_2.
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    vp = YLISTP(cur_ind0, c_m - 1);
    up = YLISTP(cur_ind0, c_m - 1) = YLIST_POP(c_m - 1);
    // TODO: Permutations.
    p1 = dd->pyarr + c_n * plist[cur_ind0];
    for (j = 0; j < n2; j++) {
      y1 = EST_XOR_YLDEC(y_in[p1[j]], vp[j]); // y1 <-- y_1[j] xor v[j].
      ADD_EST(y1, y_in[p1[j + n2]], up[j]); // y_u <-- y1 + y_2[j];
    }
  }

  // Decode u.
  polar_dec_inner(dd, c_m - 1);

#ifdef DBG
  print_list("final");
#endif

  if (dd->ret_list) {

    // Sort.
    qsort(lorder, cur_lsiz, sizeof(int), lst_comp);
#ifdef DBG
    print_list("sorted");
#endif
    for (i = 0; i < cur_lsiz; i++) {
      vgetx(lorder[i], xtmp, c_k);
      // TODO: Permutations.
      p1 = dd->pxarr + dd->c_k * plist[lorder[i]];
      xp = x_dec + i * c_k;
      for (j = 0; j < c_k; j++) xp[p1[j]] = xtmp[j];
    }

  }
  else {

    // Find the best.
    i1 = 0;
    for (i = 1; i < cur_lsiz; i++) {
      if (slist[lorder[i]] > slist[lorder[i1]]) {
        i1 = i;
      }
    }

#ifdef DBG
    printf("Best: %3.1f", slist[lorder[i1]]);
#endif

    vgetx(lorder[i1], xtmp, c_k);
    // TODO: Permutations.
    p1 = dd->pxarr + dd->c_k * plist[lorder[i1]];
    for (j = 0; j < c_k; j++) x_dec[p1[j]] = xtmp[j];

    if (s_dec != NULL) {
      (*s_dec) = slist[lorder[i1]];
    }
  }

  return 0;
}
