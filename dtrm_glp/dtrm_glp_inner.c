//=============================================================================
// RM.
// Modified G.
// Global list decoding + permutations.
// (0, m)x2, (m, m)x4.
// (m, m) add to path likelihoods.
// Optimized branching.
// List as a tree.
// Internal format can be changed.
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
#include "../common/std_defs.h"
#include "dtrm_glp_inner.h"


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

void vsetx(int listInd, xlitem *src_x, int n) {
  int i;
  xlist_item *xle = lind2xl[listInd];
  for (i = 0; i < n; i++) {
    xle->x = src_x[n - i - 1];
    xle = xle->p;
  }
}

void vpushx(int listInd, xlitem *src_x, int n) {
  int i;
  xlist_item *xle = lind2xl[listInd];
  for (i = 0; i < n; i++) {
    next_xle_ptr->x = src_x[i];
    next_xle_ptr->p = xle;
    xle = next_xle_ptr++;
  }
  lind2xl[listInd] = xle;
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

// n <= 2^15.
void code_mm(int n, xlitem *x, xlitem *y) {
  int cur_depth = 0;
  int states[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int i, n2;

  while (cur_depth >= 0) {
    if (n == 1) {
      y[0] = x[0];
      n = 2;
      cur_depth--;
      continue;
    }
    if (states[cur_depth] == 0) {
      // Calculate y_v.
      states[cur_depth++] = 1;
      n >>= 1;
      continue;
    }
    n2 = n >> 1;
    if (states[cur_depth] == 1) {
      // Calculate y_u.
      states[cur_depth++] = 2;
      x += n2;
      y += n2;
      n >>= 1;
    }
    else {
      // (states[cur_depth] == 2)
      x -= n2;
      y -= n2;
      for (i = 0; i < n2; i++) y[i] ^= y[i + n2];
      states[cur_depth--] = 0;
      n <<= 1;
    }
  } // while (cur_depth >= 0)
}

void rm1_branch(
   decoder_type *dd, // Decoder instance data.
   int m,
   int peak_lsiz
) {
  int n = 1 << m;
  int n2 = n / 2;
  int cur_lsiz_old = cur_lsiz;
  int cur_ind0, cur_ind1;
  slitem s0, s1;
  ylitem y1;
  ylitem *yp1, *yp2, *vp, *up;
  int i, j;

#ifdef DBG
  printf("rm1:    m=%d, r=%d, n=%3d, cur_lsiz=%d.\n", m, 1, n, cur_lsiz);
  print_list("a");
#endif

  // Initial branching
  for (i = 0; i < cur_lsiz_old; i++) {
    cur_ind0 = lorder[i];
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0]);
    yp1 = YLISTP(cur_ind0, m);
    yp2 = yp1 + n2;
    s0 = 0.0;
    s1 = 0.0;
    for (j = 0; j < n2; j++) {
      y1 = XOR_EST(yp1[j], yp2[j]);
      if (CHECK_EST0(y1)) {
        s0 += EST0_TO_LNP0(y1);
        s1 += EST0_TO_LNP1(y1);
      }
      else {
        s0 += EST1_TO_LNP0(y1);
        s1 += EST1_TO_LNP1(y1);
      }
    }
    if (s0 > s1) {
      PUSHZEROX(cur_ind0);
      PUSHONEX(cur_ind1);
      slist[cur_ind0] += s0; // Likelihood corresp. to 0.
      slist[cur_ind1] += s1; // Likelihood corresp. to 1.
    }
    else {
      PUSHONEX(cur_ind0);
      PUSHZEROX(cur_ind1);
      slist[cur_ind0] += s1;
      slist[cur_ind1] += s0;
    }
  }

#ifdef DBG
  print_list("b");
#endif

  // Partition and cut the list.
  if (cur_lsiz > peak_lsiz) {
    qpartition(lorder, cur_lsiz, peak_lsiz);
    while (cur_lsiz > peak_lsiz) frind[--frindp] = lorder[--cur_lsiz];
  }

  // Finalize the remaining branches
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    cur_ind0 = parent[cur_ind1];
    if (cur_ind0 >= 0) {
      memcpy(YLISTPP(cur_ind1, 0), YLISTPP(cur_ind0, 0), dd->c_m * sizeof(ylitem *));
    }
  }
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind1 = lorder[i];
    yp1 = YLISTP(cur_ind1, m);
    yp2 = yp1 + n2;
    vp = YLISTP(cur_ind1, m) = YLIST_POP(m - 1);
    up = YLISTP(cur_ind1, m - 1) = YLIST_POP(m - 1);
    if (GETX(cur_ind1) == 0) {
      VADD_EST(yp1, yp2, up, n2);
      for (j = 0; j < n2; j++) vp[j] = YLDEC0;
    }
    else {
      VADD_INV_EST(yp1, yp2, up, n2);
      for (j = 0; j < n2; j++) vp[j] = YLDEC1;
    }
  }

#ifdef DBG
  print_list("c");
#endif

}

void rm1_skip(
  decoder_type *dd, // Decoder instance data.
  int m
) {

  int n = 1 << m;
  int n2 = n / 2;
  int cur_ind0;
  slitem s1;
  ylitem y1;
  ylitem *yp1, *yp2, *vp, *up;
  int i, j;

  YLIST_RESET(m - 1);
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    yp1 = YLISTP(cur_ind0, m);
    yp2 = yp1 + n2;
    vp = YLISTP(cur_ind0, m) = YLIST_POP(m - 1);
    up = YLISTP(cur_ind0, m - 1) = YLIST_POP(m - 1);
    s1 = 0.0;
    for (j = 0; j < n2; j++) {
      y1 = XOR_EST(yp1[j], yp2[j]);
      s1 += EST_TO_LNP0(y1);
    }
    slist[cur_ind0] += s1;
    VADD_EST(yp1, yp2, up, n2);
    for (j = 0; j < n2; j++) vp[j] = YLDEC0;
  }
}

void rm11_branch(
   decoder_type *dd, // Decoder instance data.
   int peak_lsiz
) {
  int cur_lsiz_old = cur_lsiz;
  int cur_ind0, cur_ind1, cur_ind2, cur_ind3;
  ylitem y0, y1;
  ylitem *yp1, *yp3;
  xlitem xtmp2[2];
  slitem s00, s01, s10, s11;
  int i, j;

#ifdef DBG
  printf("rm11:   m=%d, r=%d, n=%3d.\n", 1, 1, 2);
  print_list("a");
#endif

  for (i = 0; i < cur_lsiz_old; i++) {

    cur_ind0 = lorder[i];
    parent[cur_ind0] = -1;
    cur_ind1 = branch(cur_ind0, lind2xl[cur_ind0]);
    cur_ind2 = branch(cur_ind0, lind2xl[cur_ind0]);
    cur_ind3 = branch(cur_ind0, lind2xl[cur_ind0]);

    yp1 = YLISTP(cur_ind0, 1);
    if (CHECK_EST0(yp1[0])) {
      y0 = yp1[0];
      xtmp[0] = 0;
    }
    else {
      y0 = INV_EST(yp1[0]);
      xtmp[0] = 1;
    }
    if (CHECK_EST0(yp1[1])) {
      y1 = yp1[1];
      xtmp[1] = 0;
    }
    else {
      y1 = INV_EST(yp1[1]);
      xtmp[1] = 1;
    }

    s00 = EST0_TO_LNP0(y0);
    s01 = EST0_TO_LNP1(y0);
    s10 = EST0_TO_LNP0(y1);
    s11 = EST0_TO_LNP1(y1);

    PUSHX(xtmp[0], cur_ind0);
    PUSHX(xtmp[1], cur_ind0);
    slist[cur_ind0] += s00 + s10;

    PUSHX(1 - xtmp[0], cur_ind1);
    PUSHX(xtmp[1], cur_ind1);
    slist[cur_ind1] += s01 + s10;

    PUSHX(xtmp[0], cur_ind2);
    PUSHX(1 - xtmp[1], cur_ind2);
    slist[cur_ind2] += s00 + s11;

    PUSHX(1 - xtmp[0], cur_ind3);
    PUSHX(1 - xtmp[1], cur_ind3);
    slist[cur_ind3] += s01 + s11;
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
    vgetx(cur_ind1, xtmp, 2);
    code_mm(2, xtmp, xtmp2);
    vsetx(cur_ind1, xtmp2, 2);
    yp3 = YLISTP(cur_ind1, 1) = YLIST_POP(1);
    yp3[0] = xtmp[0] ? YLDEC1 : YLDEC0;
    yp3[1] = xtmp[1] ? YLDEC1 : YLDEC0;
  }

#ifdef DBG
  print_list("c");
#endif

}

void rmm_branch(
   decoder_type *dd, // Decoder instance data.
   int m,
   int peak_lsiz
) {
  int n = 1 << m;
  int cur_lsiz_old = cur_lsiz;
  int cur_ind0, cur_ind1;
  xlist_item *cur_xl0;
  ylitem y1, ym, ym2, ym3;
  xlitem xtmp2[n];
  ylitem *yp1, *yp3;
  slitem s1, s2, s3;
  int i, j, i1, i2, i3;

#ifdef DBG
  printf("rmm:    m=%d, r=%d, n=%3d.\n", m, m, n);
  print_list("a");
#endif

  // Branch (m, m).
  cur_lsiz_old = cur_lsiz;

  for (i = 0; i < cur_lsiz_old; i++) {

    cur_ind0 = lorder[i];
    cur_xl0 = lind2xl[cur_ind0];
    parent[cur_ind0] = -1;
    yp1 = YLISTP(cur_ind0, m);
    i1 = 0; i2 = 1; i3 = 2;
    ym = ym2 = ym3 = STRONGEST_EST0;
    s1 = 0.0;
    for (j = 0; j < n; j++) {
      y1 = yp1[j];
      if (CHECK_EST0(y1)) {
        if (STRONGER_EST0(ym, y1)) {
          ym3 = ym2; i3 = i2;
          ym2 = ym; i2 = i1;
          ym = y1; i1 = j;
        }
        else {
          if (STRONGER_EST0(ym2, y1)) {
            ym3 = ym2; i3 = i2;
            ym2 = y1; i2 = j;
          }
          else {
            if (STRONGER_EST0(ym3, y1)) {
             ym3 = y1; i3 = j;
            }
          }
        }
        s1 += EST0_TO_LNP0(y1);
        xtmp[j] = 0;
      }
      else {
        y1 = INV_EST(y1);
        if (STRONGER_EST0(ym, y1)) {
          ym3 = ym2; i3 = i2;
          ym2 = ym; i2 = i1;
          ym = y1; i1 = j;
        }
        else {
          if (STRONGER_EST0(ym2, y1)) {
            ym3 = ym2; i3 = i2;
            ym2 = y1; i2 = j;
          }
          else {
            if (STRONGER_EST0(ym3, y1)) {
              ym3 = y1; i3 = j;
            }
          }
        }
        s1 += EST0_TO_LNP0(y1);
        xtmp[j] = 1;
      }
    }
    vpushx(cur_ind0, xtmp, n);
    slist[cur_ind0] += s1;

    cur_ind1 = branch(cur_ind0, cur_xl0);
    xtmp[i1] = 1 - xtmp[i1];
    vpushx(cur_ind1, xtmp, n);
    xtmp[i1] = 1 - xtmp[i1]; // Restore xtmp.
    slist[cur_ind1] += (s1 = EST0_ADD_LNP1_SUB_LNP0(ym));

    cur_ind1 = branch(cur_ind0, cur_xl0);
    xtmp[i2] = 1 - xtmp[i2];
    vpushx(cur_ind1, xtmp, n);
    xtmp[i2] = 1 - xtmp[i2]; // Restore xtmp.
    slist[cur_ind1] += (s2 = EST0_ADD_LNP1_SUB_LNP0(ym2));

    s3 = EST0_ADD_LNP1_SUB_LNP0(ym3);
    if (s1 + s2 > s3) {
      cur_ind1 = branch(cur_ind0, cur_xl0);
      xtmp[i1] = 1 - xtmp[i1];
      xtmp[i2] = 1 - xtmp[i2];
      vpushx(cur_ind1, xtmp, n);
      slist[cur_ind1] += s1 + s2;
    }
    else {
      cur_ind1 = branch(cur_ind0, cur_xl0);
      xtmp[i3] = 1 - xtmp[i3];
      vpushx(cur_ind1, xtmp, n);
      slist[cur_ind1] += s3;
    }
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
    vgetx(cur_ind1, xtmp, n);
    code_mm(n, xtmp, xtmp2);
    vsetx(cur_ind1, xtmp2, n);
    yp1 = YLISTP(cur_ind1, m);
    yp3 = YLISTP(cur_ind1, m) = YLIST_POP(m);
    for (j = 0; j < n; j++) {
      yp3[j] = xtmp[j] ? YLDEC1 : YLDEC0;
    }
  }

#ifdef DBG
  print_list("c");
#endif

}

void rmm_skip(
  decoder_type *dd, // Decoder instance data.
  int m
) {

  int n = 1 << m;
  int cur_ind0;
  slitem s1;
  ylitem *yp1, *yp3;
  int i, j;

  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    yp1 = YLISTP(cur_ind0, m);
    yp3 = YLISTP(cur_ind0, m) = YLIST_POP(m);
    s1 = 0.0;
    for (j = 0; j < n; j++) {
      s1 += EST_TO_LNP0(yp1[j]);
      yp3[j] = YLDEC0;
    }
    slist[cur_ind0] += s1;
  }
}

void rm_dec_inner(
  decoder_type *dd, // Decoder instance data.
  int m,
  int r
) {

  int n = 1 << m;
  int n2 = n / 2;
  int cur_ind0;
  ylitem *yp1, *yp2, *vp, *up;
  ylitem y1;
  int i, j;

  if (r == m) {
     if (node_table[node_counter] == 0) {
        rmm_skip(dd, m);
     }
     else {
        if (m == 1) {
           rm11_branch(dd, node_table[node_counter]);
        }
        else {
           rmm_branch(dd, m, node_table[node_counter]);
        }
     }
     node_counter++;
     return;
  }

  if (r == 1) {
     if (node_table[node_counter] == 0) {
        rm1_skip(dd, m);
     }
     else {
        rm1_branch(dd, m, node_table[node_counter]);
     }
     node_counter++;
  }
  else {

#ifdef DBG2
    printf("innr a: m=%d, r=%d, n=%3d.\n", m, r, n);
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

    rm_dec_inner(dd, m - 1, r - 1);

#ifdef DBG2
    printf("innr b: m=%d, r=%d, n=%3d.\n", m, r, n);
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
  }

#ifdef DBG2
  printf("innr c: m=%d, r=%d, n=%3d.\n", m, r, n);
#endif

  rm_dec_inner(dd, m - 1, r);

#ifdef DBG2
  printf("innr d: m=%d, r=%d, n=%3d.\n", m, r, n);
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
rm_dec(
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
  int c_r = dd->c_r;
  int peak_lsiz; // Peak list size.
  int flsiz; // Size of allocated list.
  ylitem *y_in;
  int n2 = dd->c_n / 2; // Half of the current length (n).
  int cur_ind0;
  ylitem *vp, *up;
  int *p1; // For permutations.
  slitem s1;
  ylitem y1;
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
  // Branch with permutations at (c_m, c_r).
  // List size is going to be dd->p_num.
  // All path metrics are 0 so far.
  for (i = 0; i < dd->p_num; i++) {
    slist[i] = 0.0;
    lorder[i] = i;
    lind2xl[i] = xlist;
    p1 = dd->pyarr + c_n * i;
    vp = YLISTP(i, c_m - 1) = YLIST_POP(c_m - 1);
    for (j = 0; j < n2; j++) vp[j] = XOR_EST(y_in[p1[j]], y_in[p1[j + n2]]);
    plist[i] = i;
  }
  cur_lsiz = frindp = dd->p_num;
  next_xle_ptr = xlist + 1;
  xlist[0].x = 0;
  xlist[0].p = NULL;

#ifdef DBG2
  printf("root 1: m=%d, r=%d, n=%3d.\n", c_m, c_r, c_n);
  print_list("root a");
#endif

  // Decode v.
  rm_dec_inner(dd, c_m - 1, c_r - 1);

#ifdef DBG2
  print_list("root b");
#endif

  // Calculate y_u = y_1 xor v + y_2.
  for (i = 0; i < cur_lsiz; i++) {
    cur_ind0 = lorder[i];
    vp = YLISTP(cur_ind0, c_m - 1);
    up = YLISTP(cur_ind0, c_m - 1) = YLIST_POP(c_m - 1);
    p1 = dd->pyarr + c_n * plist[cur_ind0];
    for (j = 0; j < n2; j++) {
      y1 = EST_XOR_YLDEC(y_in[p1[j]], vp[j]); // y1 <-- y_1[j] xor v[j].
      ADD_EST(y1, y_in[p1[j + n2]], up[j]); // y_u <-- y1 + y_2[j];
    }
  }

  // Decode u.
  rm_dec_inner(dd, c_m - 1, c_r);

  // Find the best.
  i1 = 0;
  for (i = 1; i < cur_lsiz; i++) {
    if (slist[lorder[i]] > slist[lorder[i1]]) {
      i1 = i;
    }
  }

#ifdef DBG
  print_list("final");
  printf("Best: %3.1f", slist[lorder[i1]]);
#endif

  vgetx(lorder[i1], xtmp, c_k);
  p1 = dd->pxarr + dd->c_k * plist[lorder[i1]];
  for (j = 0; j < c_k; j++) x_dec[p1[j]] = xtmp[j];

  if (s_dec != NULL) {
    if (dd->ret_s_sum > 0) {
      s1 = 0.0;
      for (i = 0; i < cur_lsiz; i++) s1 += slist[lorder[i]];
      (*s_dec) = s1;
    }
    else (*s_dec) = slist[lorder[i1]];
  }

  return 0;
}
