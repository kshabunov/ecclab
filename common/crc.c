//=============================================================================
// CRC calculation and checking.
//
// Copyright 2021 and onwards Kirill Shabunov.
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

#include <string.h>

int *calc_crc(
  const int *x,
  int x_len,
  const int *crc_poly,
  int crc_len,
  int *crc_buf
) {
  int i, j;

  memcpy(crc_buf, x, x_len * sizeof(int));
  for (i = 0; i < x_len - crc_len; i++) {
    if (crc_buf[i]) {
      for (j = 0; j < crc_len; j++) {
        crc_buf[i + j + 1] ^= crc_poly[j];
      }
    }
  }
  return crc_buf + x_len - crc_len;
}

int is_valid_crc(
  const int *x,
  int info_len,
  const int *crc_poly,
  int crc_len,
  int *buf
) {
  int s = 0, j;
  int *crc = calc_crc(x, info_len, crc_poly, crc_len, buf);
  for (j = 0; j < crc_len; j++) {
    s += x[info_len + j] ^ crc[j];
  }
  return (s == 0);
}