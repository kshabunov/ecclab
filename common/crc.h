//=============================================================================
// CRC calculation and checking header.
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

#ifndef CRC_H
#define CRC_H

/**
 * Lengths of crc_buf must be at least x_len.
 * CRC is calculated to crc_buf, starting at crc_buf + x_len - crc_len.
 * @return pointer to the first bit of the calculated CRC in crc_buf.
 */
int *calc_crc(
  const int *x,
  int x_len,
  const int *crc_poly,
  int crc_len,
  int *crc_buf
);

/**
 * x = (info_bits | crc), length is info_len + crc_len,
 * buf length must be at least info_len.
 * @return 0|1
 */
int is_valid_crc(
  const int *x,
  int info_len,
  const int *crc_poly,
  int crc_len,
  int *buf
);

#endif