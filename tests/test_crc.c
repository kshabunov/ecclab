#include <tau/tau.h>
#include "../common/crc.h"

TAU_MAIN()

TEST(face, calc_crc) {
  int x[] = {1, 0, 0, 1, 0, 0, 0, 0};
  int crc_poly[] = {0, 1, 1};
  int crc_len = 3;
  int y[8];
  int y_expected[] = {1, 0, 1, 0, 1, 1, 1, 1};
  int crc_expected[] = {1, 1, 1};
  int *crc;

  crc = calc_crc(x, 8, crc_poly, crc_len, y);

  for (int i = 0; i < crc_len; i++) {
    REQUIRE_EQ(crc[i], crc_expected[i]);
  }
  for (int i = 0; i < 8; i++) {
    REQUIRE_EQ(y[i], y_expected[i]);
  }
}

TEST(face, is_valid_crc) {
  int x_valid[] = {1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1};
  int x_invalid[] = {1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1};
  int crc[] = {0, 1, 1};
  int buf[8];

  REQUIRE_TRUE(is_valid_crc(x_valid, 8, crc, 3, buf));
  REQUIRE_FALSE(is_valid_crc(x_invalid, 8, crc, 3, buf));
}