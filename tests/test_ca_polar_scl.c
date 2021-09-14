#include "../interfaces/codec.h"
#include <tau/tau.h>

TAU_MAIN()

TEST(basic_init, cdc_init_emptyCfg) {
  void *cdc;
  int rc = cdc_init("", &cdc);
  REQUIRE_NE(rc, 0);
}

TEST(basic_init, cdc_init_noCm) {
  void *cdc;
  char config[] =
    "info_bits_mask\1""0000001101111111\1"
    "ca_polar_crc\1""1\1"
    "list_size\1 8";
  int rc = cdc_init(config, &cdc);
  REQUIRE_NE(rc, 0);
}

static char config_polar04_ca1_k8_l8[] =
  "c_m\1 4\1"
  "info_bits_mask\1""0000001101111111\1"
  "ca_polar_crc\1""1\1"
  "list_size\1 8";

TEST(basic_init, cdc_init_ok) {
  void *cdc;
  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);
}

TEST(cdc_face, cdc_getNandK) {
  void *cdc;

  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);

  int n = cdc_get_n(cdc);
  int k = cdc_get_k(cdc);

  REQUIRE_EQ(n, 16);
  REQUIRE_EQ(k, 8);
}

TEST(cdc_face, enc_bpsk) {
  void *cdc;
  int x[] = {0, 1, 1, 0, 1, 0, 1, 1};
  double y[16];
  double exp_y[] = {-1, 1, -1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1};

  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);

  enc_bpsk(cdc, x, y);

  for (int i = 0; i < 16; i++) {
    REQUIRE_EQ(y[i], exp_y[i]);
  }
}

TEST(cdc_face, dec_bpsk_noNoise) {
  void *cdc;
  int x[8];
  double y[16];
  int xd[8];

  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);

  for (int u = 0; u < 256; u++) {
    for (int j = 0; j < 8; j++) {
      x[j] = (u >> j) & 1;
    }
    enc_bpsk(cdc, x, y);
    dec_bpsk(cdc, y, xd);
    for (int i = 0; i < 8; i++) {
      REQUIRE_EQ(xd[i], x[i]);
    }
  }
}

TEST(cdc_face, dec_bpsk_withNoise_corrected) {
  void *cdc;
  int x[] = {0, 1, 1, 0, 1, 0, 1, 1};
  double y[16];
  int xd[8];

  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);

  for (int u = 0; u < 256; u++) {
    for (int j = 0; j < 8; j++) {
      x[j] = (u >> j) & 1;
    }
    enc_bpsk(cdc, x, y);
    y[u % 8] = -y[u % 8];
    dec_bpsk(cdc, y, xd);
    for (int i = 0; i < 8; i++) {
      REQUIRE_EQ(xd[i], x[i]);
    }
  }
}

TEST(cdc_face, dec_bpsk_withNoise_error) {
  void *cdc;
  int x[] = {0, 1, 1, 0, 1, 0, 1, 1};
  double y[16];
  int xd[8];

  int rc = cdc_init(config_polar04_ca1_k8_l8, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);
  enc_bpsk(cdc, x, y);
  for (int i = 0; i < 16; i++) y[i] = -y[i];

  dec_bpsk(cdc, y, xd);

  int isFailed = 0;
  for (int i = 0; i < 8; i++) {
    isFailed |= (xd[i] != x[i]);
  }
  REQUIRE_TRUE(isFailed);
}
