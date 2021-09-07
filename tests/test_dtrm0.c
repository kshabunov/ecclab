#include "../interfaces/codec.h"
#include <tau/tau.h>

TAU_MAIN()

TEST(basic_init, cdc_init_emptyCfg) {
  void *cdc;
  int rc = cdc_init("", &cdc);
  REQUIRE_NE(rc, 0);
}

TEST(basic_init, cdc_init_noRMm) {
  void *cdc;
  int rc = cdc_init("RM_r\1 2", &cdc);
  REQUIRE_NE(rc, 0);
}

TEST(basic_init, cdc_init_noRMr) {
  void *cdc;
  int rc = cdc_init("RM_m\1 4\1", &cdc);
  REQUIRE_NE(rc, 0);
}

TEST(basic_init, cdc_init_ok) {
  void *cdc;
  char config[] =
    "RM_m\1 4\1"
    "RM_r\1 2";
  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);
}

TEST(cdc_face, cdc_getNandK) {
  void *cdc;
  char config[] =
    "RM_m\1 4\1"
    "RM_r\1 2";

  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);

  int n = cdc_get_n(cdc);
  int k = cdc_get_k(cdc);

  REQUIRE_EQ(n, 16);
  REQUIRE_EQ(k, 11);
}

TEST(cdc_face, enc_bpsk) {
  void *cdc;
  char config[] =
    "RM_m\1 3\1"
    "RM_r\1 1";
  int x[] = {0, 1, 1, 0};
  double y[8];
  double exp_y[] = {-1, 1, 1, -1, -1, 1, 1, -1};

  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);

  enc_bpsk(cdc, x, y);

  for (int i = 0; i < 8; i++) {
    REQUIRE_EQ(y[i], exp_y[i]);
  }
}

TEST(cdc_face, dec_bpsk_noNoise) {
  void *cdc;
  char config[] =
    "RM_m\1 4\1"
    "RM_r\1 1";
  int x[5];
  double y[16];
  int xd[5];

  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);

  for (int u = 0; u < 32; u++) {
    for (int j = 0; j < 5; j++) {
      x[j] = (u >> j) & 1;
    }
    enc_bpsk(cdc, x, y);
    dec_bpsk(cdc, y, xd);
    for (int i = 0; i < 5; i++) {
      REQUIRE_EQ(xd[i], x[i]);
    }
  }
}

TEST(cdc_face, dec_bpsk_withNoise_corrected) {
  void *cdc;
  char config[] =
    "RM_m\1 4\1"
    "RM_r\1 1";
  int x[] = {0, 1, 1, 0, 1};
  double y[16];
  int xd[5];

  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);
  enc_bpsk(cdc, x, y);
  y[0] = -y[0];
  y[1] = -y[1];
  y[2] = -y[2];

  dec_bpsk(cdc, y, xd);

  for (int i = 0; i < 5; i++) {
    REQUIRE_EQ(xd[i], x[i]);
  }
}

TEST(cdc_face, dec_bpsk_withNoise_error) {
  void *cdc;
  char config[] =
    "RM_m\1 4\1"
    "RM_r\1 1";
  int x[] = {0, 1, 1, 0, 1};
  double y[16];
  int xd[5];

  int rc = cdc_init(config, &cdc);
  REQUIRE_EQ(rc, 0);
  cdc_set_sg(cdc, 1.0);
  enc_bpsk(cdc, x, y);
  for (int i = 0; i < 16; i++) y[i] = -y[i];

  dec_bpsk(cdc, y, xd);

  int isFailed = 0;
  for (int i = 0; i < 5; i++) {
    isFailed |= (xd[i] != x[i]);
  }
  REQUIRE_TRUE(isFailed);
}
