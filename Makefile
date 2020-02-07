CC = gcc -O3
BUILD_DIR = work

.PHONY: all
all: $(BUILD_DIR) $(BUILD_DIR)/dtrm0_bg $(BUILD_DIR)/dtrm_glp_bg

$(BUILD_DIR)/dtrm0_bg: dtrm/dtrm0.c common/srm_utils.c common/spf_par.c common/ui_txt.c simulators/sim_bg.c simulators/main_txt.c
	$(CC) -DDEC_NEEDS_SIGMA -DUSE_STDLIB_RND -o $@ $^

$(BUILD_DIR)/dtrm_glp_bg: dtrm_glp/dtrm_glp_inner.c dtrm_glp/dtrm_glp_main.c common/srm_utils.c common/spf_par.c common/ui_txt.c simulators/sim_bg.c simulators/main_txt.c
	$(CC) -DDEC_NEEDS_SIGMA -DDEC_NEEDS_CSNRN -DUSE_STDLIB_RND -o $@ $^

$(BUILD_DIR):
	mkdir $@
