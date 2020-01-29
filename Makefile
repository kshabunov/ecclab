CC = gcc -O3
BUILD_DIR = build

.PHONY: all
all: $(BUILD_DIR) $(BUILD_DIR)/dtrm_glp_gc

$(BUILD_DIR)/dtrm_glp_gc: dtrm_glp/dtrm_glp_inner.c dtrm_glp/dtrm_glp_main.c common/srm_utils.c common/spf_par.c common/ui_txt.c simulators/sim_bg.c simulators/main_txt.c
	$(CC) -DDEC_NEEDS_SIGMA -DDEC_NEEDS_CSNRN -o $@ $^

$(BUILD_DIR):
	mkdir $@
