cmake_minimum_required(VERSION 3.13)
project(ecclib)

include_directories(${CMAKE_SOURCE_DIR}/lib)

enable_testing()

set(SPF common/spf_par.c common/spf_par.h)
set(UI_TXT common/ui_txt.c interfaces/ui_utils.h)
set(SRM common/srm_utils.c common/srm_utils.h)
set(SIM_BG simulators/main_txt.c simulators/sim_bg.c  ${SPF} ${UI_TXT})
set(SIM_SRM_BG common/srm_utils.c ${SIM_BG})

add_compile_options(-DUSE_STDLIB_RND)

add_executable(dtrm0_bg dtrm/dtrm0.c ${SIM_SRM_BG})
target_compile_options(dtrm0_bg PUBLIC -DDEC_NEEDS_SIGMA)

add_executable(test_dtrm0 tests/test_dtrm0.c dtrm/dtrm0.c ${SRM} ${SPF} ${UI_TXT})
target_compile_options(test_dtrm0 PUBLIC -DDEC_NEEDS_SIGMA)
add_test(dtrm0 test_dtrm0)

add_executable(dtrm1_bg dtrm/dtrm1.c rm1_ml/rm1_ml.c ${SIM_SRM_BG})
target_compile_options(dtrm1_bg PUBLIC -DDEC_NEEDS_SIGMA)

add_executable(dtrm_glp_bg dtrm_glp/dtrm_glp_main.c dtrm_glp/dtrm_glp_inner.c ${SIM_SRM_BG})
target_compile_options(dtrm_glp_bg PUBLIC -DDEC_NEEDS_SIGMA -DDEC_NEEDS_CSNRN)

add_executable(ca_polar_scl_bg polar_scl/ca_polar_scl_main.c polar_scl/polar_scl_inner.c ${SIM_SRM_BG})
target_compile_options(ca_polar_scl_bg PUBLIC -DDEC_NEEDS_SIGMA)

INSTALL(TARGETS dtrm0_bg dtrm1_bg dtrm_glp_bg ca_polar_scl_bg DESTINATION ${CMAKE_SOURCE_DIR}/work)