set(CMAKE_CXX_COMPILER mpic++)

IRLAddExample_cpp(
  TARGET_NAME
  amrex_advector
  HDRS
  "vtk.h"
  "amrcore_advector.h"
  "amrcore_advector.tpp"
  "case_deformation_3d.h"
  "case_translation_3d.h"
  "case_rotation_3d.h"
  "cases.h"
  "reconstruction_elvira.h"
  "reconstruction_mof1.h"
  "reconstruction_lvira.h"
  "reconstruction_plicnet.h"
  "reconstruction_vf.h"
  "reconstruction_vf2.h"
  "reconstruction_ivf.h"
  "reconstruction_pu.h"
  "reconstruction_cf.h"
  "reconstruction_hybrid.h"
  # "reconstruction.h"
  "advection_remap.h"
  "advection_helpers.h"
  SRCS
  "main.cpp"
  "vtk.cpp"
  # "amrcore_advector.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/amrex_advector"
  COLLECTION
  "examples_amrex_advector"
)

IRLAddExample_cpp(
  TARGET_NAME
  amrex_advector_checkpoint_postprocess
  HDRS
  "vtk.h"
  "amrcore_advector.h"
  "amrcore_advector.tpp"
  "case_deformation_3d.h"
  "case_translation_3d.h"
  "case_rotation_3d.h"
  "cases.h"
  "reconstruction_elvira.h"
  "reconstruction_mof1.h"
  "reconstruction_lvira.h"
  "reconstruction_plicnet.h"
  "reconstruction_vf.h"
  "reconstruction_vf2.h"
  "reconstruction_ivf.h"
  "reconstruction_pu.h"
  "reconstruction_cf.h"
  "reconstruction_hybrid.h"
  # "reconstruction.h"
  "advection_remap.h"
  "advection_helpers.h"
  SRCS
  "checkpoint_postprocess.cpp"
  "vtk.cpp"
  # "amrcore_advector.cpp"
  OUTPUT_DIR
  "${IRL_EXAMPLES_OUTPUT_DIR}/amrex_advector"
  COLLECTION
  "examples_amrex_advector"
)
