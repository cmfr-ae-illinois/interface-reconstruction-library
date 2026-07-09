#ifndef ML_CLASSIFIER_C_API_H
#define ML_CLASSIFIER_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

int ml_classifier_fortran(double* vfrac, double* liq_bary);

#ifdef __cplusplus
}
#endif
#endif  // ML_CLASSIFIER_C_API_H