#include <arm_neon.h>
#include "ec.h"

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

typedef struct {
	uint64_t val[8];
} ec_naf;

ec_point_lproj ec_scalarmull_single(ec_point_laffine P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_single_lproj(ec_point_lproj P, uint64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, uint64x2x2_t k, ec_point_lproj Q, uint64x2x2_t l);

ec_point_laffine ec_scalarmull_single_endo(ec_point_laffine P, uint64x2x2_t k);

ec_naf ec_to_naf(poly64x2x2_t k);

void ec_print_naf(ec_naf k);

#endif
