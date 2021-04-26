#include <arm_neon.h>
#include "ec.h"

#ifndef EC_SCALARMULL_H
#define EC_SCALARMULL_H

ec_point_lproj ec_scalarmull_single(ec_point_lproj P, poly64x2x2_t k);

ec_point_lproj ec_scalarmull_double(ec_point_lproj P, poly64x2x2_t k, ec_point_lproj Q, poly64x2x2_t l);

#endif
