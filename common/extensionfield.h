#include <arm_neon.h>

#ifndef EXTENSIONFIELD_H
#define EXTENSIONFIELD_H

typedef struct {
	poly64x2_t p0, p1;
} ef_elem;

ef_elem ef_create_elem(poly64x2_t p0, poly64x2_t p1);

void ef_print_expr(ef_elem a);

void ef_print_expr_nl(ef_elem a);

void ef_print_hex(ef_elem a);

void ef_print_hex_nl(ef_elem a);

ef_elem ef_rand_elem();

ef_elem ef_add(ef_elem a, ef_elem b);

ef_elem ef_mull(ef_elem a, ef_elem b);

ef_elem ef_square(ef_elem a);

ef_elem ef_inv(ef_elem a);
ef_elem ef_inv2(ef_elem a);

#endif
