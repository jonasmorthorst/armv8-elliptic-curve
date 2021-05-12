#include "extensionfield_interleaved.h"

#include "basefield.h"
#include "extensionfield.h"
#include <stdio.h>

/* Field elements are now shuffled like this:
 * a = {{a[0][0], a[1][0]}, {a[0][1], a[1][1]}}
 */

void ef_intrl_print_expr(ef_intrl_elem a) {
	ef_print_expr(ef_intrl_disentangle(a));
}

void ef_intrl_print_expr_nl(ef_intrl_elem a) {
	ef_intrl_print_expr(a);
	printf("\n");
}

void ef_intrl_print_hex(ef_intrl_elem a) {
	ef_print_hex(ef_intrl_disentangle(a));
}

void ef_intrl_print_hex_nl(ef_intrl_elem a) {
	ef_intrl_print_hex(a);
	printf("\n");
}

void ef_intrl_print_unred_expr(ef_intrl_elem_unred a) {
	poly64x2x2_t a0 = ef_intrl_disentangle_unred_lower(a);
	poly64x2x2_t a1 = ef_intrl_disentangle_unred_higher(a);
	printf("(");
	bf_print_unred_expr(a1);
	printf(")u + (");
	bf_print_unred_expr(a0);
	printf(")");
}

void ef_intrl_print_unred_expr_nl(ef_intrl_elem_unred a) {
	ef_intrl_print_unred_expr(a);
	printf("\n");
}
 
ef_elem ef_intrl_disentangle(ef_intrl_elem a) {
	ef_elem res;
	poly64x2_t t = vextq_p64(a.val[0], a.val[0], 1); //t = {a[1][0], a[0][0]}
	res.val[0] = vextq_p64(t, a.val[1], 1); //{a[0][0], a[0][1]}
	res.val[1] = a.val[1]; //{a[0][1], a[1][1]}
	res.val[1][0] = t[0]; //{a[1][0], a[1][1]}
	return res;
}

//{{a[0][0], a[0][1]}, {a[1][0], a[1][1]}}
ef_intrl_elem ef_intrl_interleave(ef_elem a) {
	ef_intrl_elem res;
	poly64x2_t t = vextq_p64(a.val[0], a.val[0], 1); //t = {a[0][1], a[0][0]}
	res.val[0] = vextq_p64(t, a.val[1], 1); //{a[0][0], a[1][0]}
	res.val[1] = a.val[1]; //{a[1][0], a[1][1]}
	res.val[1][0] = t[0]; //{a[0][1], a[1][1]}
	return res;
}

poly64x2x2_t ef_intrl_disentangle_unred_lower(ef_intrl_elem_unred c) {
	return (poly64x2x2_t) {{{c.val[0][0], c.val[1][0]}, {c.val[2][0], c.val[3][0]}}};
}

poly64x2x2_t ef_intrl_disentangle_unred_higher(ef_intrl_elem_unred c) {
	return (poly64x2x2_t) {{{c.val[0][1], c.val[1][1]}, {c.val[2][1], c.val[3][1]}}};
}

ef_intrl_elem ef_intrl_rand_elem() {
	return ef_intrl_interleave(ef_rand_elem());
}
