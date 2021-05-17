.text
.global lin_pass
.global c_sel

lin_pass:
	subs x6, x2, #0
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #1
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #2
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #3
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #4
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #5
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #6
	csel x15, x1, x15, EQ
	add x1, x1, #64
	subs x6, x2, #7
	csel x15, x1, x15, EQ

	// Load table[i] into x7-x14
	ldp x7, x8, [x15]
	ldp x9, x10, [x15, #16]
	ldp x11, x12, [x15, #32]
	ldp x13, x14, [x15, #48]

	// Store table[i] at P
	stp x7, x8, [x0]
	stp x9, x10, [x0, #16]
	stp x11, x12, [x0, #32]
	stp x13, x14, [x0, #48]
	ret

c_sel:
	subs x6, x0, x1
	csel x15, x3, x2, EQ

	// Load value into x7-x14
	ldp x7, x8, [x15]
	ldp x9, x10, [x15, #16]
	ldp x11, x12, [x15, #32]
	ldp x13, x14, [x15, #48]

	// Store value at P
	stp x7, x8, [x2]
	stp x9, x10, [x2, #16]
	stp x11, x12, [x2, #32]
	stp x13, x14, [x2, #48]
	ret
