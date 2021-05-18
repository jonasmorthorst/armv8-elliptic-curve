.text
.global lin_pass
.global c_sel

lin_pass:
	// Save
	stp x19, x20, [SP, #-16]!
	stp x21, x22, [SP, #-16]!
	stp x23, x24, [SP, #-16]!
	stp x25, x26, [SP, #-16]!

	// Iteration i=0
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #0

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=1
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #1

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=2
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #2

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=3
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #3

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=4
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #4

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=5
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #5

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=6
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #6

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	add x1, x1, #64

	// Iteration i=7
	ldp x7, x8, [x1]
	ldp x9, x10, [x1, #16]
	ldp x11, x12, [x1, #32]
	ldp x13, x14, [x1, #48]

	subs x6, x2, #7

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	stp x19, x20, [x0]
	stp x21, x22, [x0, #16]
	stp x23, x24, [x0, #32]
	stp x25, x26, [x0, #48]

	// Restore
	ldp x25, x26, [SP], #16
	ldp x23, x24, [SP], #16
	ldp x21, x22, [SP], #16
	ldp x19, x20, [SP], #16
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
