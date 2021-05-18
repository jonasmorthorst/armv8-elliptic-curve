.text
.global lin_pass
.global c_sel

lin_pass:
	// Save
	stp x19, x20, [SP, #-16]!
	stp x21, x22, [SP, #-16]!
	stp x23, x24, [SP, #-16]!
	stp x25, x26, [SP, #-16]!
	stp x0, x1, [SP, #-16]!

	// Iteration i=0
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #0

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #0

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=1
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #1

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #1

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=2
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #2

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #2

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=3
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #3

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #3

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=4
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #4

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #4

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=5
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #5

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #5

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=6
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #6

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #6

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	add x2, x2, #64

	// Iteration i=7
	ldp x7, x8, [x2]
	ldp x9, x10, [x2, #16]
	ldp x11, x12, [x2, #32]
	ldp x13, x14, [x2, #48]

	cmp x3, #7

	csel x19, x7, x19, EQ
	csel x20, x8, x20, EQ
	csel x21, x9, x21, EQ
	csel x22, x10, x22, EQ
	csel x23, x11, x23, EQ
	csel x24, x12, x24, EQ
	csel x25, x13, x25, EQ
	csel x26, x14, x26, EQ

	cmp x4, #7

	csel x0, x7, x0, EQ
	csel x1, x8, x1, EQ
	csel x5, x9, x5, EQ
	csel x6, x10, x6, EQ
	csel x15, x11, x15, EQ
	csel x16, x12, x16, EQ
	csel x17, x13, x17, EQ
	csel x18, x14, x18, EQ

	// Store in output register addresses
	ldp x2, x3, [SP], #16

	stp x19, x20, [x2]
	stp x21, x22, [x2, #16]
	stp x23, x24, [x2, #32]
	stp x25, x26, [x2, #48]

	stp x0, x1, [x3]
	stp x5, x6, [x3, #16]
	stp x15, x16, [x3, #32]
	stp x17, x18, [x3, #48]

	// Restore
	ldp x25, x26, [SP], #16
	ldp x23, x24, [SP], #16
	ldp x21, x22, [SP], #16
	ldp x19, x20, [SP], #16
	ret

c_sel:
	cmp x0, x1
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
