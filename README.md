# ARMv8 GLS Elliptic Curve implementation
This library contains a GLS-254 elliptic curve implementation for ARMv8 Arch64.

## Structure
The `common` folder includes the following components for building the elliptic curve:

* [`common/setup.c`](common/setup.c): code used for initialization of necessary components.
* [`common/basefield.c`](common/basefield.c): code for the implementation of the basefield.
* [`common/extensionfield.c`](common/extensionfield.c): code for the first implementation of the extention field.  
* [`common/extensionfield_interleaved.c`](common/extensionfield_interleaved.c): code for the current implementation of the extention field using interleaved representation.
* [`common/ec.c`](common/ec.c): code for the implementation of the group law for the elliptic curve.
* [`common/ec_scalarmull.c`](common/ec_scalarmull.c): code for the implementation of scalar multiplication of points.
* [`common/utils.c`](common/utils.c): various helper methods.
* [`common/linear_pass.S`](common/linear_pass.S): external assembly implementing linear pass.


The `tests` folder includes test-code, while `benchmark` includes benchmarking-code.

## Supported platforms
This implementation is targeted ARMv8 Arch64 and heavily relies on Aarch64 ARM NEON instructions along with the ARM `crypto`.

## Running
To run the tests:
```
$ make runtests
```

To run in sandbox code:

```
$ make sandbox
```

Per default `clang` is used for compiling. To use `GCC`, add `GCC=1` flag to the command.

## Benchmarking
In order to read the cycle counter, we need to enable the cycle counter in kernel mode.

See this repository for instructions: https://github.com/zhiyisun/enable_arm_pmu

After enabling, we can run the following:

```
$ make runbench
```