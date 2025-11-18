# SST Program - Fixed Array Implementation

## Overview

This version of the SST program uses **fixed-size arrays** instead of allocatable arrays. Array dimensions are determined at **compile time** using the `kLength` parameter set via a SLURM environment variable.

### Benefits:
- **Better performance**: Fixed arrays avoid runtime allocation overhead
- **Stack allocation**: Arrays can be allocated on the stack for better cache locality
- **Compile-time optimization**: Compiler can better optimize with known array sizes
- **Simpler memory management**: No need for allocate/deallocate

## Quick Start

### 1. Set kLength Environment Variable

```bash
export kLength=144
```

### 2. Compile

```bash
./compile_with_klength.sh
```

This will create `mainSSTProgram.x` (optimized) or `mainSST_debug.x` (debug mode).

### 3. Run

```bash
mpirun -np 16 ./mainSSTProgram.x
```

## Using with SLURM

### Method 1: Set in SLURM Script

Edit `slurm_compile_and_run.sh` and set the default kLength, then submit:

```bash
sbatch slurm_compile_and_run.sh
```

### Method 2: Pass as Environment Variable

```bash
sbatch --export=kLength=200 slurm_compile_and_run.sh
```

### Method 3: Interactive SLURM Session

```bash
salloc -N 1 -n 16
export kLength=144
./compile_with_klength.sh
mpirun -np 16 ./mainSSTProgram.x
```

## Important Notes

### kLength Must Match Data Files

The compiled program **MUST** use the same `kLength` as your data files (checkpoint.bin, weightStuff.bin, etc.).

The program will check this at runtime and error out if there's a mismatch:

```
ERROR: File kLength=150 does not match compiled KLENGTH_PARAM=144
kLength mismatch - recompile with correct kLength
```

If you see this error:
1. Check what kLength your data files were created with
2. Recompile with the correct value: `export kLength=150 && ./compile_with_klength.sh`

### Finding Your Data's kLength

You can check the kLength in your data files by looking at the first integer in any .bin file:

```bash
# For checkpoint files (example with Python)
python3 -c "import struct; f=open('results/param_1.0/checkpoint.bin','rb'); print('kLength =', struct.unpack('i', f.read(4))[0])"
```

### Recompiling for Different kLength

If you need to run with a different `kLength`:

```bash
# Clean previous build
rm -f *.mod *.o mainSSTProgram.x mainSST_debug.x

# Set new kLength
export kLength=200

# Recompile
./compile_with_klength.sh
```

## How It Works

### 1. Preprocessor Defines kLength

The compilation script passes `-DkLength=$kLength` to the compiler, which defines it as a preprocessor macro.

### 2. array_dimensions Module

The `array_dimensions.F90` module converts the preprocessor macro to a Fortran parameter:

```fortran
module array_dimensions
    integer, parameter :: KLENGTH_PARAM = kLength
end module
```

### 3. Fixed Array Declarations

All modules use `KLENGTH_PARAM` for array dimensions:

```fortran
real(dp) :: E(KLENGTH_PARAM), ET(KLENGTH_PARAM), F(KLENGTH_PARAM)
real(dp) :: weight(KLENGTH_PARAM, KLENGTH_PARAM, KLENGTH_PARAM)
```

### 4. Runtime Verification

Data loader subroutines verify that file kLength matches compiled kLength:

```fortran
read(10) kLength_file
if (kLength_file /= KLENGTH_PARAM) then
    stop 'kLength mismatch - recompile with correct kLength'
end if
```

## Compilation Options

### Optimized Build (default)

```bash
MODE="optimized" ./compile_with_klength.sh
```

Creates: `mainSSTProgram.x`
Flags: `-O3 -DkLength=$kLength`

### Debug Build

Edit `compile_with_klength.sh` and change:
```bash
MODE="debug"
```

Creates: `mainSST_debug.x`
Flags: `-g -O0 -DkLength=$kLength -fcheck=all -fbacktrace`

## Files Modified

The following files were modified to use fixed arrays:

1. **array_dimensions.F90** (NEW) - Defines compile-time array dimensions
2. **mainSST.F90** - Main program arrays converted to fixed size
3. **data_loader_all.F90** - Loader subroutines updated, added kLength verification
4. **timeRoutines.F90** - Time integration arrays converted to fixed size
5. **EDQNMstratifiedModules.F90** - Module arrays converted to fixed size
6. **compile_with_klength.sh** (NEW) - Compilation script
7. **slurm_compile_and_run.sh** (NEW) - SLURM submission script

## Troubleshooting

### Error: "kLength must be defined at compile time"

You forgot to set the kLength environment variable before compiling.

**Solution:**
```bash
export kLength=144
./compile_with_klength.sh
```

### Error: "kLength mismatch - recompile with correct kLength"

The executable was compiled with a different kLength than your data files.

**Solution:** Find the correct kLength from your data files and recompile.

### Segmentation Fault

If kLength is very large, you may exceed stack size limits.

**Solution:** Increase stack size:
```bash
ulimit -s unlimited
```

### Out of Memory During Compilation

Very large kLength values may require more memory during compilation.

**Solution:** Use a login node with more memory, or reduce kLength.

## Performance Considerations

### Array Size Impact

- **Small kLength** (< 100): Fast compilation, small memory footprint
- **Medium kLength** (100-300): Reasonable compilation time, moderate memory
- **Large kLength** (> 500): Slow compilation, large memory usage

### Stack vs Heap

Arrays are declared as local variables, which may use stack memory. For very large arrays, ensure adequate stack size with `ulimit -s unlimited`.

## Example Workflow

```bash
# 1. Set up environment
export kLength=144

# 2. Compile
./compile_with_klength.sh

# 3. Verify executable exists
ls -lh mainSSTProgram.x

# 4. Run with MPI
mpirun -np 16 ./mainSSTProgram.x

# 5. For different kLength, clean and rebuild
rm -f *.mod *.o mainSSTProgram.x
export kLength=200
./compile_with_klength.sh
```

## Questions?

For issues or questions about the fixed-array implementation, check:
1. Is `kLength` environment variable set correctly?
2. Does the kLength match your data files?
3. Did compilation complete without errors?
4. Is your stack size sufficient (`ulimit -s`)?
