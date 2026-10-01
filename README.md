# triso_sim

A Monte Carlo neutron transport simulator simulating a TRISO fuel particle. Written in 
Fortran and build with [fpm](https://fpm.fortran-lang.org/).

**Status** Currently a work in progress, unfinished modules and routines that need to be
built are notated as such with comments in the files.

# Final Objectives

- parallel on a M4 macbook with gpu acceleration based on the apple metal engine
- parallel on a windows pc with an i7 and RTX3080 with NVIDIA RTX acceleration

# Repository Structure

- `src/`  - All modules containing types, subroutines, and functions.
- `app/`  - Where the main programs are located, they handle file names and some io
- `test/` - Test script programs, checking statistical accuracy or debugging

# Build

- Requires a Fortran compiler, I'm using gfortran on the macbook. I will use an NVIDIA 
based compiler for the RTX3080. Only requirement is free-form source because of the ISO C
Binding use.

- To build with fpm installed simply bash
```bash
fpm build
```

# Dependencies 

- [Philox_Fortran](https://github.com/RJaBi/Philox_Fortran) (MIT)
- Read THIRDPARTY_LICENSES.txt for more info



