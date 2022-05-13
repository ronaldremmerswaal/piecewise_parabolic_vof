# Piecewise parabolic interface calculation for geometric VOF methods

## Description
Implementation of the intersection algorithm described in https://arxiv.org/abs/2111.09627, which is based on the two-dimensional part of https://github.com/devonmpowell/r3d.
A FORTRAN wrapper is provided, where in particular the subroutine intersect_with_parabola is of interest.

## TODO
- Make tests
  - Intersection (in C)
    - Compare intersection moments to at least one analytical example
    - Verify that intersection splits into two parts, whose moments add up to the original part
    - Timing (compare to linear)
  - PPIC reconstruction (in Fortran)
    - Exact parabolic reconstruction (n, s, kappa -[intersection]> moments -> [reconstruction] n, s, kappa)
    - (Approximate) parabolic reconstruction; validation using piecewise linear approximation of the exact domain
    - Timing (compare to PLIC)
- Include computation of shift using Brent
- Include PPIC reconstruction functions: PLVIRA, PROST, PMoF
- Include linear (PLIC) module for comparison