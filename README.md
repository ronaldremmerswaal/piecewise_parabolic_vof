# Piecewise parabolic interface calculation for geometric VOF methods

## Description
Implementation of the intersection algorithm described in https://arxiv.org/abs/2111.09627, which is based on the two-dimensional part of https://github.com/devonmpowell/r3d.

## TODO
- Make tests
  - PPIC reconstruction
    - Exact reconstruction: (n, s, (kappa) -[intersection]> moments -> [reconstruction] n, s, (kappa))
    - Approximate reconstruction: validation using symmetric difference with exact domain
    - Timing (compare to PLIC)
- Include reconstruction functions: PLVIRA, PROST, LVIRA