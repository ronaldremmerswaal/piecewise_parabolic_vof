# Piecewise parabolic interface calculation for geometric VOF methods

## Description
Implementation of the intersection algorithm described in https://arxiv.org/abs/2111.09627, which is based on the two-dimensional part of https://github.com/devonmpowell/r3d.

## TODO
- Make tests
  - PPIC reconstruction
    - Approximate reconstruction: validation using symmetric difference with exact domain
    - Timing (compare to PLIC)
- Merge m_r2d_parabolic & m_reconstruction_util