module m_common
  use, intrinsic          :: ieee_arithmetic
  use, intrinsic          :: iso_fortran_env, only: real64, int64

  real(real64), parameter :: d_neg_inf = -1D100
  real(real64), parameter :: d_pos_inf = 1D100

  real(real64), parameter :: d_qnan = transfer(9221120237041090560_int64, 1._real64)

  integer, parameter      :: LIQUID_PHASE = 1, GAS_PHASE = 4
end module