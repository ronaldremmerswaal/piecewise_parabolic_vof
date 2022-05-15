module m_ppic_util
  use m_common

  implicit none


  private
  public :: cmpMoments2d_parabolic, cmpShift2d_parabolic, cmpSymmDiff2d_parabolic
contains

  function cmpMoments2d_parabolic(normal, dx, shift, kappa0, x0) result(moments)
    use m_r2d_parabolic

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), shift, kappa0
    real*8, intent(in), optional :: x0(2)
    real*8                :: moments(3)

    ! Local variables
    type(r2d_poly_f)      :: liquid
    real*8                :: x0_(2)
    
    if (shift == d_pos_inf) then
      moments(1) = product(dx)
      moments(2:3) = 0.D0
      return
    elseif (shift == d_neg_inf) then
      moments = 0.D0
      return
    endif

    
    x0_ = normal * shift
    if (present(x0)) then
      x0_ = x0_ + x0
    endif

    call init_box(liquid, [-dx/2.0, dx/2.0])
    call intersect_with_parabola(moments, liquid, normal, kappa0, x0_)
  end function

  real*8 function cmpShift2d_parabolic(normal, dx, liqVol, kappa0, relTol, moments, grad_s) result(shift)
    use m_plic_util,      only: cmpShift2d
    use m_r2d_parabolic
    use m_optimization,   only: brent

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), liqVol, kappa0
    real*8, optional, intent(in) :: relTol
    real*8, optional, intent(out) :: moments(3), grad_s(2)

    ! Local variables
    type(r2d_poly_f)      :: cell
    real*8                :: moments_(3), max_shift_plane, cellVol, shift_plane, plane_err
    real*8                :: shift_l, shift_r, err_l, err_r, normal3(3), dx3(3)
    real*8                :: relTol_, grad_s_(2)

    max_shift_plane = dot_product(abs(normal), dx)/2

    cellVol = product(dx)
    if (liqVol <= 0.0) then
      shift = d_neg_inf
      if (present(moments)) moments = 0.0
    elseif (liqVol >= cellVol) then
      shift = d_pos_inf
      if (present(moments)) moments = 0.0
    else
      call init_box(cell, [-dx/2.0, dx/2.0])

      ! Use PLIC to get a good initial bracket (one side at least)
      shift_plane = cmpShift2d(normal, dx, liqVol)
      plane_err = volume_error_function(shift_plane)

      ! Try to get a better bracket with an (educated) guess
      if (plane_err == 0.0) then
        ! iff kappa0 == 0.0
        if (present(moments)) moments = moments_
        if (present(grad_s)) grad_s = grad_s_
        return
      elseif (plane_err > 0.0) then
        ! iff kappa0 < 0.0
        shift_r = shift_plane
        err_r = plane_err

        shift_l = max_shift_plane * (kappa0 * max_shift_plane / 2. - 1.)
        err_l = -liqVol
      else
        ! iff kappa0 > 0.0
        shift_l = shift_plane
        err_l = plane_err

        shift_r = max_shift_plane * (kappa0 * max_shift_plane / 2. + 1.)
        err_r = cellVol - liqVol
      endif

      if (present(relTol)) then
        relTol_ = relTol
      else
        relTol_ = 1E-15
      endif

      ! NB the centered limits -pc_shift, pc_shift are incorrect for a parabolic interface
      shift = brent(volume_error_function, shift_l, shift_r, max_shift_plane * relTol_, 30, err_l, err_r)
    endif

    if (present(moments)) moments = moments_
    if (present(grad_s)) grad_s = grad_s_
  contains

    real*8 function volume_error_function(pc_tmp) result(err)
      implicit none

      real*8, intent(in)    :: pc_tmp

      ! Local variables
      type(r2d_poly_f)    :: liquid
      integer             :: idx

      grad_s_ = [d_qnan, d_qnan]
      call copy(to=liquid, from=cell)
      call intersect_with_parabola(moments_, liquid, normal, kappa0, normal * pc_tmp, grad_s=grad_s_)
      err = moments_(1) - liqVol
    end function
  end function


  real function cmpSymmDiff2d_parabolic(x, dx, normal, shift, kappa0, levelSet) result(sd)
    use m_r2d_parabolic
    implicit none

    real*8, intent(in)     :: x(2), dx(2), normal(2), shift, kappa0
    real*8, external       :: levelSet

    ! Local variables
    real*8                :: sd_1(3), sd_2(3)
    type(r2d_poly_f)      :: exact_gas, exact_liq


    ! Construct polygonal approximation of exact gas & liquid domains
    ! (relative to the cell centroid)
    call get_polygonal_approximation_of_exact_domain(exact_gas, x, dx, levelSet, GAS_PHASE)
    call get_polygonal_approximation_of_exact_domain(exact_liq, x, dx, levelSet, LIQUID_PHASE)

    ! Compute symmetric difference
    call intersect_with_parabola(sd_1, exact_gas, normal, kappa0, x + normal * shift)
    call intersect_with_parabola(sd_2, exact_liq, -normal, -kappa0, x + normal * shift)
    sd = sd_1(1) + sd_2(1)
  end
end module
