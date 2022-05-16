module m_reconstruction_util
  use m_common

  private
  public :: cmpMoments, cmpShift, cmpSymmDiff

  real*8, parameter     :: NORMAL_TOL = 1E-12

  interface cmpMoments
    module procedure cmpMoments2d_plane, cmpMoments2d_parabolic
  end interface

  interface cmpShift
    module procedure cmpShift2d_plane, cmpShift2d_parabolic, cmpShift2d_poly
  end interface

  interface cmpSymmDiff
    module procedure cmpSymmDiff2d_plane, cmpSymmDiff2d_parabolic
  end interface

contains

   function cmpMoments2d_plane(normal, dx, shift) result(moments)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), shift
    real*8              :: moments(3)

    ! Local variables:
    real*8              :: normal_(2), dx_(2), max_shift_plane
    integer             :: perm(2)

    normal_ = abs(normal)
    if (dx(1) * normal_(1) < dx(2) * normal_(2)) then
      dx_ = dx
      perm = [1, 2]
    else
      normal_ = [normal_(2), normal_(1)]
      dx_ = [dx(2), dx(1)]
      perm = [2, 1]
    endif

    max_shift_plane = dot_product(normal_, dx_)/2
    if (normal_(1) < NORMAL_TOL) then
      moments(1) = ((shift + max_shift_plane) / normal_(2)) * dx_(1)
      moments(2) = 0.
      moments(3) = 0.5 * (((shift + max_shift_plane) / normal_(2) - dx_(2) / 2.0)**2 - (dx_(2)/2.0)**2) * dx_(1)
    else
      moments = cmpMoments2d_pc_ordered_plane(normal_, dx_, shift + max_shift_plane, max_shift_plane)
    endif

    if (normal(perm(1)) < 0.0) moments(2) = -moments(2)
    if (normal(perm(2)) < 0.0) moments(3) = -moments(3)
    if (perm(1) /= 1) then 
      moments(2:3) = [moments(3), moments(2)]
    endif
  end function


  pure function cmpMoments2d_pc_ordered_plane(normal, dx, pc, max_shift_plane) result(moments)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), pc, max_shift_plane
    real*8              :: moments(3)

    ! Local variables:
    real*8              :: pcMod
    logical             :: largerThanHalf

    largerThanHalf = pc > max_shift_plane
    if (largerThanHalf) then
      pcMod = 2*max_shift_plane - pc
    else
      pcMod = pc
    endif

    ! We compute moment w.r.t. cell corner
    if (pcMod <= 0.0) then
      moments = 0.0
    elseif (pcMod<normal(1)*dx(1)) then
      moments(1) = pcMod**2/(2.0*normal(1)*normal(2))
      moments(2) = (pcMod**2*(2*pcMod - 3*dx(1)*normal(1)))/(12*normal(1)**2*normal(2))
      moments(3) = (pcMod**2*(2*pcMod - 3*dx(2)*normal(2)))/(12*normal(1)*normal(2)**2)
    else ! if (pcMod<=max_shift_plane) then
      moments(1) = pcMod*dx(1)/normal(2) - normal(1)*dx(1)**2/(2.0*normal(2))
      moments(2) = -(dx(1)**3*normal(1))/(12*normal(2))
      moments(3) = (dx(1)*(6*pcMod**2 + 2*(dx(1)*normal(1))**2 - 6*dx(1)*normal(1)*pcMod &
        - 6*dx(2)*normal(2)*pcMod + 3*dx(1)*dx(2)*normal(1)*normal(2)))/(12*normal(2)**2)
    endif

    if (largerThanHalf) moments(1) = product(dx) - moments(1)
  end function

  pure real*8 function cmpShift2d_plane(normal, dx, volume) result(shift)

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), volume
    
    ! Local variables:
    real*8              :: normal_(2), dx_(2), max_shift_plane

    normal_ = abs(normal)
    if (dx(1) * normal_(1) < dx(2) * normal_(2)) then
      dx_ = dx
    else
      normal_ = [normal_(2), normal_(1)]
      dx_ = [dx(2), dx(1)]
    endif

    max_shift_plane = dot_product(normal_, dx_)/2
    if (normal_(1) < NORMAL_TOL) then
      shift = normal_(2) * volume / dx_(1) - max_shift_plane
    else
      shift = cmpShift2d_pc_ordered_plane(normal_, dx_, volume, max_shift_plane) - max_shift_plane
    endif


  end function

  pure real*8 function cmpShift2d_pc_ordered_plane(normal, dx, volume, max_shift_plane) result(pc)

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), volume, max_shift_plane

    ! Local variables:
    real*8                :: max_volume_plane        
    real*8                :: v1, vMod
    logical               :: largerThanHalf

    max_volume_plane = product(dx)
    v1  = normal(1) * dx(1)**2 / (2*normal(2))

    largerThanHalf = volume > max_volume_plane / 2.0
    if (largerThanHalf) then
      vMod = max_volume_plane - volume
    else
      vMod = volume
    endif

    if (vMod <= 0.0) then
      pc = 0.0D0
    elseif (vMod < v1) then
      pc = sqrt(2*normal(1)*normal(2)*vMod)
    else ! if (vMod <= .5*max_volume_plane) then
      pc = (vMod + normal(1)*dx(1)**2/(2*normal(2))) * normal(2) / dx(1)
    endif

    if (largerThanHalf) pc = 2*max_shift_plane - pc

  end function

  real function cmpSymmDiff2d_plane(x, dx, normal, shift, levelSet) result(sd)
    use m_r2d_parabolic
    implicit none

    real*8, intent(in)     :: x(2), dx(2), normal(2), shift
    real*8, external       :: levelSet

    ! Local variables
    real*8                :: sd_1(3), sd_2(3)
    type(r2d_poly_f)      :: exact_gas, exact_liq


    ! Construct polygonal approximation of exact gas & liquid domains
    ! (relative to the cell centroid)
    exact_gas = polygonalApproximation(x, dx, levelSet, GAS_PHASE)
    exact_liq = polygonalApproximation(x, dx, levelSet, LIQUID_PHASE)
    call shift_by(exact_gas, -x)
    call shift_by(exact_liq, -x)

    ! Compute symmetric difference
    call intersect_by_plane(exact_gas, normal, shift)
    sd_1 = cmpMoments(exact_gas)
    call intersect_by_plane(exact_liq, -normal, -shift)
    sd_2 = cmpMoments(exact_liq)
    sd = sd_1(1) + sd_2(1)
  end

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

  real*8 function cmpShift2d_parabolic(normal, dx, liqVol, kappa0, relTol, moments) result(shift)
    use m_r2d_parabolic
    use m_optimization,   only: brent

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), liqVol, kappa0
    real*8, optional, intent(in) :: relTol
    real*8, optional, intent(out) :: moments(3)

    ! Local variables
    type(r2d_poly_f)      :: cell
    real*8                :: moments_(3), max_shift_plane, cellVol, shift_plane, plane_err
    real*8                :: shift_l, shift_r, err_l, err_r, normal3(3), dx3(3)
    real*8                :: relTol_

    max_shift_plane = dot_product(abs(normal), dx)/2

    cellVol = product(dx)
    if (liqVol <= 0.0) then
      shift = -100 * max_shift_plane
      if (present(moments)) moments = 0.0D0
      return
    elseif (liqVol >= cellVol) then
      shift = 100 * max_shift_plane
      if (present(moments)) moments = [cellVol, 0.0D0, 0.0D0]
      return
    endif

    call init_box(cell, [-dx/2.0, dx/2.0])

    ! Use PLIC to get a good initial bracket (one side at least)
    shift_plane = cmpShift2d_plane(normal, dx, liqVol)
    plane_err = volume_error_function(shift_plane)

    ! Try to get a better bracket with an (educated) guess
    if (plane_err == 0.0) then
      ! iff kappa0 == 0.0
      shift = shift_plane
      if (present(moments)) moments = moments_
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

    relTol_ = merge(relTol, 1D-15, present(relTol))

    shift = brent(volume_error_function, shift_l, shift_r, max_shift_plane * relTol_, 30, err_l, err_r)
    if (present(moments)) moments = moments_

  contains

    real*8 function volume_error_function(shift_tmp) result(err)
      implicit none

      real*8, intent(in)    :: shift_tmp

      ! Local variables
      type(r2d_poly_f)    :: liquid
      integer             :: idx

      call copy(to=liquid, from=cell)
      call intersect_with_parabola(moments_, liquid, normal, kappa0, normal * shift_tmp)
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
    exact_gas = polygonalApproximation(x, dx, levelSet, GAS_PHASE)
    exact_liq = polygonalApproximation(x, dx, levelSet, LIQUID_PHASE)

    ! Compute symmetric difference
    call intersect_with_parabola(sd_1, exact_gas, normal, kappa0, x + normal * shift)
    call intersect_with_parabola(sd_2, exact_liq, -normal, -kappa0, x + normal * shift)
    sd = sd_1(1) + sd_2(1)
  end

  real*8 function cmpShift2d_poly(normal, poly, liqVol, kappa0, relTol, moments) result(shift)
    use m_r2d_parabolic
    use m_optimization,   only: brent

    implicit none

    real*8, intent(in)    :: normal(2), liqVol, kappa0
    type(r2d_poly_f), intent(in) :: poly
    real*8, optional, intent(in) :: relTol
    real*8, optional, intent(out) :: moments(3)

    ! TODO support reconstruction in arbitrary cells; make function for common part with cmpShift2d_parabolic
  end function
end module