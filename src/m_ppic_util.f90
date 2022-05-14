module m_ppic_util
  use, intrinsic          :: ieee_arithmetic
  use, intrinsic          :: iso_fortran_env, only: real64, int64

  implicit none

  integer(kind=2), parameter :: LIQUID_PHASE = 1, GAS_PHASE = 4

  real(real64)            :: d_qnan = transfer(9221120237041090560_int64, 1._real64)

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

    x0_ = normal * shift
    if (present(x0)) then
      x0_ = x0_ + x0
    endif

    call init_box(liquid, [-dx/2.0, dx/2.0])
    call intersect_with_parabola(moments, liquid, normal, kappa0, x0_)
  end function

  real function cmpShift2d_parabolic(normal, dx, liqVol, kappa0, relTol, moments, grad_s) result(shift)
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
      shift = ieee_value(1.0_real64,  ieee_negative_inf)
      if (present(moments)) moments = 0.0
    elseif (liqVol >= cellVol) then
      shift = ieee_value(1.0_real64,  ieee_positive_inf)
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


  real function cmpSymmetricDifference2d_parabolic(x, dx, normal, shift, kappa0, customFun) result(sd)
    use m_r2d_parabolic
    implicit none

    real*8, intent(in)     :: x(2), dx(2), normal(2), shift, kappa0
    real*8, external, optional :: customFun

    ! Local variables
    real*8                :: sd_1(3), sd_2(3), normal_(2)
    type(r2d_poly_f)      :: exact_gas, exact_liq

    normal_ = -normal

    ! Construct polygonal approximation of exact gas & liquid domains
    ! (relative to the cell centroid)
    if (present(customFun)) then
      call get_polygonal_approximation_of_exact_domain(exact_gas, GAS_PHASE, x, dx, customFun)
      call get_polygonal_approximation_of_exact_domain(exact_liq, LIQUID_PHASE, x, dx, customFun)
    else
      call get_polygonal_approximation_of_exact_domain(exact_gas, GAS_PHASE, x, dx)
      call get_polygonal_approximation_of_exact_domain(exact_liq, LIQUID_PHASE, x, dx)
    endif

    ! Compute symmetric difference
    call intersect_with_parabola(sd_1, exact_gas, normal_, kappa0, x + normal_ * shift)
    call intersect_with_parabola(sd_2, exact_liq, -normal_, -kappa0, x + normal_ * shift)
    sd = sd_1(1) + sd_2(1)
  end

  subroutine get_polygonal_approximation_of_exact_domain(poly, phase, xc, dx, customFun)
    use m_optimization,   only: brent
    use m_r2d_parabolic
    implicit none

    type(r2d_poly_f), intent(out) :: poly
    integer*2, intent(in) :: phase
    real*8, intent(in)    :: xc(2), dx(2)
    real*8, external, optional :: customFun

    ! Local variables
    real*8                :: pos(2, R2D_MAX_VERTS), pos_skeleton(2, 8), corners(2, 4), funVals(4)
    real*8                :: x0(2), dir(2), step, tDir(2)
    integer               :: edx, vdx, vdx_first_inside, nrPos, vdx_next, nrPos_skelelton, rdx
    integer, parameter    :: VERTS_PER_SEGMENT = R2D_MAX_VERTS / 3
    logical               :: vdx_is_inside, vdx_next_is_inside, is_on_interface(8)

    corners(:,1) = xc + [-dx(1), -dx(2)]/2
    corners(:,2) = xc + [dx(1), -dx(2)]/2
    corners(:,3) = xc + [dx(1), dx(2)]/2
    corners(:,4) = xc + [-dx(1), dx(2)]/2

    ! Find out which corners of the cell are inside the domain
    vdx_first_inside = 0
    do vdx=1,4
      funVals(vdx) = interfaceFun(corners(:,vdx))
      if (funVals(vdx) >= 0.0 .and. vdx_first_inside == 0) vdx_first_inside = vdx
    enddo
    if (vdx_first_inside == 0) then
      poly%nverts = 0
      return
    endif

    ! Loop over the edges and construct the polygonal 'skeleton'
    vdx = vdx_first_inside
    vdx_is_inside = .true.
    nrPos_skelelton = 0

    do edx=1,4
      vdx_next = merge(1, vdx + 1, vdx == 4)
      vdx_next_is_inside = funVals(vdx_next) >= 0.0

      ! TODO so far we assume that an edge has at most one intersection
      if (vdx_is_inside .neqv. vdx_next_is_inside) then
        ! Find and add new position
        x0 = corners(:,vdx)
        dir = corners(:,vdx_next) - corners(:,vdx)
        step = brent(interfaceFun_step, 0.0D0, 1.0D0, 1D-15, 30, funVals(vdx), funVals(vdx_next))
        nrPos_skelelton = nrPos_skelelton + 1
        pos_skeleton(:,nrPos_skelelton) = x0 + step * dir
        is_on_interface(nrPos_skelelton) = .true.
      endif
      if (vdx_next_is_inside) then
        ! And add next node (corner)
        nrPos_skelelton = nrPos_skelelton + 1
        pos_skeleton(:,nrPos_skelelton) = corners(:,vdx_next)
        is_on_interface(nrPos_skelelton) = .false.
      endif

      vdx = vdx_next
      vdx_is_inside = vdx_next_is_inside
    enddo

    ! Now we add a refined approximation on edges that are on the interface
    nrPos = 0
    vdx = 1
    do edx=1,nrPos_skelelton
      vdx_next = merge(1, vdx + 1, vdx == nrPos_skelelton)

      ! Add (refinement of) the half open interval (pos_skeleton(:,vdx),pos_skeleton(:,vdx_next)]
      if (is_on_interface(vdx) .and. is_on_interface(vdx_next)) then

        tDir = pos_skeleton(:,vdx_next) - pos_skeleton(:,vdx)
        if (norm2(tDir) < 1E-15 * dx(1)) then
          nrPos = nrPos + 1
          pos(:,nrPos) = pos_skeleton(:,vdx_next)
        else
          ! Refine the face

          ! Make dir normal to the face
          dir = [-tDir(2), tDir(1)]
          do rdx=1,VERTS_PER_SEGMENT
            x0 = pos_skeleton(:,vdx) + rdx * tDir / VERTS_PER_SEGMENT

            ! We impose here that the radius of curvature of the interface is bounded from below by half (relative to the mesh spacing)
            step = brent(interfaceFun_step, -.5D0, .5D0, 1D-15, 30)

            nrPos = nrPos + 1
            pos(:,nrPos) = x0 + step * dir
          enddo
        endif
      else
        nrPos = nrPos + 1
        pos(:,nrPos) = pos_skeleton(:,vdx_next)
      endif

      vdx = vdx_next
    enddo
    call init_from_pos(poly, pos(:,1:nrPos))
  contains

    real*8 function interfaceFun(x) result(f)
  
      implicit none

      real*8, intent(in)  :: x(2)

      f = customFun(x)
      if (phase == GAS_PHASE) f = -f
    end

    real*8 function interfaceFun_step(step_) result(f)
      implicit none

      real*8, intent(in)    :: step_

      f = interfaceFun(x0 + step_ * dir)
    end
  end
end module
