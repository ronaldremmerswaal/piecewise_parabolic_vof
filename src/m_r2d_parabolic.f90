module m_r2d_parabolic
  use m_r2d

  implicit none

  type, bind(C) :: r2d_parabola_f
    type(r2d_rvec2_f)     :: n
    type(r2d_rvec2_f)     :: x0
    real(c_double)        :: kappa0
  end type

  interface

    subroutine r2d_clip_parabola_cmpMoments(poly, parabola, grad_s, moments, derivative, &
      compute_derivative) bind(C, name="r2d_clip_parabola_cmpMoments")
      import r2d_poly_f, r2d_parabola_f, c_double, c_bool

      implicit none

      type(r2d_poly_f), intent(inout) :: poly
      type(r2d_parabola_f), intent(in) :: parabola(1)
      real(c_double), intent(inout):: grad_s(2)
      real(c_double), intent(out)  :: moments(3), derivative(4)
      logical(c_bool), intent(in)  :: compute_derivative
    end subroutine r2d_clip_parabola_cmpMoments
  end interface

  interface cmpMoments
    module procedure cmpMoments_levelset, cmpMoments_levelset_poly, cmpMoments_parabola_poly
  end interface

  interface cmpMoments_
    module procedure cmpMoments_parabola_poly_
  end interface

  interface polygonalApproximation
    module procedure polygonalApproximation_rectangularIn, polygonalApproximation_polyIn
  end interface

  interface makeParabola
    module procedure makeParabola_x0, makeParabola_shift, makeParabola_angle
  end interface
contains

  function makeParabola_x0(normal, kappa0, x0) result(parabola)
    implicit none
    
    real*8, intent(in)    :: normal(2), kappa0, x0(2)
    type(r2d_parabola_f)  :: parabola

    parabola%n%xyz = normal
    parabola%kappa0 = kappa0
    parabola%x0%xyz = x0
  end function

  function makeParabola_shift(normal, kappa0, shift) result(parabola)
    implicit none
    
    real*8, intent(in)    :: normal(2), kappa0, shift
    type(r2d_parabola_f)  :: parabola

    parabola = makeParabola_x0(normal, kappa0, shift * normal)
  end function

  function makeParabola_angle(angle, kappa0, shift) result(parabola)
    implicit none
    
    real*8, intent(in)    :: angle, kappa0, shift
    type(r2d_parabola_f)  :: parabola

    real*8                :: normal(2)

    normal = [dcos(angle), dsin(angle)]
    parabola = makeParabola_x0(normal, kappa0, shift * normal)
  end function

  ! Removes the part of poly for which 
  !   \eta \cdot (x - x_0) + (\kappa/2) * (\tau \cdot (x - x_0))^2 > 0
  ! and returns the zeroth (M_0) and first (M_1) moments
  !
  ! IN/OUTPUT
  ! 
  ! moments01             The zeroth and first moments of intersection (so moments01(1) is the 
  !                         volume and moments01(2:3)/moments01(1) is the centroid of the
  !                         intersection volume)
  ! poly                  Polygon which is to be intersected
  ! parabola              The parabola
  ! derivative            The derivatives of the zeroth and first moments (optional):  
  !                         derivative(1)   = d M_0/ d \varphi
  !                         derivative(2:3) = d M_1/ d \varphi
  !                         derivative(4)   = d M_0/ d \kappa
  ! grad_s                The gradient of the shift s w.r.t. the normal angle (\varphi) and curvature (\kappa)
  !                         If initially NaN, then this gradient is computed such that d M_0/ d \varphi = 0
  !                         If not, then it is assumed that the gradient of the shift s is given by grad_s
  !                         (This only affects the computation of the gradient)
  function cmpMoments_parabola_poly(poly, parabola, derivative, grad_s) result(moments01)
    use m_common

    implicit none

    type(r2d_poly_f), intent(in) :: poly
    type(r2d_parabola_f), intent(in) :: parabola
    real*8, intent(inout), optional :: derivative(4), grad_s(2)
    real*8                :: moments01(3)

    type(r2d_poly_f)      :: poly_

    poly_ = copy(poly)
    moments01 = cmpMoments_parabola_poly_(poly_, parabola, derivative, grad_s)
  end function

  ! This one modifies poly itself!
  function cmpMoments_parabola_poly_(poly, parabola, derivative, grad_s) result(moments01)
    use m_common

    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    type(r2d_parabola_f), intent(in) :: parabola
    real*8, intent(inout), optional :: derivative(4), grad_s(2)
    real*8                :: moments01(3)

    ! Local variables
    type(r2d_parabola_f)  :: parabolas(1)
    real*8                :: moments01_poly(3)
    real*8                :: derivative_(4), grad_s_(2)
    logical*1             :: compute_derivative_
    integer               :: vdx

    
    ! The algorithm is implemented for kappa0 <= 0, so we must compute the moments of the complement
    ! (for which the sign of the curvature is swapped) if kappa0 > 0
    if (parabola%kappa0 > 0) then
      parabolas(1)%n%xyz = -parabola%n%xyz
      parabolas(1)%x0 = parabola%x0
      parabolas(1)%kappa0 = -parabola%kappa0
      moments01_poly = cmpMoments(poly)
    else
      parabolas(1) = parabola
    endif

    ! Setting grad_s_ to NaN ensures that grad_s_ is computed such that the derivative of the zeroth 
    ! moment w.r.t. the normal angle is zero (derivative_(1) = 0)
    grad_s_ = merge(grad_s, [d_qnan, d_qnan], present(grad_s))

    compute_derivative_ = present(derivative) .or. present(grad_s)
    call r2d_clip_parabola_cmpMoments(poly, parabolas, grad_s_, moments01, derivative_, compute_derivative_)

    if (parabola%kappa0 > 0) then
      moments01 = moments01_poly - moments01
      derivative_(1:3) = -derivative_(1:3)
    endif

    if (present(derivative)) derivative = derivative_
    if (present(grad_s)) then 
      grad_s = grad_s_
      if (parabola%kappa0 > 0) grad_s(1) = -grad_s(1)
    endif
  end function

  function cmpMoments_levelset(x, dx, levelSet, phase, verts_per_segment) result(moments)
    implicit none
    
    real*8, intent(in)    :: x(2), dx(2)
    real*8, external      :: levelSet
    real*8                :: moments(3)
    integer, intent(in), optional :: phase
    integer, intent(in), optional :: verts_per_segment

    moments = cmpMoments(polygonalApproximation_rectangularIn(x, dx, levelSet, phase, verts_per_segment))
  end function

  function cmpMoments_levelset_poly(poly, levelSet, phase, verts_per_segment) result(moments)
    implicit none
    
    type(r2d_poly_f), intent(in) :: poly
    real*8, external      :: levelSet
    real*8                :: moments(3)
    integer, intent(in), optional :: phase
    integer, intent(in), optional :: verts_per_segment

    moments = cmpMoments(polygonalApproximation_polyIn(poly, levelSet, phase, verts_per_segment))
  end function

  function polygonalApproximation_rectangularIn(x, dx, levelSet, phase, verts_per_segment) result(poly)
    use m_common
    use m_optimization,   only: brent
  
    implicit none

    real*8, intent(in)    :: x(2), dx(2)
    real*8, external      :: levelSet
    integer, intent(in), optional :: phase
    integer, intent(in), optional :: verts_per_segment
    type(r2d_poly_f)      :: poly

    poly = polygonalApproximation_polyIn(makeBox(x, dx), levelSet, phase, verts_per_segment)
  end function

  function polygonalApproximation_polyIn(cell, levelSet, phase, verts_per_segment) result(poly)
    use m_common
    use m_optimization,   only: brent
  
    implicit none

    type(r2d_poly_f), intent(in) :: cell
    real*8, external      :: levelSet
    integer, intent(in), optional :: phase
    integer, intent(in), optional :: verts_per_segment
    type(r2d_poly_f)      :: poly

    ! Local variables
    real*8                :: pos(2, R2D_MAX_VERTS), pos_skeleton(2, 2*R2D_MAX_VERTS), funVals(R2D_MAX_VERTS)
    real*8                :: x0(2), dir(2), step, tDir(2), bbox(2, 2), lengthScale
    integer               :: edx, vdx, vdx_first_inside, nrPos, vdx_next, nrPos_skelelton, rdx
    integer               :: verts_per_segment_, phase_
    logical               :: vdx_is_inside, vdx_next_is_inside, is_on_interface(2*R2D_MAX_VERTS)

    verts_per_segment_ = R2D_MAX_VERTS / 3
    if (present(verts_per_segment)) verts_per_segment_ = min(verts_per_segment_, verts_per_segment)

    phase_ = merge(phase, LIQUID_PHASE, present(phase))

    ! Find out which vertices of the cell are inside the domain
    vdx_first_inside = 0
    do vdx=1,cell%nverts
      funVals(vdx) = interfaceFun(cell%verts(vdx)%pos%xyz)
      if (funVals(vdx) < 0 .and. vdx_first_inside == 0) vdx_first_inside = vdx
    enddo
    if (vdx_first_inside == 0) then
      poly%nverts = 0
      return
    endif

    ! Loop over the edges and construct the polygonal 'skeleton'
    vdx = vdx_first_inside
    vdx_is_inside = .true.
    nrPos_skelelton = 0

    do edx=1,cell%nverts
      vdx_next = cell%verts(vdx)%pnbrs(1)+1
      vdx_next_is_inside = funVals(vdx_next) < 0.0

      ! TODO so far we assume that an edge has at most one intersection
      if (vdx_is_inside .neqv. vdx_next_is_inside) then
        ! Find and add new position
        x0 = cell%verts(vdx)%pos%xyz
        dir = cell%verts(vdx_next)%pos%xyz - x0
        step = brent(interfaceFun_step, 0.0D0, 1.0D0, 1D-15, 30, funVals(vdx), funVals(vdx_next))
        nrPos_skelelton = nrPos_skelelton + 1
        pos_skeleton(:,nrPos_skelelton) = x0 + step * dir
        is_on_interface(nrPos_skelelton) = .true.
      endif
      if (vdx_next_is_inside) then
        ! And add next node (corner)
        nrPos_skelelton = nrPos_skelelton + 1
        pos_skeleton(:,nrPos_skelelton) = cell%verts(vdx_next)%pos%xyz
        is_on_interface(nrPos_skelelton) = .false.
      endif

      vdx = vdx_next
      vdx_is_inside = vdx_next_is_inside
    enddo

    call bounding_box(bbox, cell)
    lengthScale = norm2(bbox(:,2) - bbox(:,1))

    ! Now we add a refined approximation on edges that are on the interface
    nrPos = 0
    vdx = 1
    do edx=1,nrPos_skelelton
      vdx_next = merge(1, vdx + 1, vdx == nrPos_skelelton)

      ! Add (refinement of) the half open interval (pos_skeleton(:,vdx),pos_skeleton(:,vdx_next)]
      if (.not. is_on_interface(vdx) .or. .not. is_on_interface(vdx_next)) then
        nrPos = nrPos + 1
        pos(:,nrPos) = pos_skeleton(:,vdx_next)
      else

        tDir = pos_skeleton(:,vdx_next) - pos_skeleton(:,vdx)
        if (norm2(tDir) < 1E-15 * lengthScale) then
          nrPos = nrPos + 1
          pos(:,nrPos) = pos_skeleton(:,vdx_next)
        else
          ! Refine the face

          ! Make dir normal to the face
          dir = [-tDir(2), tDir(1)]
          do rdx=1,verts_per_segment_
            x0 = pos_skeleton(:,vdx) + rdx * tDir / verts_per_segment_

            ! We impose here that the radius of curvature of the interface is bounded from below by half (relative to the mesh spacing)
            step = brent(interfaceFun_step, -.5D0, .5D0, 1D-15, 30)

            nrPos = nrPos + 1
            pos(:,nrPos) = x0 + step * dir
          enddo
        endif
      endif

      vdx = vdx_next
    enddo
    call init_from_pos(poly, pos(:,1:nrPos))
  contains

    real*8 function interfaceFun(x) result(f)
  
      implicit none

      real*8, intent(in)  :: x(2)

      f = levelSet(x)
      if (phase_ == GAS_PHASE) f = -f
    end

    real*8 function interfaceFun_step(step_) result(f)
      implicit none

      real*8, intent(in)    :: step_

      f = interfaceFun(x0 + step_ * dir)
    end
  end

end module
