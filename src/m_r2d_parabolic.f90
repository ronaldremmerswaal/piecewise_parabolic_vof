module m_r2d_parabolic
  use, intrinsic :: iso_c_binding

  implicit none

  integer, parameter      :: R2D_MAX_VERTS = 256

  type, bind(C) :: r2d_rvec2_f
    real(c_double)        :: xyz(2)
  end type

  type, bind(C) :: r2d_plane_f
    type(r2d_rvec2_f)     :: n
    real(c_double)        :: d
  end type

  type, bind(C) :: r2d_parabola_f
    type(r2d_rvec2_f)     :: n
    real(c_double)        :: kappa0
  end type

  type, bind(C) :: r2d_vertex_f
    integer(c_int32_t)    :: pnbrs(2)
    type(r2d_rvec2_f)     :: pos
  end type

  type, bind(C) :: r2d_poly_f
    type(r2d_vertex_f)    :: verts(R2D_MAX_VERTS)
    integer(c_int32_t)    :: nverts
  end type

  interface
    subroutine r2d_init_box_f(poly, rbounds) bind(C, name="r2d_init_box")
      import r2d_poly_f, r2d_rvec2_f

      implicit none

      type(r2d_poly_f), intent(out) :: poly
      type(r2d_rvec2_f), intent(in) :: rbounds(2)
    end subroutine r2d_init_box_f
  end interface

  interface
    subroutine r2d_init_poly_f(poly, vertices, numverts) bind(C, name="r2d_init_poly")
      import r2d_poly_f, r2d_rvec2_f, c_int32_t

      implicit none

      type(r2d_poly_f), intent(out) :: poly
      type(r2d_rvec2_f), intent(in) :: vertices(numverts)
      integer(c_int32_t), intent(in):: numverts
    end subroutine r2d_init_poly_f
  end interface

  interface
    subroutine r2d_print(poly) bind(C, name="r2d_print")
      import r2d_poly_f

      implicit none

      type(r2d_poly_f), intent(in):: poly
    end subroutine r2d_print
  end interface

  interface
    subroutine r2d_volume(poly, volume) bind(C, name="r2d_volume")
      import r2d_poly_f, c_double, c_int32_t

      implicit none

      type(r2d_poly_f), intent(in):: poly
      real(c_double), intent(out) :: volume
    end subroutine r2d_volume
  end interface

  interface
    subroutine r2d_reduce(poly, moments, order) bind(C, name="r2d_reduce")
      import r2d_poly_f, c_double, c_int32_t

      implicit none

      type(r2d_poly_f), intent(in):: poly
      real(c_double), intent(out) :: moments(3)
      integer(c_int32_t), intent(in) :: order
    end subroutine r2d_reduce
  end interface

  interface
    subroutine r2d_clip(poly, planes, nplanes) bind(C, name="r2d_clip")
      import r2d_poly_f, r2d_plane_f, c_int32_t

      implicit none

      type(r2d_poly_f), intent(inout) :: poly
      type(r2d_plane_f), intent(in)   :: planes(nplanes)
      integer(c_int32_t), intent(in)  :: nplanes
    end subroutine r2d_clip
  end interface

  interface
    subroutine r2d_split_ptr(polys, npolys, plane, out_pos, out_neg) bind(C, name="r2d_split_ptr")
      import r2d_poly_f, r2d_plane_f, c_int32_t

      implicit none

      type(r2d_poly_f), intent(inout) :: polys(npolys)
      integer(c_int32_t), intent(in)  :: npolys
      type(r2d_plane_f), intent(in)   :: plane
      type(r2d_poly_f), intent(out)   :: out_pos(npolys)
      type(r2d_poly_f), intent(out)   :: out_neg(npolys)
    end subroutine r2d_split_ptr
  end interface

  interface

    subroutine r2d_clip_parabola_moments_01(poly, parabola, grad_s, moments, derivative, &
      compute_derivative) bind(C, name="r2d_clip_parabola_moments_01")
      import r2d_poly_f, r2d_parabola_f, c_double, c_bool

      implicit none

      type(r2d_poly_f), intent(inout) :: poly
      type(r2d_parabola_f), intent(in) :: parabola(1)
      real(c_double), intent(inout):: grad_s(2)
      real(c_double), intent(out)  :: moments(3), derivative(4)
      logical(c_bool), intent(in)  :: compute_derivative
    end subroutine r2d_clip_parabola_moments_01
  end interface
contains

  subroutine bounding_box(box, poly)
    implicit none

    type(r2d_poly_f), intent(in) :: poly
    real*8, intent(out)   :: box(2, 2)

    ! Local variables
    integer               :: v, dim

    if (poly%nverts == 0) then
      box = 0.0
      return
    endif

    box(:, 1) = poly%verts(1)%pos%xyz
    box(:, 2) = poly%verts(1)%pos%xyz
    do v=2,poly%nverts
      do dim=1,2
        box(dim, 1) = min(box(dim, 1), poly%verts(v)%pos%xyz(dim))
        box(dim, 2) = max(box(dim, 2), poly%verts(v)%pos%xyz(dim))
      enddo
    enddo
  end subroutine

  function moments_01(poly) result(mom)
    implicit none

    type(r2d_poly_f), intent(in) :: poly
    real*8                :: mom(3)

    call r2d_reduce(poly, mom, %val(1))
  end function

  subroutine intersect_with_cell(poly, x0, dx)
    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    real*8, intent(in)    :: x0(2), dx(2)

    ! Local variables
    type(r2d_plane_f)     :: planes(2, 2)
    integer               :: dim, side
    real*8                :: sgn

    do dim=1,2
    do side=1,2
      sgn = 2 * side - 3.0D0

      planes(dim, side)%n%xyz = 0.0D0
      planes(dim, side)%n%xyz(dim) = sgn
      planes(dim, side)%d = -sgn * x0(dim) + dx(dim) / 2
    enddo
    enddo

    call r2d_clip(poly, planes, %val(4))
  end subroutine

  ! Remove part on outward normal side of plane defined by η ⋅ x = s, or η ⋅ x - s > 0
  subroutine intersect_by_plane(poly, normal, shift)
    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    real*8, intent(in)    :: normal(2), shift

    type(r2d_plane_f)     :: plane(1)

    plane(1)%n%xyz = -normal
    plane(1)%d = shift

    ! r2d convention: remove part for which η . x + s < 0
    call r2d_clip(poly, plane, %val(1))
  end subroutine

  subroutine split_by_plane(polys, normal, shift, out_pos, out_neg)
    implicit none

    type(r2d_poly_f), intent(inout) :: polys(1:)
    real*8, intent(in)    :: normal(2), shift
    type(r2d_poly_f), intent(out) :: out_pos(1:), out_neg(1:)

    type(r2d_plane_f)     :: plane

    plane%n%xyz = -normal
    plane%d = shift

    ! r2d convention: remove part for which η . x + s < 0
    call r2d_split_ptr(polys, %val(size(polys, 1)), plane, out_pos, out_neg)
  end subroutine

  ! Removes the part of poly for which 
  !   \eta \cdot (x - x_0) + (\kappa/2) * (\tau \cdot (x - x_0))^2 > 0
  ! and returns the zeroth (M_0) and first (M_1) moments
  !
  ! IN/OUTPUT
  ! 
  ! moments01             The zeroth and first moments of intersection (so moments01(1) is the 
  !                         volume and moments01(2:3)/moments01(1) is the centroid of the
  !                         intersection volume)
  ! poly                  Polygon which is to be intersected (is modified!)
  ! normal                The normal \eta of the interface at x = x_0 (\eta = [cos(\varphi), sin(\varphi)])
  ! kappa0                The curvature \kappa of the interface at x = x_0
  ! x0                    The position x_0 on the interface
  ! derivative            The derivatives of the zeroth and first moments (optional):  
  !                         derivative(1)   = d M_0/ d \varphi
  !                         derivative(2:3) = d M_1/ d \varphi
  !                         derivative(4)   = d M_0/ d \kappa
  ! grad_s                The gradient of the shift s w.r.t. the normal angle (\varphi) and curvature (\kappa)
  !                         If initially NaN, then this gradient is computed such that d M_0/ d \varphi = 0
  !                         If not, then it is assumed that the gradient of the shift s is given by grad_s
  !                         (This only affects the computation of the gradient)
  subroutine intersect_with_parabola(moments01, poly, normal, kappa0, x0, derivative, grad_s)
    use m_common

    implicit none

    real*8, intent(out)   :: moments01(3)
    type(r2d_poly_f), intent(inout) :: poly
    real*8, intent(in)    :: normal(2), kappa0, x0(2)
    real*8, intent(inout), optional :: derivative(4), grad_s(2)

    ! Local variables
    type(r2d_parabola_f)  :: parabola(1)
    real*8                :: moments01_poly(3)
    real*8                :: derivative_(4), grad_s_(2)
    logical*1             :: compute_derivative_
    integer               :: vdx

    ! The algorithm is implemented assuming that x0 = 0, so we must shift the positions of the polygon beforehand
    do vdx=1,poly%nverts
      poly%verts(vdx)%pos%xyz = poly%verts(vdx)%pos%xyz - x0
    enddo

    ! The algorithm is implemented for kappa0 <= 0, so we must compute the moments of the complement
    ! (for which the sign of the curvature is swapped) if kappa0 > 0
    if (kappa0 > 0) then
      parabola(1)%n%xyz = -normal
      parabola(1)%kappa0 = -kappa0
      moments01_poly = moments_01(poly)
    else
      parabola(1)%n%xyz = normal
      parabola(1)%kappa0 = kappa0
    endif

    grad_s_ = merge(grad_s, [d_qnan, d_qnan], present(grad_s))

    compute_derivative_ = present(derivative) .or. present(grad_s)
    call r2d_clip_parabola_moments_01(poly, parabola, grad_s_, moments01, derivative_, compute_derivative_)

    if (kappa0 > 0) then
      moments01 = moments01_poly - moments01
      derivative_(1:3) = -derivative_(1:3)
    endif

    ! Shift the position back relative to x0
    moments01(2:3) = moments01(2:3) + moments01(1) * x0

    if (present(derivative)) then 
      derivative = derivative_
      derivative(2:3) = derivative(2:3) + derivative(1) * x0
    endif
    if (present(grad_s)) grad_s = grad_s_
  end subroutine

  subroutine copy(to, from)
    implicit none

    type(r2d_poly_f), intent(in) :: from
    type(r2d_poly_f), intent(out) :: to

    integer               :: v

    to%nverts = from%nverts
    do v=1,to%nverts
      to%verts(v)%pnbrs = from%verts(v)%pnbrs
      to%verts(v)%pos%xyz = from%verts(v)%pos%xyz
    enddo

  end subroutine

  subroutine init_box(poly, rbounds)
    implicit none

    type(r2d_poly_f), intent(out) :: poly
    real*8, intent(in)    :: rbounds(2, 2)

    ! Local variables
    type(r2d_rvec2_f)     :: rbounds_vec(2)

    rbounds_vec(1)%xyz = rbounds(:, 1)
    rbounds_vec(2)%xyz = rbounds(:, 2)

    call r2d_init_box_f(poly, rbounds_vec)

  end subroutine

  subroutine init_from_pos(poly, pos)
    implicit none

    type(r2d_poly_f), intent(out) :: poly
    real*8, intent(in)    :: pos(1:, 1:)

    ! Local variables
    type(r2d_rvec2_f)     :: vertices(size(pos, 2))
    integer               :: vdx, nverts

    nverts = size(pos, 2)
    if (nverts > R2D_MAX_VERTS .or. size(pos, 1) /= 2) return

    do vdx=1,nverts
      vertices(vdx)%xyz = pos(:, vdx)
    enddo
    call r2d_init_poly_f(poly, vertices, %val(nverts))
  end subroutine


  subroutine print(poly)
    implicit none

    type(r2d_poly_f), intent(in) :: poly

    call r2d_print(poly)
  end

  function reference_moments(phase, x, dx, levelSet, verts_per_segment) result(moments)
    implicit none
    
    integer*2, intent(in) :: phase
    real*8, intent(in)    :: x(2), dx(2)
    real*8, external      :: levelSet
    real*8                :: moments(3)
    integer, intent(in), optional :: verts_per_segment

    ! Local variables
    type(r2d_poly_f)      :: poly
    
    call get_polygonal_approximation_of_exact_domain(poly, phase, x, dx, levelSet, verts_per_segment)
    moments = moments_01(poly)
  end function

  subroutine get_polygonal_approximation_of_exact_domain(poly, phase, x, dx, levelSet, verts_per_segment)
    use m_common
    use m_optimization,   only: brent
  
    implicit none

    type(r2d_poly_f), intent(out) :: poly
    integer*2, intent(in) :: phase
    real*8, intent(in)    :: x(2), dx(2)
    real*8, external      :: levelSet
    integer, intent(in), optional :: verts_per_segment

    ! Local variables
    real*8                :: pos(2, R2D_MAX_VERTS), pos_skeleton(2, 8), corners(2, 4), funVals(4)
    real*8                :: x0(2), dir(2), step, tDir(2)
    integer               :: edx, vdx, vdx_first_inside, nrPos, vdx_next, nrPos_skelelton, rdx
    integer               :: verts_per_segment_
    logical               :: vdx_is_inside, vdx_next_is_inside, is_on_interface(8)

    verts_per_segment_ = R2D_MAX_VERTS / 3
    if (present(verts_per_segment)) verts_per_segment_ = min(verts_per_segment_, verts_per_segment)

    corners(:,1) = x + [-dx(1), -dx(2)]/2
    corners(:,2) = x + [dx(1), -dx(2)]/2
    corners(:,3) = x + [dx(1), dx(2)]/2
    corners(:,4) = x + [-dx(1), dx(2)]/2

    ! Find out which corners of the cell are inside the domain
    vdx_first_inside = 0
    do vdx=1,4
      funVals(vdx) = interfaceFun(corners(:,vdx))
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

    do edx=1,4
      vdx_next = merge(1, vdx + 1, vdx == 4)
      vdx_next_is_inside = funVals(vdx_next) < 0.0

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
      if (.not. is_on_interface(vdx) .or. .not. is_on_interface(vdx_next)) then
        nrPos = nrPos + 1
        pos(:,nrPos) = pos_skeleton(:,vdx_next)
      else

        tDir = pos_skeleton(:,vdx_next) - pos_skeleton(:,vdx)
        if (norm2(tDir) < 1E-15 * dx(1)) then
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
      if (phase == GAS_PHASE) f = -f
    end

    real*8 function interfaceFun_step(step_) result(f)
      implicit none

      real*8, intent(in)    :: step_

      f = interfaceFun(x0 + step_ * dir)
    end
  end

end module
