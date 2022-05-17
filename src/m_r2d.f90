module m_r2d
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


  interface cmpMoments
    module procedure reduce_1, cmpMoments_poly_plane
  end interface

  interface cmpMoments_
    module procedure cmpMoments_poly_plane_
  end interface

  interface makeBox
    module procedure makeBox_bounds, makeBox_xdx
  end interface
contains

  function makePlane(normal, shift) result(plane)
    implicit none
    
    real*8                :: normal(2), shift
    type(r2d_plane_f)     :: plane

    plane%n%xyz = normal
    plane%d = shift
  end function

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

  function reduce_1(poly) result(mom)
    implicit none

    type(r2d_poly_f), intent(in) :: poly
    real*8                :: mom(3)

    call r2d_reduce(poly, mom, %val(1))
  end function

  function cmpMoments_poly_plane(poly, plane) result(mom)
    implicit none

    type(r2d_poly_f), intent(in) :: poly
    type(r2d_plane_f), intent(in) :: plane
    real*8                :: mom(3)

    ! Local variables
    type(r2d_poly_f)      :: poly_

    poly_ = copy(poly)
    call intersect(poly_, plane) 
    mom = reduce_1(poly_)
  end function

  function cmpMoments_poly_plane_(poly, plane) result(mom)
    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    type(r2d_plane_f), intent(in) :: plane
    real*8                :: mom(3)

    call intersect(poly, plane)
    mom = reduce_1(poly)
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
  subroutine intersect(poly, plane)
    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    type(r2d_plane_f), intent(in)  :: plane

    type(r2d_plane_f)     :: planes(1)

    planes(1)%n%xyz = -plane%n%xyz
    planes(1)%d = plane%d

    ! r2d convention: remove part for which η . x + s < 0
    call r2d_clip(poly, planes, %val(1))
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

  function copy(from) result(to)
    implicit none

    type(r2d_poly_f), intent(in) :: from
    type(r2d_poly_f)      :: to

    integer               :: v

    to%nverts = from%nverts
    do v=1,to%nverts
      to%verts(v)%pnbrs = from%verts(v)%pnbrs
      to%verts(v)%pos%xyz = from%verts(v)%pos%xyz
    enddo

  end function

  function makeBox_bounds(rbounds) result(poly)
    implicit none

    real*8, intent(in)    :: rbounds(2, 2)
    type(r2d_poly_f)      :: poly

    ! Local variables
    type(r2d_rvec2_f)     :: rbounds_vec(2)

    rbounds_vec(1)%xyz = rbounds(:, 1) 
    rbounds_vec(2)%xyz = rbounds(:, 2)

    call r2d_init_box_f(poly, rbounds_vec)
  end function

  function makeBox_xdx(x, dx) result(poly)
    implicit none
    
    real*8, intent(in)    :: x(2), dx(2)
    type(r2d_poly_f)      :: poly

        ! Local variables
    type(r2d_rvec2_f)     :: rbounds_vec(2)

    rbounds_vec(1)%xyz = x - dx/2
    rbounds_vec(2)%xyz = x + dx/2

    call r2d_init_box_f(poly, rbounds_vec)
  end function  

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

  subroutine shift_by(poly, dpos)
    implicit none

    type(r2d_poly_f), intent(inout) :: poly
    real*8, intent(in)    :: dpos(2)

    ! Local variables
    integer               :: vdx

    do vdx=1,poly%nverts
      poly%verts(vdx)%pos%xyz = poly%verts(vdx)%pos%xyz + dpos
    enddo
  end subroutine


  subroutine print(poly)
    implicit none

    type(r2d_poly_f), intent(in) :: poly

    call r2d_print(poly)
  end
end module
