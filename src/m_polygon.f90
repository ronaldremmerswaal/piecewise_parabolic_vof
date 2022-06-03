module m_polygon
  integer, parameter      :: MAX_NR_VERTS = 512

  type tPolygon
    real*8                :: verts(2,MAX_NR_VERTS)
    integer               :: nverts
  end type

  type tPlane
    real*8                :: normal(2), shift
  end type

  interface makeBox
    module procedure makeBox_bounds, makeBox_dx
  end interface
contains

  real*8 function cmpVolume(poly) result(vol)
    implicit none
    
    type(tPolygon), intent(in) :: poly

    ! Local variables
    integer               :: vdx, ndx

    vol = 0
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      vol = vol + poly%verts(1,vdx)*poly%verts(2,ndx) - poly%verts(2,vdx)*poly%verts(1,ndx)
    enddo

    vol = vol/2
  end function

  subroutine intersect(iPoly, poly, plane)
    implicit none
    
    type(tPolygon), intent(in) :: poly
    type(tPlane), intent(in) :: plane

    type(tPolygon), intent(out)        :: iPoly

    call copy(out=iPoly, in=poly)
    call intersect_(iPoly, plane)
  end subroutine

  subroutine intersect_(poly, plane)
    implicit none

    type(tPolygon), intent(inout) :: poly
    type(tPlane), intent(in) :: plane

    ! Local variables
    real*8                :: dist(MAX_NR_VERTS), buffer(2,MAX_NR_VERTS), coeff
    integer               :: edx, vdx, ndx, inside_count, first_inside, new_count

    ! Check for each vertex if it under or above the plane
    inside_count = 0
    first_inside = 0
    do vdx=1,poly%nverts
      dist(vdx) = plane%normal(1)*poly%verts(1,vdx) + plane%normal(2)*poly%verts(2,vdx) - plane%shift
      if (dist(vdx) <= 0) then 
        inside_count = inside_count + 1
        if (first_inside == 0) first_inside = vdx
      endif
    enddo
    if (inside_count==0) then
      poly%nverts = 0
      return
    elseif (inside_count==poly%nverts) then
      return
    endif

    ! Store old vertices
    buffer(:,1:poly%nverts) = poly%verts(:,1:poly%nverts)

    ! Loop over edges of old polygon: insert new vertices and keep some old vertices
    new_count = 0
    vdx = first_inside
    do edx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      if (dist(vdx)<=0 .neqv. dist(ndx)<=0) then
        ! Insert new
        new_count = new_count + 1
        coeff = abs(dist(ndx) / (dist(ndx) - dist(vdx)))

        poly%verts(:,new_count) = coeff * buffer(:,vdx) + (1 - coeff) * buffer(:,ndx)
      endif

      if (dist(ndx)<=0) then
        ! Keep old
        new_count = new_count + 1
        poly%verts(:,new_count) = buffer(:,ndx)
      endif

      vdx = ndx
    enddo

    poly%nverts = new_count
  end subroutine

  subroutine makeBox_bounds(poly, rbounds)
    implicit none

    real*8, intent(in)    :: rbounds(2, 2)
    type(tPolygon), intent(out) :: poly

    poly%verts(:,1) = [rbounds(1,1), rbounds(2,1)]
    poly%verts(:,2) = [rbounds(1,2), rbounds(2,1)]
    poly%verts(:,3) = [rbounds(1,2), rbounds(2,2)]
    poly%verts(:,4) = [rbounds(1,1), rbounds(2,2)]

    poly%nverts = 4

  end subroutine

  subroutine makeBox_dx(poly, dx)
    implicit none

    real*8, intent(in)    :: dx(2)
    type(tPolygon), intent(out) :: poly

    poly%verts(:,1) = [-dx(1), -dx(2)]
    poly%verts(:,2) = [dx(1), -dx(2)]
    poly%verts(:,3) = [dx(1), dx(2)]
    poly%verts(:,4) = [-dx(1), dx(2)]

    poly%nverts = 4

  end subroutine

  ! NOTE: returning poly by function is much more expensive than returning via subroutine!
  subroutine copy(out, in)
    implicit none
    
    type(tPolygon), intent(in) :: in
    type(tPolygon), intent(out):: out

    out%nverts = in%nverts
    out%verts(:,out%nverts) = in%verts(:,out%nverts)
  end subroutine
end