module m_polygon
  integer, parameter      :: MAX_NR_VERTS = 2**6
  integer, parameter      :: MAX_MONOMIAL = 5

  type tParabola
    real*8                :: normal(2), shift
    real*8                :: kappa0 = 0.
  end type

  type tPolygon
    real*8                :: verts(2,MAX_NR_VERTS)
    integer               :: nverts = 0

    logical               :: intersected = .false. ! NB intersected by a parabola
    logical               :: on_parabola(MAX_NR_VERTS)
    type(tParabola)       :: parabola

    real*8                :: x_tau(MAX_NR_VERTS)
    real*8                :: x_eta(MAX_NR_VERTS)
    
    integer               :: avail_monomial = -1
    real*8                :: x_tau_power(MAX_NR_VERTS)
    real*8                :: monomials(0:MAX_MONOMIAL)
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

  function cmpFirstMoment(poly) result(mom)
    implicit none
    
    type(tPolygon), intent(in) :: poly
    real*8                :: mom(2), areaTerm

    ! Local variables
    integer               :: vdx, ndx

    mom = 0
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      areaTerm = poly%verts(1,vdx)*poly%verts(2,ndx) - poly%verts(2,vdx)*poly%verts(1,ndx)
      mom = mom + areaTerm * (poly%verts(:,vdx) + poly%verts(:,ndx))
    enddo

    mom = mom/6
  end function

  function cmpMoments(poly) result(mom)
    implicit none
    
    type(tPolygon), intent(in) :: poly
    real*8                :: mom(3), areaTerm

    ! Local variables
    integer               :: vdx, ndx

    mom = 0
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      areaTerm = poly%verts(1,vdx)*poly%verts(2,ndx) - poly%verts(2,vdx)*poly%verts(1,ndx)
      mom(1) = mom(1) + areaTerm
      mom(2:3) = mom(2:3) + areaTerm * (poly%verts(:,vdx) + poly%verts(:,ndx))
    enddo

    mom(1) = mom(1)/2
    mom(2:3) = mom(2:3)/6
  end function

  subroutine intersect(iPoly, poly, parabola)
    implicit none
    
    type(tPolygon), intent(in) :: poly
    type(tParabola), intent(in) :: parabola

    type(tPolygon), intent(out)        :: iPoly

    call copy(out=iPoly, in=poly)
    call intersect_(iPoly, parabola)
  end subroutine

  subroutine intersect_(poly, parabola)
    use m_poly_roots,     only: real_roots_2

    implicit none

    type(tPolygon), intent(inout) :: poly
    type(tParabola), intent(in) :: parabola

    ! Local variables
    real*8                :: dist(MAX_NR_VERTS), buffer(2,MAX_NR_VERTS), coeff
    integer               :: edx, vdx, ndx, inside_count, first_inside, new_count
    logical               :: is_parabolic
    real*8                :: x_eta, x_tau

    is_parabolic = parabola%kappa0 /= 0.0

    ! Check for each vertex if it under or above the parabola
    inside_count = 0
    first_inside = 0
    do vdx=1,poly%nverts
      x_eta = parabola%normal(1)*poly%verts(1,vdx) + parabola%normal(2)*poly%verts(2,vdx)
      dist(vdx) = x_eta - parabola%shift
      if (is_parabolic) then
        x_tau = -parabola%normal(2)*poly%verts(1,vdx) + parabola%normal(1)*poly%verts(2,vdx)
        dist(vdx) = dist(vdx) + (parabola%kappa0/2) * x_tau**2
      endif

      if (dist(vdx) <= 0) then 
        inside_count = inside_count + 1
        if (first_inside == 0) first_inside = vdx
      endif
    enddo
    if (.not. is_parabolic) then
      if (inside_count==0) then
        poly%nverts = 0
        return
      elseif (inside_count==poly%nverts) then
        return
      endif
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
        ! Insert new vertex (edge bisection)
        new_count = new_count + 1
        if (.not. is_parabolic) then
          coeff = abs(dist(ndx) / (dist(ndx) - dist(vdx)))
        else

        endif

        poly%verts(:,new_count) = coeff * buffer(:,vdx) + (1 - coeff) * buffer(:,ndx)
      elseif (is_parabolic) then
        ! Consider if the edge is trisected
        if (dist(vdx)<=0) then !  .and. dist(ndx)<=0

        else ! if (dist(vdx)>0 .and. dist(ndx)>0) then

        endif
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

  subroutine compute_momonial(poly, nr)
    implicit none
    
    type(tPolygon), intent(inout) :: poly
    integer, intent(in)   :: nr

    ! Local variables
    integer               :: mdx, vdx, ndx, old_nr
    
    old_nr = poly%avail_monomial
    if (nr <= old_nr) return
    if (nr > MAX_MONOMIAL) then
      print*, 'ERROR in compute_momonial: too many monomials requested'
    endif

    if (old_nr<0) old_nr = -1

    ! We compute integral of x_tau^mdx over edges which are parabola
    do mdx=old_nr+1,nr
      if (mdx==0) then
        do vdx=1,poly%nverts
          if (poly%on_parabola(vdx)) poly%x_tau_power(vdx) = poly%x_tau(vdx)
        enddo
      else
        do vdx=1,poly%nverts
          if (poly%on_parabola(vdx)) poly%x_tau_power(vdx) = poly%x_tau_power(vdx) * poly%x_tau(vdx)
        enddo
      endif

      poly%monomials(mdx) = 0
      do vdx=1,poly%nverts
        ndx = vdx + 1
        if (vdx==poly%nverts) ndx = 1

        if (poly%on_parabola(vdx) .and. poly%on_parabola(ndx)) then
          poly%monomials(mdx) = poly%monomials(mdx) + poly%x_tau_power(ndx) - poly%x_tau_power(vdx)
        endif
      enddo
    enddo

    poly%avail_monomial = nr
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
    out%verts(:,1:out%nverts) = in%verts(:,1:out%nverts)
    
    out%intersected = in%intersected
    if (in%intersected) then
      out%parabola%normal = in%parabola%normal
      out%parabola%shift = in%parabola%shift
      out%parabola%kappa0 = in%parabola%kappa0

      out%on_parabola(1:out%nverts) = in%on_parabola(1:out%nverts)

      out%x_tau(1:out%nverts) = in%x_tau(1:out%nverts)
      out%x_eta(1:out%nverts) = in%x_eta(1:out%nverts)
      
      out%avail_monomial = in%avail_monomial
      if (out%avail_monomial>=0) out%x_tau_power(1:out%nverts) = in%x_tau_power(1:out%nverts)
      out%monomials = in%monomials
    endif

  end subroutine
end