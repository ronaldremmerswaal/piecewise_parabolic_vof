module m_polygon
  integer, parameter      :: MAX_NR_VERTS = 2**7
  integer, parameter      :: MAX_NR_PARA_EDGES = 2**2
  integer, parameter      :: MAX_MONOMIAL = 5

  type tParabola
    real*8                :: normal(2), shift
    real*8                :: kappa0 = 0
  end type

  type tPolygon
    ! Polygon is stored as a list of positions
    real*8                :: verts(2,MAX_NR_VERTS)
    integer               :: nverts = 0

    ! When a polygon is intersected with a parabola we store the corresponding parabola
    ! which is needed for computing moments/derivatives
    logical               :: intersected = .false.
    logical               :: parabolic = .false.
    type(tParabola)       :: parabola

    ! For each vertex we store whether or not it was part of the original polygon
    logical               :: on_parabola(MAX_NR_VERTS)

    ! We compute the complement if the parabola has kappa0 > 0, in which case the moments
    ! of the original polygon must be stored
    logical               :: complement = .false. 
    real*8                :: original_moments(3)

    ! The tangential and normal coordinates of the vertices which lie on the parabola
    ! (needed for computation of moments/derivatives)
    real*8                :: x_tau(2,MAX_NR_PARA_EDGES)
    real*8                :: x_eta(2,MAX_NR_PARA_EDGES)
    
    ! Integrated powers of x_tau, which are needed for computation of moments/derivatives
    integer               :: avail_monomial = -1
    real*8                :: x_tau_power(2,MAX_NR_PARA_EDGES)
    real*8                :: monomials(0:MAX_MONOMIAL,MAX_NR_PARA_EDGES)
    real*8                :: monomials_sum(0:MAX_MONOMIAL)
  end type


  interface makeBox
    module procedure makeBox_dx, makeBox_xdx
  end interface

  interface cmpVolume
    module procedure cmpVolume_poly
  end interface

  interface cmpMoments
    module procedure cmpMoments_poly
  end interface

  interface makeParabola
    module procedure makeParabola_poly, makeParabola_angle_poly
  end interface

  interface makePlane
    module procedure makePlane_def, makePlane_angle
  end interface
contains
  function makePlane_def(normal, shift) result(plane)
    implicit none
    
    real*8, intent(in)    :: normal(2), shift
    type(tParabola)       :: plane

    plane%normal = normal
    plane%shift = shift
  end function

  function makePlane_angle(angle, shift) result(plane)
    implicit none
    
    real*8, intent(in)    :: angle, shift
    type(tParabola)       :: plane

    plane%normal = [dcos(angle), dsin(angle)]
    plane%shift = shift
  end function

  function makeParabola_poly(normal, kappa0, shift) result(parabola)
    implicit none
    
    real*8, intent(in)    :: normal(2), kappa0, shift
    type(tParabola)       :: parabola

    parabola%normal = normal
    parabola%kappa0 = kappa0
    parabola%shift = shift
  end function

  function makeParabola_angle_poly(angle, kappa0, shift) result(parabola)
    implicit none
    
    real*8, intent(in)    :: angle, kappa0, shift
    type(tParabola)       :: parabola

    parabola%normal = [dcos(angle), dsin(angle)]
    parabola%kappa0 = kappa0
    parabola%shift = shift
  end function

  function complement(parabola) result(c_parabola)
    implicit none
    
    type(tParabola), intent(in) :: parabola
    type(tParabola)       :: c_parabola

    c_parabola%normal = -parabola%normal
    c_parabola%kappa0 = -parabola%kappa0
    c_parabola%shift = -parabola%shift
  end function
    

  real*8 function cmpVolume_poly(poly) result(vol)
    implicit none
    
    type(tPolygon), intent(inout) :: poly ! inout because volume correction may update monomials and taupower

    ! Local variables
    integer               :: vdx, ndx

    vol = 0
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      vol = vol + poly%verts(1,vdx)*poly%verts(2,ndx) - poly%verts(2,vdx)*poly%verts(1,ndx)
    enddo

    vol = vol/2

    if (poly%intersected .and. poly%parabolic) then
      vol = vol + parabola_volume_correction(poly)
      if (poly%complement) then
        vol = poly%original_moments(1) - vol
      endif
    endif
  end function

  function cmpMoments_poly(poly) result(mom)
    implicit none
    
    type(tPolygon), intent(inout) :: poly
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

    if (poly%intersected .and. poly%parabolic) then
      mom = mom + parabola_moments_correction(poly)
      if (poly%complement) then
        mom = poly%original_moments - mom
      endif
    endif
  end function

  subroutine cmpMoments_poly_SUB(mom, poly)
    implicit none
    
    type(tPolygon), intent(inout) :: poly
    real*8, intent(out)           :: mom(3)

    ! Local variables
    integer               :: vdx, ndx
    real*8                :: areaTerm

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

    if (poly%intersected .and. poly%parabolic) then
      mom = mom + parabola_moments_correction(poly)
      if (poly%complement) then
        mom = poly%original_moments - mom
      endif
    endif
  end subroutine

  real*8 function cmpDerivative_volAngle(poly, shiftAngleDerivative) result(der)
    implicit none

    type(tPolygon), intent(inout) :: poly
    real*8, intent(in), optional :: shiftAngleDerivative

    ! Local variables
    real*8                :: shiftAngleDerivative_    

    if (.not. poly%intersected) then
      ! print*, 'ERROR: derivative cannot be computed because polygon was not yet intersected'
      der = 0
      return
    endif
    
    if (.not. present(shiftAngleDerivative)) then 
      der = 0
      ! Force zero volume, so derivative is always zero
      return
      ! shiftAngleDerivative_ = cmpDerivative_shiftAngle(poly)
    endif

    shiftAngleDerivative_ = shiftAngleDerivative
    if (poly%complement) shiftAngleDerivative_ = -shiftAngleDerivative_
    
    if (poly%parabolic) then
      call compute_momonial(poly, 3)
    else
      call compute_momonial(poly, 1)
    endif

    der = poly%monomials_sum(0) * shiftAngleDerivative_ - poly%monomials_sum(1)
    if (poly%parabolic) then
      der = der + poly%monomials_sum(1) * poly%parabola%shift * poly%parabola%kappa0 - &
        (poly%monomials_sum(3) * poly%parabola%kappa0**2) / 2
    endif
    
    if (poly%complement) der = -der
  end function
  
  real*8 function cmpDerivative_volKappa(poly, shiftKappaDerivative) result(der)
    implicit none

    type(tPolygon), intent(inout) :: poly
    real*8, intent(in), optional :: shiftKappaDerivative

    if (.not. poly%intersected) then
      der = 0
      ! print*, 'ERROR: derivative cannot be computed because polygon was not yet intersected'
      return
    endif
    
    if (.not. present(shiftKappaDerivative)) then 
      ! Force zero volume
      der = 0
      return
      ! shiftKappaDerivative = cmpDerivative_shiftKappa(poly)
    endif
    
    call compute_momonial(poly, 2)

    der = poly%monomials_sum(0) * shiftKappaDerivative - poly%monomials_sum(2)/2
    
  end function
  
  function cmpDerivative_firstMomentAngle(poly, shiftAngleDerivative) result(der)
    implicit none

    type(tPolygon), intent(inout) :: poly
    real*8, intent(in), optional :: shiftAngleDerivative
    real*8                :: der(2)

    ! Local variables
    real*8                :: shiftAngleDerivative_, der1_eta, der1_tau

    if (.not. poly%intersected) then
      der = 0
      ! print*, 'ERROR: derivative cannot be computed because polygon was not yet intersected'
      return
    endif
    
    if (present(shiftAngleDerivative)) then 
      shiftAngleDerivative_ = shiftAngleDerivative
    else
      ! Force zero volume
      shiftAngleDerivative_ = cmpDerivative_shiftAngle(poly)
    endif
    
    if (poly%parabolic) then
      call compute_momonial(poly, 5)
    else
      call compute_momonial(poly, 2)
    endif

    if (poly%complement) shiftAngleDerivative_ = -shiftAngleDerivative_

    ! we write the derivative in terms of its normal and tangential components
	  ! the normal component is given by
	  ! 	int [grad_s - τ + κ τ(s - κ/2 τ^2)][s - κ/2 τ^2] dτ
    der1_eta = poly%parabola%shift * shiftAngleDerivative_ * poly%monomials_sum(0) - poly%parabola%shift * poly%monomials_sum(1)
    if (poly%parabolic) then
      der1_eta = der1_eta + poly%monomials_sum(3) * poly%parabola%kappa0 / 2 - &
        poly%monomials_sum(2) * shiftAngleDerivative_ * poly%parabola%kappa0 / 2 + &
        poly%parabola%kappa0 * (poly%monomials_sum(1) * poly%parabola%shift**2 + &
        poly%monomials_sum(5) * poly%parabola%kappa0**2 / 4 - &
        poly%parabola%shift * poly%parabola%kappa0 * poly%monomials_sum(3))
    endif

    ! and the tangential component
    ! 	int [grad_s - τ + κ τ(s - κ/2 τ^2)]τ dτ
    der1_tau = shiftAngleDerivative_ * poly%monomials_sum(1) - poly%monomials_sum(2)
    if (poly%parabolic) then
      der1_tau = der1_tau + poly%parabola%kappa0 * poly%parabola%shift * poly%monomials_sum(2) -&
        poly%monomials_sum(4) * poly%parabola%kappa0**2 / 2
    endif
    
    der(1) = poly%parabola%normal(1) * der1_eta - poly%parabola%normal(2) * der1_tau
    der(2) = poly%parabola%normal(2) * der1_eta + poly%parabola%normal(1) * der1_tau

    if (poly%complement) der = -der
  end function

    ! The derivative of the shift w.r.t. to the normal angle is given by
	! 	-int [(s - κ/2 τ^2) κ τ - τ] dτ / int 1 dτ
  real*8 function cmpDerivative_shiftAngle(poly) result(der)
    implicit none

    type(tPolygon), intent(inout) :: poly

    if (.not. poly%intersected) then
      der = 0
      ! print*, 'ERROR: derivative cannot be computed because polygon was not yet intersected'
      return
    endif

    if (poly%parabolic) then
      call compute_momonial(poly, 3)
    else
      call compute_momonial(poly, 1)
    endif

    if (poly%monomials_sum(0)==0) then
      der = 0
      return
    endif

    der = poly%monomials_sum(1)
    if (poly%parabolic) then
      der = der - poly%monomials_sum(1) * poly%parabola%shift * poly%parabola%kappa0 &
        + (poly%monomials_sum(3) * poly%parabola%kappa0**2) / 2
    endif
    der = der / poly%monomials_sum(0)

    if (poly%complement) der = -der
  end function

  ! The derivative of the shift w.r.t. to the curvature is given by
	! 	int [τ^2/2] dτ / int 1 dτ
  real*8 function cmpDerivative_shiftKappa(poly) result(der)
    implicit none

    type(tPolygon), intent(inout) :: poly

    if (.not. poly%intersected) then
      der = 0
      ! print*, 'ERROR: derivative cannot be computed because polygon was not yet intersected'
      return
    endif

    call compute_momonial(poly, 2)

    if (poly%monomials_sum(0)==0) then
      der = 0
      return
    endif

    der = (poly%monomials_sum(2)/2) / poly%monomials_sum(0)
  end function

  subroutine intersect(poly, parabola)

    implicit none

    type(tPolygon), intent(inout) :: poly
    type(tParabola), intent(in) :: parabola

    ! Local variables
    real*8                :: dist(MAX_NR_VERTS), buffer(2,MAX_NR_VERTS), coeffs(2)
    integer               :: edx, vdx, ndx, inside_count, first_inside, new_count, nr_coeffs, tdx, prev_idx, next_idx
    logical               :: is_parabolic, edge_is_bisected, edge_could_be_trisected, new_vertex(2)
    real*8                :: x_eta, x_tau

    if (poly%intersected .and. poly%parabolic) then
      print*, 'ERROR: a polygon cannot be intersected if it was previously intersected by a parabola'
      return
    endif

    is_parabolic = parabola%kappa0 /= 0
    if (is_parabolic .and. parabola%kappa0 > 0) then 
      poly%original_moments = cmpMoments(poly)
      poly%complement = .true.
    endif

    ! Check for each vertex if it under or above the parabola
    inside_count = 0
    first_inside = 0
    do vdx=1,poly%nverts
      x_eta = dot_product(poly%verts(:,vdx), parabola%normal) - parabola%shift
      dist(vdx) = x_eta
      if (is_parabolic) then
        x_tau = dot_rotate(poly%verts(:,vdx), parabola%normal)
        dist(vdx) = dist(vdx) + (parabola%kappa0/2) * x_tau**2
        if (poly%complement) dist(vdx) = -dist(vdx)
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
    nr_coeffs = 0
    vdx = merge(first_inside, 1, first_inside>0)
    do edx=1,poly%nverts
      ndx = vdx + 1
      if (ndx > poly%nverts) ndx = 1

      new_vertex = .false.
      edge_is_bisected = dist(vdx)<=0 .neqv. dist(ndx)<=0
      edge_could_be_trisected = .false.
      if (is_parabolic) then 
        if (.not. edge_is_bisected) then
          edge_could_be_trisected = dist(vdx) <= 0
        endif
        if (edge_is_bisected .or. edge_could_be_trisected) then
          nr_coeffs = parabola_line_intersection(coeffs, parabola, buffer(:,vdx), buffer(:,ndx))
        endif
      endif

      if (edge_is_bisected) then
        ! Insert new vertex (edge bisection)
        if (.not. is_parabolic) then
          coeffs(1) = abs(dist(vdx) / (dist(ndx) - dist(vdx)))
        else
          if (abs(coeffs(2) - .5) < abs(coeffs(1) - .5)) then 
            coeffs(1) = coeffs(2)
          endif
        endif
        new_vertex(1) = .true.
      elseif (edge_could_be_trisected) then
        if (nr_coeffs==2) then
          new_vertex = .true.
        else
          if (coeffs(2) < coeffs(1)) coeffs = coeffs([2, 1])
          if (abs(coeffs(1)) < 1D-14 .and. abs(coeffs(2)-1) < 1D-14) then
            prev_idx = merge(vdx - 1, poly%nverts, vdx>1)
            next_idx = merge(ndx + 1, 1, ndx<poly%nverts)

            new_vertex(1) = dist(prev_idx)<=0 .eqv. dist(vdx)<=0
            new_vertex(2) = dist(next_idx)<=0 .eqv. dist(ndx)<=0
          endif
        endif


      endif

      do tdx=1,2
        if (new_vertex(tdx)) then
          new_count = new_count + 1
          poly%on_parabola(new_count) = .true.
          poly%verts(:,new_count) = (1 - coeffs(tdx)) * buffer(:,vdx) + coeffs(tdx) * buffer(:,ndx)
        endif
      enddo

      if (dist(ndx)<=0) then
        ! Keep old
        new_count = new_count + 1
        poly%on_parabola(new_count) = .false.
        poly%verts(:,new_count) = buffer(:,ndx)
      endif

      vdx = ndx
    enddo

    poly%nverts = new_count

    poly%parabolic = is_parabolic
    poly%intersected = .true.
    if (.not. poly%complement) then
      poly%parabola = parabola
    else
      poly%parabola = complement(parabola)
    endif
   
  end subroutine

  subroutine compute_momonial(poly, nr)
    implicit none
    
    type(tPolygon), intent(inout) :: poly
    integer, intent(in)   :: nr

    ! Local variables
    integer               :: mdx, vdx, ndx, old_nr, edx
    
    old_nr = poly%avail_monomial
    if (nr <= old_nr) return
    if (nr > MAX_MONOMIAL) then
      print*, 'ERROR in compute_momonial: too many monomials requested'
    endif

    if (old_nr<0) old_nr = -1

    ! We compute integral of x_tau^mdx over edges which are parabola
    do mdx=old_nr+1,nr
      edx = 0 ! Parabolic edge index

      poly%monomials_sum(mdx) = 0
      do vdx=1,poly%nverts
        ndx = vdx + 1
        if (vdx==poly%nverts) ndx = 1
        
        if (poly%on_parabola(vdx) .and. poly%on_parabola(ndx)) then
          edx = edx + 1

          if (mdx==0) then
            poly%x_tau(1,edx) = dot_rotate(poly%verts(:,vdx), poly%parabola%normal)
            poly%x_tau(2,edx) = dot_rotate(poly%verts(:,ndx), poly%parabola%normal)

            poly%x_eta(1,edx) = dot_product(poly%verts(:,vdx), poly%parabola%normal) - poly%parabola%shift
            poly%x_eta(2,edx) = dot_product(poly%verts(:,ndx), poly%parabola%normal) - poly%parabola%shift

            poly%x_tau_power(1,edx) = poly%x_tau(1,edx)
            poly%x_tau_power(2,edx) = poly%x_tau(2,edx)
          else
            poly%x_tau_power(1,edx) = poly%x_tau_power(1,edx) * poly%x_tau(1,edx)
            poly%x_tau_power(2,edx) = poly%x_tau_power(2,edx) * poly%x_tau(2,edx)
          endif
          
          poly%monomials(mdx,edx) = (poly%x_tau_power(2,edx) - poly%x_tau_power(1,edx)) / (mdx+1)
          poly%monomials_sum(mdx) = poly%monomials_sum(mdx) + poly%monomials(mdx,edx)
        endif
      enddo
    enddo

    poly%avail_monomial = nr
  end subroutine

  real*8 function parabola_volume_correction(poly) result(vol)
    implicit none
    
    type(tPolygon), intent(inout) :: poly

    ! Local variables
    real*8                :: coeff(2), dtau
    integer               :: edx, vdx, ndx

    call compute_momonial(poly, 2)

    vol = 0
    edx = 0 ! Parabolic edge index
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (vdx==poly%nverts) ndx = 1
      
      if (poly%on_parabola(vdx) .and. poly%on_parabola(ndx)) then
        edx = edx + 1

        dtau = poly%x_tau(2,edx) - poly%x_tau(1,edx)
        if (dtau==0) cycle

        ! in the local coordinates the polygon face is given by
        ! x_η = c_1 * x_τ + c_2
        coeff(1) = (poly%x_eta(2,edx) - poly%x_eta(1,edx)) / dtau
        coeff(2) = (poly%x_eta(2,edx) + poly%x_eta(1,edx)) / 2 - coeff(1) * (poly%x_tau(1,edx) + poly%x_tau(2,edx)) / 2
    
        vol = vol - (poly%parabola%kappa0/2) * poly%monomials(2,edx) - &
          coeff(1) * poly%monomials(1,edx) - coeff(2) * poly%monomials(0,edx)
      endif
    enddo
  end function

  function parabola_moments_correction(poly) result(moments)
    implicit none
    
    type(tPolygon), target, intent(inout) :: poly
    real*8                :: moments(3)

    ! Local variables
    real*8                :: coeff(2), vol_corr, mom_corr(2), dtau
    integer               :: edx, vdx, ndx
    type(tParabola), pointer :: parabola

    call compute_momonial(poly, 4)

    parabola => poly%parabola

    moments = 0
    edx = 0 ! Parabolic edge index
    do vdx=1,poly%nverts
      ndx = vdx + 1
      if (vdx==poly%nverts) ndx = 1
      
      if (poly%on_parabola(vdx) .and. poly%on_parabola(ndx)) then
        edx = edx + 1

        ! in the local coordinates the polygon face is given by
        ! x_η = c_1 * x_τ + c_2
        dtau = poly%x_tau(2,edx) - poly%x_tau(1,edx)
        if (dtau == 0) cycle
        
        coeff(1) = (poly%x_eta(2,edx) - poly%x_eta(1,edx)) / dtau
        coeff(2) = (poly%x_eta(2,edx) + poly%x_eta(1,edx)) / 2 - coeff(1) * (poly%x_tau(1,edx) + poly%x_tau(2,edx)) / 2
    
        vol_corr = -(parabola%kappa0/2) * poly%monomials(2,edx) - &
          coeff(1) * poly%monomials(1,edx) - coeff(2) * poly%monomials(0,edx)

        ! Corrections to first moment in parabola coordinates
        mom_corr(1) = -(parabola%kappa0/2) * poly%monomials(3,edx) - (coeff(1) * poly%monomials(2,edx) + &
          coeff(2) * poly%monomials(1,edx))
        mom_corr(2) = parabola%kappa0**2 * poly%monomials(4,edx)/8 - (coeff(1)**2 * poly%monomials(2,edx) + &
          2 * coeff(1) * coeff(2) * poly%monomials(1,edx) + coeff(2)**2 * poly%monomials(0,edx))/2 

        moments(1) = moments(1) + vol_corr
        moments(2) = moments(2) + parabola%normal(1) * (mom_corr(2) + parabola%shift * vol_corr) &
          - parabola%normal(2) * mom_corr(1);
        moments(3) = moments(3) + parabola%normal(2) * (mom_corr(2) + parabola%shift * vol_corr) &
          + parabola%normal(1) * mom_corr(1);
      endif
    enddo
  end function

  ! Given a parabola and a line connecting the points pos1, pos2; find the intersection
  integer function parabola_line_intersection(roots, parabola, pos1, pos2) result(nr_roots)
    use m_poly_roots,     only: real_roots_2
    implicit none
    
    type(tParabola), intent(in) :: parabola
    real*8, intent(in)    :: pos1(2), pos2(2)
    real*8, intent(out)   :: roots(2)

    ! Local variables
    real*8                :: coeff(3)            
    logical               :: root_is_good(2)
    integer               :: rdx

    ! We parametrize the line as: l(t) = pos1 + t * (pos2 - pos1),
    ! and solve for t
    coeff(1) = (parabola%kappa0/2) * dot_relative_rotate(pos2, parabola%normal, pos1)**2
    coeff(2) = parabola%kappa0 * dot_rotate(pos1, parabola%normal) * dot_relative_rotate(pos2, parabola%normal, pos1) &
      + dot_relative(pos2, parabola%normal, pos1)
    coeff(3) = dot_product(pos1, parabola%normal) - parabola%shift + (parabola%kappa0/2) * dot_rotate(pos1, parabola%normal)**2

    call real_roots_2(roots, coeff)

    nr_roots = 0
    do rdx=1,2
      root_is_good(rdx) = .not. isnan(roots(rdx)) .and. roots(rdx) >= 0 .and. roots(rdx) <= 1
      if (root_is_good(rdx)) nr_roots = nr_roots + 1
    enddo
  
    if (root_is_good(2) .and. .not. root_is_good(1)) then 
      roots = roots([2, 1])
    endif
  end function

  pure real*8 function dot_relative_rotate(va, vb, vr) result(drr)
    implicit none
    
    real*8, intent(in)    :: va(2), vb(2), vr(2)

    drr = (va(2) - vr(2))*vb(1) - (va(1) - vr(1))*vb(2)
  end function

  pure real*8 function dot_rotate(va, vb) result(dr)
    implicit none
    
    real*8, intent(in)    :: va(2), vb(2)

    dr = va(2)*vb(1) - va(1)*vb(2)
  end function

  pure real*8 function dot_relative(va, vb, vr) result(dr)
    implicit none
    
    real*8, intent(in)    :: va(2), vb(2), vr(2)

    dr = (va(1) - vr(1))*vb(1) + (va(2) - vr(2))*vb(2)
  end function

  ! subroutine makeBox_bounds(poly, rbounds)
  !   implicit none

  !   real*8, intent(in)    :: rbounds(2, 2)
  !   type(tPolygon), intent(out) :: poly

  !   poly%verts(:,1) = [rbounds(1,1), rbounds(2,1)]
  !   poly%verts(:,2) = [rbounds(1,2), rbounds(2,1)]
  !   poly%verts(:,3) = [rbounds(1,2), rbounds(2,2)]
  !   poly%verts(:,4) = [rbounds(1,1), rbounds(2,2)]

  !   poly%nverts = 4

  ! end subroutine

  subroutine makeBox_dx(poly, dx)
    implicit none

    real*8, intent(in)    :: dx(2)
    type(tPolygon), intent(out) :: poly

    poly%verts(:,1) = [-dx(1)/2, -dx(2)/2]
    poly%verts(:,2) = [dx(1)/2, -dx(2)/2]
    poly%verts(:,3) = [dx(1)/2, dx(2)/2]
    poly%verts(:,4) = [-dx(1)/2, dx(2)/2]

    poly%nverts = 4

  end subroutine

  subroutine makeBox_xdx(poly, x, dx)
    implicit none
    
    real*8, intent(in)    :: x(2), dx(2)
    type(tPolygon), intent(out) :: poly

    poly%verts(:,1) = x + [-dx(1)/2, -dx(2)/2]
    poly%verts(:,2) = x + [dx(1)/2, -dx(2)/2]
    poly%verts(:,3) = x + [dx(1)/2, dx(2)/2]
    poly%verts(:,4) = x + [-dx(1)/2, dx(2)/2]

    poly%nverts = 4

  end subroutine  

  subroutine init(poly, pos)
    implicit none

    type(tPolygon), intent(out) :: poly
    real*8, intent(in)    :: pos(1:, 1:)

    ! Local variables
    integer               :: nverts

    nverts = size(pos, 2)
    if (nverts > MAX_NR_VERTS .or. size(pos, 1) /= 2) return

    poly%verts(1:nverts,:) = pos(1:nverts,:)
    poly%nverts = nverts
  end subroutine

  subroutine split(polys, parabola, out_pos, out_neg)
    implicit none

    type(tPolygon), intent(in) :: polys(1:)
    type(tParabola), intent(in) :: parabola
    type(tPolygon), intent(out) :: out_pos(1:), out_neg(1:)

    ! Local variables
    integer               :: vdx
    type(tParabola)       :: c_parabola

    c_parabola = complement(parabola)

    do vdx=1,size(polys, 1)
      call copy(out=out_pos(vdx), in=polys(vdx))
      call intersect(out_pos(vdx), parabola)
    
      call copy(out=out_neg(vdx), in=polys(vdx))
      call intersect(out_neg(vdx), c_parabola)
    enddo
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

      out%parabolic = in%parabolic
      if (out%parabolic) then
        out%on_parabola(1:out%nverts) = in%on_parabola(1:out%nverts)

        out%complement = in%complement
        if (in%complement) out%original_moments = in%original_moments
  
        out%x_tau = in%x_tau
        out%x_eta = in%x_eta
        
        out%avail_monomial = in%avail_monomial
        if (out%avail_monomial>=0) then 
          out%x_tau_power = in%x_tau_power
          out%monomials = in%monomials
          out%monomials_sum = in%monomials_sum
        endif
      endif

    endif

  end subroutine
end