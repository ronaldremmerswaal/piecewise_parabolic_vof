module m_recon_util
  use m_common

  private
  public :: cmpVolume, cmpMoments, cmpInterfaceMoments, cmpShift, cmpSymmDiff, makeParabola, cmpSymmDiffVolume

  real*8, parameter     :: NORMAL_TOL = 1E-12
  real*8, parameter     :: DEFAULT_SHIFT_TOL = 1E-15
  integer, parameter    :: MAX_SHIFT_ITERS = 30

  interface cmpMoments
    module procedure cmpMoments2d_plane, cmpMoments_dx
  end interface

  interface cmpVolume
    module procedure cmpVolume2d_plane, cmpVolume_dx
  end interface

  interface cmpShift
    module procedure cmpShift2d_plane, cmpShift2d_parabolic, cmpShift2d_poly
  end interface

  interface makeParabola
    module procedure makeParabola_noshift
  end interface

  interface cmpSymmDiffVolume
    module procedure cmpSymmDiffVolume2d_dxIn, cmpSymmDiffVolume2d_polyIn
  end interface

  interface cmpInterfaceMoments
    module procedure cmpInterfaceMoments_plane
  end interface

contains

   function cmpVolume2d_plane(normal, dx, shift) result(volume)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), shift
    real*8              :: volume

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
      volume = ((shift + max_shift_plane) / normal_(2)) * dx_(1)
    else
      volume = cmpVolume2d_pc_ordered_plane(normal_, dx_, shift + max_shift_plane, max_shift_plane)
    endif

  end function

   subroutine cmpMoments2d_plane(moments, normal, dx, shift)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), shift
    real*8, intent(out) :: moments(3)

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
      call cmpMoments2d_pc_ordered_plane(moments, normal_, dx_, shift + max_shift_plane, max_shift_plane)
    endif

    if (normal(perm(1)) < 0.0) moments(2) = -moments(2)
    if (normal(perm(2)) < 0.0) moments(3) = -moments(3)
    if (perm(1) /= 1) then 
      moments(2:3) = [moments(3), moments(2)]
    endif
  end subroutine

  pure function cmpVolume2d_pc_ordered_plane(normal, dx, pc, max_shift_plane) result(volume)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), pc, max_shift_plane
    real*8              :: volume

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
      volume = 0.0
    elseif (pcMod<normal(1)*dx(1)) then
      volume = pcMod**2/(2.0*normal(1)*normal(2))
    else ! if (pcMod<=max_shift_plane) then
      volume = pcMod*dx(1)/normal(2) - normal(1)*dx(1)**2/(2.0*normal(2))
    endif

    if (largerThanHalf) volume = product(dx) - volume
  end function

  subroutine cmpMoments2d_pc_ordered_plane(moments, normal, dx, pc, max_shift_plane)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), pc, max_shift_plane
    real*8, intent(out) :: moments(3)

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
  end subroutine

  real*8 function cmpVolume_dx(dx, parabola) result(vol)
    use m_polygon
    implicit none
    
    real*8, intent(in)    :: dx(2)
    type(tParabola), intent(in) :: parabola

    ! Local variables
    type(tPolygon)        :: poly

    call makeBox(poly, dx)
    call intersect(poly, parabola)
    vol = cmpVolume(poly)

  end function

  subroutine cmpMoments_dx(mom, dx, parabola)
    use m_polygon
    implicit none
    
    real*8, intent(in)    :: dx(2)
    type(tParabola), intent(in) :: parabola
    real*8, intent(out)   :: mom(3)

    ! Local variables
    type(tPolygon)        :: poly

    call makeBox(poly, dx)
    call intersect(poly, parabola)
    call cmpMoments(mom, poly)

  end subroutine


  function cmpInterfaceMoments_plane(normal, dx, shift) result(moments)
    implicit none

    real*8, intent(in)  :: normal(2), dx(2), shift
    real*8              :: moments(3)

    ! Local variables:
    real*8              :: normal_(2), dx_(2), max_shift_plane, pMoment1(2)
    integer             :: perm(2)


    moments = 0
    normal_ = abs(normal)
    max_shift_plane = dot_product(normal_, dx)/2
    if (abs(shift) >= max_shift_plane) return

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
      moments(1) = dx_(1)
      pMoment1(1) = 0
      pMoment1(2) = shift * moments(1)
      if (normal(perm(2)) < 0.0) pMoment1(2) = -pMoment1(2)
    else
      moments(1) = cmpInterfaceMoment0_pc_ordered_plane(normal_, dx_, shift + max_shift_plane, max_shift_plane)
      if (moments(1)==0) then
        pMoment1 = 0
      else
        pMoment1 = cmpInterfaceMoment1_pc_ordered_plane(normal_, dx_, shift + max_shift_plane, max_shift_plane)
        if (normal(perm(1)) < 0.0) pMoment1(1) = -pMoment1(1)
        if (normal(perm(2)) < 0.0) pMoment1(2) = -pMoment1(2)
      endif
    endif

    moments(1+perm(1)) = pMoment1(1); moments(1+perm(2)) = pMoment1(2)
  end function

  pure function cmpInterfaceMoment0_pc_ordered_plane(normal, dx, pc, max_shift_plane) result(area)
    implicit none

    ! Input/output arguments:
    real*8, intent(in)  :: normal(2), dx(2), pc, max_shift_plane
    real*8              :: area

    ! Local variables:
    real*8              :: pcMod
    logical             :: largerThanHalf

    largerThanHalf = pc > max_shift_plane
    if (largerThanHalf) then
      pcMod = 2*max_shift_plane - pc
    else
      pcMod = pc
    endif

    if (pcMod <= 0.0) then
      area = 0
    elseif (pcMod<normal(1)*dx(1)) then
      area = sqrt(pcMod**2/normal(1)**2 + pcMod**2/normal(2)**2)
    else ! if (pcMod<=max_shift_plane) then
      area = sqrt((dx(1)**2*(normal(1)**2 + normal(2)**2))/normal(2)**2)
    endif
  end function

  pure function cmpInterfaceMoment1_pc_ordered_plane(normal, dx, pc, max_shift_plane) result(moment1)
    implicit none

    ! Input/output arguments:
    real*8,intent(in)   :: normal(2), dx(2), pc, max_shift_plane
    real*8              :: moment1(2)

    ! Local variables:
    real*8              :: pcMod
    logical             :: largerThanHalf

    largerThanHalf = pc > max_shift_plane
    if (largerThanHalf) then
      pcMod = 2*max_shift_plane - pc
    else
      pcMod = pc
    endif

    ! We compute moment1 w.r.t. cell corner
    if (pcMod <= 0.0) then
      moment1 = 0
    elseif (pcMod<normal(1)*dx(1)) then
      moment1(1) = -(sqrt(pcMod**2/normal(1)**2 + pcMod**2/normal(2)**2)*(dx(1) - pcMod/normal(1)))/2
      moment1(2) = -(sqrt(pcMod**2/normal(1)**2 + pcMod**2/normal(2)**2)*(dx(2) - pcMod/normal(2)))/2
    else ! if (pcMod<=max_shift_plane) then
      moment1(1) = 0
      moment1(2) = -((dx(1)*normal(1) - 2*pcMod + dx(2)*normal(2))*sqrt((dx(1)**2*(normal(1)**2 +& 
        normal(2)**2))/normal(2)**2))/(2*normal(2))
    endif

    if (largerThanHalf) moment1 = -moment1
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

  real*8 function cmpSymmDiffVolume2d_dxIn(x, dx, parabola, levelSet) result(sd)
    use m_polygon
    implicit none

    real*8, intent(in)    :: x(2), dx(2)
    type(tParabola), intent(in) :: parabola
    procedure(levelset_fun) :: levelSet
    
    type(tPolygon)        :: cell

    call makeBox(cell, x, dx)
    sd = cmpSymmDiffVolume(cell, parabola, levelset, x0=x)
  end

  real*8 function cmpSymmDiffVolume2d_polyIn(cell, parabola, levelSet, x0) result(sd)
    use m_polygon
    implicit none

    type(tPolygon), intent(in) :: cell
    type(tParabola), intent(in) :: parabola
    procedure(levelset_fun) :: levelSet
    real*8, optional, intent(in) :: x0(2)

    ! Local variables
    type(tPolygon)        :: exact_gas, exact_liq


    ! Construct polygonal approximation of exact gas & liquid domains
    ! (relative to the cell centroid)
    call polyApprox(exact_gas, cell, levelSet, GAS_PHASE)
    call polyApprox(exact_liq, cell, levelSet, LIQUID_PHASE)
    if (present(x0)) then
      call shift_by(exact_gas, -x0)
      call shift_by(exact_liq, -x0)  
    endif

    ! Compute symmetric difference
    call intersect(exact_gas, parabola)
    call intersect(exact_liq, complement(parabola))
    
    sd = cmpVolume(exact_gas) + cmpVolume(exact_liq)
  end

  real*8 function cmpShift2d_parabolic(normal, dx, liqVol, kappa0, relTol, volume, intersected, shift0) result(shift)
    use m_polygon
    use m_optimization,   only: brent

    implicit none

    real*8, intent(in)    :: normal(2), dx(2), liqVol, kappa0
    real*8, optional, intent(in) :: relTol, shift0
    real*8, optional, intent(out) :: volume
    type(tPolygon), optional, intent(out) :: intersected

    ! Local variables
    type(tPolygon)        :: cell
    real*8                :: volume_, cellVol, shift0_, err0
    real*8                :: max_shift_plane_eta, max_shift_plane_tau
    real*8                :: shift_l, shift_r, err_l, err_r
    real*8                :: relTol_
    type(tParabola)       :: parabola
    ! integer, save         :: nr_shifts = 0, nr_evals = 0

    max_shift_plane_eta = dot_product(abs(normal), dx)/2
    max_shift_plane_tau = dot_product(abs([normal(2), normal(1)]), dx)/2

    cellVol = product(dx)
    if (liqVol <= 0) then
      if (kappa0 > 0) then
        shift = -max_shift_plane_eta
      else
        shift = kappa0 * max_shift_plane_tau**2 / 2. - max_shift_plane_eta
      endif
      if (present(volume)) volume = 0.0D0
      if (present(intersected)) intersected%nverts = 0
      return
    elseif (liqVol >= cellVol) then
      if (kappa0 > 0) then
        shift = kappa0 * max_shift_plane_tau**2 / 2. + max_shift_plane_eta
      else
        shift = max_shift_plane_eta
      endif
      if (present(volume)) volume = cellVol
      if (present(intersected)) call makeBox(intersected, dx)
      return
    endif

    shift0_ = merge(shift0, d_qnan, present(shift0))

    if (isnan(shift0_)) then
      ! Use PLIC to get a good initial bracket (one side at least)
      shift0_ = cmpShift2d_plane(normal, dx, liqVol)
    endif
    call makeParabola(parabola, normal, kappa0, shift0_)
    err0 = volume_error_function(shift0_)

    relTol_ = merge(relTol, DEFAULT_SHIFT_TOL, present(relTol))

    ! Try to get a better bracket with an (educated) guess
    if (abs(err0) < cellVol * relTol_) then
      ! iff kappa0 == 0.0
      shift = shift0_
      if (present(volume)) volume = volume_
      if (present(intersected)) call copy(out=intersected, in=cell)
      return
    endif

    if (err0 > 0.0) then
      ! iff kappa0 < 0.0
      shift_r = shift0_
      err_r = err0

      shift_l = -abs(kappa0) * max_shift_plane_tau**2 / 2 - max_shift_plane_eta
      err_l = -liqVol
    else
      ! iff kappa0 > 0.0
      shift_l = shift0_
      err_l = err0

      shift_r = abs(kappa0) * max_shift_plane_tau**2 / 2 + max_shift_plane_eta
      err_r = cellVol - liqVol
    endif

    shift = brent(volume_error_function, shift_l, shift_r, max_shift_plane_eta * relTol_, MAX_SHIFT_ITERS, err_l, err_r)
    if (present(volume)) volume = volume_
    if (present(intersected)) call copy(out=intersected, in=cell)
  contains

    real*8 function volume_error_function(shift_tmp) result(err)
      implicit none

      real*8, intent(in)  :: shift_tmp  

      call makeBox(cell, dx)

      parabola%shift = shift_tmp
      call intersect(cell, parabola)
      volume_ = cmpVolume(cell)

      err = volume_ - liqVol
    end function
  end function


  subroutine makeParabola_noshift(parabola, normal, kappa0, dx, liqVol)
    use m_polygon
    implicit none
    
    real*8, intent(in)    :: normal(2), kappa0, dx(2), liqVol
    type(tParabola), intent(out) :: parabola

    parabola%normal = normal
    parabola%kappa0 = kappa0
    parabola%shift = cmpShift(normal, dx, liqVol, kappa0)
  end subroutine

  real*8 function cmpShift2d_poly(normal, cell, liqVol, kappa0, x0, relTol, volume, intersected) result(shift)
    use m_polygon
    use m_optimization,   only: brent

    implicit none

    real*8, intent(in)    :: normal(2), liqVol, kappa0
    type(tPolygon), intent(inout) :: cell
    real*8, optional, intent(in) :: x0(2), relTol
    real*8, optional, intent(out) :: volume
    type(tPolygon), optional, intent(out) :: intersected

    ! Local variables
    real*8                :: volume_, cellVolume, shift0, err0
    real*8                :: shift_l, shift_r, err_l, err_r, eta_dist, tau_dist_sq
    real*8                :: min_eta_dist, max_eta_dist, max_tau_dist_sq
    real*8                :: relTol_, tau(2), pos(2)
    integer               :: vdx
    type(tPolygon)        :: cell_copy
    type(tParabola)       :: parabola

    relTol_ = merge(relTol, DEFAULT_SHIFT_TOL, present(relTol))

    min_eta_dist = 0
    max_eta_dist = 0
    max_tau_dist_sq = 0
    
    tau = [-normal(2), normal(1)]
    do vdx=1,cell%nverts
      pos = cell%verts(:,vdx)
      if (present(x0)) pos = pos - x0
      eta_dist = dot_product(normal, pos)
      tau_dist_sq = dot_product(tau, pos)**2
      if (vdx == 1) then
        min_eta_dist = eta_dist
        max_eta_dist = eta_dist
        max_tau_dist_sq = tau_dist_sq
      else
        if (eta_dist < min_eta_dist) min_eta_dist = eta_dist
        if (eta_dist > max_eta_dist) max_eta_dist = eta_dist
        if (tau_dist_sq > max_tau_dist_sq) max_tau_dist_sq = tau_dist_sq
      endif
    enddo

    cellVolume = cmpVolume(cell)
    if (liqVol <= 0.0) then
      if (kappa0 > 0) then
        shift = min_eta_dist
      else
        shift = kappa0 * max_tau_dist_sq / 2 + min_eta_dist
      endif
      if (present(volume)) volume = 0.0D0
      if (present(intersected)) intersected%nverts = 0
      return
    elseif (liqVol >= cellVolume) then
      if (kappa0 > 0) then
        shift = kappa0 * max_tau_dist_sq / 2 + max_eta_dist
      else
        shift = max_eta_dist
      endif
      if (present(volume)) volume = cellVolume
      if (present(intersected)) call copy(out=intersected, in=cell)
      return
    endif
    
    shift0 = (min_eta_dist + max_eta_dist) / 2
    call makeParabola(parabola, normal, kappa0, shift0)
    err0 = volume_error_function(shift0)
    
    if (err0 == 0) then
      shift = shift0
      if (present(volume)) volume = volume_
      if (present(intersected)) call copy(out=intersected, in=cell_copy)
      return
    elseif (err0 > 0) then
      shift_r = shift0
      err_r = err0

      if (kappa0 > 0) then
        shift_l = min_eta_dist
      else
        shift_l = kappa0 * max_tau_dist_sq / 2 + min_eta_dist
      endif

      err_l = -liqVol
    else
      shift_l = shift0
      err_l = err0

      if (kappa0 > 0) then
        shift_r = kappa0 * max_tau_dist_sq / 2 + max_eta_dist
      else
        shift_r = max_eta_dist
      endif
      err_r = cellVolume - liqVol
    endif

    shift = brent(volume_error_function, shift_l, shift_r, &
      (max_eta_dist - min_eta_dist) * relTol_, MAX_SHIFT_ITERS, err_l, err_r)
    if (present(volume)) volume = volume_
    if (present(intersected)) call copy(out=intersected, in=cell_copy)

    contains

      real*8 function volume_error_function(shift_tmp) result(err)
        implicit none

        real*8, intent(in)    :: shift_tmp

        call copy(out=cell_copy, in=cell)
        parabola%shift = shift_tmp
        call intersect(cell_copy, parabola)
        volume_ = cmpVolume(cell_copy)!, x0=x0)
        err = volume_ - liqVol
      end function
  end function
end module