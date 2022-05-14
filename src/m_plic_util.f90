module m_plic_util
  real*8, parameter     :: NORMAL_TOL = 1E-12
contains

   function cmpMoments2d(normal, dx, shift) result(moments)
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
      moments = cmpMoments2d_pc_ordered(normal_, dx_, shift + max_shift_plane, max_shift_plane)
    endif

      if (normal(perm(1)) < 0.0) moments(2) = -moments(2)
      if (normal(perm(2)) < 0.0) moments(3) = -moments(3)
      if (perm(1) /= 1) then 
        moments(2:3) = [moments(3), moments(2)]
      endif
  end function


  pure function cmpMoments2d_pc_ordered(normal, dx, pc, max_shift_plane) result(moments)
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

  pure real*8 function cmpShift2d(normal, dx, volume) result(shift)

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

    if (normal_(1) < NORMAL_TOL) then
      shift = volume / product(dx_) - max_shift_plane
    else
      max_shift_plane = dot_product(normal_, dx_)/2
      shift = cmpShift2d_pc_ordered(normal_, dx_, volume, max_shift_plane) - max_shift_plane
    endif


  end function

  pure real*8 function cmpShift2d_pc_ordered(normal, dx, volume, max_shift_plane) result(pc)

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
end module