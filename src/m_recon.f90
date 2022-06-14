module m_recon

  private
  public :: mofNormal, pmofNormal, lviraNormal, plviraNormal, prostNormal

  interface pmofNormal
    module procedure pmofNormal_poly, pmofNormal_rect
  end interface

contains

  function lviraNormal(refVolumes, dxs, verbose, errTol) result(normal)
    use m_optimization
    use m_recon_util

    implicit none

    real*8, intent(in)    :: refVolumes(-1:1,-1:1), dxs(-1:1,2)
    logical, intent(in), optional :: verbose
    real*8, intent(in), optional :: errTol
    real*8                :: normal(2)

        ! Local variables:
    real*8                :: errTol_, lviraAngle
    logical               :: verbose_

    verbose_ = merge(verbose, .false., present(verbose))
    errTol_ = merge(errTol, 1D-8, present(errTol))

    lviraAngle = lvira_angle_guess(refVolumes, dxs)
    lviraAngle = brent_min(dcost, lviraAngle, errTol_, 25, verbose_, maxStep=0.5D0)
    normal = [dcos(lviraAngle), dsin(lviraAngle)]
  contains

    real*8 function dcost(angle, err) result(derr)

      implicit none

      real*8, intent(in)  :: angle
      real*8, intent(out), optional :: err

      real*8              :: err_

      err_ = lvira_error(refVolumes, angle, dxs, derivative=derr)
      if (present(err)) err = err_
    end function
  end function

  real*8 function lvira_error(refVolumes, angle, dxs, derivative) result(err)
    use m_common
    use m_recon_util
    use m_r2d_parabolic

    implicit none

    real*8, intent(in)    :: refVolumes(-1:1,-1:1), dxs(-1:1,2), angle
    real*8, intent(out), optional :: derivative

    ! Local variables
    real*8                :: normal(2), shift, volume, derivative_local, err_local, absNormal(2)
    real*8                :: xc_neighbour(2), dx_neighbour(2), cellVol_neighbour
    real*8                :: shift_neighbour, shift_derivative, tangent(2), iFaceMoments(3), shift_max
    integer               :: i, j

    normal = [dcos(angle), dsin(angle)]
    absNormal = abs(normal)
    shift = cmpShift(normal, dxs(0,:), refVolumes(0,0))
    
    err = 0
    if (present(derivative)) then      
      derivative = 0
      ! In general the derivative of the volume is given by
      !   d M_0(c^l) /d \theta = M_0(I_c) d shift / d \theta - tangent \cdot M_1(I_c),
      ! where M_0(c^l) is the liquid volume inside the corresponding cell, \theta is the normal angle
      ! M_0(I_c) is the interface area inside c and M_1(I_c) is the first moment of the interface inside c

      tangent = [-normal(2), normal(1)]
      iFaceMoments = cmpInterfaceMoments(normal, dxs(0,:), shift)

      ! In the central cell (i=j=0) we impose conservation of volume: d M_0(c^l) /d \theta = 0
      ! which determines d shift / d \theta
      shift_derivative = dot_product(tangent, iFaceMoments(2:3)) / iFaceMoments(1)
    endif
    do j=-1,1
    do i=-1,1
      if (i==0 .and. j==0) cycle

      dx_neighbour(1) = dxs(i, 1)
      dx_neighbour(2) = dxs(j, 2)
      cellVol_neighbour = product(dx_neighbour)

      xc_neighbour(1) = i * (dxs(0,1) + dx_neighbour(1))/2
      xc_neighbour(2) = j * (dxs(0,2) + dx_neighbour(2))/2

      shift_neighbour = shift - dot_product(xc_neighbour, normal)

      shift_max = dot_product(dx_neighbour, absNormal)/2
      derivative_local = 0
      if (shift_neighbour < -shift_max) then
        volume = 0
      elseif (shift_neighbour > shift_max) then
        volume = product(dx_neighbour)
      else
        if (present(derivative)) then
          iFaceMoments = cmpInterfaceMoments(normal, dx_neighbour, shift_neighbour)
          iFaceMoments(2:3) = iFaceMoments(2:3) + iFaceMoments(1) * xc_neighbour
          derivative_local = iFaceMoments(1) * shift_derivative - dot_product(tangent, iFaceMoments(2:3))
        endif
        volume = cmpVolume(normal, dx_neighbour, shift_neighbour)
      endif
      err_local = (volume - refVolumes(i,j)) / cellVol_neighbour
      
      err = err + err_local**2
      if (present(derivative)) then
        derivative = derivative + 2 * err_local * derivative_local / cellVol_neighbour
      endif
    enddo
    enddo

  end function

  function plviraNormal(refVolumes, kappa0, dxs, verbose, errTol) result(normal)
    use m_optimization
    use m_recon_util

    implicit none

    real*8, intent(in)    :: refVolumes(-1:1,-1:1), kappa0, dxs(-1:1,2)
    logical, intent(in), optional :: verbose
    real*8, intent(in), optional :: errTol
    real*8                :: normal(2)

        ! Local variables:
    real*8                :: errTol_, plviraAngle
    logical               :: verbose_

    verbose_ = merge(verbose, .false., present(verbose))
    errTol_ = merge(errTol, 1D-8, present(errTol))

    plviraAngle = lvira_angle_guess(refVolumes, dxs)
    plviraAngle = brent_min(dcost, plviraAngle, errTol_, 25, verbose_, maxStep=0.5D0)
    normal = [dcos(plviraAngle), dsin(plviraAngle)]
  contains

    real*8 function dcost(angle, err) result(derr)

      implicit none

      real*8, intent(in)  :: angle
      real*8, intent(out), optional :: err

      real*8              :: err_, derivatives(2)

      err_ = plvira_error(refVolumes, angle, kappa0, dxs, derivatives)
      derr = derivatives(1)
      if (present(err)) err = err_
    end function
  end function

  function prostNormal(refVolumes, kappa0, dxs, verbose, errTol) result(normal)
    use m_optimization
    use m_recon_util

    implicit none

    real*8, intent(in)    :: refVolumes(-1:1,-1:1), dxs(-1:1,2)
    logical, intent(in), optional :: verbose
    real*8, intent(in), optional :: errTol
    real*8, intent(inout) :: kappa0
    real*8                :: normal(2)

        ! Local variables:
    real*8                :: prostAngle, XSOL(2), lengthScale
    type(optimOpts)       :: opts = optimOpts(maxFEval=100)
    type(lsOpts)          :: ls = lsOpts(decreaseCondition=1E-3)
    type(optimInfo)       :: info

    opts%verbose = merge(verbose, .false., present(verbose))
    opts%errTol = merge(errTol, 1D-8, present(errTol))

    prostAngle = lvira_angle_guess(refVolumes, dxs)
    lengthScale = norm2(dxs(0,:))
    XSOL = [prostAngle, lengthScale * kappa0]

    call optimize(XSOL, cost, opts, info, ls_opts=ls, fun_and_grad=cost_and_grad)
    prostAngle = XSOL(1)
    kappa0 = XSOL(2) / lengthScale

    normal = [dcos(prostAngle), dsin(prostAngle)]

  contains
    real*8 function cost(X) result(err)

      implicit none

      real*8, intent(in)  :: X(2)

      ! Local variables
      real*8              :: angle_, kappa0_

      angle_ = X(1)
      kappa0_ = X(2) / lengthScale

      err = plvira_error(refVolumes, angle_, kappa0_, dxs)
    end function

    real*8 function cost_and_grad(grad, X) result(err)

      implicit none

      real*8, intent(out) :: grad(2)
      real*8, intent(in)  :: X(2)

      ! Local variables
      real*8              :: angle_, kappa0_

      angle_ = X(1)
      kappa0_ = X(2) / lengthScale

      err = plvira_error(refVolumes, angle_, kappa0_, dxs, derivatives=grad)
      grad(2) = grad(2) / lengthScale

    end function
  end function

  real*8 function lvira_angle_guess(refVolumes, dxs) result(angle)
    implicit none
    
    real*8, intent(in)    :: refVolumes(-1:1,-1:1), dxs(-1:1,2)

    ! Local variables 
    real*8                :: normal(2)

    normal(1) = (refVolumes(-1,0)/(dxs(-1,1)*dxs(0,2)) - refVolumes(1,0)/(dxs(1,1)*dxs(0,2)))/&
      (dxs(-1,1) + 2*dxs(0,1) + dxs(1,1))
    normal(2) = (refVolumes(0,-1)/(dxs(-1,2)*dxs(0,1)) - refVolumes(0,1)/(dxs(1,2)*dxs(0,1)))/&
      (dxs(-1,2) + 2*dxs(0,2) + dxs(1,2))
    angle = datan2(normal(2), normal(1))
  end function  

  real*8 function plvira_error(refVolumes, angle, kappa0, dxs, derivatives) result(err)
    use m_common
    use m_recon_util
    use m_polygon

    implicit none

    real*8, intent(in)    :: refVolumes(-1:1,-1:1), kappa0, dxs(-1:1,2), angle
    real*8, intent(out), optional :: derivatives(2)   ! w.r.t. angle and curvature respectively

    ! Local variables
    type(tPolygon)        :: cell
    type(tParabola)       :: parabola
    real*8                :: normal(2), shift, volume, derivatives_local(2), err_local
    real*8                :: xc_neighbour(2), dx_neighbour(2), cellVol_neighbour, grad_s(2)
    integer               :: i, j

    normal = [dcos(angle), dsin(angle)]
    shift = cmpShift(normal, dxs(0,:), refVolumes(0,0), kappa0, intersected=cell)

    parabola = makeParabola(normal, kappa0, shift)
    
    err = 0
    if (present(derivatives)) then 
      derivatives = 0
       
      ! Compute derivatives which ensure that the centred volume is conserved
      grad_s(1) = cmpDerivative_shiftAngle(cell)
      grad_s(2) = cmpDerivative_shiftKappa(cell)

    endif
    do j=-1,1
    do i=-1,1
      if (i==0 .and. j==0) cycle

      dx_neighbour(1) = dxs(i, 1)
      dx_neighbour(2) = dxs(j, 2)
      cellVol_neighbour = product(dx_neighbour)

      xc_neighbour(1) = i * (dxs(0,1) + dx_neighbour(1))/2
      xc_neighbour(2) = j * (dxs(0,2) + dx_neighbour(2))/2

      call makeBox(cell, xc_neighbour, dx_neighbour)
      call intersect(cell, parabola)

      volume = cmpVolume(cell)
      err_local = (volume - refVolumes(i,j)) / cellVol_neighbour

      err = err + err_local**2
      if (present(derivatives) .and. volume > 0 .and. volume < cellVol_neighbour) then
        derivatives_local(1) = cmpDerivative_volAngle(cell, shiftAngleDerivative=grad_s(1))
        derivatives_local(2) = cmpDerivative_volKappa(cell, shiftKappaDerivative=grad_s(2))

        derivatives(1) = derivatives(1) + 2 * err_local * derivatives_local(1) / cellVol_neighbour
        derivatives(2) = derivatives(2) + 2 * err_local * derivatives_local(2) / cellVol_neighbour
      endif
    enddo
    enddo
  end function

  function pmofNormal_rect(refMoments, kappa0, dx, verbose, errTol) result(normal)
    use m_optimization
    use m_recon_util

    implicit none

    real*8, intent(in)    :: refMoments(3), kappa0, dx(2)
    logical, optional     :: verbose
    real*8, intent(in), optional :: errTol
    real*8                :: normal(2)

    ! Local variables:
    real*8                :: cost_fun_scaling, centNorm, mofAngle
    real*8                :: cellVol, mofMoments_(3), errTol_
    logical               :: verbose_

    verbose_ = merge(verbose, .false., present(verbose))
    errTol_ = merge(errTol, 1D-8, present(errTol))

    cellVol = product(dx)
    cost_fun_scaling = cellVol**1.5D0

    ! Initial guess based on the reference centroid
    centNorm = norm2(refMoments(2:3))
    if (centNorm > 0.0) then
      mofAngle = datan2(-refMoments(3), -refMoments(2))
    else
      mofAngle = 0
    endif
    
    mofAngle =  brent_min(dcost, mofAngle, errTol_, 25, verbose_, maxStep=0.5D0)

    normal = [dcos(mofAngle), dsin(mofAngle)]
  contains 

    real*8 function dcost(angle, err) result(derr)
      use m_recon_util
      use m_polygon
  
      implicit none

      real*8, intent(in)  :: angle
      real*8, intent(out), optional :: err

      ! Local variables:
      real*8              :: difference(2), derivative(2), normal_(2), shift, err_
      type(tPolygon)      :: poly

      normal_ = [dcos(angle), dsin(angle)]

      shift = cmpShift(normal_, dx, refMoments(1), kappa0, intersected=poly)
      
      mofMoments_ = cmpMoments(poly)
      derivative = cmpDerivative_firstMomentAngle(poly)
      
      difference = (mofMoments_(2:3) - refMoments(2:3)) / cost_fun_scaling
      derivative = derivative / cost_fun_scaling
      
      err_ = sum(difference**2)
      derr = dot_product(derivative, difference)*2
      if (present(err)) err = err_
    end function

  end function

  function pmofNormal_poly(refMoments, kappa0, cell, x0, verbose, errTol) result(normal)
    use m_optimization
    use m_recon_util
    use m_polygon

    implicit none

    real*8, intent(in)    :: refMoments(3), kappa0
    type(tPolygon)        :: cell
    logical, optional     :: verbose
    real*8, intent(in), optional :: x0(2), errTol
    real*8                :: normal(2)

    ! Local variables:
    real*8                :: cost_fun_scaling, centNorm, mofAngle
    real*8                :: cellMoments(3), mofMoments_(3), errTol_
    logical               :: verbose_

    verbose_ = merge(verbose, .false., present(verbose))
    errTol_ = merge(errTol, 1D-8, present(errTol))

    cellMoments = cmpMoments(cell)
    cost_fun_scaling = cellMoments(1)**1.5D0

    ! Initial guess based on the reference centroid
    centNorm = norm2(refMoments(2:3))
    if (centNorm > 0.0) then
      mofAngle = datan2(-refMoments(3), -refMoments(2))
    else
      mofAngle = 0
    endif
    
    mofAngle =  brent_min(dcost, mofAngle, errTol_, 25, verbose_, maxStep=0.5D0)

    normal = [dcos(mofAngle), dsin(mofAngle)]

  contains 

    real*8 function dcost(angle, err) result(derr)
      use m_recon_util
      use m_r2d_parabolic
  
      implicit none

      real*8, intent(in)  :: angle
      real*8, intent(out), optional :: err

      ! Local variables:
      real*8              :: difference(2), derivative(2), normal_(2), shift, err_
      type(tPolygon)      :: intersected

      normal_ = [dcos(angle), dsin(angle)]

      shift = cmpShift(normal_, cell, refMoments(1), kappa0, intersected=intersected) !x0=x0)
      
      mofMoments_ = cmpMoments(intersected)

      derivative = cmpDerivative_firstMomentAngle(intersected)
      
      difference = (mofMoments_(2:3) - refMoments(2:3)) / cost_fun_scaling
      derivative = derivative / cost_fun_scaling
      
      err_ = sum(difference**2)
      derr = dot_product(derivative, difference)*2

      if (present(err)) err = err_
    end function

  end function


  !NB The 2D MoF optimization problem is solved exactly using a direct method
  ! See also "Moment-of-fluid analytic reconstruction on 2D Cartesian grids", JCP 2017
  function mofNormal(refMoments, dx, mofMoments) result(normal)

    implicit none

    real*8, intent(in)    :: refMoments(3), dx(2)
    real*8                :: normal(2)
    real*8, intent(out), optional :: mofMoments(3)

    ! Local variables
    integer               :: selectedCases(3), cdx
    real*8                :: tmpNormal(2), tmpCentroid(2), tmpError, mofCentroid_(2), centroidError_
    real*8                :: refFrac, refCentroid(2)
    logical               :: largerThanHalf

    normal = 0.0
    centroidError_ = norm2(dx) * 43  ! Larger than diameter of cell

    refFrac = refMoments(1) / product(dx)
    refCentroid = refMoments(2:3) / refMoments(1)
    
    largerThanHalf = refFrac > 0.5
    if (largerThanHalf) then
      refCentroid = -refCentroid * refFrac
      refFrac = 1.0D0 - refFrac
      refCentroid = refCentroid / refFrac
    endif

    ! There are 8 possible liquid configurations for which the optimization problem restricted to this configuration can be solved analytically
    ! Based on the centroid we may consider only 3 out of those 8 cases
    if (refCentroid(1) <= 0 .and. refCentroid(2) <= 0) then
      selectedCases = [1, 2, 8]
    elseif (refCentroid(1) > 0 .and. refCentroid(2) <= 0) then
      selectedCases = [2, 3, 4]
    elseif (refCentroid(1) > 0 .and. refCentroid(2) > 0) then
      selectedCases = [4, 5, 6]
    else
      selectedCases = [6, 7, 8]
    endif

    ! For each selected case we solve the optimization problem; thus yielding
    ! 3 candidate normal/pc. We then select the one which results in the
    ! smallest centroid error.
    do cdx=1,3
      ! All of the cases can be reduce to only 2 cases by appropriately
      ! transforming the problem
      select case(selectedCases(cdx))
      case(1)
        call solve_2dMoF_hyperbola(refFrac, refCentroid, dx, tmpNormal, tmpCentroid, tmpError)
      case(2)
        call solve_2dMoF_parabola(refFrac, refCentroid, dx, tmpNormal, tmpCentroid, tmpError)
      case(3)
        call solve_2dMoF_hyperbola(refFrac, [refCentroid(2), -refCentroid(1)], &
          [dx(2), dx(1)], tmpNormal, tmpCentroid, tmpError)
        tmpNormal = [-tmpNormal(2), tmpNormal(1)]
        tmpCentroid = [-tmpCentroid(2), tmpCentroid(1)]
      case(4)
        call solve_2dMoF_parabola(refFrac, [refCentroid(2), -refCentroid(1)], & 
          [dx(2), dx(1)], tmpNormal, tmpCentroid, tmpError)
        tmpNormal = [-tmpNormal(2), tmpNormal(1)]
        tmpCentroid = [-tmpCentroid(2), tmpCentroid(1)]
      case(5)
        call solve_2dMoF_hyperbola(refFrac, -refCentroid, dx, tmpNormal, tmpCentroid, tmpError)
        tmpNormal = -tmpNormal
        tmpCentroid = -tmpCentroid
      case(6)
        call solve_2dMoF_parabola(refFrac, -refCentroid, dx, tmpNormal, tmpCentroid, tmpError)
        tmpNormal = -tmpNormal
        tmpCentroid = -tmpCentroid
      case(7)
        call solve_2dMoF_hyperbola(refFrac, [-refCentroid(2), refCentroid(1)], [dx(2), dx(1)], &
          tmpNormal, tmpCentroid, tmpError)
        tmpNormal = [tmpNormal(2), -tmpNormal(1)]
        tmpCentroid = [tmpCentroid(2), -tmpCentroid(1)]
      case(8)
        call solve_2dMoF_parabola(refFrac, [-refCentroid(2), refCentroid(1)], [dx(2), dx(1)], &
          tmpNormal, tmpCentroid, tmpError)
        tmpNormal = [tmpNormal(2), -tmpNormal(1)]
        tmpCentroid = [tmpCentroid(2), -tmpCentroid(1)]
      end select
      if (tmpError < centroidError_) then
        normal = tmpNormal
        centroidError_ = tmpError
        mofCentroid_ = tmpCentroid
      endif
    end do

    if (largerThanHalf) then
      normal = -normal
      mofCentroid_ = -mofCentroid_ * refFrac / (1.0D0 - refFrac)
    endif

    if (present(mofMoments)) mofMoments = [refMoments(1), refMoments(1) * mofCentroid_]

  end function mofNormal

  subroutine solve_2dMoF_hyperbola(refFrac, refCentroid, dx, mofNormal, mofCentroid, mofError)
    use m_poly_roots, only: polynomial_roots_deg4

    implicit none

    real*8, intent(in)    :: refFrac, refCentroid(2), dx(2)
    real*8, intent(out)   :: mofNormal(2), mofCentroid(2), mofError

    ! Local variables
    integer               :: rdx
    real*8                :: tmpCentroid(2), tmpError, scaledAndShiftedCx
    real*8                :: minVal, maxVal, coeffs(5), aspectRatio
    real*8                :: roots_real(4), roots_imag(4)

    aspectRatio = dx(2) / dx(1)

    ! We solve for (centroid_x / dx(1) + 1/2) by nondimensionalising everything by dx(1)
    coeffs = [1.0D0, -(refCentroid(1) / dx(1) + 0.5D0), 0.0D0, 2 * refFrac * aspectRatio * &
      (refCentroid(2) / dx(1) + aspectRatio / 2) / 9, -((2.0D0 / 9.0D0) * refFrac * aspectRatio)**2]
    call polynomial_roots_deg4(coeffs, roots_real, roots_imag)

    minVal = 2 * refFrac / 3
    maxVal = 1.0D0 / 3.0D0

    ! Consider at most 4 roots
    scaledAndShiftedCx = 0
    mofError = norm2(dx) * 42  ! Larger than diameter of cell
    do rdx=1,4
      if (.not. isnan(roots_real(rdx)) .and. roots_imag(rdx) == 0 .and. &
        roots_real(rdx) >= minVal .and. roots_real(rdx) <= maxVal) then

        tmpCentroid = [dx(1) * (roots_real(rdx) - 0.5D0), dx(2) * ((2 * refFrac / (9 * (roots_real(rdx)))) - 0.5D0)]
        tmpError = norm2(tmpCentroid - refCentroid)
        if (tmpError < mofError) then
          scaledAndShiftedCx = roots_real(rdx)
          mofError = tmpError
          mofCentroid = tmpCentroid
        endif
      endif
    enddo
    if (mofError < norm2(dx) * 42) then
      mofNormal = [refFrac * aspectRatio, (9.0D0 / 2.0D0) * (scaledAndShiftedCx**2)]
      mofNormal = mofNormal / norm2(mofNormal)
    else
      mofNormal = 0
      mofCentroid = 0
    endif
  end subroutine

  subroutine solve_2dMoF_parabola(refFrac, refCentroid, dx, mofNormal, mofCentroid, mofError)
    use m_poly_roots, only: polynomial_roots_deg3
    implicit none

    real*8, intent(in)    :: refFrac, refCentroid(2), dx(2)
    real*8, intent(out)   :: mofNormal(2), mofCentroid(2), mofError

    ! Local variables
    integer               :: rdx
    real*8                :: tmpCentroid(2), tmpError, scaledCx
    real*8                :: minVal, maxVal, aspectRatio, coeffs(4)
    real*8                :: roots_real(3), roots_imag(3)

    aspectRatio = dx(2) / dx(1)

    ! We solve for (centroid_x / dx(1)) (this yields a dimensionless cubic equation)
    coeffs = [72 * (refFrac * aspectRatio)**2, 0.0D0, 1.0D0 + 6 * (refFrac * aspectRatio)**2 -& 
      12 * aspectRatio * refFrac * (refCentroid(2) / dx(1) + aspectRatio / 2), -refCentroid(1) / dx(1)]
    call polynomial_roots_deg3(coeffs, roots_real, roots_imag)

    minVal = -1.0D0 / 6.0D0
    maxVal = 1.0D0 / 6.0D0

    ! Consider at most 3 roots
    mofError = norm2(dx) * 42  ! Larger than diameter of cell
    scaledCx = 0.0

    do rdx=1,3
      if (.not. isnan(roots_real(rdx)) .and. roots_imag(rdx) == 0 .and. &
        roots_real(rdx) >= minVal .and. roots_real(rdx) <= maxVal) then
        tmpCentroid = [dx(1) * roots_real(rdx), dx(2) * (refFrac - 1) / 2 + &
          6 * refFrac * dx(2) * roots_real(rdx)**2]
        tmpError = norm2(tmpCentroid - refCentroid)
        if (tmpError < mofError) then
          scaledCx = roots_real(rdx)
          mofError = tmpError
          mofCentroid = tmpCentroid
        endif
      endif
    enddo

    if (mofError < norm2(dx) * 42) then
      mofNormal = [-aspectRatio * scaledCx * refFrac, 1.0D0 / 12]
      mofNormal = mofNormal / norm2(mofNormal)
    else
      mofNormal = 0
      mofCentroid = 0
    endif
  end subroutine
end module