module m_reconstruction

  private
  public :: plic_normal_mof2d, ppic_normal_pmof2d

contains

  function ppic_normal_pmof2d(refMoments, kappa0, dx, verbose, errTol) result(mofNormal)
    use m_optimization
    use m_ppic_util

    implicit none

    real*8, intent(in)    :: refMoments(3), kappa0, dx(2)
    real*8                :: mofNormal(2)
    logical, optional     :: verbose
    real*8, intent(in), optional :: errTol

    ! Local variables:
    type(optimOpts)       :: LBFGS_OPTIONS = optimOpts(fTol=1D-12)  ! NOTE: defaults are used as set in m_optimization
    type(lsOpts)          :: MT_OPTIONS = lsOpts()
    type(optimInfo)       :: info

    real*8                :: cost_fun_scaling, mofMoments_(3), centNorm, mofAngle(1), tmp(1)
    real*8                :: refMoments_(3), cellVol
    logical               :: largerThanHalf

    if (present(verbose)) LBFGS_OPTIONS%verbose = verbose
    if (present(errTol))  LBFGS_OPTIONS%errTol = errTol

    cellVol = product(dx)
    largerThanHalf = refMoments(1) > cellVol/2
    if (.not. largerThanHalf) then
      refMoments_ = refMoments
    else
      refMoments_(1) = cellVol - refMoments(1)
      refMoments_(2:3) = -refMoments(2:3)
    endif

    cost_fun_scaling = cellVol**1.5D0

    ! Initial guess based on the reference centroid
    centNorm = norm2(refMoments_(2:3))
    if (centNorm > 0.0) then
      mofNormal = -refMoments_(2:3) / centNorm
    else
      mofNormal = [1.0, 0.0]
    endif
    
    mofAngle = datan2(mofNormal(2), mofNormal(1))

    call optimize(mofAngle, cost, LBFGS_OPTIONS, info, fun_and_grad=cost_fun_and_grad, ls_opts=MT_OPTIONS)

    mofNormal = [dcos(mofAngle(1)), dsin(mofAngle(1))]
    if (largerThanHalf) then
      mofNormal = -mofNormal
    endif

  contains 

    real*8 function cost(angle) result(err)
      use m_ppic_util
      implicit none

      real*8, intent(in)    :: angle(1)
      real*8                :: difference(2)

      ! Local variables
      real*8                :: normal(2), shift

      normal = [cos(angle(1)), sin(angle(1))]
      shift = cmpShift2d_parabolic(normal, dx, refMoments_(1), kappa0, moments=mofMoments_)

      difference = (mofMoments_(2:3) - refMoments_(2:3)) / cost_fun_scaling

      err = norm2(difference)
    end function

    real*8 function cost_fun_and_grad(grad, angle) result(err)
      use m_ppic_util
      use m_r2d_parabolic
      implicit none

      real*8, intent(out) :: grad(1)
      real*8, intent(in)  :: angle(1)

      ! Local variables:
      type(r2d_poly_f)    :: poly
      real*8              :: difference(2), derivative(4), normal(2), shift

      normal = [cos(angle(1)), sin(angle(1))]
      shift = cmpShift2d_parabolic(normal, dx, refMoments_(1), kappa0)

      call init_box(poly, [-dx/2, dx/2])
      call intersect_with_parabola(mofMoments_, poly, normal, kappa0, normal * shift, derivative)
      difference = (mofMoments_(2:3) - refMoments_(2:3)) / cost_fun_scaling
      derivative(2:3) = derivative(2:3) / cost_fun_scaling

      err = norm2(difference)
      grad = dot_product(derivative(2:3), difference) / err
    end function

  end function


  !NB The 2D MoF optimization problem is solved exactly using a direct method
  ! See also "Moment-of-fluid analytic reconstruction on 2D Cartesian grids", JCP 2017
  function plic_normal_mof2d(refMoments, dx, mofMoments) result(mofNormal)

    implicit none

    real*8, intent(in)    :: refMoments(3), dx(2)
    real*8                :: mofNormal(2)
    real*8, intent(out), optional :: mofMoments(3)

    ! Local variables
    integer               :: selectedCases(3), cdx
    real*8                :: tmpNormal(2), tmpCentroid(2), tmpError, mofCentroid_(2), centroidError_
    real*8                :: refFrac, refCentroid(2)
    logical               :: largerThanHalf

    mofNormal = 0.0
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
        mofNormal = tmpNormal
        centroidError_ = tmpError
        mofCentroid_ = tmpCentroid
      endif
    end do

    if (largerThanHalf) then
      mofNormal = -mofNormal
      mofCentroid_ = -mofCentroid_ * refFrac / (1.0D0 - refFrac)
    endif

    if (present(mofMoments)) mofMoments = [refMoments(1), refMoments(1) * mofCentroid_]

  end function plic_normal_mof2d

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