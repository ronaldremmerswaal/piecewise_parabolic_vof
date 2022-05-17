program reconstruction_demo
  ! Computation of moments / shift
  use m_reconstruction_util, only: cmpMoments, cmpShift, cmpSymmDiff

  ! Tools for polygon intersection
  use m_r2d_parabolic,    only: cmpMoments

  ! Reconstruction methods
  use m_reconstruction,   only: mofNormal, pmofNormal

  implicit none

  ! The zeroth (refMoments(1)) and first moment (refMoments(2:3)) are stored in refMoments
  real*8                  :: refMoments(3)

  ! The centroid and mesh widths of the rectangular control volume
  real*8                  :: xc(2), dx(2)

  ! The PLIC/PPIC reconstruction is defined using the normal, shift and curvature:
  !   normal \cdot (x - x_c) - shift + (kappa0/2) * ((tangent \cdot (x - x_c)))^2
  real*8                  :: normal(2), shift, kappa0

  ! We store the error in the moments as well as the moments of the symmetric difference
  real*8                  :: errMoments(3), errSD(3)

  xc = [1/sqrt(2.), 1/sqrt(2.)]
  dx = [0.3, 0.3]

  ! Compute the reference zeroth and first moment (equivalent to the volume fraction and centroid)
  refMoments = cmpMoments(xc, dx, exact_interface)

  ! Note that the moments should be relative to xc:
  refMoments(2:3) = refMoments(2:3) - xc * refMoments(1)

  ! First we try a linear reconstruction, for which kappa0 = 0
  write(*, '(A)') 'The MOF method'
  normal = mofNormal(refMoments, dx)

  ! The shift is computed by enforcing volume conservation
  shift = cmpShift(normal, dx, refMoments(1))

  ! ... hence we expect the zeroth moment have zero error
  errMoments = abs(cmpMoments(normal, dx, shift) - refMoments)

  ! The symmetric difference between the exact and approximate liquid domain can be approximated as follows
  errSD = cmpSymmDiff(xc, dx, normal, shift, exact_interface)

  write(*, '(A,1PD9.3,A,1PD9.3,A)') '... yields an interface normal (', normal(1), ', ', normal(2), ')'
  write(*, '(A,1PD9.3,A,1PD9.3,A,1PD9.3,A)') '... and zeroth and first moment error given by ', errMoments(1), &
    ' and (', errMoments(2), ', ', errMoments(3), '), respectively'
  write(*, '(A,1PD9.3)') '... whereas the area of the symmetric difference is given by ', errSD(1)
  

  ! We repeat the same steps for the PMOF method, for which we use an exact curvature
  kappa0 = 1.
  normal = pmofNormal(refMoments, kappa0, dx)
  shift = cmpShift(normal, dx, refMoments(1), kappa0)
  errMoments = abs(cmpMoments(normal, dx, shift, kappa0) - refMoments)
  errSD = cmpSymmDiff(xc, dx, normal, shift, kappa0, exact_interface)

  write(*, '(A)') ''
  write(*, '(A)') 'The PMOF method'
  write(*, '(A,1PD9.3,A,1PD9.3,A)') '... yields an interface normal (', normal(1), ', ', normal(2), ')'
  write(*, '(A,1PD9.3,A,1PD9.3,A,1PD9.3,A)') '... and zeroth and first moment error given by ', errMoments(1), &
    ' and (', errMoments(2), ', ', errMoments(3), '), respectively'
  write(*, '(A,1PD9.3)') '... whereas the area of the symmetric difference is given by ', errSD(1)

  ! Note that both methods conserve volume, but PMOF yields a better approximation;
  ! both in terms of the error of the first moment (or centroid) as well as the area
  ! of the symmetric difference
contains

  real*8 function exact_interface(x) result(ans)
  implicit none

  real*8, intent(in)    :: x(2)

    ! The exact interface is defined as the unit circle
    ans = norm2(x) - 1.0
  end
end