program reconstruction_demo

  write(*, '(A)') 'In a rectilinear control volume:'
  call reconstruction_rect

  write(*, '(A)') ''
  write(*, '(A)') 'In a nonrectilinear polygonal control volume:'
  call reconstruction_nonrect

contains

  subroutine reconstruction_rect
    ! Computation of moments / shift
    use m_recon_util,     only: cmpMoments, cmpShift, cmpSymmDiffVolume, makeParabola

    ! Tools for polygon intersection
    use m_polygon,        only: cmpMoments, makePlane, makeParabola, tParabola, tPolygon, polyApprox

    ! Reconstruction methods
    use m_recon,          only: mofNormal, pmofNormal

    implicit none

    ! The zeroth (refMoments(1)) and first moment (refMoments(2:3)) are stored in refMoments
    real*8                :: refMoments(3)

    ! The centroid and mesh widths of the rectangular control volume
    real*8                :: xc(2), dx(2)

    ! The PLIC/PPIC reconstruction is defined using the normal, shift and curvature:
    !   normal \cdot (x - x_c) - shift + (kappa0/2) * ((tangent \cdot (x - x_c)))^2
    real*8                :: normal(2), shift, kappa0

    ! We store the error in the moments as well as the moments of the symmetric difference
    real*8                :: errMoments(3), errSD(3)

    ! The reconstructed parabola
    type(tParabola)       :: parabola

    ! Temporary polygon
    type(tPolygon)        :: poly

    xc = [1/sqrt(2.), 1/sqrt(2.)]
    dx = [0.3, 0.3]

    ! Compute the reference zeroth and first moment (equivalent to the volume fraction and centroid)
    call polyApprox(poly, xc, dx, exact_interface)
    refMoments = cmpMoments(poly)

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
    errSD = cmpSymmDiffVolume(xc, dx, makePlane(normal, shift), exact_interface)

    write(*, '(A,1PD9.3,A,1PD9.3,A)') '... yields an interface normal (', normal(1), ', ', normal(2), ')'
    write(*, '(A,1PD9.3,A,1PD9.3,A,1PD9.3,A)') '... and zeroth and first moment error given by ', errMoments(1), &
      ' and (', errMoments(2), ', ', errMoments(3), '), respectively'
    write(*, '(A,1PD9.3)') '... whereas the area of the symmetric difference is given by ', errSD(1)
    

    ! We repeat the same steps for the PMOF method, for which we use an exact curvature
    kappa0 = 1.
    normal = pmofNormal(refMoments, kappa0, dx)
    parabola = makeParabola(normal, kappa0, dx, refMoments(1))
    errMoments = abs(cmpMoments(dx, parabola) - refMoments)
    errSD = cmpSymmDiffVolume(xc, dx, parabola, exact_interface)

    write(*, '(A)') ''
    write(*, '(A)') 'The PMOF method'
    write(*, '(A,1PD9.3,A,1PD9.3,A)') '... yields an interface normal (', normal(1), ', ', normal(2), ')'
    write(*, '(A,1PD9.3,A,1PD9.3,A,1PD9.3,A)') '... and zeroth and first moment error given by ', errMoments(1), &
      ' and (', errMoments(2), ', ', errMoments(3), '), respectively'
    write(*, '(A,1PD9.3)') '... whereas the area of the symmetric difference is given by ', errSD(1)

    ! Note that both methods conserve volume, but PMOF yields a better approximation;
    ! both in terms of the error of the first moment (or centroid) as well as the area
    ! of the symmetric difference
  end subroutine

  subroutine reconstruction_nonrect
    use m_recon_util, only: cmpMoments, cmpShift, cmpSymmDiff, makeParabola
    use m_r2d_parabolic,    only: cmpMoments, r2d_poly_f, init_from_pos, r2d_parabola_f
    use m_recon,   only: mofNormal, pmofNormal

    implicit none

    real*8                  :: refMoments(3)

    ! The control volume is now a polygon
    type(r2d_poly_f)        :: cell

    real*8                  :: pos(2,32), xc(2), pacmanRadius
    integer                 :: vdx, count
    logical                 :: midPointAdded

    real*8                  :: normal(2), kappa0, angle, pi
    real*8                  :: errMoments(3), errSD(3)

    type(r2d_parabola_f)    :: parabola

    pi = 4 * datan(1.0D0)
    pacmanRadius = 0.3D0
    ! Move inside such that interface intersects Pacman's mouth
    xc = [1/dsqrt(2.0D0), 1/dsqrt(2.0D0)] * (1.0D0 - pacmanRadius/2)

    ! We reconstruct inside a polygonal approximation of Pacman
    count = 1
    midPointAdded = .false.
    do vdx=1,size(pos,2)
      angle = 2 * pi * (vdx - 1.D0) / size(pos,2) - pi
      pos(1,count) = xc(1)
      pos(2,count) = xc(2)
      if (angle > pi/2 .or. angle < 0) then
        pos(1,count) = pos(1,count) + pacmanRadius * dcos(angle)
        pos(2,count) = pos(2,count) + pacmanRadius * dsin(angle)
        count = count + 1
      elseif (.not. midPointAdded) then
        count = count + 1
        midPointAdded = .true.
      endif
    enddo
    call init_from_pos(cell, pos(:,1:count))

    ! Other than the initialisation of the control volume, the rest is the same as before
    refMoments = cmpMoments(cell, exact_interface, verts_per_segment=1000)

    ! We pass the optional argoment x0=xc to tell the intersection functions that our
    ! parabola is defined relative to xc, rather than 0, which is the default
    kappa0 = 1.
    print*, 'TODO using m_polygon'
    ! normal = pmofNormal(refMoments, kappa0, cell, x0=xc)
    parabola = makeParabola(normal, kappa0, cell, refMoments(1), x0=xc)
    errMoments = abs(cmpMoments(cell, parabola, x0=xc) - refMoments)
    errSD = cmpSymmDiff(cell, parabola, exact_interface, x0=xc)

    write(*, '(A)') ''
    write(*, '(A)') 'The PMOF method can reconstruct inside non-rectilinear (and non-convex) cells as well, and'
    write(*, '(A,1PD9.3,A,1PD9.3,A)') '... yields an interface normal (', normal(1), ', ', normal(2), ')'
    write(*, '(A,1PD9.3,A,1PD9.3,A,1PD9.3,A)') '... and zeroth and first moment error given by ', errMoments(1), &
      ' and (', errMoments(2), ', ', errMoments(3), '), respectively'
    write(*, '(A,1PD9.3)') '... whereas the area of the symmetric difference is given by ', errSD(1)
  end subroutine


  real*8 function exact_interface(x) result(ans)
    implicit none

    real*8, intent(in)    :: x(2)

    ! The exact interface is defined as the unit circle
    ans = norm2(x) - 1.0
  end
end