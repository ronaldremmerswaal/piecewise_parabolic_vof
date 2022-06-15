module m_benchmark
  real*8, parameter         :: IS_INTERFACE_TOLERANCE = 1E-8
  real*8, parameter         :: RADIUS = 0.25

contains
  subroutine run_recon_benchmark()
    use m_recon_util
    use m_polygon
    use m_recon

    implicit none

    integer, parameter    :: NR_METHODS = 7, N = 2**11
    integer, parameter    :: LVIRA_IDX = 1, LVIRAP_IDX = 2, PLVIRA_IDX = 3, PROST_IDX = 4, MOF_IDX = 5, PMOF_IDX = 6, PMOFP_IDX = 7
    character*6, parameter :: METHODS(NR_METHODS) = ["LVIRA ", "LVIRAP", "PLVIRA", "PROST ", "MOF   ", "PMOF  ", "PMOFP "]
    ! LVIRAP refers to LVIRA using the PLVIRA implementation

    real*8                :: pi, Lx, Ly, xc(2), dx(2), dxs(-1:1,2)
    real*8                :: refMoments(3), normal(2), tmp, reconstruction_time(NR_METHODS)
    real*8                :: refVolumes(-1:1,-1:1), reconMoments(3), errTol, kappa0
    real*8                :: volume_grid(N,N), momx_grid(N,N), momy_grid(N,N), kappa0_exact(N,N)
    integer               :: i, j, mdx, nr_interface_cells
    type(tParabola)       :: parabola
    type(tPolygon)        :: cell

    pi = 4 * datan(1.0D0)
    
    Lx = 1.333333
    Ly = 0.966666

    dx = [Lx, Ly] / N
    dxs(:,1) = dx(1)
    dxs(:,2) = dx(2)

    ! Precompute reference moments
    nr_interface_cells = 0
    do j=1,N
    do i=1,N
      xc(1) = (i-0.5) * dx(1)
      xc(2) = (j-0.5) * dx(2)

      call polyApprox(cell, xc, dx, levelSet)
      refMoments = cmpMoments(cell)
      volume_grid(i,j) = refMoments(1)
      momx_grid(i,j) = refMoments(2) - refMoments(1) * xc(1)
      momy_grid(i,j) = refMoments(3) - refMoments(1) * xc(2)

      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        kappa0_exact(i,j) = levelSet_curvature(xc, minval(dx) / 5)
        nr_interface_cells = nr_interface_cells + 1
      endif
    enddo
    enddo

    errTol = min(1D-2, (minval(dx) / max(Lx, Ly))**2)

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refVolumes = volume_grid(i-1:i+1,j-1:j+1)

        normal = lviraNormal(refVolumes, dxs, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(LVIRA_IDX))
    reconstruction_time(LVIRA_IDX) = reconstruction_time(LVIRA_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refVolumes = volume_grid(i-1:i+1,j-1:j+1)
        kappa0 = 0
        normal = plviraNormal(refVolumes, kappa0, dxs, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(LVIRAP_IDX))
    reconstruction_time(LVIRAP_IDX) = reconstruction_time(LVIRAP_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refVolumes = volume_grid(i-1:i+1,j-1:j+1)
        kappa0 = kappa0_exact(i,j)

        normal = plviraNormal(refVolumes, kappa0, dxs, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(PLVIRA_IDX))
    reconstruction_time(PLVIRA_IDX) = reconstruction_time(PLVIRA_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refVolumes = volume_grid(i-1:i+1,j-1:j+1)
        kappa0 = kappa0_exact(i,j)

        normal = prostNormal(refVolumes, kappa0, dxs, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(PROST_IDX))
    reconstruction_time(PROST_IDX) = reconstruction_time(PROST_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refMoments(2) = momx_grid(i,j)
        refMoments(3) = momy_grid(i,j)

        normal = mofNormal(refMoments, dx)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(MOF_IDX))
    reconstruction_time(MOF_IDX) = reconstruction_time(MOF_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refMoments(2) = momx_grid(i,j)
        refMoments(3) = momy_grid(i,j)
        kappa0 = kappa0_exact(i,j)

        normal = pmofNormal(refMoments, kappa0, dx, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(PMOF_IDX))
    reconstruction_time(PMOF_IDX) = reconstruction_time(PMOF_IDX) - tmp

    call cpu_time(tmp)
    do j=2,N-1
    do i=2,N-1
      refMoments(1) = volume_grid(i,j)
      
      if (refMoments(1) > IS_INTERFACE_TOLERANCE .and. product(dx) - refMoments(1) > IS_INTERFACE_TOLERANCE) then
        refMoments(2) = momx_grid(i,j)
        refMoments(3) = momy_grid(i,j)
        kappa0 = kappa0_exact(i,j)

        call makeBox(cell, dx)
        normal = pmofNormal(refMoments, kappa0, cell, errTol=errTol)
      endif   
    enddo
    enddo
    call cpu_time(reconstruction_time(PMOFP_IDX))
    reconstruction_time(PMOFP_IDX) = reconstruction_time(PMOFP_IDX) - tmp

    reconstruction_time = reconstruction_time / nr_interface_cells

    print*, 'RECONSTRUCTION TIMING:'
    do mdx=1,NR_METHODS
      print*, '... ', METHODS(mdx), ' = ', reconstruction_time(mdx)
    enddo
  end subroutine

  real*8 function levelSet(x) result(ans)

    implicit none
    
    real*8, intent(in)  :: x(2)
    
    ! Local variables
    real*8, parameter   :: PERT_REL_AMP = 0.1, PERT_FREQ = 5., X0(2) = [0.5, 0.5], PERT_SHIFT = 0.1
    real*8              :: angle

    angle = datan2(x(2) - X0(2), x(1) - X0(1))
    ans = norm2(x - X0) - RADIUS * (1 + PERT_REL_AMP * dcos(PERT_SHIFT + angle * PERT_FREQ))

  end function

  real*8 function levelSet_curvature(x, h) result(kappa)
    use m_optimization

    implicit none

    real*8, intent(in)  :: x(2), h

    ! Local variables
    real*8, parameter   :: DSTEP = 1D-6
    real*8              :: lsNormal(2), lsTangent(2), lhfVals(-1:1)
    integer             :: ldx

    lsNormal = [levelSet(x + [1D-6,.0D0]) - levelSet(x - [1D-6,.0D0]), &
                levelSet(x + [.0D0,1D-6]) - levelSet(x - [.0D0,1D-6])]
    lsNormal = lsNormal / norm2(lsNormal)
    lsTangent = [-lsNormal(2), lsNormal(1)]

    ! Construct LHF in lsNormal direction
    do ldx=-1,1
      lhfVals(ldx) = brent(rootfun, -RADIUS/5, RADIUS/5, 1D-12, 52)
    enddo

    kappa = -((lhfVals(1) - 2 * lhfVals(0) + lhfVals(-1)) / h**2) / (sqrt(1 + ((lhfVals(1) - lhfVals(-1))/(2*h))**2)**3)
    
  contains
    real*8 function rootfun(s) result(err)
      implicit none
      
      real*8, intent(in)  :: s

      err = levelSet(x + lsNormal * s + lsTangent * ldx * h)
    end function
  end function
end module

program reconstruction_demo
  use m_benchmark

  call run_recon_benchmark()
end program