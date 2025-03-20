recursive subroutine gdist(coord1, coord2, n1, n2, proj, rescale, diag, nthreads, mdist) bind(C, name="gdist")
    use iso_c_binding, only: c_int, c_double
    use omp_lib
    implicit none
    integer(c_int), value :: n1, n2
    integer(c_int), value :: proj, rescale, diag
    integer(c_int), value :: nthreads
    real(c_double), intent(in) :: coord1(n1,2)
    real(c_double), intent(in) :: coord2(n2,2)
    real(c_double), intent(out) :: mdist

    ! Local variables
    integer(c_int) :: i, j, ndist
    real(c_double) :: diffx, diffy, dist, mdist1, mdist2, sdist
    real(c_double) :: a, c_val
    real(c_double), parameter :: radius = 6372.8_c_double

    ! Determine number of distances to compute
    if (diag /= 0) then
       ndist = sum((/(i, i=1,n1, 1)/))
    else
       ndist = n1 * n2
    end if

    mdist = 0.0_c_double

    if (proj /= 0) then  ! Euclidean distance
      if (diag /= 0) then
      !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
        do i = 1, n1
          sdist = 0.0_c_double
          do j = i, n2
            diffx = coord1(i,1)-coord2(j,1)
            diffy = coord1(i,2)-coord2(j,2)
            dist = sqrt(diffx*diffx + diffy*diffy)
            sdist = sdist + dist
          end do
          mdist = mdist + sdist / real(ndist, c_double)
        end do
      !$omp end parallel do
      else
      !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
        do i = 1, n1
          sdist = 0.0_c_double
          do j = 1, n2
            diffx = coord1(i,1)-coord2(j,1)
            diffy = coord1(i,2)-coord2(j,2)
            dist = sqrt(diffx*diffx + diffy*diffy)
            sdist = sdist + dist
          end do
          mdist = mdist + sdist / real(ndist, c_double)
        end do
      !$omp end parallel do
      end if
    else  ! Great-circle distance
      if (diag /= 0) then
      !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist, a, c_val) reduction(+:mdist)
        do i = 1, n1
          sdist = 0.0_c_double
          do j = i, n2
            diffx = coord1(i,1)-coord2(j,1)
            diffy = coord1(i,2)-coord2(j,2)
            a = (sin(diffx/2.0_c_double))**2 + cos(coord1(i,1))*cos(coord2(j,1))*(sin(diffy/2.0_c_double))**2
            c_val = 2.0_c_double * asin(sqrt(a))
            dist = radius * c_val
            sdist = sdist + dist
          end do
          mdist = mdist + sdist / real(ndist, c_double)
        end do
      !$omp end parallel do
      else
      !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist, a, c_val) reduction(+:mdist)
        do i = 1, n1
          sdist = 0.0_c_double
          do j = 1, n2
            diffx = coord1(i,1)-coord2(j,1)
            diffy = coord1(i,2)-coord2(j,2)
            a = (sin(diffx/2.0_c_double))**2 + cos(coord1(i,1))*cos(coord2(j,1))*(sin(diffy/2.0_c_double))**2
            c_val = 2.0_c_double * asin(sqrt(a))
            dist = radius * c_val
            sdist = sdist + dist
          end do
          mdist = mdist + sdist / real(ndist, c_double)
        end do
      !$omp end parallel do
      end if
    end if

    if (rescale /= 0) then
      call gdist(coord1, coord1, n1, n1, proj, 0, 1, nthreads, mdist1)
      call gdist(coord2, coord2, n2, n2, proj, 0, 1, nthreads, mdist2)
      mdist = mdist - 0.5_c_double * (mdist1 + mdist2)
    end if

  end subroutine gdist