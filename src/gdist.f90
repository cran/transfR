!For compiling the Fortran file :
!R CMD SHLIB gdist.f90

recursive subroutine gdist(coord1,coord2,n1,n2,proj,rescale,diag,nthreads,mdist)
    use omp_lib

    implicit none

    integer, intent(in) :: n1, n2
    double precision, intent(in)  :: coord1(n1,2)
    double precision, intent(in)  :: coord2(n2,2)
    logical, intent(in)  :: proj, rescale, diag
    integer, intent(inout) :: nthreads
    double precision, intent(out) :: mdist


    double precision :: dist, diffx, diffy, a, c
    double precision :: mdist1, mdist2, sdist
    real, parameter :: radius = 6372.8
    integer :: i,j,ndist

!    if(size(coord1,1).LT.size(coord2,1)) then
!        call gdist(coord2,coord1,n2,n1,proj,rescale,diag,mdist)
!        return
!    end if

    if(diag) then
        ndist=sum((/(i, i=1,n1, 1)/))
    else
        ndist=n1*n2
    end if
!    if(ndist.GE.1e9) WRITE(*,'(A)') 'More than 1e9 distances to compute. Precision of the results not reliable.'

    mdist = 0.
    if(proj) then !euclidian distance
        if(diag) then
            !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
            do i=1,n1
                sdist = 0.
                do j=i,n2
                    diffx = coord1(i,1)-coord2(j,1)
                    diffy = coord1(i,2)-coord2(j,2)
                    dist = sqrt(diffx*diffx+diffy*diffy)
                    sdist = sdist+dist
!                    mdist = mdist+dist/ndist
                end do
                mdist = mdist+sdist/ndist
            end do
            !$omp end parallel do
        else
            !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
            do i=1,n1
                sdist = 0.
                do j=1,n2
                    diffx = coord1(i,1)-coord2(j,1)
                    diffy = coord1(i,2)-coord2(j,2)
                    dist = sqrt(diffx*diffx+diffy*diffy)
                    sdist = sdist+dist
!                    mdist = mdist+dist/ndist
                end do
                mdist = mdist+sdist/ndist
            end do
            !$omp end parallel do
        end if
    else !great-circle distance
        if(diag) then
            !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
            do i=1,n1
                sdist = 0.
                do j=i,n2
                    diffx = coord1(i,1)-coord2(j,1)
                    diffy = coord1(i,2)-coord2(j,2)
                    a = (sin(diffx/2))**2 + cos(coord1(i,1))*cos(coord2(j,1))*(sin(diffy/2))**2
                    c = 2*asin(sqrt(a))
                    dist = radius*c
                    sdist = sdist+dist
                    !mdist = mdist+dist/ndist
                end do
                mdist = mdist+sdist/ndist
            end do
            !$omp end parallel do
        else
            !$omp parallel do num_threads(nthreads) private(i, j, diffx, diffy, dist, sdist) reduction(+:mdist)
            do i=1,n1
                sdist = 0.
                do j=1,n2
                    diffx = coord1(i,1)-coord2(j,1)
                    diffy = coord1(i,2)-coord2(j,2)
                    a = (sin(diffx/2))**2 + cos(coord1(i,1))*cos(coord2(j,1))*(sin(diffy/2))**2
                    c = 2*asin(sqrt(a))
                    dist = radius*c
                    sdist = sdist+dist
!                    mdist = mdist+dist/ndist
                end do
                mdist = mdist+sdist/ndist
           end do
           !$omp end parallel do
        end if
    end if

    if(rescale) then
        call gdist(coord1,coord1,n1,n1,proj,.FALSE.,.TRUE.,nthreads,mdist1)
        call gdist(coord2,coord2,n2,n2,proj,.FALSE.,.TRUE.,nthreads,mdist2)
!        write(*,'(2I6,3F20.1)') n1, n2, mdist, mdist1, mdist2
        mdist = mdist - 0.5*(mdist1+mdist2)
    end if

end subroutine
