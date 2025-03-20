module similarity_module
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic
  implicit none

  abstract interface
    function objective_fun(sim, obs, n) result(val)
      use iso_c_binding, only: c_int, c_double
      implicit none
      integer(c_int), intent(in) :: n
      real(c_double), intent(in) :: sim(n), obs(n)
      real(c_double) :: val
    end function objective_fun
  end interface

  procedure(objective_fun), pointer :: FUN => null()

contains

  subroutine set_objective_function(crit)
    implicit none
    integer(c_int), intent(in) :: crit

    select case (crit)
      case(1)
        FUN => kge
      case(2)
        FUN => rmse
      case(3)
        FUN => invkge
      case(4)
        FUN => invrmse
      case default
        FUN => invkge
    end select

  end subroutine set_objective_function

  subroutine similarity(Rn, nt, nb, crit, nthreads, sim_matrix) bind(C, name="similarity")
    use, intrinsic :: ieee_arithmetic
    use omp_lib
    implicit none
    integer(c_int), value :: nt, nb, crit, nthreads
    real(c_double), intent(in) :: Rn(nt,nb)
    real(c_double), intent(out) :: sim_matrix(nb,nb)

    integer(c_int) :: i, j
    real(c_double) :: val, nan_val
    nan_val = ieee_value(0.0_c_double, ieee_quiet_nan)

    call set_objective_function(crit)
    if (.not. associated(FUN)) then
      sim_matrix = nan_val
      return
    end if
    
    select case (crit) 
      case(2, 4) ! if symetric, like RMSE
        !$omp parallel do num_threads(nthreads) default(shared) private(i, j) schedule(static)
        do i = 1, nb
          do j = i, nb
            sim_matrix(i,j) = FUN(Rn(:,i), Rn(:,j), nt)
            sim_matrix(j,i) = sim_matrix(i,j)
          end do
        end do
        !$omp end parallel do
      case default ! if not symetric, like KGE
        !$omp parallel do num_threads(nthreads) default(shared) private(i, j) schedule(static)
        do i = 1, nb
          do j = 1, nb
            sim_matrix(i,j) = FUN(Rn(:,i), Rn(:,j), nt)
          end do
        end do
        !$omp end parallel do
    end select

  end subroutine similarity

  function kge(sim, obs, n) result(kge_val)
    use, intrinsic :: ieee_arithmetic
    implicit none
    integer, intent(in) :: n
    real(c_double), dimension(n), intent(in) :: sim, obs
    real(c_double) :: kge_val
    integer :: i, count
    real(c_double) :: sum_obs, sum_sim, sum_obs2, sum_sim2, sum_obs_sim
    real(c_double) :: mean_obs, mean_sim, std_obs, std_sim, r, beta, gamma, diff
    real(c_double) :: nan_val, eps
    eps = 1.0e-12_c_double

    count = 0
    sum_obs = 0.0_c_double
    sum_sim = 0.0_c_double
    sum_obs2 = 0.0_c_double
    sum_sim2 = 0.0_c_double
    sum_obs_sim = 0.0_c_double
    nan_val = ieee_value(0.0_c_double, ieee_quiet_nan)

    ! Loop to filter NaN values and accumulate sums
    do i = 1, n
      if (.not. ieee_is_nan(obs(i)) .and. .not. ieee_is_nan(sim(i))) then
        count = count + 1
        sum_obs = sum_obs + obs(i)
        sum_sim = sum_sim + sim(i)
        sum_obs2 = sum_obs2 + obs(i) * obs(i)
        sum_sim2 = sum_sim2 + sim(i) * sim(i)
        sum_obs_sim = sum_obs_sim + obs(i) * sim(i)
      end if
    end do

    ! Return NA if the common period is not long enough
    if (count <= 1) then
      kge_val = nan_val
      return
    end if

    ! Compute means
    mean_obs = sum_obs / count
    mean_sim = sum_sim / count

    ! Check to avoid division by zero
    if (abs(mean_obs) < eps) then
      kge_val = nan_val
      return
    end if

    ! Compute standard deviations (using Welford's correction)
    std_obs = sqrt((sum_obs2 - sum_obs * mean_obs) / (count - 1))
    std_sim = sqrt((sum_sim2 - sum_sim * mean_sim) / (count - 1))

    ! Check standard deviations
    if (std_obs < eps .or. std_sim < eps) then
      kge_val = nan_val
      return
    end if

    ! Compute correlation
    r = ((sum_obs_sim - count * mean_obs * mean_sim) / (count - 1)) / (std_obs * std_sim)

    ! Compute KGE components
    beta = mean_sim / mean_obs
    gamma = (std_sim / mean_sim) / (std_obs / mean_obs)

    ! Final KGE computation
    diff = (r - 1.0_c_double)**2 + (beta - 1.0_c_double)**2 + (gamma - 1.0_c_double)**2
    kge_val = 1.0_c_double - sqrt(diff)
    if (abs(kge_val - 1.0_c_double) < eps) kge_val = 1.0_c_double
  end function kge
  
  function rmse(sim, obs, n) result(res)
    integer, intent(in) :: n
    real(c_double), dimension(n), intent(in) :: sim, obs
    real(c_double) :: res
    integer :: i, count
    real(c_double) :: sum_sq, diff
    real(c_double) :: nan_val, eps
    nan_val = ieee_value(0.0_c_double, ieee_quiet_nan)
    eps = 1.0e-12_c_double

    count = 0
    sum_sq = 0.0_c_double
    do i = 1, n
      if (.not. ieee_is_nan(obs(i)) .and. .not. ieee_is_nan(sim(i))) then
        count = count + 1
        diff = sim(i) - obs(i)
        sum_sq = sum_sq + diff**2
      end if
    end do
    if (count == 0) then
      res = nan_val
      return
    end if
    res = sqrt(sum_sq / count)
    if (res < eps) res = 0.0_c_double
    
  end function rmse

  function invrmse(sim, obs, n) result(res)
    integer, intent(in) :: n
    real(c_double), dimension(n), intent(in) :: sim, obs
    real(c_double) :: rmse_val, res
    real(c_double) :: inf_val, eps
    inf_val = ieee_value(0.0_c_double, ieee_positive_inf)
    eps = 1.0e-12_c_double
    
    rmse_val = rmse(obs, sim, n)
    if (rmse_val < eps) then
      res = inf_val
    else
      res = 1.0_c_double / rmse_val
    end if

  end function invrmse

  function invkge(sim, obs, n) result(res)
    integer, intent(in) :: n
    real(c_double), dimension(n), intent(in) :: sim, obs
    real(c_double) :: kge_val, res
    real(c_double) :: inf_val, eps
    inf_val = ieee_value(0.0_c_double, ieee_positive_inf)
    eps = 1.0e-12_c_double

    kge_val = kge(obs, sim, n)
    if (abs(kge_val - 1.0_c_double) < eps) then
      res = inf_val
    else
      res = 1.0_c_double / (1-kge_val)
    end if

  end function invkge

end module similarity_module

