subroutine convolution(rn, uh, nrn, nuh, q) bind(C, name="convolution")
  use iso_c_binding, only: c_int, c_double
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer(c_int), value :: nrn, nuh
  real(c_double), intent(in)  :: rn(nrn) 
  real(c_double), intent(in)  :: uh(nuh)
  real(c_double), intent(out) :: q(nrn+nuh)
  real(c_double) :: nan_val
  integer :: t

  nan_val = ieee_value(0.0_c_double, ieee_quiet_nan)
  q = 0.0_c_double

  do t = 1, nrn
    if (ieee_is_nan(rn(t))) then
      q(t:t+nuh-1) = nan_val
    else if (rn(t) > 0.0_c_double) then
      q(t:t+nuh-1) = q(t:t+nuh-1) + rn(t) * uh(1:nuh)
    end if
  end do

end subroutine convolution
