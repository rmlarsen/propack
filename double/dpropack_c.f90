! C-callable wrappers for dlansvd and dlansvd_irl.
! Each wrapper allocates workspace internally and translates a C function
! pointer + void* userdata pair into the Fortran APROD interface expected
! by the original drivers.

module dpropack_c_mod
  use iso_c_binding
  implicit none

  abstract interface
     subroutine c_aprod_d(transa, m, n, x, y, userdata) bind(C)
       import :: c_char, c_int, c_double, c_ptr
       character(c_char), value :: transa
       integer(c_int),    value :: m, n
       real(c_double), intent(in)  :: x(*)
       real(c_double), intent(out) :: y(*)
       type(c_ptr),       value :: userdata
     end subroutine
  end interface

end module dpropack_c_mod


subroutine propack_dlansvd(jobu, jobv, m, n, k, kmax, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_dlansvd")
  use iso_c_binding
  use dpropack_c_mod
  implicit none

  character(c_char), value :: jobu, jobv
  integer(c_int),    value :: m, n, k, kmax, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  real(c_double), intent(inout) :: U(ldu,*)
  real(c_double), intent(out)   :: Sigma(*)
  real(c_double), intent(out)   :: bnd(*)
  real(c_double), intent(inout) :: V(ldv,*)
  real(c_double), value         :: tolin
  real(c_double), intent(in)    :: doption(*)
  integer(c_int), intent(in)    :: ioption(*)
  integer(c_int), intent(out)   :: info

  ! Local variables
  integer :: lwork, liwork, nb, mn_max
  real(c_double), allocatable :: work(:)
  integer,        allocatable :: iwork(:)
  real(c_double) :: dparm(1)
  integer        :: iparm(1)

  procedure(c_aprod_d), pointer :: fp
  type(c_ptr) :: ud

  ! Save callback info for host association
  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lwork  = m + n + 9*kmax + 2*kmax**2 + 4 + max(m+n, 4*kmax+4)
     liwork = 2*kmax + 1
  else
     lwork  = m + n + 9*kmax + 5*kmax**2 + 4 + max(3*kmax**2+4*kmax+4, nb*mn_max)
     liwork = 8*kmax
  end if

  allocate(work(lwork), iwork(liwork))
  work  = 0.0d0
  iwork = 0

  call dlansvd(jobu, jobv, m, n, k, kmax, aprod_adapter, &
       U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, iwork, liwork, doption, ioption, info, dparm, iparm)

  deallocate(work, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, dparm, iparm)
    character*1    :: transa
    integer        :: m, n, iparm(*)
    real(c_double) :: x(*), y(*), dparm(*)
    call fp(transa, m, n, x, y, ud)
  end subroutine

end subroutine propack_dlansvd


subroutine propack_dlansvd_irl(which, jobu, jobv, m, n, &
     dim_in, p, nwanted, maxiter, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_dlansvd_irl")
  use iso_c_binding
  use dpropack_c_mod
  implicit none

  character(c_char), value :: which, jobu, jobv
  integer(c_int),    value :: m, n, dim_in, p, nwanted, maxiter, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  real(c_double), intent(inout) :: U(ldu,*)
  real(c_double), intent(out)   :: Sigma(*)
  real(c_double), intent(out)   :: bnd(*)
  real(c_double), intent(inout) :: V(ldv,*)
  real(c_double), value         :: tolin
  real(c_double), intent(in)    :: doption(*)
  integer(c_int), intent(in)    :: ioption(*)
  integer(c_int), intent(out)   :: info

  ! Local variables
  integer :: dim, lwork, liwork, nb, mn_max
  real(c_double), allocatable :: work(:)
  integer,        allocatable :: iwork(:)
  real(c_double) :: dparm(1)
  integer        :: iparm(1)

  procedure(c_aprod_d), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  ! Copy dim so the driver can modify it internally
  dim = dim_in

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lwork  = m + n + 10*dim + 2*dim**2 + 5 + max(mn_max, 4*dim+4)
     liwork = 2*dim + 1
  else
     lwork  = m + n + 10*dim + 5*dim**2 + 4 + max(3*dim**2+4*dim+4, nb*mn_max)
     liwork = 8*dim
  end if

  allocate(work(lwork), iwork(liwork))
  work  = 0.0d0
  iwork = 0

  call dlansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter, &
       aprod_adapter, U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, iwork, liwork, doption, ioption, info, dparm, iparm)

  deallocate(work, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, dparm, iparm)
    character*1    :: transa
    integer        :: m, n, iparm(*)
    real(c_double) :: x(*), y(*), dparm(*)
    call fp(transa, m, n, x, y, ud)
  end subroutine

end subroutine propack_dlansvd_irl
