! C-callable wrappers for slansvd and slansvd_irl.

module spropack_c_mod
  use iso_c_binding
  implicit none

  abstract interface
     subroutine c_aprod_s(transa, m, n, x, y, userdata) bind(C)
       import :: c_char, c_int, c_float, c_ptr
       character(c_char), value :: transa
       integer(c_int),    value :: m, n
       real(c_float), intent(in)  :: x(*)
       real(c_float), intent(out) :: y(*)
       type(c_ptr),       value :: userdata
     end subroutine
  end interface

end module spropack_c_mod


subroutine propack_slansvd(jobu, jobv, m, n, k, kmax, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_slansvd")
  use iso_c_binding
  use spropack_c_mod
  implicit none

  character(c_char), value :: jobu, jobv
  integer(c_int),    value :: m, n, k, kmax, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  real(c_float), intent(inout) :: U(ldu,*)
  real(c_float), intent(out)   :: Sigma(*)
  real(c_float), intent(out)   :: bnd(*)
  real(c_float), intent(inout) :: V(ldv,*)
  real(c_float), value         :: tolin
  real(c_float), intent(in)    :: doption(*)
  integer(c_int), intent(in)   :: ioption(*)
  integer(c_int), intent(out)  :: info

  integer :: lwork, liwork, nb, mn_max
  real(c_float), allocatable :: work(:)
  integer,       allocatable :: iwork(:)
  real(c_float) :: dparm(1)
  integer       :: iparm(1)

  procedure(c_aprod_s), pointer :: fp
  type(c_ptr) :: ud

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
  work  = 0.0
  iwork = 0

  call slansvd(jobu, jobv, m, n, k, kmax, aprod_adapter, &
       U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, iwork, liwork, doption, ioption, info, dparm, iparm)

  deallocate(work, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, dparm, iparm)
    use iso_c_binding, only: c_char
    character*1   :: transa
    integer       :: m, n, iparm(*)
    real(c_float) :: x(*), y(*), dparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_slansvd


subroutine propack_slansvd_irl(which, jobu, jobv, m, n, &
     dim_in, p, nwanted, maxiter, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_slansvd_irl")
  use iso_c_binding
  use spropack_c_mod
  implicit none

  character(c_char), value :: which, jobu, jobv
  integer(c_int),    value :: m, n, dim_in, p, nwanted, maxiter, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  real(c_float), intent(inout) :: U(ldu,*)
  real(c_float), intent(out)   :: Sigma(*)
  real(c_float), intent(out)   :: bnd(*)
  real(c_float), intent(inout) :: V(ldv,*)
  real(c_float), value         :: tolin
  real(c_float), intent(in)    :: doption(*)
  integer(c_int), intent(in)   :: ioption(*)
  integer(c_int), intent(out)  :: info

  integer :: dim, lwork, liwork, nb, mn_max
  real(c_float), allocatable :: work(:)
  integer,       allocatable :: iwork(:)
  real(c_float) :: dparm(1)
  integer       :: iparm(1)

  procedure(c_aprod_s), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

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
  work  = 0.0
  iwork = 0

  call slansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter, &
       aprod_adapter, U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, iwork, liwork, doption, ioption, info, dparm, iparm)

  deallocate(work, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, dparm, iparm)
    use iso_c_binding, only: c_char
    character*1   :: transa
    integer       :: m, n, iparm(*)
    real(c_float) :: x(*), y(*), dparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_slansvd_irl
