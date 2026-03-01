! C-callable wrappers for clansvd and clansvd_irl.

module cpropack_c_mod
  use iso_c_binding
  implicit none

  abstract interface
     subroutine c_aprod_c(transa, m, n, x, y, userdata) bind(C)
       import :: c_char, c_int, c_float_complex, c_ptr
       character(c_char),      value :: transa
       integer(c_int),         value :: m, n
       complex(c_float_complex), intent(in)  :: x(*)
       complex(c_float_complex), intent(out) :: y(*)
       type(c_ptr),            value :: userdata
     end subroutine
  end interface

end module cpropack_c_mod


subroutine propack_clansvd(jobu, jobv, m, n, k, kmax, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, soption, ioption, info) &
     bind(C, name="propack_clansvd")
  use iso_c_binding
  use cpropack_c_mod
  implicit none

  character(c_char), value :: jobu, jobv
  integer(c_int),    value :: m, n, k, kmax, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  complex(c_float_complex), intent(inout) :: U(ldu,*)
  real(c_float), intent(out)   :: Sigma(*)
  real(c_float), intent(out)   :: bnd(*)
  complex(c_float_complex), intent(inout) :: V(ldv,*)
  real(c_float), value         :: tolin
  real(c_float), intent(in)    :: soption(*)
  integer(c_int), intent(in)   :: ioption(*)
  integer(c_int), intent(out)  :: info

  integer :: lswork, lcwork, liwork, nb, mn_max
  real(c_float),              allocatable :: swork(:)
  complex(c_float_complex),   allocatable :: cwork(:)
  integer,                    allocatable :: iwork(:)
  complex(c_float_complex) :: cparm(1)
  integer                  :: iparm(1)

  procedure(c_aprod_c), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lswork = m + n + 9*kmax + 2*kmax**2 + 4 + max(m+n, 4*kmax+4)
     lcwork = m + n
     liwork = 2*kmax + 1
  else
     lswork = m + n + 9*kmax + 5*kmax**2 + 4 + max(3*kmax**2+4*kmax+4, nb*mn_max)
     lcwork = m + n + nb*mn_max
     liwork = 8*kmax
  end if

  allocate(swork(lswork), cwork(lcwork), iwork(liwork))
  swork = 0.0
  cwork = (0.0, 0.0)
  iwork = 0

  call clansvd(jobu, jobv, m, n, k, kmax, aprod_adapter, &
       U, ldu, Sigma, bnd, V, ldv, tolin, &
       swork, lswork, cwork, lcwork, iwork, liwork, &
       soption, ioption, info, cparm, iparm)

  deallocate(swork, cwork, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, cparm, iparm)
    use iso_c_binding, only: c_char
    character*1               :: transa
    integer                   :: m, n, iparm(*)
    complex(c_float_complex)  :: x(*), y(*), cparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_clansvd


subroutine propack_clansvd_irl(which, jobu, jobv, m, n, &
     dim_in, p, nwanted, maxiter, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, soption, ioption, info) &
     bind(C, name="propack_clansvd_irl")
  use iso_c_binding
  use cpropack_c_mod
  implicit none

  character(c_char), value :: which, jobu, jobv
  integer(c_int),    value :: m, n, dim_in, p, nwanted, maxiter, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  complex(c_float_complex), intent(inout) :: U(ldu,*)
  real(c_float), intent(out)   :: Sigma(*)
  real(c_float), intent(out)   :: bnd(*)
  complex(c_float_complex), intent(inout) :: V(ldv,*)
  real(c_float), value         :: tolin
  real(c_float), intent(in)    :: soption(*)
  integer(c_int), intent(in)   :: ioption(*)
  integer(c_int), intent(out)  :: info

  integer :: dim, lswork, lcwork, liwork, nb, mn_max
  real(c_float),              allocatable :: swork(:)
  complex(c_float_complex),   allocatable :: cwork(:)
  integer,                    allocatable :: iwork(:)
  complex(c_float_complex) :: cparm(1)
  integer                  :: iparm(1)

  procedure(c_aprod_c), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  dim = dim_in

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lswork = m + n + 10*dim + 2*dim**2 + 5 + max(mn_max, 4*dim+4)
     lcwork = m + n
     liwork = 2*dim + 1
  else
     lswork = m + n + 10*dim + 5*dim**2 + 4 + max(3*dim**2+4*dim+4, nb*mn_max)
     lcwork = m + n + nb*mn_max
     liwork = 8*dim
  end if

  allocate(swork(lswork), cwork(lcwork), iwork(liwork))
  swork = 0.0
  cwork = (0.0, 0.0)
  iwork = 0

  call clansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter, &
       aprod_adapter, U, ldu, Sigma, bnd, V, ldv, tolin, &
       swork, lswork, cwork, lcwork, iwork, liwork, &
       soption, ioption, info, cparm, iparm)

  deallocate(swork, cwork, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, cparm, iparm)
    use iso_c_binding, only: c_char
    character*1               :: transa
    integer                   :: m, n, iparm(*)
    complex(c_float_complex)  :: x(*), y(*), cparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_clansvd_irl
