! C-callable wrappers for zlansvd and zlansvd_irl.

module zpropack_c_mod
  use iso_c_binding
  implicit none

  abstract interface
     subroutine c_aprod_z(transa, m, n, x, y, userdata) bind(C)
       import :: c_char, c_int, c_double_complex, c_ptr
       character(c_char),       value :: transa
       integer(c_int),          value :: m, n
       complex(c_double_complex), intent(in)  :: x(*)
       complex(c_double_complex), intent(out) :: y(*)
       type(c_ptr),             value :: userdata
     end subroutine
  end interface

end module zpropack_c_mod


subroutine propack_zlansvd(jobu, jobv, m, n, k, kmax, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_zlansvd")
  use iso_c_binding
  use zpropack_c_mod
  implicit none

  character(c_char), value :: jobu, jobv
  integer(c_int),    value :: m, n, k, kmax, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  complex(c_double_complex), intent(inout) :: U(ldu,*)
  real(c_double), intent(out)   :: Sigma(*)
  real(c_double), intent(out)   :: bnd(*)
  complex(c_double_complex), intent(inout) :: V(ldv,*)
  real(c_double), value         :: tolin
  real(c_double), intent(in)    :: doption(*)
  integer(c_int), intent(in)    :: ioption(*)
  integer(c_int), intent(out)   :: info

  integer :: lwork, lzwork, liwork, nb, mn_max
  real(c_double),              allocatable :: work(:)
  complex(c_double_complex),   allocatable :: zwork(:)
  integer,                     allocatable :: iwork(:)
  complex(c_double_complex) :: zparm(1)
  integer                   :: iparm(1)

  procedure(c_aprod_z), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lwork  = m + n + 9*kmax + 2*kmax**2 + 4 + max(m+n, 4*kmax+4)
     lzwork = m + n
     liwork = 2*kmax + 1
  else
     lwork  = m + n + 9*kmax + 5*kmax**2 + 4 + max(3*kmax**2+4*kmax+4, nb*mn_max)
     lzwork = m + n + nb*mn_max
     liwork = 8*kmax
  end if

  allocate(work(lwork), zwork(lzwork), iwork(liwork))
  work  = 0.0d0
  zwork = (0.0d0, 0.0d0)
  iwork = 0

  call zlansvd(jobu, jobv, m, n, k, kmax, aprod_adapter, &
       U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, zwork, lzwork, iwork, liwork, &
       doption, ioption, info, zparm, iparm)

  deallocate(work, zwork, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, zparm, iparm)
    use iso_c_binding, only: c_char
    character*1                :: transa
    integer                    :: m, n, iparm(*)
    complex(c_double_complex)  :: x(*), y(*), zparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_zlansvd


subroutine propack_zlansvd_irl(which, jobu, jobv, m, n, &
     dim_in, p, nwanted, maxiter, &
     c_aprod, userdata, U, ldu, Sigma, bnd, V, ldv, &
     tolin, doption, ioption, info) &
     bind(C, name="propack_zlansvd_irl")
  use iso_c_binding
  use zpropack_c_mod
  implicit none

  character(c_char), value :: which, jobu, jobv
  integer(c_int),    value :: m, n, dim_in, p, nwanted, maxiter, ldu, ldv
  type(c_funptr),    value :: c_aprod
  type(c_ptr),       value :: userdata
  complex(c_double_complex), intent(inout) :: U(ldu,*)
  real(c_double), intent(out)   :: Sigma(*)
  real(c_double), intent(out)   :: bnd(*)
  complex(c_double_complex), intent(inout) :: V(ldv,*)
  real(c_double), value         :: tolin
  real(c_double), intent(in)    :: doption(*)
  integer(c_int), intent(in)    :: ioption(*)
  integer(c_int), intent(out)   :: info

  integer :: dim, lwork, lzwork, liwork, nb, mn_max
  real(c_double),              allocatable :: work(:)
  complex(c_double_complex),   allocatable :: zwork(:)
  integer,                     allocatable :: iwork(:)
  complex(c_double_complex) :: zparm(1)
  integer                   :: iparm(1)

  procedure(c_aprod_z), pointer :: fp
  type(c_ptr) :: ud

  call c_f_procpointer(c_aprod, fp)
  ud = userdata

  dim = dim_in

  nb = 32
  mn_max = max(m, n)

  if (jobu == 'N' .and. jobv == 'N') then
     lwork  = m + n + 10*dim + 2*dim**2 + 5 + max(mn_max, 4*dim+4)
     lzwork = m + n
     liwork = 2*dim + 1
  else
     lwork  = m + n + 10*dim + 5*dim**2 + 4 + max(3*dim**2+4*dim+4, nb*mn_max)
     lzwork = m + n + nb*mn_max
     liwork = 8*dim
  end if

  allocate(work(lwork), zwork(lzwork), iwork(liwork))
  work  = 0.0d0
  zwork = (0.0d0, 0.0d0)
  iwork = 0

  call zlansvd_irl(which, jobu, jobv, m, n, dim, p, nwanted, maxiter, &
       aprod_adapter, U, ldu, Sigma, bnd, V, ldv, tolin, &
       work, lwork, zwork, lzwork, iwork, liwork, &
       doption, ioption, info, zparm, iparm)

  deallocate(work, zwork, iwork)

contains

  subroutine aprod_adapter(transa, m, n, x, y, zparm, iparm)
    use iso_c_binding, only: c_char
    character*1                :: transa
    integer                    :: m, n, iparm(*)
    complex(c_double_complex)  :: x(*), y(*), zparm(*)
    character(c_char) :: transa_c
    transa_c = transa
    call fp(transa_c, m, n, x, y, ud)
  end subroutine

end subroutine propack_zlansvd_irl
