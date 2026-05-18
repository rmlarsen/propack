! Global Lanczos / SVD counters and timers shared across all PROPACK
! precision variants. Replaces the legacy include 'stat.h' COMMON block.
! Variable names are unchanged from stat.h so callers only need a `use
! stat_mod` instead of the include.
!
! (C) Rasmus Munk Larsen, Stanford University, 2000

module stat_mod
   implicit none
   public

   ! Counters
   integer :: nopx, nreorth, ndot, nrestart, nbsvd, nlandim, nsing

   ! Timers
   real :: tmvopx, tgetu0, tlanbpro, tbsvd, tlansvd, tritzvec, &
           trestart, treorth
end module stat_mod
