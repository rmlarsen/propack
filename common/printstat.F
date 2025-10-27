c
c     (C) Rasmus Munk Larsen, Stanford University, 2000
c
      subroutine clearstat
      implicit none
      include '../common/stat.h'
c     Reset global counters
      nopx = 0
      nreorth = 0
      ndot = 0
      nbsvd = 0
      nrestart = 0
      nlandim = 0
      nsing = 0
c     Reset global timers
      tmvopx = 0
      tgetu0 = 0
      tlanbpro = 0
      tbsvd = 0
      tlansvd = 0
      tritzvec = 0
      trestart = 0
      treorth = 0
      end

      subroutine getstat(counters, timers)
      implicit none
      integer counters(*)
      real timers(*)
      include '../common/stat.h'
c     Return global counter values.
      counters(1) = nopx
      counters(2) = nreorth
      counters(3) = ndot
      counters(4) = nbsvd
      counters(5) = nrestart
      counters(6) = nlandim
      counters(7) = nsing
c     Return global timer values.
      timers(1) = tmvopx
      timers(2) = tgetu0
      timers(3) = tlanbpro
      timers(4) = tbsvd
      timers(5) = tlansvd
      timers(6) = tritzvec
      timers(7) = trestart
      timers(8) = treorth
      end

      subroutine printstat
      implicit none
      include '../common/stat.h'

      print *,'+---------------------------------------------------',
     c     '--------+'
      print *,'Dimension of Lanczos basis                  = ',nlandim
      print *,'Number of singular values requested         = ',nsing
      print *,'Number of restarts                          = ',nrestart
      print *,'Number of matrix-vector multiplications     = ',nopx
      print *,'Number of reorthogonalizations              = ',nreorth
      print *,'Number of inner products in reorth.         = ',ndot
      print *,'Number of bidiagonal SVDs calculated        = ',nbsvd
      print *
      print *,'Total Time spent in LANBPRO                 = ',tlanbpro
      print *,'  Time spent doing matrix-vector multiply   = ',tmvopx
      print *,'  Time spent generating starting vectors    = ',tgetu0
      print *,'  Time spent reorthogonalizing              = ',treorth
      print *,'Time spent computing bidiagonal SVDs        = ',tbsvd
      print *,'Time spent doing implicit restarts          = ',trestart
      print *,'Time spent computing Ritz vectors           = ',tritzvec
      print *
      print *,'Total Time spent in LANSVD                  = ',tlansvd
      print *,'+----------------------------------------------------',
     c     '-------+'
      end


