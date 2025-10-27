c
c     (C) Rasmus Munk Larsen, Stanford University, 2000
c

      integer nopx, nreorth, ndot, nrestart, nbsvd, nlandim, nsing
      real   tmvopx, tgetu0, tlanbpro, tbsvd, tlansvd, tritzvec,
     c       trestart, treorth

      common/timing/ nopx, nreorth, ndot, nrestart, nbsvd, nlandim,
     c     nsing, tmvopx, tgetu0, tlanbpro, tbsvd, tlansvd, tritzvec,
     c     trestart, treorth

      external second
