c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c
      subroutine dsafescal(n,alpha,x)
c
c     Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n
      double precision alpha, x(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one
      parameter(one = 1.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,info

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dscal,dlascl

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch
      external dlamch

      if (abs(alpha).ge.dlamch('s')) then
         call pdscal(n,one/alpha, x, 1)
      else
         call dlascl('General',i,i,alpha,one,n,1,x,n,info)
      endif

      end
