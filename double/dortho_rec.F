      subroutine dset_mu(k,mu,index,val)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer k,index(*)
      double precision mu(*), val

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
         p = index(i)
         q = index(i+1)
         do j=p,q
            mu(j) = val
         enddo
         i = i+2
      enddo
      end
c
c**********************************************************************
c
      subroutine dcompute_int(mu,j,delta,eta,index)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer j,index(*)
      double precision mu(*)
      double precision delta,eta

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,k,s,ip

      if (delta.lt.eta) then
         return
      endif

      ip = 0
      index(1) = 0
      i = 0
      outer: do while(i.lt.j)
c        find the next mu(k), k>i where abs(mu(k)) > delta
         do k=i+1,j
            if (abs(mu(k)).gt.delta) exit
         enddo
         if (k.gt.j) exit outer
c        find smallest s<=k such that abs(mu(s..k)) >= eta
         do s=k,max(i,1),-1
            if (abs(mu(s)).lt.eta) exit
         enddo
         ip = ip+1
         index(ip) = s+1
         do i=s+1,j
            if (abs(mu(i)).lt.eta) exit
         enddo
         ip = ip+1
         index(ip) = i-1
      enddo outer
      ip = ip+1
      index(ip) = j+1
      end
c
c**********************************************************************
c
      subroutine dupdate_mu(mumax,mu,nu,j,alpha,beta,anorm,
     c     eps1_aprod,eps1_atprod)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer j
      double precision mumax,eps1_aprod,eps1_atprod,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one
      parameter(one = 1d0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d
      integer k

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlapy2
      external dlapy2

      if (j.eq.1) then
         d = eps1_aprod*dlapy2(alpha(j), beta(j)) +
     c        eps1_atprod*alpha(1) +
     c        max(eps1_aprod,eps1_atprod)*anorm
         mu(1) = max(eps1_aprod,eps1_atprod)/beta(1)
         mumax = abs(mu(1))
      else
         mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
         d = eps1_aprod*dlapy2(alpha(j), beta(j)) +
     c        eps1_atprod*alpha(1) +
     c        max(eps1_aprod,eps1_atprod)*anorm
         mu(1) = (mu(1) + dsign(d,mu(1))) / beta(j)
         mumax = abs(mu(1))
         do k=2,j-1
            mu(k) = alpha(k)*nu(k) +beta(k-1)*nu(k-1)-alpha(j)*mu(k)
            d = eps1_aprod*dlapy2(alpha(j), beta(j)) +
     c           eps1_atprod*dlapy2(alpha(k), beta(k-1)) +
     c           max(eps1_aprod,eps1_atprod)*anorm
            mu(k) = (mu(k) + dsign(d,mu(k))) / beta(j)
            mumax = max(mumax,abs(mu(k)))
         enddo
         mu(j) = beta(j-1)*nu(j-1)
         d = eps1_aprod*dlapy2(alpha(j), beta(j)) +
     c        eps1_atprod*dlapy2(alpha(j), beta(j-1)) +
     c        max(eps1_aprod,eps1_atprod)*anorm
         mu(j) = (mu(j) + dsign(d,mu(j))) / beta(j)
         mumax = max(mumax,abs(mu(j)))
      endif
      mu(j+1) = one
      end
c
c**********************************************************************
c

      subroutine dupdate_nu(numax,mu,nu,j,alpha,beta,anorm,
     c     eps1_aprod,eps1_atprod)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer j
      double precision numax,eps1_aprod,eps1_atprod,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1d0, zero = 0d0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d
      integer k

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlapy2
      external dlapy2

      if (j.gt.1) then
         numax = zero
         do k=1,j-1
            nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
            d = eps1_atprod*dlapy2(alpha(k),beta(k)) +
     c           eps1_aprod*dlapy2(alpha(j),beta(j-1)) +
     c           max(eps1_aprod,eps1_atprod)*anorm
            nu(k) = (nu(k) + dsign(d,nu(k))) / alpha(j)
            numax = max(numax,abs(nu(k)))
         enddo
         nu(j) = one
      endif
      end
