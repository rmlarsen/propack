c
c     PROPACK benchmark - complex single precision (complex)
c
c     Benchmarks clanbpro, clansvd, and clansvd_irl across varying
c     numbers of requested singular values.  All parameters are
c     hardcoded for reproducibility (no stdin input required).
c
      program benchmark

      implicit none

      include 'matvec.h'

c     %------------%
c     | Parameters |
c     %------------%
      integer lswrk,lcwrk,liwrk,maxiter,nreps,nneig
      parameter(lswrk=mmax+nmax+14*kmax+8*kmax**2+32*mmax+9)
      parameter(lcwrk=mmax+nmax+32*mmax)
      parameter(liwrk = 8*kmax)
      parameter(maxiter=1000)
      parameter(nreps=5,nneig=4)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,m,n,neig,k,k0,dim,p,ioption(10)
      integer iparm(1),info,ierr,counters(7)
      integer neig_list(nneig),iwork(liwrk)
      complex U(mmax,kmax+1),V(nmax,kmax)
      complex cwork(lcwrk),cparm(1)
      real soption(10),Sigma(kmax),swork(lswrk)
      real tol,bnd(kmax),B(kmax,2),rnorm,eps
      real times(nreps),timers(8),t1,t2
c     Save large arrays in common block to avoid stack overflow
      common/cbenchmark/U,V,B,Sigma,bnd,swork,cwork,iwork

c     %--------------------%
c     | External Functions |
c     %--------------------%
      real slamch
      external slamch, matvecmult

      data neig_list/10, 50, 100, 200/

c-------------------- Here begins executable code ---------------------

      matform = 1
      tol = 1e-6

      write (*,*) '=== PROPACK Complex8 Benchmark ==='
      write (*,*) 'Reading matrix from file mhd1280b.cua'
      call readmatrix('mhd1280b.cua',m,n)
      write (*,*) 'Dimensions: m=',m,' n=',n
      if (m.gt.mmax .or. n.gt.nmax) then
         stop 'Matrix exceeded maximum size defined in matvec.h.'
      endif

c
c     Set options
c
      eps = slamch('e')
      ioption(1) = 0
      ioption(2) = 1
      soption(1) = sqrt(eps)
      soption(2) = eps**(3.0/4.0)
      soption(3) = 0.0
      soption(4) = 0.002

c
c     Loop over neig values
c
      do i=1,nneig
         neig = neig_list(i)
         write (*,*)
         write (*,'(a,i4)') ' *** neig = ', neig
         write (*,'(a)')
     c        ' -----------------------------------------------'

c        --- Benchmark clanbpro ---
         do j=1,nreps
            call czero(m,U,1)
            k = min(neig+10, kmax, min(m,n))
            k0 = 0
            call clearstat
            call cpu_time(t1)
            call clanbpro(m,n,k0,k,matvecmult,U,mmax,V,nmax,
     c           B,kmax,rnorm,soption,ioption,swork,cwork,iwork,
     c           cparm,iparm,ierr)
            call cpu_time(t2)
            times(j) = t2 - t1
         enddo
         call getstat(counters, timers)
         call report_times('clanbpro', neig, nreps, times,
     c        counters, timers)

c        --- Benchmark clansvd ---
         do j=1,nreps
            call czero(m,U,1)
            k = neig
            call clearstat
            call cpu_time(t1)
            call clansvd('n','n',m,n,k,kmax,matvecmult,
     c           U,mmax,Sigma,bnd,V,nmax,tol,swork,lswrk,
     c           cwork,lcwrk,iwork,liwrk,soption,ioption,
     c           info,cparm,iparm)
            call cpu_time(t2)
            times(j) = t2 - t1
         enddo
         call getstat(counters, timers)
         call report_times('clansvd', neig, nreps, times,
     c        counters, timers)

c        --- Benchmark clansvd_irl ---
         dim = min(2*neig+10, kmax, min(m,n))
         p = min(dim-neig, neig)
         if (p.gt.0) then
            do j=1,nreps
               call czero(m,U,1)
               dim = min(2*neig+10, kmax, min(m,n))
               p = min(dim-neig, neig)
               k = neig
               call clearstat
               call cpu_time(t1)
               call clansvd_irl('l','n','n',m,n,dim,p,k,
     c              maxiter,matvecmult,U,mmax,Sigma,bnd,V,nmax,
     c              tol,swork,lswrk,cwork,lcwrk,iwork,liwrk,
     c              soption,ioption,info,cparm,iparm)
               call cpu_time(t2)
               times(j) = t2 - t1
            enddo
            call getstat(counters, timers)
            call report_times('clansvd_irl', neig, nreps, times,
     c           counters, timers)
         else
            write (*,*) 'Skipping clansvd_irl: p<=0 for neig=',neig
         endif

      enddo

      write (*,*)
      write (*,*) '=== Benchmark complete ==='
      end


c-----------------------------------------------------------------------
c     report_times: Print timing statistics and PROPACK internal stats
c-----------------------------------------------------------------------
      subroutine report_times(drvname, neig, nreps, times,
     c     counters, timers)
      implicit none
      character*(*) drvname
      integer neig, nreps
      real times(nreps), timers(8)
      integer counters(7)

      integer i, j
      real tmin, tmax, tmean, tmedian, tsum, tmp

c     Insertion sort
      do i=2,nreps
         tmp = times(i)
         j = i - 1
         do while (j.ge.1 .and. times(j).gt.tmp)
            times(j+1) = times(j)
            j = j - 1
         enddo
         times(j+1) = tmp
      enddo

c     Compute statistics
      tmin = times(1)
      tmax = times(nreps)
      tsum = 0.0
      do i=1,nreps
         tsum = tsum + times(i)
      enddo
      tmean = tsum / real(nreps)
      if (mod(nreps,2).eq.1) then
         tmedian = times((nreps+1)/2)
      else
         tmedian = (times(nreps/2) + times(nreps/2+1)) / 2.0
      endif

c     Print results
      write (*,'(a,a,a,i4,a)') '  Driver: ', drvname,
     c     '  neig=', neig, ':'
      write (*,'(a,1p,e12.4,a,e12.4,a,e12.4,a,e12.4)')
     c     '    Time(s) min=', tmin, '  median=', tmedian,
     c     '  mean=', tmean, '  max=', tmax
      write (*,'(a,i6,a,i6,a,i8,a,i6)')
     c     '    nopx=', counters(1), '  nreorth=', counters(2),
     c     '  ndot=', counters(3), '  nbsvd=', counters(4)
      write (*,'(a,i4,a,i6,a,i4)')
     c     '    nrestart=', counters(5),
     c     '  nlandim=', counters(6),
     c     '  nsing=', counters(7)
      write (*,'(a,1p,e10.3,a,e10.3,a,e10.3,a,e10.3)')
     c     '    tmvopx=', timers(1), ' tgetu0=', timers(2),
     c     ' tlanbpro=', timers(3), ' tbsvd=', timers(4)
      write (*,'(a,1p,e10.3,a,e10.3,a,e10.3,a,e10.3)')
     c     '    tlansvd=', timers(5), ' tritzvec=', timers(6),
     c     ' trestart=', timers(7), ' treorth=', timers(8)
      end
