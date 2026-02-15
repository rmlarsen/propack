c
c     PROPACK scaling benchmark - complex8 (generated matrices)
c
c     Self-contained benchmark using programmatically generated matrices
c     for systematic scaling studies.  No external data files needed.
c
      program benchmark_gen

      implicit none

c     %------------%
c     | Parameters |
c     %------------%
      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)
      integer lswrk,lcwrk,liwrk,maxiter,nreps
      parameter(lswrk=gmmax+gnmax+14*gkmax+8*gkmax**2+32*gmmax+9)
      parameter(lcwrk=gmmax+gnmax+32*gmmax)
      parameter(liwrk = 8*gkmax)
      parameter(maxiter=1000)
      parameter(nreps=3)

c     %----------------------------%
c     | Generated matrix storage   |
c     %----------------------------%
      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

c     %-----------------%
c     | Local variables |
c     %-----------------%
      complex U(gmmax,gkmax+1),V(gnmax,gkmax)
      complex cwork(lcwrk)
      real Sigma(gkmax),swork(lswrk)
      real bnd(gkmax),B(gkmax,2)
      integer iwork(liwrk)
      common/cbenchgen/U,V,B,Sigma,bnd,swork,cwork,iwork

c-------------------- Here begins executable code ---------------------

      write (*,*) '=== PROPACK Complex8 Scaling Benchmark ==='

c
c     Random sparse configurations (nnz_per_col=20, seed=12345)
c
      call bench_config(1000,1000, 50,20,12345)
      call bench_config(2000,2000, 50,20,12345)
      call bench_config(3000,3000, 50,20,12345)
      call bench_config(3000,3000,100,20,12345)
      call bench_config(5000,5000, 50,20,12345)
      call bench_config(5000,5000,100,20,12345)
      call bench_config(5000,2000, 50,20,12345)

c
c     Diagonal configurations (alpha=1.0)
c
      call bench_diag_config(2000,2000, 50,1.0)
      call bench_diag_config(2000,2000,100,1.0)
      call bench_diag_config(5000,5000, 50,1.0)
      call bench_diag_config(5000,5000,100,1.0)

      write (*,*)
      write (*,*) '=== Benchmark complete ==='
      end


c-----------------------------------------------------------------------
c     bench_config: Generate random sparse matrix and run benchmarks
c-----------------------------------------------------------------------
      subroutine bench_config(m,n,neig,nnz_per_col,seed)
      implicit none
      integer m,n,neig,nnz_per_col,seed

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      gmattype = 1
      call gen_random_sparse(m,n,nnz_per_col,seed)

      write (*,*)
      write (*,'(a,i5,a,i5,a,i8,a)')
     c     ' Generated random sparse matrix (',m,'x',n,
     c     ', nnz=',gnnz,')'
      write (*,'(a,i5,a,i5)') ' Dimensions: m=',m,'  n=',n

      call bench_drivers(m,n,neig)
      end


c-----------------------------------------------------------------------
c     bench_diag_config: Generate diagonal matrix and run benchmarks
c-----------------------------------------------------------------------
      subroutine bench_diag_config(m,n,neig,alpha)
      implicit none
      integer m,n,neig
      real alpha

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      gmattype = 2
      call gen_diagonal(m,n,alpha)

      write (*,*)
      write (*,'(a,i5,a,i5,a,f4.1,a)')
     c     ' Generated diagonal matrix (',m,'x',n,
     c     ', alpha=',alpha,')'
      write (*,'(a,i5,a,i5)') ' Dimensions: m=',m,'  n=',n

      call bench_drivers(m,n,neig)
      end


c-----------------------------------------------------------------------
c     bench_drivers: Run all 3 PROPACK drivers for given configuration
c-----------------------------------------------------------------------
      subroutine bench_drivers(m,n,neig)
      implicit none
      integer m,n,neig

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)
      integer lswrk,lcwrk,liwrk,maxiter,nreps
      parameter(lswrk=gmmax+gnmax+14*gkmax+8*gkmax**2+32*gmmax+9)
      parameter(lcwrk=gmmax+gnmax+32*gmmax)
      parameter(liwrk = 8*gkmax)
      parameter(maxiter=1000)
      parameter(nreps=3)

      integer j,k,k0,dim,p,ioption(10)
      integer iparm(1),info,ierr,counters(7)
      integer iwork(liwrk)
      complex U(gmmax,gkmax+1),V(gnmax,gkmax)
      complex cwork(lcwrk),cparm(1)
      real soption(10),Sigma(gkmax),swork(lswrk)
      real tol,bnd(gkmax),B(gkmax,2),rnorm,eps
      real times(nreps),timers(8),t1,t2
      common/cbenchgen/U,V,B,Sigma,bnd,swork,cwork,iwork

      real slamch
      external slamch, matvecmult

      tol = 1e-6

      write (*,'(a,i4)') ' *** neig = ', neig
      write (*,'(a)')
     c     ' -----------------------------------------------'

c     Set options
      eps = slamch('e')
      ioption(1) = 0
      ioption(2) = 1
      soption(1) = sqrt(eps)
      soption(2) = eps**(3.0/4.0)
      soption(3) = 0.0
      soption(4) = 0.002

c     --- Benchmark clanbpro ---
      do j=1,nreps
         call czero(m,U,1)
         k = min(neig+10, gkmax, min(m,n))
         k0 = 0
         call clearstat
         call cpu_time(t1)
         call clanbpro(m,n,k0,k,matvecmult,U,gmmax,V,gnmax,
     c        B,gkmax,rnorm,soption,ioption,swork,cwork,iwork,
     c        cparm,iparm,ierr)
         call cpu_time(t2)
         times(j) = t2 - t1
      enddo
      call getstat(counters, timers)
      call report_times('clanbpro', neig, nreps, times,
     c     counters, timers)

c     --- Benchmark clansvd ---
      do j=1,nreps
         call czero(m,U,1)
         k = neig
         call clearstat
         call cpu_time(t1)
         call clansvd('n','n',m,n,k,gkmax,matvecmult,
     c        U,gmmax,Sigma,bnd,V,gnmax,tol,swork,lswrk,
     c        cwork,lcwrk,iwork,liwrk,soption,ioption,
     c        info,cparm,iparm)
         call cpu_time(t2)
         times(j) = t2 - t1
      enddo
      call getstat(counters, timers)
      call report_times('clansvd', neig, nreps, times,
     c     counters, timers)

c     --- Benchmark clansvd_irl ---
      dim = min(2*neig+10, gkmax, min(m,n))
      p = min(dim-neig, neig)
      if (p.gt.0) then
         do j=1,nreps
            call czero(m,U,1)
            dim = min(2*neig+10, gkmax, min(m,n))
            p = min(dim-neig, neig)
            k = neig
            call clearstat
            call cpu_time(t1)
            call clansvd_irl('l','n','n',m,n,dim,p,k,
     c           maxiter,matvecmult,U,gmmax,Sigma,bnd,V,gnmax,
     c           tol,swork,lswrk,cwork,lcwrk,iwork,liwrk,
     c           soption,ioption,info,cparm,iparm)
            call cpu_time(t2)
            times(j) = t2 - t1
         enddo
         call getstat(counters, timers)
         call report_times('clansvd_irl', neig, nreps, times,
     c        counters, timers)
      else
         write (*,*) 'Skipping clansvd_irl: p<=0 for neig=',neig
      endif

      end


c-----------------------------------------------------------------------
c     matvecmult: Dispatch to sparse or diagonal matvec
c-----------------------------------------------------------------------
      subroutine matvecmult(tression,m,n,x,y,dparm,iparm)
      implicit none
      character*(*) tression
      integer m,n
      complex x(*),y(*),dparm(*)
      integer iparm(*)

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      if (gmattype.eq.1) then
         call sparse_mv(tression,m,n,x,y)
      else
         call diag_mv(tression,m,n,x,y)
      endif
      end


c-----------------------------------------------------------------------
c     sparse_mv: CCS sparse matrix-vector multiply (complex8)
c-----------------------------------------------------------------------
      subroutine sparse_mv(trans,m,n,x,y)
      implicit none
      character*(*) trans
      integer m,n

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      complex x(*),y(*)
      integer i,j,k

      if (trans.eq.'n' .or. trans.eq.'N') then
c        y = A*x
         do i=1,m
            y(i) = cmplx(0.0,0.0)
         enddo
         do j=1,n
            do k=gcolptr(j),gcolptr(j+1)-1
               i = growind(k)
               y(i) = y(i) + gvalues(k)*x(j)
            enddo
         enddo
      else
c        y = A^H*x  (conjugate transpose)
         do j=1,n
            y(j) = cmplx(0.0,0.0)
         enddo
         do j=1,n
            do k=gcolptr(j),gcolptr(j+1)-1
               i = growind(k)
               y(j) = y(j) + conjg(gvalues(k))*x(i)
            enddo
         enddo
      endif
      end


c-----------------------------------------------------------------------
c     diag_mv: Diagonal matrix-vector multiply (real diagonal, complex vectors)
c-----------------------------------------------------------------------
      subroutine diag_mv(trans,m,n,x,y)
      implicit none
      character*(*) trans
      integer m,n

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      complex x(*),y(*)
      integer i,mn

      mn = min(m,n)
      if (trans.eq.'n' .or. trans.eq.'N') then
         do i=1,mn
            y(i) = cmplx(gdiag(i),0.0)*x(i)
         enddo
         do i=mn+1,m
            y(i) = cmplx(0.0,0.0)
         enddo
      else
c        A^H = A^T for real diagonal
         do i=1,mn
            y(i) = cmplx(gdiag(i),0.0)*x(i)
         enddo
         do i=mn+1,n
            y(i) = cmplx(0.0,0.0)
         enddo
      endif
      end


c-----------------------------------------------------------------------
c     gen_random_sparse: Fill CCS arrays with random complex sparse matrix
c-----------------------------------------------------------------------
      subroutine gen_random_sparse(m,n,nnz_per_col,seed)
      implicit none
      integer m,n,nnz_per_col,seed

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      integer j,k,idx,iseed
      integer lcg_next
      real rval,ival

      gm = m
      gn = n
      iseed = seed
      idx = 1

      do j=1,n
         gcolptr(j) = idx
         do k=1,nnz_per_col
            iseed = lcg_next(iseed)
            growind(idx) = mod(abs(iseed), m) + 1
            iseed = lcg_next(iseed)
            rval = 2.0 * real(iseed) / 2147483647.0 - 1.0
            iseed = lcg_next(iseed)
            ival = 2.0 * real(iseed) / 2147483647.0 - 1.0
            gvalues(idx) = cmplx(rval, ival)
            idx = idx + 1
         enddo
      enddo
      gcolptr(n+1) = idx
      gnnz = idx - 1
      end


c-----------------------------------------------------------------------
c     gen_diagonal: Fill diagonal array with prescribed decay
c-----------------------------------------------------------------------
      subroutine gen_diagonal(m,n,alpha)
      implicit none
      integer m,n
      real alpha

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      complex gvalues(gnnzmax)
      real gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      integer i,mn
      real rmn

      gm = m
      gn = n
      gnnz = 0
      mn = min(m,n)
      rmn = real(mn)

      do i=1,mn
         gdiag(i) = (rmn / real(i))**alpha
      enddo
      end


c-----------------------------------------------------------------------
c     lcg_next: Park-Miller LCG for portable reproducibility
c-----------------------------------------------------------------------
      function lcg_next(seed)
      implicit none
      integer seed, lcg_next
      integer*8 s8
      s8 = seed
      s8 = mod(s8 * 16807, 2147483647_8)
      seed = int(s8)
      lcg_next = seed
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
