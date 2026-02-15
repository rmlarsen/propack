c
c     PROPACK scaling benchmark - double precision (generated matrices)
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
      integer lwrk,liwrk,maxiter,nreps
      parameter(lwrk=gmmax+gnmax+14*gkmax+8*gkmax**2+32*gmmax+9)
      parameter(liwrk = 8*gkmax)
      parameter(maxiter=1000)
      parameter(nreps=3)

c     %----------------------------%
c     | Generated matrix storage   |
c     %----------------------------%
c     CCS sparse: colptr, rowind, values
c     Diagonal: diag array
c     gmattype: 1=sparse, 2=diagonal
      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision U(gmmax,gkmax+1),V(gnmax,gkmax),work(lwrk)
      double precision Sigma(gkmax)
      double precision bnd(gkmax),B(gkmax,2)
      integer iwork(liwrk)
c     Save large arrays in common block to avoid stack overflow
      common/cbenchgen/U,V,B,Sigma,bnd,work,iwork

c-------------------- Here begins executable code ---------------------

      write (*,*) '=== PROPACK Double Precision Scaling Benchmark ==='

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
      call bench_diag_config(2000,2000, 50,1.0d0)
      call bench_diag_config(2000,2000,100,1.0d0)
      call bench_diag_config(5000,5000, 50,1.0d0)
      call bench_diag_config(5000,5000,100,1.0d0)

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
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
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
      double precision alpha

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
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
      integer lwrk,liwrk,maxiter,nreps
      parameter(lwrk=gmmax+gnmax+14*gkmax+8*gkmax**2+32*gmmax+9)
      parameter(liwrk = 8*gkmax)
      parameter(maxiter=1000)
      parameter(nreps=3)

      integer j,ii,k,k0,dim,p,ioption(10)
      integer iparm(1),info,ierr,counters(7)
      integer iwork(liwrk)
      double precision U(gmmax,gkmax+1),V(gnmax,gkmax),work(lwrk)
      double precision doption(10),Sigma(gkmax),dparm(1)
      double precision tol,bnd(gkmax),B(gkmax,2),rnorm,eps
      real times(nreps),timers(8),t1,t2
      common/cbenchgen/U,V,B,Sigma,bnd,work,iwork

      double precision dlamch
      external dlamch, matvecmult

      tol = 1d-12

      write (*,'(a,i4)') ' *** neig = ', neig
      write (*,'(a)')
     c     ' -----------------------------------------------'

c     Set options
      eps = dlamch('e')
      ioption(1) = 0
      ioption(2) = 1
      doption(1) = sqrt(eps)
      doption(2) = eps**(3.0/4.0)
      doption(3) = 0d0
      doption(4) = 0.002

c     --- Benchmark dlanbpro ---
      do j=1,nreps
         do ii=1,m
            U(ii,1) = 0d0
         enddo
         k = min(neig+10, gkmax, min(m,n))
         k0 = 0
         call clearstat
         call cpu_time(t1)
         call dlanbpro(m,n,k0,k,matvecmult,U,gmmax,V,gnmax,
     c        B,gkmax,rnorm,doption,ioption,work,iwork,
     c        dparm,iparm,ierr)
         call cpu_time(t2)
         times(j) = t2 - t1
      enddo
      call getstat(counters, timers)
      call report_times('dlanbpro', neig, nreps, times,
     c     counters, timers)

c     --- Benchmark dlansvd ---
      do j=1,nreps
         do ii=1,m
            U(ii,1) = 0d0
         enddo
         k = neig
         call clearstat
         call cpu_time(t1)
         call dlansvd('n','n',m,n,k,gkmax,matvecmult,
     c        U,gmmax,Sigma,bnd,V,gnmax,tol,work,lwrk,
     c        iwork,liwrk,doption,ioption,info,dparm,iparm)
         call cpu_time(t2)
         times(j) = t2 - t1
      enddo
      call getstat(counters, timers)
      call report_times('dlansvd', neig, nreps, times,
     c     counters, timers)

c     --- Benchmark dlansvd_irl ---
      dim = min(2*neig+10, gkmax, min(m,n))
      p = min(dim-neig, neig)
      if (p.gt.0) then
         do j=1,nreps
            do ii=1,m
               U(ii,1) = 0d0
            enddo
            dim = min(2*neig+10, gkmax, min(m,n))
            p = min(dim-neig, neig)
            k = neig
            call clearstat
            call cpu_time(t1)
            call dlansvd_irl('l','n','n',m,n,dim,p,k,
     c           maxiter,matvecmult,U,gmmax,Sigma,bnd,V,gnmax,
     c           tol,work,lwrk,iwork,liwrk,doption,ioption,
     c           info,dparm,iparm)
            call cpu_time(t2)
            times(j) = t2 - t1
         enddo
         call getstat(counters, timers)
         call report_times('dlansvd_irl', neig, nreps, times,
     c        counters, timers)
      else
         write (*,*) 'Skipping dlansvd_irl: p<=0 for neig=',neig
      endif

      end


c-----------------------------------------------------------------------
c     matvecmult: Dispatch to sparse or diagonal matvec
c-----------------------------------------------------------------------
      subroutine matvecmult(tression,m,n,x,y,dparm,iparm)
      implicit none
      character*(*) tression
      integer m,n
      double precision x(*),y(*),dparm(*)
      integer iparm(*)

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      if (gmattype.eq.1) then
         call sparse_mv(tression,m,n,x,y)
      else
         call diag_mv(tression,m,n,x,y)
      endif
      end


c-----------------------------------------------------------------------
c     sparse_mv: CCS sparse matrix-vector multiply
c-----------------------------------------------------------------------
      subroutine sparse_mv(trans,m,n,x,y)
      implicit none
      character*(*) trans
      integer m,n

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      double precision x(*),y(*)
      integer i,j,k

      if (trans.eq.'n' .or. trans.eq.'N') then
c        y = A*x
         do i=1,m
            y(i) = 0d0
         enddo
         do j=1,n
            do k=gcolptr(j),gcolptr(j+1)-1
               i = growind(k)
               y(i) = y(i) + gvalues(k)*x(j)
            enddo
         enddo
      else
c        y = A^T*x
         do j=1,n
            y(j) = 0d0
         enddo
         do j=1,n
            do k=gcolptr(j),gcolptr(j+1)-1
               i = growind(k)
               y(j) = y(j) + gvalues(k)*x(i)
            enddo
         enddo
      endif
      end


c-----------------------------------------------------------------------
c     diag_mv: Diagonal matrix-vector multiply
c-----------------------------------------------------------------------
      subroutine diag_mv(trans,m,n,x,y)
      implicit none
      character*(*) trans
      integer m,n

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      double precision x(*),y(*)
      integer i,mn

      mn = min(m,n)
      if (trans.eq.'n' .or. trans.eq.'N') then
c        y = A*x  (m output)
         do i=1,mn
            y(i) = gdiag(i)*x(i)
         enddo
         do i=mn+1,m
            y(i) = 0d0
         enddo
      else
c        y = A^T*x  (n output)
         do i=1,mn
            y(i) = gdiag(i)*x(i)
         enddo
         do i=mn+1,n
            y(i) = 0d0
         enddo
      endif
      end


c-----------------------------------------------------------------------
c     gen_random_sparse: Fill CCS arrays with random sparse matrix
c-----------------------------------------------------------------------
      subroutine gen_random_sparse(m,n,nnz_per_col,seed)
      implicit none
      integer m,n,nnz_per_col,seed

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      integer j,k,idx,iseed
      integer lcg_next

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
            gvalues(idx) = 2.0d0 * dble(iseed) / 2147483647.0d0 - 1.0d0
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
      double precision alpha

      integer gmmax,gnmax,gkmax,gnnzmax
      parameter(gmmax=5000,gnmax=5000,gkmax=300,gnnzmax=200000)

      integer gcolptr(gnmax+1),growind(gnnzmax)
      double precision gvalues(gnnzmax)
      double precision gdiag(gnmax)
      integer gm,gn,gnnz,gmattype
      common/genmat/gvalues,gdiag,gcolptr,growind,gm,gn,gnnz,gmattype

      integer i,mn
      double precision dmn

      gm = m
      gn = n
      gnnz = 0
      mn = min(m,n)
      dmn = dble(mn)

      do i=1,mn
         gdiag(i) = (dmn / dble(i))**alpha
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
