      subroutine clansvd_irl(which,jobu,jobv,m,n,dim,p,nwanted,maxiter,
     c     aprod,U,ldu,Sigma,bnd,V,ldv,tolin,swork,lswork,cwork,lcwrk,
     c     iwork,liwork,soption,ioption,info,cparm,iparm)

c     CLANSVD_IRL: Compute the leading or trailing singular triplets of a
c     sparse matrix A by implicitly restarted Lanczos bidiagonalization
c     with partial reorthogonalization.
c
c     Parameters:
c
c     WHICH: CHARACTER*1. Decides which singular triplets to compute.
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values.
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values.
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     DIM:  INTEGER. Dimension of the Krylov subspace. DIM <= MIN(M,N).
c     P:    INTEGER. Number of shift per restart.
c     NWANTED: INTEGER. Number of desired singular triplets.
c              NWANTED should be at most DIM-P-1.
c     MAXITER: INTEGER. Maximum number of restarts.
c     APROD: Subroutine defining the linear operator A.
c            APROD should be of the form:
c
c           SUBROUTINE APROD(TRANSA,M,N,X,Y,CPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           COMPLEX X(*),Y(*),CPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'C' then the function should compute the matrix-vector
c           product Y = A^H * X, where A^H is the conjugate transpose (adjoint)
c           of A.
c           The arrays IPARM and CPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors
c               corresponding to the K largest or smallest (depending on the
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     SIGMA(K): REAL array. On return Sigma contains approximation
c               to the K largest or smallest (depending on the
c               value of WHICH) singular values of A.
c     BND(K)  : REAL array. Error estimates on the computed
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,KMAX): COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors
c               corresponding to the K largest or smallest (depending on the
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: REAL. Desired accuracy of computed singular values.
c            SIGMA(K) is considered converged when
c            BND(K) <= MAX(TOLIN*SIGMA(K), 16*EPS*||A||)
c     SWORK(LSWORK): REAL array. Workspace of dimension LSWORK.
c     LSWORK: INTEGER. Dimension of SWORK.
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LSWORK should be at least
c            M + N + 10*KMAX + 2*KMAX**2 + 5 + MAX(M+N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LSWORK should be at least
c            M + N + 10*KMAX + 5*KMAX**2 + 4 +
c            MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>1 is a block
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3
c            operation.
c     CWORK: COMPLEX array of dimension LCWORK.
c     LCWORK: INTEGER. Dimension of CWORK.
c             If JOBU.EQ.'N' and JOBV.EQ.'N' then LCWORK should be at least
c             M + N.
c             If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LCWORK should be at least
c             M + N + NB*MAX(M,N), where NB>0 is a block size, which determines
c             how large a fraction of the work in setting up the singular
c             vectors is done using fast BLAS-3 operations.
c     IWORK: INTEGER array. Integer workspace of dimension LIWORK.
c     LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
c             JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
c     DOPTION: REAL array. Parameters for LANBPRO.
c        soption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        soption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        soption(3) = anorm. Estimate of || A ||.
c        soption(4) = min relgap. Smallest relgap allowed between any shift
c                     the smallest requested Ritz value.
c
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical Gram-Schmidt. IF CGS.EQ.0 then
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c
c     INFO: INTEGER.
c         INFO = 0: The NWANTED largest or smallest (depending on the value of
c                   WHICH) singular triplets were computed succesfully.
c         INFO = J>0, J<NWANTED: An invariant subspace of dimension J was found.
c         INFO = J<0, |J|<NWANTED: Only J singular triplets converge before
c                    exceeding MAXITER iterations.
c     CPARM: COMPLEX array. Array used for passing data to the APROD
c         function.
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include '../common/stat.h'
      character*1 which,jobu,jobv
      integer m,n,p,nwanted,maxiter,ldu,ldv,iter,liwork,lcwrk
      integer iwork(liwork),lswork,info,ioption(*)
      complex U(ldu,*),V(ldv,*), cwork(lcwrk),cparm(*)
      real Sigma(*),bnd(*),swork(lswork)
      real tolin,soption(*)
      integer iparm(*)
      external aprod

c     %------------%
c     | Parameters |
c     %------------%
      real one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer k,i,ibnd,iwrk,ierr,ip,iq,nconv,lwrk,kold,dim,first,last
      integer ialpha,ibeta,ialpha1,ibeta1,ishift,nshft,lapinfo,idx
      real eps,eps34,epsn2,sfmin,anorm,rnorm,tol
      real shift, relgap, deflate_tol
      real t0,t1,t2,t3
      integer st,cnk,wst,wcnk, tid, nt

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external szero,izero,scopy,sbdsqr

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      real slamch,pscnrm2
      external pscnrm2,lsame
      external slamch
#ifdef _OPENMP
      integer omp_get_thread_num, omp_get_num_threads
      external omp_get_thread_num, omp_get_num_threads
#endif

c-------------------- Here begins executable code ---------------------

c     %-------------%
c     | Start timer |
c     %-------------%
      call cpu_time(t0)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = slamch('e')
      eps34 = eps**(3.0/4.0)
      epsn2 = sqrt(real(max(m,n)))*eps/2.0
      sfmin = slamch('s')


c     %--------------------------------%
c     | Guard against absurd arguments |
c     %--------------------------------%
      dim = min(dim,n+1,m+1)
      k = dim-p
      tol = min(one,max(16.0*eps,tolin))
      anorm = zero

c     %------------------------------%
c     | Set pointers into swork array |
c     %------------------------------%
      ibnd = 1
      ialpha = ibnd + dim+1
      ibeta = ialpha + dim
      ialpha1 = ibeta + dim
      ibeta1 = ialpha1 + dim
      ishift = ibeta1 + dim
      ip = ishift + dim
      iq = ip + (dim+1)**2
      iwrk = iq + dim**2
      lwrk = lswork-iwrk+1
      call szero(8*dim + 3 + 2*dim**2,swork,1)

c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
      rnorm = pscnrm2(m,U(1,1),1)
      iter = 0
      info = 0
      nconv = 0
      kold = 0
c     %------------------------------%
c     | Iterate until convergence... |
c     %------------------------------%
      do while (nconv.lt.nwanted .and. iter.lt.maxiter)

c     %---------------------------------------------------%
c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
c     %---------------------------------------------------%

         call clanbpro(m, n, kold, dim, aprod, U, ldu, V, ldv,
     c        swork(ialpha),dim,rnorm,soption(1),ioption(1),
     c        swork(iwrk), cwork, iwork, cparm, iparm, ierr)
         kold = k
c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%
         call cpu_time(t2)
         call scopy(dim, swork(ialpha),1,swork(ialpha1),1)
         call scopy(dim, swork(ibeta),1,swork(ibeta1),1)
         call szero(dim+1,swork(ibnd),1)
         call sbdqr((dim.eq.min(m,n)),'N',dim,swork(ialpha1),
     c        swork(ibeta1),swork(ibnd+dim-1),swork(ibnd+dim),
     c        swork(ip),dim+1)
         call sbdsqr('u',dim,0,1,0,swork(ialpha1),swork(ibeta1),swork,1,
     c        swork(ibnd),1,swork,1,swork(iwrk),lapinfo)
         call  cpu_time(t3)
         tbsvd = tbsvd + (t3-t2)
         nbsvd = nbsvd + 1

         if (dim.gt.5 .and. lsame(which,'l')) then
            anorm = swork(ialpha1)
         else
            anorm = max(anorm,swork(ialpha1))
         endif
         do i=1,dim
            swork(ibnd+i-1) = abs(rnorm*swork(ibnd+i-1))
         enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%
         if (lsame(which,'s')) then
            call srefinebounds(dim,dim,swork(ialpha1),swork(ibnd),
     c           epsn2*anorm,eps34)
         else
            call srefinebounds(min(m,n),dim,swork(ialpha1),swork(ibnd),
     c           epsn2*anorm,eps34)
         endif

c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
         if (lsame(which,'s')) then
            first = max(1,dim-nwanted+1)
            last = dim
         else
            first = 1
            last = min(nwanted,dim)
         endif
         nconv = 0
         do i=first,last
            idx = i-first+1
            sigma(idx) = swork(ialpha1+i-1)
            bnd(idx) = swork(ibnd+i-1)
            if (bnd(idx).le.max(tol*sigma(idx), 16*eps*anorm)) then
               nconv = nconv + 1
            endif
         enddo

c     %-----------------------------------------------%
c     | Test if an invariant subspace have been found |
c     %-----------------------------------------------%
         if (ierr.lt.0) then
            if (dim.lt.k) then
               info = dim
            endif
            goto 50
         endif
         if (nconv.lt.nwanted) then

c     %---------------------------------------------------------------------%
c     | Implicit restart:                                                   |
c     |   Apply  shifts mu_1, mu_2,...,mu_p to "back up" the          |
c     |   bidiagonalization (dim-k) steps to                                |
c     |                                                                     |
c     |      A V_{k}^{+} = U_{k+1}^{+} B_{k}^{+}                            |
c     | corresponding to starting vector                                    |
c     |      u_1^{+} = \prod_{i=1}^{p} (A A^T - mu_i^2) u_1             |
c     |                                                                     |
c     | We use exact shifts mu_i for which the relative gap between mu_i    |
c     | and the lower bound on the k'th Ritzvalue is larger than soption(4) |
c     %---------------------------------------------------------------------%
            call cpu_time(t2)
            call szero(dim-k,swork(ishift),1)
            nshft = 0
            if (lsame(which,'s')) then
               do i=1,dim-k
                  relgap = swork(ialpha1+i-1) -
     c                 swork(ialpha1+dim-nwanted-1)
                  if (relgap.gt.soption(4)*swork(ialpha1+dim-nwanted-1))
     c                 then
                     swork(ishift + nshft) = swork(ialpha1+i-1)
                  else
                     swork(ishift + nshft) = swork(ialpha1)
                  endif
                  nshft = nshft + 1
               enddo
            else
               do i=k+1,dim
                  relgap = swork(ialpha1+k-1) - swork(ialpha1+i-1)
                  if (relgap.gt.soption(4)*swork(ialpha1+k-1)) then
                     swork(ishift + nshft) = swork(ialpha1+i-1)
                  else
                     swork(ishift + nshft) = zero
                  endif
                  nshft = nshft + 1
               enddo
            endif

c     %----------------------------------------------------%
c     | Apply dim-k steps of the implicitly shifted QL     |
c     | algorithm and accumulate rotations such that       |
c     |   B_{dim}^{+} = P * B_{dim} * Q^T                  |
c     | Each unreduced sub-matrix is treated independently |
c     %----------------------------------------------------%
            call szero((dim+1)*(dim+1),swork(ip),1)
            call szero(dim*dim,swork(iq),1)
            do i=1,dim+1
               swork(ip+(i-1)*(dim+2)) = one
            enddo
            do i=1,dim
               swork(iq+(i-1)*(dim+1)) = one
            enddo

c     Note: The implicit LQ steps applied in sbsvdstep called below are not
c     stable if the off-diagonal elements of B*B^T (which equal
c     alpha(i)*beta(i)) are near zero. Here we split the lower bidiagonal
c     matrix into independent blocks based on whether
c     alpha(i)*beta(i) < eps*||B||^2. See Parlett SEP section 7.11.            
            deflate_tol = eps*anorm*anorm
            do i=dim,k+1,-1
               shift = swork(ishift+dim-i)
               first = 1
               do while (first.lt.dim)
                  last = first
c     Find the end of the next unreduced sub-matrix.
                  do while(last.lt.dim .and.
     c                 abs(swork(ibeta+last-1)*swork(ialpha+last-1))
     c                 .gt. deflate_tol)
                     last = last + 1
                  enddo
                  if (last.gt.first) then
                     call sbsvdstep('y','y',dim+1,dim,
     c                    min(dim,last)-first+1,
     c                    shift,
     c                    swork(ialpha+first-1),
     c                    swork(ibeta+first-1),
     c                    swork(ip+(first-1)*(dim+1)),dim+1,
     c                    swork(iq+(first-1)*dim), dim)
                  endif
                  first = last+1
               enddo
            enddo

c     %---------------------------------------------------%
c     | Compute first k+1 left and first k right updated  |
c     | Lanczos vectors                                   |
c     |   U_{dim+1}^{+} = U_{dim+1} * P(:,1:k+1)          |
c     |   V_{dim}^{+} = V_{dim} * Q(:,1:k)                |
c     %---------------------------------------------------%
c$OMP PARALLEL private(tid,nt,cnk,st,wcnk,wst)
#ifdef _OPENMP
            tid = omp_get_thread_num()
            nt = omp_get_num_threads()
#else
            tid = 0
            nt = 1
#endif
            wcnk = lcwrk/nt
            wst = tid*wcnk+1
            cnk = m/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               wcnk = lcwrk-wst+1
               cnk = m-st+1
            endif
            call csgemm_ovwr_left('n',cnk,k+1,dim+1,U(st,1),ldu,
     c           swork(ip),dim+1,cwork(wst),wcnk)
            cnk = n/nt
            st = tid*cnk+1
            if (tid.eq.nt-1) then
               cnk = n-st+1
            endif
            call csgemm_ovwr_left('n',cnk,k,dim,V(st,1),ldv,swork(iq),
     c           dim,cwork(wst),wcnk)
c$OMP END PARALLEL
            rnorm = swork(ibeta+k-1)
            call cpu_time(t3)
            trestart = trestart + (t3-t2)
            nrestart = nrestart + 1
         endif
         iter = iter + 1
      enddo

 50   if (lsame(jobu,'y') .or. lsame(jobv,'y')) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
         call scopy(dim, swork(ialpha),1,swork(ialpha1),1)
         call scopy(dim, swork(ibeta),1,swork(ibeta1),1)
         lwrk = lwrk + dim**2 + (dim+1)**2
         call critzvec(which, jobu,jobv,m,n,nconv,dim,swork(ialpha1),
     c        swork(ibeta1),U,ldu,V,ldv,swork(ip),
     c        lwrk,cwork,lcwrk,iwork)
      endif

      if (nconv.lt.nwanted) then
         info = -nconv
      endif

      nsing = nconv
      nlandim = dim
      call cpu_time(t1)
      tlansvd = t1-t0
      end
