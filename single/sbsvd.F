c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      subroutine sbsvdstep(jobu,jobv,m,n,k,sigma,D,E,U,ldu,V,ldv)
c
c     Perform one implicit LQ SVD sweep with shift SIGMA.
c

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      character*1 jobu,jobv
      integer m,n,k,ldu,ldv
      real D(*),E(*),U(ldu,*),V(ldv,*),sigma

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i
      real c,s,x,y,r
      logical dou,dov

      logical lsame
      external lsame,slartg,srot

c-------------------- Here begins executable code ---------------------

      if (k.le.1) return

      dou = lsame(jobu,'y')
      dov = lsame(jobv,'y')

c     Compute the initial rotation based on B*B^T-sigma^2
      x = D(1)*D(1) - sigma*sigma
      y = E(1)*D(1)

c     Chase the "bulge" down the lower bidiagonal with Givens rotations.
c     Below 'y' is the "bulge" and 'x' is the element used to eliminate it.
      do i=1,k-1
        if (i.gt.1) then
           call slartg(x,y,c,s,E(i-1))
        else
           call slartg(x,y,c,s,r)
        endif
        x = c*D(i) + s*E(i)
        E(i) = -s*D(i) + c*E(i)
        D(i) = x
        y = s*D(i+1)
        D(i+1) = c*D(i+1)

        if (dou .and. m.gt.0) then
           call srot(m,U(1,i),1,U(1,i+1),1,c,s)
        endif

        call slartg(x,y,c,s,D(i))
        x = c*E(i) + s*D(i+1)
        D(i+1) = -s*E(i) + c*D(i+1)
        E(i) = x
        y = s*E(i+1)
        E(i+1) = c*E(i+1)

        if (dov .and. n.gt.0) then
           call srot(n,V(1,i),1,V(1,i+1),1,c,s)
        endif
      enddo
      call slartg(x,y,c,s,E(k-1))
      x = c*D(k) + s*E(k)
      E(k) = -s*D(k) + c*E(k)
      D(k) = x
      if (dou .and. m.gt.0) then
         call srot(m,U(1,k),1,U(1,k+1),1,c,s)
      endif
      return
      end




      subroutine sbdqr(ignorelast, jobq, n, D, E, c1, c2, Qt, ldq)
      implicit none

c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix
c with diagonal elements d(1)...d(n) and first subdiagonal elements
c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
c If ignorelast.eq..true. then e(n) is assumed to be zero.
c
c If jobq=='Y' then on return Qt contains Q^T.

c     %------------%
c     | Parameters |
c     %------------%
      character*1 jobq
      logical ignorelast
      integer n,ldq
      real D(*),E(*),c1,c2,Qt(ldq,*)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j
      real cs,sn,r

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      logical lsame
      external lsame

c-------------------- Here begins executable code ---------------------

      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         call slaset('A',n+1,n+1,0.0,1.0,Qt,ldq)
      endif
      do i=1,n-1
         call slartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      enddo
      if (.not.ignorelast) then
         call slartg(d(n),e(n),cs,sn,r)
         d(n) = r
         e(n) = 0.0
         c1 = sn
         c2 = cs
         if (lsame(jobq,'Y')) then
            do j=1,n
               Qt(n+1,j) = -sn*Qt(n,j)
               Qt(n,j) = cs*Qt(n,j)
            enddo
            Qt(n,n+1) = sn
            Qt(n+1,n+1) = cs
         endif
      endif
      end




      subroutine srefinebounds(n,k,theta,bound,eps34)
c
c     Refine Lanczos error bounds using the gap theorem.
c     Two-pass algorithm: cluster detection followed by gap refinement.
c
c     Pass 1 identifies clusters of Ritz values within eps34*|theta|
c     of each other. For each cluster, the combined residual and
c     external gap are used to compute a cluster-level bound.
c     Cluster members' bounds are set to max(current, cluster bound)
c     -- never zero, avoiding false convergence of doppelganger Ritz
c     values when orthogonality is lost.
c
c     Pass 2 applies Davis-Kahan gap refinement to individual values.
c
c     Input arguments:
c              n:     smallest dimension of original matrix
c              k:     number of Ritz values to refine
c              theta: array of Ritz values (sorted descending)
c              bound: array of unrefined error bounds
c              eps34: machine epsilon to the power 3/4.

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n,k
      real theta(*), bound(*), eps34

c     %-----------------%
c     | Local variables |
c     %-----------------%
      real gap, r_comb, delta_ext, b_cluster, gl, gr
      integer i,l,cstart,cend

c     %------------------------------------%
c     | External Functions and Subroutines |
c     %------------------------------------%
      real slapy2
      external slapy2

c-------------------- Here begins executable code ---------------------
      if (k.le.1) return

c     %-----------------------------------------------------------%
c     | Pass 1: Cluster detection and cluster-level gap theorem   |
c     %-----------------------------------------------------------%
      i = 1
      do while (i.le.k)
c        Find end of cluster starting at i
         cstart = i
         cend = i
         do while (cend.lt.k)
            if (abs(theta(cend)-theta(cend+1)) .lt.
     c           eps34*abs(theta(cend))) then
               cend = cend + 1
            else
               goto 5
            endif
         enddo
 5       continue
         if (cend.gt.cstart) then
c           Multi-element cluster [cstart..cend]
c           Combined residual: r_comb = sqrt(sum(bound(l)^2))
            r_comb = 0.0
            do l=cstart,cend
               r_comb = slapy2(r_comb, bound(l))
            enddo
c           External gap: distance to nearest theta outside cluster
            gl = huge(1.0)
            gr = huge(1.0)
            if (cstart.gt.1) then
               gl = abs(theta(cstart-1) - theta(cstart))
            endif
            if (cend.lt.k) then
               gr = abs(theta(cend) - theta(cend+1))
            endif
            delta_ext = min(gl, gr)
c           Cluster-level bound via gap theorem
            if (delta_ext.gt.0.0) then
               b_cluster = r_comb*r_comb / delta_ext
            else
               b_cluster = r_comb
            endif
c           Set cluster members' bounds (never zero)
            do l=cstart,cend
               bound(l) = max(bound(l), b_cluster)
            enddo
         endif
         i = cend + 1
      enddo

c     %-----------------------------------------------------------%
c     | Pass 2: Davis-Kahan gap refinement                        |
c     %-----------------------------------------------------------%
      do i=1,k
         if (i.lt.k .or. k.eq.n) then
c
c     We cannot compute a reliable value for the gap of the last
c     Ritz value unless we know it is an approximation to the
c     smallest singular value (k.eq.n). In this case we can take the
c     distance to the next bigger one as the gap, which can really
c     save us from getting stuck on matrices with a single isolated tiny
c     singular value.
c
            if (i.le.1) then
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1))
            else if (i.eq.n) then
               gap = abs(theta(i-1)-theta(i))-max(bound(i-1),bound(i))
            else
               gap = abs(theta(i)-theta(i+1))-max(bound(i),bound(i+1))
               gap = min(gap,abs(theta(i-1) - theta(i)) -
     c              max(bound(i-1),bound(i)))
            endif
            if (gap.gt.bound(i)) then
               bound(i) = bound(i) * (bound(i)/gap)
            endif
         endif
      enddo
      end

