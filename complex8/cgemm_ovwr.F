c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c


      subroutine csgemm_ovwr_left(transb,m,n,k,A,lda,B,ldb,cwork,lcwork)
c
c     compute  A <- A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,lcwork
      complex A(lda,*),cwork(lcwork)
      real B(ldb,*)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (lcwork.lt.n) stop 'Too little workspace in CSGEMM_OVWR_LEFT'
      if (n.gt.k) stop 'Result does not fit in A in CSGEMM_OVWR_LEFT'
      blocksize = int(lcwork/n)
      i = 1
      do i=1,m-blocksize+1,blocksize
         call csgemm(transb,blocksize,n,k, A(i,1),lda,
     c        B,ldb,cwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = cwork(j*blocksize+1+l)
            enddo
         enddo
      enddo
      blocksize = m-i+1
      if (blocksize.le.0) return
      call csgemm(transb,blocksize,n,k,A(i,1),lda,
     c     B,ldb,cwork,blocksize)
      do j=0,n-1
         do l=0,blocksize-1
            A(i+l,j+1) = cwork(j*blocksize+1+l)
         enddo
      enddo
      return
      end

      subroutine  csgemm(transb,m,n,k,A,lda,B,ldb,C,ldc)
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex*8 A(lda,*), C(ldc,*)
      real B(ldb,*)
      external sgemm
      call sgemm('n',transb,2*m,n,k,1.0,A,2*lda,B,ldb,0.0,C,2*ldc)
      end

      subroutine  csgemm_naive(transb,m,n,k,A,lda,B,ldb,C,ldc)
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex A(lda,*), C(ldc,*)
      real B(ldb,*)
      integer i,j,l
      logical lsame
      external lsame

      do j=1,n
         do i=1,m
            C(i,j) = cmplx(0.0, 0.0)
         enddo
      enddo
      if (lsame('T',transb)) then
         do l=1,k
            do j=1,n
               do i=1,m
                C(i,j) = A(i,l)*B(j,l) + C(i,j)
               enddo
            enddo
         enddo
      else
         do l=1,k
            do j=1,n
               do i=1,m
                  C(i,j) = A(i,l)*B(l,j) + C(i,j)
               enddo
            enddo
         enddo
      endif
      end
