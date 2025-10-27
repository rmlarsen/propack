c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c
      subroutine zdgemm_ovwr_left(transb,m,n,k,A,lda,B,ldb,zwork,lzwork)
c
c     compute  A <- A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,lzwork
      complex*16 A(lda,*),zwork(lzwork)
      double precision B(ldb,*)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (lzwork.lt.n) stop 'Too little workspace in ZDGEMM_OVWR_LEFT'
      if (n.gt.k) stop 'Result does not fit in A in ZDGEMM_OVWR_LEFT'
      blocksize = int(lzwork/n)
      i = 1
      do i=1,m-blocksize+1,blocksize
         call zdgemm(transb,blocksize,n,k, A(i,1),lda,
     c        B,ldb,zwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = zwork(j*blocksize+1+l)
            enddo
         enddo
      enddo
      blocksize = m-i+1
      if (blocksize.le.0) return
      call zdgemm(transb,blocksize,n,k,A(i,1),lda,
     c     B,ldb,zwork,blocksize)
      do j=0,n-1
         do l=0,blocksize-1
            A(i+l,j+1) = zwork(j*blocksize+1+l)
         enddo
      enddo
      return
      end

      subroutine  zdgemm(transb,m,n,k,A,lda,B,ldb,C,ldc)
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex*16 A(lda,*), C(ldc,*)
      double precision B(ldb,*)
      external dgemm
      call dgemm('n',transb,2*m,n,k,1d0,A,2*lda,B,ldb,0d0,C,2*ldc)
      end

      subroutine  zdgemm_naive(transb,m,n,k,A,lda,B,ldb,C,ldc)
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldc
      complex*16 A(lda,*), C(ldc,*)
      double precision B(ldb,*)
      integer i,j,l
      logical lsame
      external lsame

      do i=1,m
         do j=1,n
            C(i,j) = dcmplx(0d0,0d0)
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
