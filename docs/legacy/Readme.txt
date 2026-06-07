1. OVERVIEW OF THE PACKAGE

Below is a list giving a brief description of the files found in the
package.
  
============================================================================

  Readme        :  This description.

Main routines
-------------
  lanbpro.m     :  Lanczos bidiagonalization with partial 
                   reorthogonalization (BPRO).
  lansvd.m      :  Sparse SVD routine built on BPRO.
  lanpro.m      :  Lanczos tridiagonalization with partial 
                   reorthogonalization (PRO).
  laneig.m      :  Sparse symmetric eigenvalue routine built on PRO.
                   Similar to the "lanso" subroutine in the LANSO package.


Computational routines
----------------------
  refinebounds.m:  Refines the error bounds based on the gap structure of
                   the Ritz values.
  compute_int.m :  Routine for determining which Lanczos vectors to purge
                   in a reorthogonalization.
  reorth.m      :  Iterated (modified or classical) Gram-Schmidt 
                   reorthogonalization.  
  reorth.mex<...> :  MEX file using a FORTRAN 77 version of reorth.m.	
  tqlb.m        :  Computes the eigenvalues and the top and bottom elements
                   of the eigenvectors for a symmetric tridiagonal matrix.
  tqlb.mex<...> :  The corresponding MEX files call a modified version of
                   the EISPACK routine TQL1 (see below).
  bdsqr.m       :  Computes the singular values and the bottom elements
                   of the left singular vectors for a lower bidiagonal matrix.
  bdsqr.mex<...>:  The corresponding MEX files call the LAPACK routine DBDSQR.

Test scripts and test problems
------------------------------
  test.m        :  Main script for testing sparse SVD calculations.
  testtqlb.m    :  Script for comparing different versions of TQLB.
  Afunc.m       :  Function that computes Y <- A*X.
  Atfunc.m      :  Function that computes Y <- A'*X.
  AtAfunc.m     :  Function that computes Y <- A'*A*X.
  Cfunc.m       :  Function that computes 
                     Y <- [A*X(m+1:m+n,:); A'*X(1:m,:)] 
                   i.e.
                     Y = [ 0  A ]  * X
                         [ A' 0 ].
  helio.mat     :  Testproblem from helioseismology. It is a full 212 x 100
                   matrix with condition number 2.9*10^12.

  Harwell-Boeing/ :
                   Directory with various testmatrices from the 
                   Harwell-Boeing collection (in Matrix Market format).
                   The directory also contains Matlab function for
                   querying, reading and writing matrices stored in
                   Matrix Market file format.


============================================================================

The routine TQLB is used for computing the eigenvalues and top and
bottom elements of the eigenvectors of the Lanczos tridiagonal, which
enter into the error bounds. TQLB is based on an optimized version of
the EISPACK routine TQL1 (see below), and for the sake of speed
FORTRAN versions of TQLB, the reorthogonalization routine REORTH and
the LAPACK routine DBDSQR have been included as MEX files. The current
versions of the MEX files have been compiled for SGI64, SUN and Linux
x86 machines. Purely Matlab based versions are included for use on
other platforms, but these can be rather slow in comparison.


2. IMPLEMENTATION DETAILS

2.1 Modifications to tqlb relative to the version used in the lanso
package:

In TQLB I have added a little extra code to compute the bottom
elements of the normalized eigenvectors.  In addition I have replaced
the calls to PYTHAG in the inner loop of the QL iteration with in-line
code according to the suggestions in

  A. Cline and J. Meyering, "Converting EISPACK to run efficiently on
  a vector processor", Tech. Memo., Pleasant Valley Software,
  Austin TX, 1989.

(see their code in http://www.netlib.org/eispack/3090vf/double/tql1.f). 

Apart from decreasing the execution time by a factor of 3 compared to
the original version, it also causes the small eigenvalues to be
computed more accurately.  As a test I computed the eigenvalues of the
standard testmatrix:

    (  2  -1           )
    ( -1   2 -1        )
T = (     -1  2  .     )   ,
    (         .  . -1  )
    (           -1  2  )

where the true eigenvalues are

  lambda_i = 4 * cos(pi/2 * i / n+1)^2,  i = 1,2,...,n .

Below are the results (T=time in seconds, E=max relative error)
obtained on an SGI Origin 200 with a 180MHz MIPS R10000 CPU, 32Kbytes
primary, 1Mb secondary cache:
 
        Original TQLB       Modified TQLB       MATLAB EIG (eigenvalues only)
-----------------------------------------------------------------------------
   n      T        E           T        E           T       E   
  100    0.020   1.4e-12      0.007   3.1e-13      0.024  6.2e-14 
 1000    1.68    1.2e-9       0.55    2.9e-10      1.98   1.5e-10
10000  154.1     1.9e-6      52.2     8.2e-8      181.2   3.0e-8 



============================================================================== 

Rasmus Munk Larsen, Stanford University, 2000
