This folder contains a Matlab interface for the c-port by Naoaki Okazaki 
(http://www.chokkan.org/software/liblbfgs/) of the lbfgs code
by Jorge Nocedal (http://www.ece.northwestern.edu/~nocedal/lbfgs.html),
with the addition of the Orthant-Wise Limited-memory Quasi-Newton algorithm
due to Galen Andrew and Jianfeng Gao. See references below.

To install:

1. Download liblbfgs
2. Make necessary modification in make.m such that the files "lbfgs.c" and "lbfgs.h" (and optionally files for SSE optimziations)
   can be found. This should normally just be to change the "library_path"-variable.
3. Run make.m in Matlab

If liblbfgs is compiled outside Matlab, then the installation procedure looks as follows:

1. Download liblbfgs
2. Follow instructions to compile liblbfgs
3. Make necessary modification in make.m such that the library can be found
4. Run make.m in Matlab

(Note that it is currently assumed that lbfgsfloatval_t coincide with matlabs double.)


Run example.m for an example.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Tested with liblbfgs v.1.9 and Matlab R2010a and R2008a, 32-bit, on Windows Vista
and Matlab R2011a on Windows 7 64-bit and MSVC++ 2010 express with Windows SDK 7.1

Update history:

7/12 2011:
Modified the file "make.m" so as to make more compilation options available

21/2 2011: 
Changed the name of the file "make.m" to "build.m"
Changed in lbfgs.c so that the solver always returns the current value of x at exit
Changed "max_iterations" to "MaxIter" to comply with fmincon synatax


References:

The L-BFGS algorithm is described in:
    - Jorge Nocedal.
      Updating Quasi-Newton Matrices with Limited Storage.
      <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
    - Dong C. Liu and Jorge Nocedal.
      On the limited memory BFGS method for large scale optimization.
      <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

The line search algorithms are described in:
    - John E. Dennis and Robert B. Schnabel.
      <i>Numerical Methods for Unconstrained Optimization and Nonlinear
      Equations</i>, Englewood Cliffs, 1983.
    - Jorge J. More and David J. Thuente.
      Line search algorithm with guaranteed sufficient decrease.
      <i>ACM Transactions on Mathematical Software (TOMS)</i>, Vol. 20, No. 3,
      pp. 286-307, 1994.

The Orthant-Wise Limited-memory Quasi-Newton (OWL-QN) method is presented in:
    - Galen Andrew and Jianfeng Gao.
      Scalable training of L1-regularized log-linear models.
      In <i>Proceedings of the 24th International Conference on Machine
      Learning (ICML 2007)</i>, pp. 33-40, 2007.


TODO: Make sure input validation is not duplicated
TODO: Improve documentation
TODO: lbfgs returns the zero vector if the minimzation fails for some reason,
      should (whenever appropriate) return the current value of x
TODO: Make sure we don't get into trouble by mixing double and lbfgsfloatval_t
TODO: Return gradient at x
TODO: Display function should be written in Matlab
