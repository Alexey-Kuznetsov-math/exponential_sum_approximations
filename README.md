# exponential_sum_approximations
Fortran code for the paper "Approximating functions on R^+ by exponential sums", https://doi.org/10.1016/j.cam.2026.117756

This code was tested on Ubuntu with the gfortran compiler. 

HOW TO RUN THE CODE: 

1) Download the MPFUN2020 arbitrary precision package from https://www.davidhbailey.com/dhbsoftware/

2) Copy required files from the MPFUN2020 package
 
mpfuna.f90 mpfunb.f90 mpfunc.f90 mpfund.f90 mpfune.f90 
mpfunf.f90 mpfung2.f90 mpfunh2.f90 mpmodule.f90 mpmask13.f90 second.f90

into the same directory as the following project files: 

Makefile
main_Gaussian.f90
main_hockey_stick.f90
polynomials_module.f90

3) Open a terminal, navigate to the directory containing all the files, and execute the following commands:

make mp
make
make run

The first command compiles all MPFUN2020 modules. The second compiles the program main_hockey_stick.f90. The third command runs the program.  
 
The program should produce the results for the 30-term exponential sum approximation to the hockey stick function. The computed coefficients c_j and lambda_j will be saved (in quadruple precision) in the files 

c_30.txt
lambda_30.txt

ADDITIONAL COMMENTS: 

a) The MPFUN2020 modules need to be compiled only once. If you modify parameters in main_hockey_stick.f90, you 
can recompile and run with

make && make run

b) To run the program main_Gaussian.f90 (which computes exponential sum approximations to the Gaussian function), uncomment line #9 in the Makefile (and comment out line #8). 

To run this program in parallel, ensure that the library "libomp-dev" is installed, and you also need to uncomment line #6 in the Makefile (and comment out line #4) before compiling and running this program.  

APPROXIMATING OTHER FUNCTIONS: 

To find exponential sum approximations for a different function f (other than hockey stick or Gaussian), edit the following functions in the main program: 

	function f(x) 
	
	function f_derivatives(n) 
	
If the Laplace transform of f is known in closed form, enter it in the function
	
	function Laplace_transform(z) 
	
and set 

	explicit_Laplace_transform=.true.
	
If the Laplace transform of f is not known in closed form, set 

	explicit_Laplace_transform=.false. 
	
In this case, you do not need to edit *Laplace_transform* function, as it will not be used by the algorithm.	

