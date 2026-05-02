	module polynomials_module
!#######################################################################
	use mpmodule
	implicit none		
	integer, parameter		:: Nmax=150 ! the maximum possible degree of a polynomial
	integer, parameter		:: precision_level=100, nwds=int(precision_level/mpdpw+2)
	type(mp_real)			:: zero, one, eps, pi
	type(mp_complex)		:: one_c, zero_c, ii
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! a polynomial is represented in the form a(0)+a(1)x+a(2)x^2+...+a(N-1)x^(N-1)+a(N)x^N
! the coefficients of polynomials are stored as the first N+1 elements of the vector a(0:Nmax) (of type mp_complex)
! the remaining elements of this vector are set to be zero
	type polynomial
		integer			:: degree
		type(mp_complex)	:: a(0:Nmax)
	end type polynomial
	type real_polynomial
		integer			:: degree
		type(mp_real)		:: a(0:Nmax)
	end type real_polynomial	
! a rational function is represented by two complex polynomials -- numerator and denominator
	type rational 
		type(polynomial)	:: num, den
	end type rational
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	interface assignment (=)
		module procedure p_eq_mpc, R_eq_p, R_eq_mpc, mpc_eq_int		
	end interface
	interface operator (+)
		module procedure add_p_p, add_p_mpc, add_mpc_p, add_R_R, add_p_R, &
				 add_R_p, add_R_mpc, add_mpc_R
	end interface	
	interface operator (-)
		module procedure minus_p, minus_R, sub_p_p, sub_p_mpc, sub_mpc_p, &
		sub_R_R, sub_R_p, sub_p_R, sub_R_mpc, sub_mpc_R
	end interface		
	interface operator (*)
		module procedure mul_p_p, mul_p_mpc, mul_mpc_p, mul_R_R, mul_R_p, mul_p_R, &
		mul_R_mpc, mul_mpc_R, mul_int_mpc
	end interface	
	interface operator (/)
		module procedure div_p_p, div_p_mpc, div_mpc_p, div_R_R, div_R_p, div_p_R, &
		div_R_mpc, div_mpc_R
	end interface	
	interface mp2q
		module procedure mp_real_to_quad, mp_complex_to_quad
	end interface
	interface q2mp
		module procedure quad_to_mp_real, quad_to_mp_complex
	end interface			
!#######################################################################
	contains
!#######################################################################
	subroutine initialize_constants() 
		! eps is a small number, used mostly for comparing if mp_real values are equal to zero
		eps=mpreal('0.1',nwds)**(floor(0.75*precision_level)) 
		zero=mpreal('0.0',nwds) ! mp_real value of zero
		one=mpreal('1.0',nwds)  ! mp_real value of one
		pi=mppi(nwds)		! mp_real value of pi
		ii=mpcmplx(zero,one,nwds)  ! mp_complex value of imaginary unit
		zero_c=mpcmplx(zero,zero,nwds) !mp_complex value of zero
		one_c=mpcmplx(one,zero,nwds)   !mp_complex value of one
	end subroutine initialize_constants
!#######################################################################
	function monomial(n) result(p) ! this function creates a polynomial p(z)=z^n
	integer, intent(in):: n; type(polynomial):: p
	p%a=zero_c; p%a(n)=one; p%degree=n
	end function monomial 	
!#######################################################################
	function find_degree(p) result(d)
	! this function finds degree of a polynomial, by finding the last nonzero coefficient a(j)
	implicit none
	type(polynomial), intent(in):: p; integer:: d
	d=Nmax
	do while ((abs(p%a(d))<eps).and.(d>0)) 
		d=d-1
	end do
	end function find_degree  													
!#######################################################################
!assignment subroutines
	subroutine p_eq_mpc(p,c) ! polynomial = constant
	type(polynomial), intent(out):: p; type(mp_complex), intent(in):: c
	p%a=zero_c; p%a(0)=c; p%degree=0; end subroutine p_eq_mpc

	subroutine R_eq_p(R,p) ! rational function  = polynomial
	type(rational), intent(out):: R; type(polynomial), intent(in):: p
	R%num=p; R%den=one_c; end subroutine R_eq_p

	subroutine R_eq_mpc(R,c) ! rational function  = constant
	type(rational), intent(out):: R; type(mp_complex), intent(in):: c
	R%num=c; R%den=one_c; end subroutine R_eq_mpc
	
	subroutine mpc_eq_int(a,c) ! mp_complex number = integer
	type(mp_complex), intent(out):: a; integer, intent(in):: c
	a=c+zero; end subroutine mpc_eq_int
!#######################################################################
!addition subroutines
	function add_p_p(p,q) result(r)  ! add two polynomials
	implicit none
	type(polynomial), intent(in):: p, q; type(polynomial):: r; integer:: i
	r%a=zero_c; r%degree=max(p%degree,q%degree)
	do i=0,r%degree
		r%a(i)=p%a(i)+q%a(i) 
	end do
	r%degree=find_degree(r); end function add_p_p
	
	function add_p_mpc(p,c) result(r) ! add a polynomial and a constant
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r
	r=p; r%a(0)=r%a(0)+c; end function add_p_mpc
	
	function add_mpc_p(c,p) result(r) ! add a constant and a polynomial
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r
	r=p; r%a(0)=r%a(0)+c; end function add_mpc_p	

	function add_R_R(R1,R2) result(Q) ! add two rational functions
	type(rational), intent(in):: R1, R2; type(rational):: Q
	Q%den=R1%den*R2%den; Q%num=R1%num*R2%den+R1%den*R2%num; end function add_R_R
	
	function add_R_p(R,p) result(Q) ! add a rational function and a polynomial
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%den=R%den; Q%num=R%num+R%den*p; end function add_R_p

	function add_p_R(p,R) result(Q) ! add a polynomial and a rational function
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%den=R%den; Q%num=R%num+R%den*p; end function add_p_R	

	function add_R_mpc(R,a) result(Q) ! add a rational function and a constant
	type(rational), intent(in):: R; type(mp_complex), intent(in):: a; type(rational):: Q
	if (abs(a)>eps) then
		Q%den=R%den; Q%num=R%num+R%den*a; 
	else
		Q=R
	end if
	end function add_R_mpc

	function add_mpc_R(a,R) result(Q) ! add a constant and a rational function
	type(rational), intent(in):: R; type(mp_complex), intent(in):: a; type(rational):: Q
	if (abs(a)>eps) then
		Q%den=R%den; Q%num=R%num+R%den*a; 
	else
		Q=R
	end if
	end function add_mpc_R		
!#######################################################################
!subtraction subroutines 
	function minus_p(p) result(q) ! compute negative of a polynomial
	type(polynomial), intent(in):: p; type(polynomial):: q; integer :: i
	q=p
	do i=0,p%degree	
		q%a(i)=-p%a(i)
	end do
	end function minus_p

	function minus_R(R) result(Q) ! compute negative of a rational function
	type(rational), intent(in):: R; type(rational):: Q
	Q%num=-R%num; Q%den=R%den; end function minus_R
	
	function sub_p_p(p,q) result(r) ! compute polynomial minus a polynomial
	implicit none
	type(polynomial), intent(in):: p, q; type(polynomial):: r; integer:: i
	r%a=zero_c; r%degree=max(p%degree,q%degree)
	do i=0,r%degree	
		r%a(i)=p%a(i)-q%a(i)
	end do
	r%degree=find_degree(r); end function sub_p_p
	
	function sub_p_mpc(p,c) result(r) ! compute polynomial minus a constant
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r
	r=p; r%a(0)=r%a(0)-c; end function sub_p_mpc
	
	function sub_mpc_p(c,p) result(r) ! compute a constant minus a polynomial 
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r
	r=-p; r%a(0)=r%a(0)+c; end function sub_mpc_p		
	
	function sub_R_R(R1,R2) result(Q) ! compute a rational function minus a rational function 
	type(rational), intent(in):: R1, R2; type(rational):: Q
	Q%den=R1%den*R2%den; Q%num=R1%num*R2%den-R1%den*R2%num; end function sub_R_R	

	function sub_R_p(R,p) result(Q) ! compute a rational function minus a polynomial
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%den=R%den; Q%num=R%num-R%den*p; end function sub_R_p
	
	function sub_p_R(p,R) result(Q) ! compute a polynomial minus a rational function
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%den=R%den; Q%num=R%den*p-R%num; end function sub_p_R
	
	function sub_R_mpc(R,c) result(Q) ! compute a rational function minus a constant
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	if (abs(c)>eps) then
		Q%den=R%den; Q%num=R%num-R%den*c; 
	else
		Q=R
	end if
	end function sub_R_mpc
	
	function sub_mpc_R(c,R) result(Q) ! compute a constant minus a rational function
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	if (abs(c)>eps) then
		Q%den=R%den; Q%num=R%den*c-R%num; 
	else
		Q=R
	end if
	end function sub_mpc_R
!#######################################################################
!multiplication subroutines
	function mul_p_p(p,q) result(r)  ! multiply two polynomials
	implicit none
	type(polynomial), intent(in):: p, q; type(polynomial):: r; integer :: n, k
	if (p%degree+q%degree>Nmax) then
		print *, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		print *, 'ERROR in function multiply_polynomials: max degree exceeded!'
		print *, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		stop
	end if	
	r%a=zero_c; r%degree=p%degree+q%degree
	do n=0,r%degree	
		do k=max(0,n-q%degree),min(n,p%degree)
			r%a(n)=r%a(n)+p%a(k)*q%a(n-k)
		end do
	end do
	end function mul_p_p

	function mul_p_mpc(p,c) result(r) ! multiply a polynomial by a constant
	implicit none
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r; integer :: i
	r=p
	do i=0,r%degree	
		r%a(i)=c*r%a(i)
	end do
	end function mul_p_mpc
	
	function mul_mpc_p(c,p) result(r) ! mutliply a constant by a polynomial
	implicit none
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(polynomial):: r; integer :: i
	r=p
	do i=0,r%degree	
		r%a(i)=c*r%a(i)
	end do
	end function mul_mpc_p
	
	function mul_R_R(R1,R2) result(Q) ! multiply two rational functions
	type(rational), intent(in):: R1, R2; type(rational):: Q
	Q%num=R1%num*R2%num; Q%den=R1%den*R2%den; end function mul_R_R
	
	function mul_R_p(R,p) result(Q)   ! multiply a rational function by a polynomial
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%num=R%num*p; Q%den=R%den; end function mul_R_p
	
	function mul_p_R(p,R) result(Q) ! mutiply a polynomial by a rational function
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%num=R%num*p; Q%den=R%den; end function mul_p_R
	
	function mul_R_mpc(R,c) result(Q) ! multipoy a rational function by a constant
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	Q%num=R%num*c; Q%den=R%den; end function mul_R_mpc	
	
	function mul_mpc_R(c,R) result(Q) ! multiply a constant by a rational function
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	Q%num=R%num*c; Q%den=R%den; end function mul_mpc_R
		
	function mul_int_mpc(a,b) result(c) ! multiply an integer by an mp_complex number
	integer, intent(in):: a; type(mp_complex), intent(in):: b; type(mp_complex):: c 
	c=(a*one)*b; end function mul_int_mpc			
!#######################################################################
!division subroutines
	function div_p_mpc(q,c) result(r) ! divide polynomial by a constant
	type(polynomial), intent(in):: q; type(mp_complex), intent(in):: c; type(polynomial):: r
	r=q*(one/c); end function div_p_mpc		

	function div_p_p(p,q) result(R)  ! divide a polynomial by a polynomial (the result is a rational function)
	type(polynomial), intent(in):: p, q; type(rational):: R
	R%num=p; R%den=q; end function div_p_p
	
	function div_mpc_p(c,p) result(R) ! divide a constant by a polynomial (the result is a rational function)
	type(polynomial), intent(in):: p; type(mp_complex), intent(in):: c; type(rational):: R	
	R%num=c; R%den=p; end function div_mpc_p
			
	function div_R_R(R1,R2) result(Q) ! divide two rational functions
	type(rational), intent(in):: R1, R2; type(rational):: Q
	Q%num=R1%num*R2%den; Q%den=R1%den*R2%num; end function div_R_R

	function div_p_R(p,R) result(Q) ! divide a polynomial by a rational function
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%num=p*R%den; Q%den=R%num; end function div_p_R
	
	function div_R_p(R,p) result(Q) ! divide a rational function by a polynomial
	type(rational), intent(in):: R; type(polynomial), intent(in):: p; type(rational):: Q
	Q%num=R%num; Q%den=p*R%den; end function div_R_p
	
	function div_mpc_R(c,R) result(Q) ! divide a constant by a rational function
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	Q%num=c*R%den; Q%den=R%num; end function div_mpc_R
	
	function div_R_mpc(R,c) result(Q) ! divide a rational function by a constant
	type(rational), intent(in):: R; type(mp_complex), intent(in):: c; type(rational):: Q
	Q%num=R%num; Q%den=c*R%den; end function div_R_mpc								
!#######################################################################
	function max_element(p) result(max_a) 
	! find the maximum coefficient of a polynomial
	implicit none
	type(polynomial), intent(in)	:: p
	type(mp_real)			:: max_a
	integer 			:: i
	max_a=zero
	do i=0,Nmax
		if (abs(p%a(i))>max_a) max_a=abs(p%a(i))
	end do
	end function max_element
!#######################################################################
	function is_zero(p) result(t_f)
	! check if the polynomial has all zero coefficients
	implicit none
	type(polynomial), intent(in):: p; logical:: t_f; integer:: i
	t_f=.true.
	do i=0,Nmax
		if (abs(p%a(i))>eps) then 
			t_f=.false.; exit
		end if
	end do
	end function is_zero
!#######################################################################
	function is_non_zero(p) result(t_f)
	! check if the polynomial has some non-zero coefficients
	implicit none
	type(polynomial), intent(in):: p; logical:: t_f; integer:: i
	if (max_element(p)>eps) then 
		t_f=.true.
	else
		t_f=.false.
	end if
	end function is_non_zero
!#######################################################################
	function find_roots(p) result(z)
	implicit none
!this function computes all roots of the real polynomial a(0)+a(1)x+a(2)x^2+...+a(N-1)x^(N-1)+a(N)x^N to precision eps_Newton
!we use Ehrlich-Aberth method combined with Newton's method
	type(real_polynomial), intent(in)	:: p
	type(mp_complex)			:: z(0:Nmax), z0, z1, f0, f_prime, A, B, dz
	type(mp_real)				:: C, theta, r, eps_Newton, eps1, r0, r1
	integer 				:: k, i, j, l, ind, d, root_found(1:Nmax)
	real					:: U, V
	integer, parameter			:: max_iteration=200 !the maximum number of iterations for Ehrlich–Aberth algorithm
	integer, parameter			:: max_iteration_Newton=15 ! the maximum number of iterations for Newton's method	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	eps_Newton=mpreal('0.1',nwds)**(floor(0.5*precision_level))
	eps1=mpreal('0.1',nwds)**7
	z=zero_c; d=p%degree
	call random_seed()
	r0=mpreal('0.1',nwds); r1=mpreal('10.0',nwds)
	do k=1,d !choose initial values of z(k) randomly from the disk |z|<R_max 
		call random_number(U)
		call random_number(V)
		theta=quad_to_mp_real(V+0.0q0)*2*pi ! the argument of z(k) -- chosen from Uniform[0,2*pi] distribution
		r=(r0+(r1-r0)*quad_to_mp_real(U+0.0q0)) ! the absolute value of z(k) -- chosen from Uniform[r0,r1] distribution
		z(k)=r*exp(ii*theta)
	end do
	root_found=0; i=0
	do while ((sum(root_found)<d).and.(i<max_iteration))	
	do k=1,d ! loop through all the roots -- do this in parallel
	if (root_found(k).eq.0) then ! if root number k wasn't found yet -- do one step of Ehrlich–Aberth algorithm
		call poly_prime_eval(p,z(k),f0,f_prime)
		A=f0/f_prime; B=one
		do j=1,d
			if (j.ne.k) B=B-A/(z(k)-z(j))				
		end do
		dz=A/B; z(k)=z(k)-dz
		if (abs(dz)<eps1) then ! if the next approximation of root #k obtained by Ehrlich–Aberth algorithm is 
		!within eps1 from the previous approximation, then use Newton's method to locate this root to desired precision eps_Newton
			z0=z(k); dz=one; j=0
			do while ((abs(dz)>eps_Newton).and.(j<max_iteration_Newton))			
				call poly_prime_eval(p,z0,f0,f_prime) 
				dz=f0/f_prime; z0=z0-dz; j=j+1
			end do			
			if (j<max_iteration_Newton) then !if Newton's method succeeded -- record this root 
				z(k)=z0; root_found(k)=1
				if (abs(aimag(z0))>eps1) then ! if z0 is complex then z1=conjg(z0) is also a root of p
					z1=conjg(z0)
					C=10000*one
					do l=1,d ! find the root closes to z1 and set it to z1
						if (abs(z(l)-z1)<C) then 
							C=abs(z(l)-z1)
							ind=l
						end if	
					end do	
					z(ind)=z1
					root_found(ind)=1			
				end if
			end if
		end if
	end if
	end do
	i=i+1;  !update the approximations to roots of p(z) and go to the next iteration of Ehrlich–Aberth algorithm
	end do
	if (i<max_iteration) then !if all the roots were found in less than max_iteration number
		z(Nmax)=d
!		write(*,'(A,I3,A)') '  function find_roots: all roots z_k of p(z) successfully found after',i,&
!			' iterations of Ehrlich–Aberth algorithm'
		call sort(z,d) !sort the roots in the increasing order of real part
!		C=zero 
!		do k=1,d !find the max value of |p(z(k))|
!			C=max(C,abs(poly_eval(p,z(k))))
!		end do
!		write(*, '(A)', advance='no')  '  max value of |p(z_k)| is ' 
!		call mpwrite(6,40,20,C)	
	else 
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		print *,'ERROR in the function "find_roots": could not find all the roots of p(z)'
		print *,'try increasing precision_level or the number of iterations of iterations of Ehrlich–Aberth algorithm'
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		stop 
	end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end function find_roots
!#######################################################################
	subroutine partial_fraction_decomposition(R,w,z) 
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! this function takes a rational function R(z)=P(z)/Q(z) as an input
! the output is the partial fraction decomposition R(x)=sum_{i=1}^M w_i/(x-z_i)
! here z_1, z_2, ..., z_M are the roots of the polynomial Q(z) (we assume that all roots are simple)
! the output is: the vector of zeros z(i), 1<=i<=M (stored as elements 1 to M of the vector z(i), i=0,1,...,Nmax)
! 		 the vector of coefficients w(i), 1<=i<=M (stored as elements 1 to M of the vector w(i), i=0,1,...,Nmax)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	type(rational), intent(in)	:: R
	type(mp_complex), intent(inout)	:: w(0:Nmax), z(0:Nmax)	
	type(real_polynomial)		:: p, q
	type(mp_real)			:: eps1
	type(mp_complex)		:: f0, f_prime
	integer				:: i, N1, N2, is_non_zero
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	N1=R%num%degree
	N2=R%den%degree
	if (N1>=N2) then
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		print *,'ERROR in the function partial_fraction_decomposition'
		print *,'the degree of numerator of R is not smaller than the degree of denominator'
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'	
		stop 
	end if
	eps1=(one/10)**(floor(precision_level*0.25))	
	is_non_zero=0
	do i=0,N1
		if (abs(aimag(R%num%a(i)))>eps1) then
			is_non_zero=1
			exit
		end if		
	end do
	if (is_non_zero.eq.1) then 
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		print *,'ERROR in the function partial_fraction_decomposition'
		print *,'the numerator of R is not a real polynomial'
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'	
		stop 
	
	end if 
	is_non_zero=0
	do i=0,N2
		if (abs(aimag(R%den%a(i)))>eps1) then
			is_non_zero=1
			exit	
		end if
	end do
	if (is_non_zero.eq.1) then 
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		print *,'ERROR in the function partial_fraction_decomposition'
		print *,'the denominator of R is not a real polynomial'
		print *,'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'	
		stop 
	end if 	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	do i=0,Nmax
		p%a(i)=mpreal(R%num%a(i),nwds)
		q%a(i)=mpreal(R%den%a(i),nwds)	
	end do	
	p%degree=N1
	q%degree=N2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	z=find_roots(q)	
	w=zero_c
	w(Nmax)=N2
	do i=1,N2
		call poly_prime_eval(q,z(i),f0,f_prime) 
		w(i)=poly_eval(p,z(i))/f_prime
	end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end subroutine partial_fraction_decomposition	
!#######################################################################
	function poly_eval(p,z) result(f0)
	! evaluate p(z)  (here p is a polynomial)
	implicit none
	type(real_polynomial), intent(in):: p
	type (mp_complex), intent(in):: z
	type (mp_complex)::f0
	integer:: i, N
	f0=p%a(p%degree)
	do i=(p%degree-1),0,(-1)
		f0=p%a(i)+z*f0
	end do
	end function poly_eval
!#######################################################################
	subroutine poly_prime_eval(p,z,f0,f_prime)
	! evaluate p(z) and p'(z) 
	implicit none
	type(real_polynomial), intent(in):: p
	type (mp_complex), intent(in):: z
	type (mp_complex), intent(out):: f0, f_prime
	integer:: i
	f0=p%a(p%degree); f_prime=zero
	do i=(p%degree-1),0,(-1)
		f_prime=f0+z*f_prime; f0=p%a(i)+z*f0
	end do
	end subroutine poly_prime_eval	
!#######################################################################	
	function mp_complex_to_quad(x_mp) result(x_16)
	! convert mp_complex to complex(kind=16)
	implicit none
	type(mp_complex),intent(in):: x_mp; complex(kind=16):: x_16
	x_16=mp_real_to_quad(mpreal(x_mp,nwds))+(0.0q0,1.0q0)*mp_real_to_quad(aimag(x_mp))		
	end function mp_complex_to_quad
!#######################################################################
	function mp_real_to_quad(x_mp) result(x_16)
	! convert mp_real to real(kind=16)
	implicit none
	type(mp_real),intent(in):: x_mp; type(mp_real):: z; real(kind=16):: x_16; integer:: i, k, j(1:5), power_of_ten
	z=x_mp; i=nint(log(abs(dble(z)))/log(10.0d0))
	z=z/(one*10)**i; power_of_ten=10**7
	do k=1,5
		j(k)=int(dble(z*power_of_ten)); z=z*power_of_ten-j(k)
	end do
	x_16=0.0q0
	do k=1,5
		x_16=x_16+j(k)/10.0q0**(7*k)
	end do
	x_16=x_16*10.0q0**i	
	end function mp_real_to_quad
!#######################################################################
	function quad_to_mp_real(x_16) result(x_mp)
	! convert real(kind=16) to mp_real
	implicit none
	real(kind=16),intent(in):: x_16; type(mp_real):: x_mp; real(kind=16):: z; integer:: i, k, j(1:5), power_of_ten
	z=x_16; i=nint(log(abs(dble(z)))/log(10.0d0))
	z=z/10.0q0**i; power_of_ten=10**7
	do k=1,5
		j(k)=floor(z*power_of_ten); z=z*power_of_ten-j(k)
	end do
	x_mp=zero
	do k=1,5
		x_mp=x_mp+j(k)/(one*10)**(7*k)
	end do
	x_mp=x_mp*(one*10)**i	
	end function quad_to_mp_real
!#######################################################################
	function quad_to_mp_complex(z_16) result(z_mp)
	! convert complex(kind=16) to mp_complex
	implicit none
	complex(kind=16),intent(in):: z_16; type(mp_complex):: z_mp
	z_mp=quad_to_mp_real(z_16%re)+ii*quad_to_mp_real(z_16%im)	
	end function quad_to_mp_complex	 
!#######################################################################	
	subroutine sort(z,n) 
	! sort elements z_i, i=1,...,n in the increasing order of real part
	implicit none
	integer, intent(in)			:: n
	type (mp_complex), intent(inout) 	:: z(0:Nmax)
	type (mp_complex)			:: z0
	integer					:: i
	i=1
   	do while (i<=n)
	   	if (i==1)  i=i+1	   	
	       	if (mpreal(z(i),nwds) > mpreal(z(i-1),nwds)-eps) then
			i=i+1
	       	else
	       		z0=z(i); z(i)=z(i-1); z(i-1)=z0; i=i-1
		end if
	end do
	end subroutine sort			
!#######################################################################	
	end module polynomials_module
