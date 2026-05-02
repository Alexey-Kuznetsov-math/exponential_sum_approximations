	program main_hockey_stick
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	use mpmodule
	use polynomials_module
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	integer, parameter	:: M=30	
	integer, parameter	:: n_infinity=4
	real, parameter		:: A=0.0
	real, parameter		:: B=83.0	
	logical, parameter	:: explicit_Laplace_transform=.true.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
	integer, parameter	:: p=2*M-n_infinity 
	integer, parameter	:: Nt=1000 !the number of discretization points for computing F(z_j) via double exponential quadrature
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	type(mp_complex)	:: W, z(1:p)
	type(mp_real)		:: v
	complex(kind=16)	:: c(1:M), lambda(1:M)
	integer			:: j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (precision_level>mpipl) then
		print *, 'Increase default precision in module MPFUNF'
		stop
	end if	
	if ((2*M+2)>Nmax) then
		print *,'Increase parameter Nmax in polynomials_module'
		stop
	end if		
	call initialize_constants()
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	W=q2mp(A+0.0q0)+ii*q2mp(B+0.0q0)
	do j=1,p
		v=abs(j-(p+one)/2)/((p-one)/2)
		if (2*j>p+1) then
			z(j)=v*W
		else
			z(j)=v*conjg(W)
		end if
	end do		
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	call compute_exp_sum_coefficients(z,c,lambda)
	open(1,file='c_30.txt')
	open(2,file='lambda_30.txt')
		do j=1,M
			write(1,'(2ES45.35)') real(c(j)), aimag(c(j))
			write(2,'(2ES45.35)') real(lambda(j)), aimag(lambda(j))
		end do
	close(1)
	close(2)	
!#######################################################################
contains
!#######################################################################
	function f(x) result(fx)
	implicit none
	type (mp_real), intent(in)		:: x
	type (mp_real)				:: fx	
	fx=max(one-x,zero)
	end function f
!#######################################################################
	function f_derivatives(n) result(f1)
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	integer, intent(in)			:: n
	type (mp_complex)			:: f1(0:n)
	integer                                 :: j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!this function returns f1(j)=f^{(j)}(0) for j=0,1,...,n
!where f(x) is the function defined above
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
	f1=zero_c
  	f1(0)=one  	
  	if (n>0) then
  		f1(1)=-one
  	end if	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end function f_derivatives	
!#######################################################################
	function Laplace_transform(z) result(fz)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! if the Laplace transform of f(x) is known explicitly, it should be entered here
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	implicit none
	type (mp_complex), intent(in)		:: z
	type (mp_complex)			:: fz	
	if (abs(z)>eps) then
		fz=(z+exp(-z)-one)/z**2
	else
		fz=one/2
	end if		
	end function Laplace_transform	
!#######################################################################	
subroutine compute_exp_sum_coefficients(z,c_16,lambda_16)
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(in)	:: z(1:p)
	complex (kind=16), intent(inout):: c_16(1:M), lambda_16(1:M)
	type (mp_complex)		:: coeffs(1:p), coeffs_infty(0:(n_infinity+2)), c_mp(0:Nmax), &
		lambda_mp(0:Nmax), a(1:p), a_inf(0:(n_infinity+2))		
	type (rational)			:: R
	integer				:: i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	call generate_coefficients(z,a,a_inf)
	call find_CF_coeffs(z,a,a_inf,coeffs,coeffs_infty) 
	call write_CF_as_rational_function(z,coeffs,coeffs_infty,R)	
	call partial_fraction_decomposition(R,c_mp,lambda_mp) 	
	do i=1,M
		c_16(i)=mp2q(c_mp(i))
		lambda_16(i)=-mp2q(lambda_mp(i))
	end do	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end subroutine compute_exp_sum_coefficients
!#######################################################################
	subroutine generate_coefficients(z,a,a_inf)
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(in)		:: z(1:p)
	type (mp_complex), intent(inout)	:: a(1:p), a_inf(0:(n_infinity+2))	
	integer                                 :: i, j, k, k_max
	type (mp_real)				:: h, t(-Nt:Nt), mu(-Nt:Nt)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!we compute a(i)=F(z_i)=int_0^{\infty} exp(-x z_i) f(x) d x for i=1,2,...,p
!where F(z)=\int_0^{\infty} exp(-zx)f(x) d x is the Laplace transform of f(x)
!these integrals are computed via double exponential quadrature 
!the coefficients a_inf(i) are computed as a_inf(0)=0 and a_inf(i)=f^{i-1}(0) for i=1,2,...,n_infinity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (explicit_Laplace_transform.eqv..false.) then
		h=mpreal('6.0',nwds)/Nt
		k_max=Nt
		!$omp parallel sections
		!$omp section
		do k=-Nt,0  !precompute the nodes and weights of the double exponential quadrature 		
			t(k)=exp(k*h-exp(-k*h)) 
			mu(k)=f(t(k))*t(k)*(one+exp(-k*h))*h 
		end do 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		!$omp section	
		do k=1,Nt  !precompute the nodes and weights of the double exponential quadrature 		
			t(k)=exp(k*h-exp(-k*h)) 
			mu(k)=f(t(k))*t(k)*(one+exp(-k*h))*h 
			if (abs(mu(k))<eps) then
				k_max=k
				exit
			end if	
		end do  	 
		!$omp end parallel sections 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		a=zero_c
		!$omp parallel do private(k) shared(mu,z)
		do i=1,(p+1)/2
			do k=-Nt,k_max	
				a(i)=a(i)+mu(k)*exp(-z(i)*t(k))		
		  	end do	
		end do
		!$omp end parallel do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
	else 
		do i=1,(p+1)/2
			a(i)=Laplace_transform(z(i))
		end do
	end if	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~		
	do i=1,(p/2)   ! we have z(i)=conjg(z(p+1-i)), thus c(p+1-i)=conjg(c(i))
			a(p+1-i)=conjg(a(i))			
	end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	a_inf=zero_c
	a_inf(1:n_infinity)=f_derivatives(n_infinity-1)
	a_inf(n_infinity+2)=n_infinity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	end subroutine generate_coefficients					
!#######################################################################
	subroutine find_CF_coeffs(z,a,a_inf,coeffs,coeffs_infty) 
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(in)		:: z(1:p)
	type (mp_complex), intent(inout)	:: a(1:p), a_inf(0:(n_infinity+2)), coeffs(1:p), coeffs_infty(0:(n_infinity+2))
	integer                                 :: n, Na, i, j
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	do i=1,p
		coeffs(i)=a(i)	
		do j=i+1,p
			a(j)=(a(i)/a(j)-one)/(z(j)-z(i))
		end do
		call update_coeffs_at_infinity_from_a_finite_point(a_inf,a(i),z(i))
	end do
	coeffs_infty=compute_coeffs_at_infinity(a_inf)	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end subroutine find_CF_coeffs
!#######################################################################
	subroutine write_CF_as_rational_function(z,d,c,R) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!this function takes the coefficients c(0:n2) of the continued fraction at infinity
!h(z)=c(0)+c(1)/(z+c(2)*z/(z+c(3)*z/(.... +c(n2-1)z/(z+c(n2)))))
!and the coefficients d(1:n1) and and returns the rational function 
!R(z)=P(z)/Q(z)=d(1)/(1+d(2)*(z-w(1))/(1+d(3)*(z-w(2))/(....(1+d(n1)*(z-w(n1-1))/(1+(z-w(n1))*h(z)))))))
!here n2=c(n_infinity+2) and n1=p
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	type(mp_complex), intent(in)	:: z(1:p), d(1:p), c(0:(n_infinity+2))
	type(rational), intent(inout)	:: R
	type(rational)			:: R1
	type(polynomial)		:: w
	type(mp_complex)		:: b
	integer				:: i, j, n2, m1, m2, L
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	n2=anint(abs(c(n_infinity+2)))
	w=monomial(1)
	if (n2.eq.0) then
		R=c(0)
	elseif (n2.eq.1) then
		R=c(0)+c(1)/w
	elseif (n2.eq.2) then
		R=c(0)+c(1)/(w+c(2))		
	else 	! compute R1(w)=c(0)+c(1)*w/(1+c(2)*w/(1+c(3)*w/...+c(n2-1)/(1+c(n2)*w))) and then set R(w)=R1(1/w)
		R1=one_c+w*c(n2)
		do i=(n2-1),2,(-1)
			R1=one_c+c(i)*w/R1
		end do
		R1=c(0)+c(1)*w/R1
		m1=R1%num%degree
		m2=R1%den%degree
		L=max(m1,m2)
		R=zero_c
		do j=0,L
			R%num%a(j)=R1%num%a(L-j)
			R%den%a(j)=R1%den%a(L-j)
		end do
		R%num%degree=find_degree(R%num)
		R%den%degree=find_degree(R%den)
	end if		
	R=one_c+(w-z(p))*R
	do i=p,2,-1
		R=one_c+d(i)*(w-z(i-1))/R
	end do
	R=d(1)/R
	b=R%den%a(0)
	R%num=R%num/b
	R%den=R%den/b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end subroutine write_CF_as_rational_function	
!#######################################################################
	function compute_coeffs_at_infinity(a) result(c)
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(in)		:: a(0:(n_infinity+2))
	type (mp_complex)			:: c(0:(n_infinity+2)), b(0:(n_infinity+2))
	integer                                 :: n, Na
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!this function takes as iput the coefficients a0,a1,a2,... of
!the expansion f(z)~a0+a1/z+a2/z^2+... (as z-> infty)
!and returns the coefficients c0=a0,c1, c2, ... of the continued fraction
!f(z)=c0+c1/(z+c2*z/(z+c3*z/(.... (z+cN))))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Na=anint(abs(a(n_infinity+2)))
	c=zero_c
	c(n_infinity+2)=a(n_infinity+2)
	if (Na.eq.0) then
		c(0)=a(0)
		c(n_infinity+2)=zero
	elseif (Na.eq.1) then	
		c(0)=a(0)	
		c(1)=a(1)
		c(n_infinity+2)=one
	else
		c(0)=a(0)	
		c(1)=a(1)
		b=a
		b(0)=zero		
		do n=2,Na-1
			call update_coeffs_at_infinity(b)
			c(n)=b(1)
		end do	
		c(Na)=-b(2)/b(1)		
	end if
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	end function compute_coeffs_at_infinity	
!#######################################################################
	subroutine update_coeffs_at_infinity_from_a_finite_point(a,c,z1) 
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(inout)	:: a(0:(n_infinity+2))
	type (mp_complex), intent(in)		:: c, z1
	type (mp_complex)			:: b(0:(n_infinity+2)), u
	integer                                 :: n, Na, i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!this function take as input the coefficients a_i of f(z)~a0+a1/z+a2/z^2+...
!where only the coefficients with i=0,1,...,a(Na) are used 
!and computes the coefficients b_i (for i=0,1,...,Na) of g~b0+b1/z+b2/z^2+...
!where f(z)=c/(1+(z-z1) g(z))
!the element b(n_infinity+2) is an integer, that represents the number of relevant coefficients  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Na=anint(abs(a(n_infinity+2)))
	b=zero_c
	if (abs(a(0))>eps) then !if a(0) is not zero
		b(1)=c/a(0)-one 
		do n=1,Na	
			u=a(n)	
			do i=1,n
				u=u+b(i)*(a(n+1-i)-z1*a(n-i))
			end do	
			b(n+1)=-u/a(0)
		end do
		b(n_infinity+2)=a(n_infinity+2)+one
	else	!if a(0) is zero
		b(0)=c/a(1)
		do n=1,Na-1
			u=a(n)	
			do i=0,(n-1)
				u=u+b(i)*(a(n+1-i)-z1*a(n-i))
			end do	
			b(n)=-u/a(1)
		end do
		b(n_infinity+2)=a(n_infinity+2)-one
	end if
	a=b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end subroutine update_coeffs_at_infinity_from_a_finite_point
!#######################################################################
	subroutine update_coeffs_at_infinity(a)
	implicit none
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
	type (mp_complex), intent(inout)	:: a(0:(n_infinity+2))
	type (mp_complex)			:: u, b(0:(n_infinity+2))
	integer                                 :: n, Na, i
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!this function take as iput the coefficients a_i of f(z)~a0+a1/z+a2/z^2+...
!where only the coefficients with i=0,1,...,a(Na) are used (the rest are zeros)
!and returns the coefficients of g~b0+b1/z+b2/z^2+...
!where f(z)=a0+a1/(z + z g(z))
!this function returns coefficients b_i for i=0,1,...,Nb=Na-1
!(the remaining values are filled by zeros and we set b(n_infinity+2)=Nb)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Na=anint(abs(a(n_infinity+2)))
	b=zero_c
	b(n_infinity+2)=a(n_infinity+2)-one
	do n=1,Na-1
		u=a(n+1)
		do i=1,n-1
			u=u+b(i)*a(n+1-i)
		end do
		b(n)=-u/a(1)
	end do
	a=b
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	end subroutine update_coeffs_at_infinity		
!#######################################################################
	end program main_hockey_stick
