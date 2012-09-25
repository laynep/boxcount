module box_count
  use types, only : dp
  use sorters
  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Written by Layne Price, University of Auckland, May 2012
!This is adapted from a number of articles written online about box-counting, as well as
!using some Numerical Recipe techniques for fitting.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This is a series of subroutines which will implement a simple box counting method.  The main boxcount subroutine outputs an array, temp (nX2), which expresses the {box size, box number} n times.  This is then plotted using linear least squares regression on the values which are deemeed to be in the linear regime.  The fractal dimension is then calculated as the slope of this line.

!The subroutines and functions are as follows:
!1.	fractal_dimn(counted,dimn,dev,chi2) SUBROUTINE: calcs fractal dimension on table 				counted AFTER box counting has been done with SUBR boxcount.
!2.	simlp(points,sad,alpha,beta,d,iter,inext,ifault) SUBROUTINE: Linear fit by 				minimizing absolute deviation.  From INTERNET.
!3.	leastsqr(points,a,b,siga,sigb,chi2) SUBROUTINE: Linear least squares reg from 						NumRec.
!4.	boxcount(list,epsilon_0,scaling,minim,temp) SUBROUTINE: Sorts the list and performs 					boxcounting.
!5.	uniq_elts(list)	FUNCTION: Calcs the uniq elemnts of a fully SORTED integer list.
!6.	init_random_seed()	SUBROUTINE: sets random seed different dependent on time.  					From the GNU Fortran webpage.
!7.	findmin(array,i)	FUNCT: Finds min of array in ith column.
!8.	findmax(array,i)	FUNCT: Finds max of array in ith column.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

interface boxcount
  module procedure boxcount_dp
  module procedure bound_boxcount
end interface


contains

!This subroutine calculates the fractal dimension after box counting has been performed.  It selects the points in the linear regime by finding the max and min values of log(N) and determining the "midpoint" of the dataset with respect to the log(N) direction.  It then systematically throws away datapoints from the end, plots it, and gives out the chi**2 value, until it's satisfied.

subroutine fractal_dimn(counted,dimn,dev,chi2)
	implicit none

	real(dp), dimension(:,:), intent(in) :: counted   !will be n x 2 array.
	real(dp), intent(out) :: dimn, dev, chi2
	real(dp) :: sad
	real(dp), dimension(size(counted,1)) :: log_n, log_d
	real(dp), dimension(:,:), allocatable :: linear
	real(dp) :: a, b, maxim, minim, delta, low, high, siga,sigb, averagey, averagex
	integer :: i, j,k,l, n, size1, averagex_numb, dropr, dropl,  numbr, numbl, counter
	logical :: fit
	integer :: iter, ifault
	real(dp), dimension(:), allocatable :: d
	integer, dimension(:), allocatable :: inext

	n = size(counted,1) !also size of log_n, log_d
	
	!calc the logs.
	do i=1, n
		log_d(i) = -1_dp*log(counted(i,1))
		log_n(i) = log(counted(i,2))
	end do

	!find max/min values.  
	maxim = log_n(n)
	minim = 0_dp
	averagey = (maxim-minim)*.8_dp !position at upper 80% of range.

doj:	do j=1,n-1
		if(log_n(j)<=averagey .and. log_n(j+1)>averagey) then
			averagex = log_d(j+1)
			averagex_numb = j+1
			exit doj
		end if
	end do doj

	dropr=0
	dropl=0

	print*, "type		", "chi**2		", "dimn		", "dev"

do1:	do l=0,n-5
		size1 = n - dropr - dropl
		allocate(linear(size1,2))
		
		do k=1,size1
			linear(k,1)=log_d(k+dropl)
			linear(k,2)=log_n(k+dropl)
		end do

		!fit least squares to linear area.
		call leastsqr(linear,a,b,siga,sigb,chi2)
		if(chi2<20) print*, "leastsqr: ",chi2, b, sigb

		if(chi2 < .01) then
			fit = .true.
			exit do1
		else
			fit = .false.
		end if
		deallocate(linear)

		!determine how many pts to left and right.
		numbr = n -averagex_numb - dropr
		numbl = averagex_numb - dropl -1
		!drop the outlying points until get good chi**2 linear fit.
		if(numbr<=numbl) then  
			dropl = dropl+1
		else
			dropr = dropr+1
		end if

	end do do1

	print*,"number of right drops ",dropr
	print*,"number of left drops ",dropl
	print*,"number of points used ",(n-dropr - dropl)

	if(fit .eqv. .false.) then
		print*,"error: couldn't find a fit within specified chi**2."
	end if


	!fractal dimension = slope (b).
	dimn = b
	dev = sigb

end subroutine fractal_dimn

!***********************************************************************************
!Subroutine to fit a line (alpha + beta*x + error) by minimizing absolute deviation.  From http://jblevins.org/mirror/amiller/as132.f90 .

!Input are the arrays x,y of length n.  Output is sad~sum of absolute deviations (like chi**2), alpha, beta.

subroutine simlp(points, sad, alpha, beta, d, iter, inext, ifault)
 
! code converted using to_f90 by alan miller
! date: 2003-06-19  time: 09:37:54

!     algorithm as 132  appl. statist. (1978) vol.27, no.3

!     simlp:   fit  y = alpha + beta.x + error

	implicit none

real(dp), dimension(:,:), intent(in) :: points
integer :: n
real(dp), dimension(size(points,1)) :: x, y
!integer, intent(in)   :: n
!real, intent(in out)  :: x(n)
!real, intent(in)      :: y(n)
real(dp), intent(out)     :: sad
real(dp), intent(out)     :: alpha
real(dp), intent(out)     :: beta
real(dp), dimension(:), intent(out), allocatable :: d
integer, intent(out)  :: iter
integer, dimension(:), allocatable, intent(out)  :: inext
integer, intent(out)  :: ifault

real(dp), parameter  :: acu = 1.0e-06_dp, big = 1.0e19_dp, half = 0.5_dp, zero = 0.0_dp,   &
                    one = 1.0_dp, two = 2.0_dp
integer  :: i, ibas1, ibas2, iflag, iin, iout, isave, j, i_1
real(dp)     :: a1, a2, aaa, aaaa, ahalf, aone, bbb, bbbb, ddd, det,   &
            ratio, rho, subt, sum, t, test, tot1, tot2, y1, y2, zzz

!change the table to two vectors of length n.
n = size(points,1)
allocate(d(n))
allocate(inext(n))
do i_1=1,n
	x(i_1)=points(i_1,1)
	y(i_1)=points(i_1,2)
end do



!     initial settings

ifault = 0
iter = 0
ahalf = half + acu
aone = ahalf + ahalf

!     determine initial basis

d(1) = zero
y1 = y(1)
ibas1 = 1
a1 = x(1)
do  i = 2, n
  if (abs(a1-x(i)) >= acu) then
    a2 = x(i)
    ibas2 = i
    y2 = y(i)
    go to 20
  end if
end do
ifault = 1
return

!     calculate initial beta value

20 det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det

!     calculate initial d-vector

do  i = 2, n
  ddd = y(i) - (aaaa+bbbb*x(i))
  d(i) = sign(one,ddd)
end do
tot1 = one
tot2 = x(ibas2)
d(ibas2) = -one
do  i = 2, n
  tot1 = tot1 + d(i)
  tot2 = tot2 + d(i) * x(i)
end do
t = (a2*tot1-tot2) * det
if (abs(t) >= aone) then
  det = -det
  go to 70
end if

!     main iterative loop begins

50 t = (tot2-a1*tot1) * det
if (abs(t) < aone) go to 130
iflag = 2
iout = ibas2
x(iout) = a1
aaa = a1
bbb = a2
go to 80

60 t = (tot2-a2*tot1) * det
if (abs(t) < aone) go to 130

70 iflag = 1
bbb = a1
aaa = a2
iout = ibas1

80 rho = sign(one,t)
t = half * abs(t)
det = det * rho

!     perform partial sort of ratios

inext(ibas1) = ibas2
ratio = big
sum = ahalf
do  i = 1, n
  ddd = (x(i)-aaa) * det
  if (ddd*d(i) > acu) then
    test = (y(i)-aaaa-bbbb*x(i)) / ddd
    if (test < ratio) then
      j = ibas1
      sum = sum + abs(ddd)
      90 isave = abs(inext(j))
      if (test < d(isave)) then
        if (sum >= t) then
          subt = abs((x(isave)-aaa)*det)
          if (sum-subt >= t) then
            sum = sum - subt
            d(isave) = sign(1,inext(j))
            inext(j) = inext(isave)
            go to 90
          end if
        end if
        100 j = isave
        isave = abs(inext(j))
        if (test < d(isave)) go to 100
      end if
      inext(i) = inext(j)
      inext(j) = sign(i,int(d(i)))
      d(i) = test
      if (sum >= t) then
        iin = abs(inext(ibas1))
        ratio = d(iin)
      end if
    end if
  end if
end do

!     update basic indicators

iin = abs(inext(ibas1))
j = iin
120 isave = abs(inext(j))
if (isave /= ibas2) then
  zzz = sign(1,inext(j))
  tot1 = tot1 - zzz - zzz
  tot2 = tot2 - two * zzz * x(isave)
  d(isave) = -zzz
  j = isave
  go to 120
end if
zzz = sign(1,inext(ibas1))
tot1 = tot1 - rho - zzz
tot2 = tot2 - rho * bbb - zzz * x(iin)
d(iout) = -rho
iter = iter + 1
if (iflag /= 1) then
  x(ibas2) = a2
  ibas2 = iin
  d(ibas2) = -one
  a2 = x(iin)
  y2 = y(iin)
  det = one / (a1-a2)
  aaaa = (a1*y2-a2*y1) * det
  bbbb = (y1-y2) * det
  go to 60
end if
ibas1 = iin
a1 = x(iin)
d(ibas1) = zero
y1 = y(iin)
det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det
go to 50

!     calculate optimal sum of absolute deviations

130 sad = zero
do  i = 1, n
  d(i) = y(i) - (aaaa+bbbb*x(i))
  sad = sad + abs(d(i))
end do
alpha = aaaa
beta = bbbb

return
end subroutine simlp


!****************************************************************************

!Linear least-squares regression on the linear regime.  Inputs: points(z,2)~the (x,y) points to plot starting with z=1--needs to be constrained to linear regime prior to call to this subroutine.  Outputs: a,b~(a+b*x) and sig's~main standard deviation.  From Numerical Recipes with unavailable standard deviations.

subroutine leastsqr(points,a,b,siga,sigb,chi2)
	implicit none

	real(dp), dimension(:,:), intent(in) :: points  !an n x 2 array of n 2-d points (x,y).
	real(dp), intent(out) :: a, b, siga, sigb, chi2
	real(dp) :: sx, sy, st2, sxoss, t, sigdat, q, ss
	integer :: i, j, k, n
	
	sx = 0_dp
	sy = 0_dp
	st2 = 0_dp
	b = 0_dp

	n = size(points,1)

	do i=1,n
		sx = sx + points(i,1)
		sy = sy + points(i,2)
	end do
	
	ss = dble(n)
	sxoss = sx/ss

	do j=1,n
		t = points(j,1) - sxoss
		st2 = st2 +(t*t)
		b = b + (t*points(j,2))
	end do

	!calculate a, b, and probable uncertainties siga, sigb.

	b=b/st2
	a=(sy-(sx*b))/ss
	siga = sqrt((1_dp + (sx*sx)/(ss*st2))/ss)
	sigb = sqrt(1_dp/st2)

	chi2 = 0_dp

	do k=1,n
		chi2 = chi2 + (points(k,2)-a-(b*points(k,1)))**2_dp
	end do
	q=1
	sigdat = sqrt(chi2/(n-2))
	siga = siga*sigdat
	sigb = sigb*sigdat

end subroutine leastsqr


!***************************************************************************************

!Performs the boxcount on an n x m dimensional list of n points (x_1,...,x_m).  It first sorts this data, then rescales it according to the grid size.  It then finds how many boxes it takes to cover the set by doing the CEILING count below and identifying the unique elements.

!Output is the array temp(x,2) of points (delta,N) according to the input scaling.

!NOTE: the list put into boxcount requires it to be sorted from low to high in first column.

subroutine boxcount_dp(list,minim,temp)
	implicit none

	real(dp), dimension(:,:), intent(inout) :: list
	integer, dimension(:,:), allocatable :: ceillist
	real(dp), intent(in) :: minim
	real(dp) :: scaling, m, p, d, epsilon_0
	real(dp), dimension(size(list,2)) :: maxelt, minelt, totalrange
	integer :: i,i_1,j,k,l, boxnumb
	integer :: nmax
	real(dp), dimension(:,:), intent(out), allocatable :: temp
	real(dp), dimension(:), allocatable :: delta
	logical :: same

	!sort the table's first column.
	call heapsort(list)

	m = dble(size(list,1))
	d= dble(size(list,2))
	p = 200_dp	!number of desired steps in "linear" regime.
	scaling = m**(-1_dp/(p*d))

	do j=1,size(list,2)
		maxelt(j)=findmax(list,j)
		minelt(j)=findmin(list,j)
		totalrange(j)=maxelt(j)-minelt(j)
	end do

	do k=1,size(list,2)

	!renormalize

	!$OMP PARALLEL SHARED(LIST,MINELT,K)
	!$OMP DO SCHEDULE(STATIC)
		do l=1,size(list,1)
			list(l,k)=list(l,k)-minelt(k)
			list(l,k)=list(l,k)/totalrange(k)
		end do
	!$OMP END DO
	!$OMP END PARALLEL
	end do

	epsilon_0=4_dp
	nmax=abs(ceiling(log(minim/epsilon_0)/log(scaling)))

	allocate(delta(nmax))
	do i_1=1,nmax	
	delta(i_1) = epsilon_0*(scaling**(i_1-1_dp))
	end do
	

	allocate(temp(nmax,2))


	print*,"box dimension	","   number of boxes"
	
	!$OMP PARALLEL DEFAULT(NONE)&
	!$OMP& SHARED(TOTALRANGE,TEMP,LIST,SCALING,NMAX,DELTA)&
	!$OMP& PRIVATE(CEILLIST,BOXNUMB)
	!$OMP DO SCHEDULE(STATIC)
do1:	do i=1,nmax
		if(delta(i)>1) then
			temp(i,1) = delta(i)
			temp(i,2) = 1
			print*,temp(i,1), temp(i,2)
			cycle
		end if
		
		
		allocate(ceillist(size(list,1),size(list,2)))
		ceillist = ceiling(list/delta(i))

				
	!!!!!very important for function uniq_elts.!!!!!	
	!sorts the rest of the table.
		if(size(list,2)>1)then
			call heapsorttotal(ceillist)
		end if

		!count all unique elements.  
		boxnumb = uniq_elts(ceillist)

		deallocate(ceillist)

		temp(i,1) = delta(i)
		temp(i,2) = dble(boxnumb)

		print*,temp(i,1), temp(i,2)
	end do do1
	!$OMP END DO
	!$OMP END PARALLEL

	deallocate(delta)

		
end subroutine boxcount_dp

!***************************************************************************************

!Performs the boxcount on an n x m dimensional list of n points (x_1,...,x_m).  It first sorts this data, then rescales it according to the grid size.  It then finds how many boxes it takes to cover the set by doing the CEILING count below and identifying the unique elements.

!Output is the array temp(x,2) of points (delta,N) according to the input scaling.

!NOTE: the list put into boxcount requires it to be sorted from low to high in first column.

subroutine bound_boxcount(success,fail,minim,temp)
	implicit none

	real(dp), dimension(:,:), intent(inout) :: success, fail
	integer, dimension(:,:), allocatable :: ceilsuccess, ceilfail, joint
	real(dp), intent(in) :: minim
	real(dp) :: scaling, m, p, d, epsilon_0
	real(dp), dimension(size(success,2)) :: maxelt, minelt, totalrange
	integer :: i,i_1,j,k,l,l_1,k_1, boxnumb
	integer :: nmax
	real(dp), dimension(:,:), intent(out), allocatable :: temp
	real(dp), dimension(:), allocatable :: delta
	logical :: same

	!sort the table's first column.
	call heapsort(success)
	call heapsort(fail)

	m = dble(size(success,1)+size(fail,1))
	d = dble(size(success,2)+size(fail,2))
	p = 200_dp	!number of desired steps in "linear" regime.
	scaling = m**(-1_dp/(p*d))

	do j=1,size(success,2)
		maxelt(j)=max(findmax(success,j),findmax(fail,j))
		minelt(j)=min(findmin(success,j),findmin(fail,j))
		totalrange(j)=maxelt(j)-minelt(j)
	end do

	!renormalize
	do l=1,size(success,1)
		do k=1,size(success,2)
			success(l,k)=success(l,k)-minelt(k)
			success(l,k)=success(l,k)/totalrange(k)
		end do
	end do
	do l_1=1,size(fail,1)
		do k_1=1,size(fail,2)
			fail(l_1,k_1)=fail(l_1,k_1)-minelt(k_1)
			fail(l_1,k_1)=fail(l_1,k_1)/totalrange(k_1)
		end do
	end do

	epsilon_0=4_dp
	nmax=abs(ceiling(log(minim/epsilon_0)/log(scaling)))

	allocate(delta(nmax))

	do i_1=1,nmax	
	delta(i_1) = epsilon_0*(scaling**(i_1-1_dp))
	end do
	

	allocate(temp(nmax,2))


	print*,"Box Dimension	","   Number of boxes"
	
	!$OMP PARALLEL DEFAULT(NONE)&
	!$OMP& SHARED(totalrange,temp,success,fail,scaling,nmax,delta)&
	!$OMP& PRIVATE(ceilsuccess,ceilfail,joint,boxnumb)
	!$OMP DO SCHEDULE(STATIC)
do1:	do i=1,nmax
		if(delta(i)>1) then
			temp(i,1) = delta(i)
			temp(i,2) = 1
			print*,temp(i,1), temp(i,2)
			cycle
		end if

		
		allocate(ceilsuccess(size(success,1),size(success,2)))
		allocate(ceilfail(size(fail,1),size(fail,2)))
		ceilsuccess = ceiling(success/delta(i))
		ceilfail = ceiling(fail/delta(i))

				
	!!!!!VERY IMPORTANT FOR FUNCTION UNIQ_ELTS.!!!!!	
	!Sorts the rest of the table.
		if(size(success,2)>1)then
			call heapsorttotal(ceilsuccess)
		end if
		if(size(fail,2)>1)then
			call heapsorttotal(ceilfail)
		end if

		call intersect(ceilsuccess,ceilfail,joint)

		if(size(joint,1)==0) then
			temp(i,1) = delta(i)
			temp(i,2) = 0_dp
			deallocate(ceilsuccess,ceilfail)
			cycle do1
		end if
		
		!count all unique elements.  
		boxnumb = uniq_elts(joint)

		deallocate(ceilsuccess,ceilfail)

		temp(i,1) = delta(i)
		temp(i,2) = dble(boxnumb)

		print*,temp(i,1), temp(i,2)
	end do do1
	!$omp end do
	!$omp end parallel

	deallocate(delta)

		
end subroutine bound_boxcount



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Count the number of unique elements of an ordered integer table.  

!Requires list to be sorted by HEAPSORTTOTAL from SORTERS module.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function uniq_elts(list)
	implicit none

	integer, dimension(:,:), intent(inout) :: list
	integer :: j
	logical :: same

	uniq_elts=1

	if(size(list,1)==1) return

	
do2:	do j=2,size(list,1)
		same = equalcond(size(list,2),j,j-1,list)
		if(same .eqv. .false.) then
			uniq_elts = uniq_elts+1
		end if
	end do do2

end function uniq_elts

!*********************************************************
!SUBROUTINE to find the intersection of two integer tables, table1 and table2.  Gives the result as joint, which will be the intersection of the two tables.

!NOTE: this requires that the table1 and table2 be totally sorted, ascending.  Else, could use technique as in boundary_finder.

subroutine intersect(table1,table2,joint)
	implicit none

	integer, dimension(:,:), intent(in) :: table1, table2
	integer, dimension(:,:), intent(out), allocatable :: joint
	logical, dimension(size(table1,1)) :: indexer
	integer :: i_1, i_2, j_1, j_2, j_3, counter, start
	logical :: same, ending

	indexer = .false.
	counter = 0
	start = 1
	same = .false.
	ending = .false.

doi1:	do i_1=1,size(table1,1)
	doj3:	do j_3=1,size(table1,2)
			if(i_1>1 .and. table1(i_1,j_3)==table1(i_1-1,j_3)) then
				same = .true.
			else
				same = .false.
				exit doj3
			end if
		end do doj3
	
		if(i_1>1 .and. same.eqv..true.) then
			indexer(i_1) = indexer(i_1-1)
			if(indexer(i_1).eqv..true.) then
				counter = counter + 1
			end if
			cycle
		end if

	doj2:	do j_2=start,size(table2,1)
		doj1:	do j_1=1,size(table2,2)
				if(table1(i_1,j_1)==table2(j_2,j_1)) then
					indexer(i_1) = .true.
				end if
				if(table1(i_1,j_1)>table2(j_2,j_1)) then
					indexer(i_1) = .false.
					start = j_2
					exit doj1
				end if
				if (table1(i_1,j_1)<table2(j_2,j_1)) then
					!start = j_2 -1
					indexer(i_1) = .false.
					ending = .true.
					exit doj1
				end if
			end do doj1
			if(indexer(i_1).eqv..true.) then
				counter = counter + 1
				!start = j_2
				exit doj2
			end if
			if(ending .eqv. .true.) then
				ending = .false.
				exit doj2
			end if
		end do doj2
	end do doi1

	allocate(joint(counter,size(table1,2)))

	if (counter == 0) return

	counter = 0
doi2:	do i_2=1,size(table1,1)
		if(indexer(i_2).eqv. .true.) then
			counter = counter + 1
			joint(counter,:)=table1(i_2,:)
		end if
	end do doi2


end subroutine intersect



!********************************************************
!Function to find minimum of a real NxN array in the ith column.

real(dp) function findmin(array,i)
	implicit none

	real(dp),dimension(:,:) :: array
	integer :: i
	real(dp), dimension(size(array,1)) :: col
	integer :: j

	do j=1,size(array,1)
		col(j)=array(j,i)
	end do

	findmin = minval(col)

end function findmin

!********************************************************
!Function to find maximum of a real NxN array in the ith column.

real(dp) function findmax(array,i)
	implicit none

	real(dp),dimension(:,:) :: array
	integer :: i
	real(dp), dimension(size(array,1)) :: col
	integer :: j

	do j=1,size(array,1)
		col(j)=array(j,i)
	end do

	findmax = maxval(col)

end function findmax


!********************************************************
!Subroutine to implement the correlation dimension technique to find the fractal dimension of an object.  This involves finding the pointwise mass function and averaging to achieve the correlation integral.  Then the correlation dimension is found by taking the limit as r -> 0 in logC/logr.  

!This subroutine will output the table temp(N,2) which catalogues the logC and logr in the columns.

!This procedure is taken from Theiler public.lanl.gov/jt/Papers/est-fractal-dim.pdf





end module box_count


