MODULE box_count
USE sorters

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




CONTAINS

!This subroutine calculates the fractal dimension after box counting has been performed.  It selects the points in the linear regime by finding the max and min values of log(N) and determining the "midpoint" of the dataset with respect to the log(N) direction.  It then systematically throws away datapoints from the end, plots it, and gives out the chi**2 value, until it's satisfied.

SUBROUTINE fractal_dimn(counted,dimn,dev,chi2)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: counted   !Will be n x 2 array.
	DOUBLE PRECISION, INTENT(OUT) :: dimn, dev, chi2
	DOUBLE PRECISION :: sad
	DOUBLE PRECISION, DIMENSION(SIZE(counted,1)) :: log_N, log_d
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: linear
	DOUBLE PRECISION :: a, b, maxim, minim, delta, low, high, siga,sigb, averagey, averagex
	INTEGER :: i, j,k,l, n, size1, averagex_numb, dropR, dropL,  numbR, numbL, counter
	LOGICAL :: fit
	INTEGER :: iter, ifault
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: d
	INTEGER, DIMENSION(:), ALLOCATABLE :: inext

	n = SIZE(counted,1) !Also size of log_N, log_d
	
	!Calc the logs.
	DO i=1, n
		log_d(i) = -1D0*LOG(counted(i,1))
		log_N(i) = LOG(counted(i,2))
	END DO

	!Find max/min values.  
	maxim = log_N(n)
	minim = 0D0
	averagey = (maxim-minim)*.8D0 !Position at upper 80% of range.

doj:	DO j=1,n-1
		IF(log_N(j)<=averagey .AND. log_N(j+1)>averagey) THEN
			averagex = log_d(j+1)
			averagex_numb = j+1
			EXIT doj
		END IF
	END DO doj

	dropR=0
	dropL=0

	PRINT*, "TYPE		", "chi**2		", "DIMN		", "DEV"

do1:	DO l=0,n-5
		size1 = n - dropR - dropL
		ALLOCATE(linear(size1,2))
		
		DO k=1,size1
			linear(k,1)=log_d(k+dropL)
			linear(k,2)=log_N(k+dropL)
		END DO

		!Fit least squares to linear area.
		CALL leastsqr(linear,a,b,siga,sigb,chi2)
		IF(chi2<20) PRINT*, "LEASTSQR: ",chi2, b, sigb

		IF(chi2 < .01) THEN
			fit = .TRUE.
			EXIT do1
		ELSE
			fit = .FALSE.
		END IF
		DEALLOCATE(linear)

		!Determine how many pts to Left and Right.
		numbR = n -averagex_numb - dropR
		numbL = averagex_numb - dropL -1
		!Drop the outlying points until get good chi**2 linear fit.
		IF(numbR<=numbL) THEN  
			dropL = dropL+1
		ELSE
			dropR = dropR+1
		END IF

	END DO do1

	PRINT*,"Number of right drops ",dropR
	PRINT*,"Number of left drops ",dropL
	PRINT*,"Number of points used ",(n-dropR - dropL)

	IF(fit .EQV. .FALSE.) THEN
		PRINT*,"ERROR: Couldn't find a fit within specified chi**2."
	END IF


	!Fractal dimension = slope (b).
	dimn = b
	dev = sigb

END SUBROUTINE fractal_dimn

!***********************************************************************************
!Subroutine to fit a line (alpha + beta*x + error) by minimizing absolute deviation.  From http://jblevins.org/mirror/amiller/as132.f90 .

!Input are the arrays x,y of length n.  Output is sad~sum of absolute deviations (like chi**2), alpha, beta.

SUBROUTINE simlp(points, sad, alpha, beta, d, iter, inext, ifault)
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-19  Time: 09:37:54

!     ALGORITHM AS 132  APPL. STATIST. (1978) VOL.27, NO.3

!     SIMLP:   Fit  Y = ALPHA + BETA.X + error

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: points
INTEGER :: n
DOUBLE PRECISION, DIMENSION(SIZE(points,1)) :: x, y
!INTEGER, INTENT(IN)   :: n
!REAL, INTENT(IN OUT)  :: x(n)
!REAL, INTENT(IN)      :: y(n)
DOUBLE PRECISION, INTENT(OUT)     :: sad
DOUBLE PRECISION, INTENT(OUT)     :: alpha
DOUBLE PRECISION, INTENT(OUT)     :: beta
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT), ALLOCATABLE :: d
INTEGER, INTENT(OUT)  :: iter
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: inext
INTEGER, INTENT(OUT)  :: ifault

DOUBLE PRECISION, PARAMETER  :: acu = 1.0D-06, big = 1.0D19, half = 0.5D0, zero = 0.0D0,   &
                    one = 1.0D0, two = 2.0D0
INTEGER  :: i, ibas1, ibas2, iflag, iin, iout, isave, j, i_1
DOUBLE PRECISION     :: a1, a2, aaa, aaaa, ahalf, aone, bbb, bbbb, ddd, det,   &
            ratio, rho, subt, sum, t, test, tot1, tot2, y1, y2, zzz

!Change the table to two vectors of length n.
n = SIZE(points,1)
ALLOCATE(d(n))
ALLOCATE(inext(n))
DO i_1=1,n
	x(i_1)=points(i_1,1)
	y(i_1)=points(i_1,2)
END DO



!     Initial settings

ifault = 0
iter = 0
ahalf = half + acu
aone = ahalf + ahalf

!     Determine initial basis

d(1) = zero
y1 = y(1)
ibas1 = 1
a1 = x(1)
DO  i = 2, n
  IF (ABS(a1-x(i)) >= acu) THEN
    a2 = x(i)
    ibas2 = i
    y2 = y(i)
    GO TO 20
  END IF
END DO
ifault = 1
RETURN

!     Calculate initial beta value

20 det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det

!     Calculate initial D-vector

DO  i = 2, n
  ddd = y(i) - (aaaa+bbbb*x(i))
  d(i) = SIGN(one,ddd)
END DO
tot1 = one
tot2 = x(ibas2)
d(ibas2) = -one
DO  i = 2, n
  tot1 = tot1 + d(i)
  tot2 = tot2 + d(i) * x(i)
END DO
t = (a2*tot1-tot2) * det
IF (ABS(t) >= aone) THEN
  det = -det
  GO TO 70
END IF

!     Main iterative loop begins

50 t = (tot2-a1*tot1) * det
IF (ABS(t) < aone) GO TO 130
iflag = 2
iout = ibas2
x(iout) = a1
aaa = a1
bbb = a2
GO TO 80

60 t = (tot2-a2*tot1) * det
IF (ABS(t) < aone) GO TO 130

70 iflag = 1
bbb = a1
aaa = a2
iout = ibas1

80 rho = SIGN(one,t)
t = half * ABS(t)
det = det * rho

!     Perform partial sort of ratios

inext(ibas1) = ibas2
ratio = big
sum = ahalf
DO  i = 1, n
  ddd = (x(i)-aaa) * det
  IF (ddd*d(i) > acu) THEN
    test = (y(i)-aaaa-bbbb*x(i)) / ddd
    IF (test < ratio) THEN
      j = ibas1
      sum = sum + ABS(ddd)
      90 isave = ABS(inext(j))
      IF (test < d(isave)) THEN
        IF (sum >= t) THEN
          subt = ABS((x(isave)-aaa)*det)
          IF (sum-subt >= t) THEN
            sum = sum - subt
            d(isave) = SIGN(1,inext(j))
            inext(j) = inext(isave)
            GO TO 90
          END IF
        END IF
        100 j = isave
        isave = ABS(inext(j))
        IF (test < d(isave)) GO TO 100
      END IF
      inext(i) = inext(j)
      inext(j) = SIGN(i,INT(d(i)))
      d(i) = test
      IF (sum >= t) THEN
        iin = ABS(inext(ibas1))
        ratio = d(iin)
      END IF
    END IF
  END IF
END DO

!     Update basic indicators

iin = ABS(inext(ibas1))
j = iin
120 isave = ABS(inext(j))
IF (isave /= ibas2) THEN
  zzz = SIGN(1,inext(j))
  tot1 = tot1 - zzz - zzz
  tot2 = tot2 - two * zzz * x(isave)
  d(isave) = -zzz
  j = isave
  GO TO 120
END IF
zzz = SIGN(1,inext(ibas1))
tot1 = tot1 - rho - zzz
tot2 = tot2 - rho * bbb - zzz * x(iin)
d(iout) = -rho
iter = iter + 1
IF (iflag /= 1) THEN
  x(ibas2) = a2
  ibas2 = iin
  d(ibas2) = -one
  a2 = x(iin)
  y2 = y(iin)
  det = one / (a1-a2)
  aaaa = (a1*y2-a2*y1) * det
  bbbb = (y1-y2) * det
  GO TO 60
END IF
ibas1 = iin
a1 = x(iin)
d(ibas1) = zero
y1 = y(iin)
det = one / (a2-a1)
aaaa = (a2*y1-a1*y2) * det
bbbb = (y2-y1) * det
GO TO 50

!     Calculate optimal sum of absolute deviations

130 sad = zero
DO  i = 1, n
  d(i) = y(i) - (aaaa+bbbb*x(i))
  sad = sad + ABS(d(i))
END DO
alpha = aaaa
beta = bbbb

RETURN
END SUBROUTINE simlp


!****************************************************************************

!Linear least-squares regression on the linear regime.  Inputs: points(z,2)~the (x,y) points to plot starting with z=1--needs to be constrained to linear regime prior to call to this subroutine.  Outputs: a,b~(a+b*x) and sig's~main standard deviation.  From Numerical Recipes with unavailable standard deviations.

SUBROUTINE leastsqr(points,a,b,siga,sigb,chi2)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: points  !An n x 2 array of n 2-D points (x,y).
	DOUBLE PRECISION, INTENT(OUT) :: a, b, siga, sigb, chi2
	DOUBLE PRECISION :: sx, sy, st2, sxoss, t, sigdat, q, ss
	INTEGER :: i, j, k, n
	
	sx = 0D0
	sy = 0D0
	st2 = 0D0
	b = 0D0

	n = SIZE(points,1)

	DO i=1,n
		sx = sx + points(i,1)
		sy = sy + points(i,2)
	END DO
	
	ss = DBLE(n)
	sxoss = sx/ss

	DO j=1,n
		t = points(j,1) - sxoss
		st2 = st2 +(t*t)
		b = b + (t*points(j,2))
	END DO

	!Calculate a, b, and probable uncertainties siga, sigb.

	b=b/st2
	a=(sy-(sx*b))/ss
	siga = SQRT((1D0 + (sx*sx)/(ss*st2))/ss)
	sigb = SQRT(1D0/st2)

	chi2 = 0D0

	DO k=1,n
		chi2 = chi2 + (points(k,2)-a-(b*points(k,1)))**2D0
	END DO
	q=1
	sigdat = SQRT(chi2/(n-2))
	siga = siga*sigdat
	sigb = sigb*sigdat

END SUBROUTINE leastsqr


!***************************************************************************************

!Performs the boxcount on an n x m dimensional list of n points (x_1,...,x_m).  It first sorts this data, then rescales it according to the grid size.  It then finds how many boxes it takes to cover the set by doing the CEILING count below and identifying the unique elements.

!Output is the array temp(x,2) of points (delta,N) according to the input scaling.

!NOTE: the list put into boxcount requires it to be sorted from low to high in first column.

SUBROUTINE boxcount(list,minim,temp)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: list
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: ceillist
	DOUBLE PRECISION, INTENT(IN) :: minim
	DOUBLE PRECISION :: scaling, m, p, d, epsilon_0
	DOUBLE PRECISION, DIMENSION(SIZE(list,2)) :: maxelt, minelt, totalrange
	INTEGER :: i,i_1,j,k,l, boxnumb
	INTEGER :: nmax
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE :: temp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: delta
	LOGICAL :: same

	!Sort the table's first column.
	CALL heapsort(list)

	m = DBLE(SIZE(list,1))
	d= DBLE(SIZE(list,2))
	p = 200D0	!Number of desired steps in "linear" regime.
	scaling = m**(-1D0/(p*d))

	DO j=1,SIZE(list,2)
		maxelt(j)=findmax(list,j)
		minelt(j)=findmin(list,j)
		totalrange(j)=maxelt(j)-minelt(j)
	END DO

	DO k=1,SIZE(list,2)

	!Renormalize

	!$OMP PARALLEL SHARED(list,minelt,k)
	!$OMP DO SCHEDULE(STATIC)
		DO l=1,SIZE(list,1)
			list(l,k)=list(l,k)-minelt(k)
			list(l,k)=list(l,k)/totalrange(k)
		END DO
	!$OMP END DO
	!$OMP END PARALLEL
	END DO

	epsilon_0=4D0
	nmax=ABS(CEILING(LOG(minim/epsilon_0)/LOG(scaling)))

	ALLOCATE(delta(nmax))
	DO i_1=1,nmax	
	delta(i_1) = epsilon_0*(scaling**(i_1-1D0))
	END DO
	

	ALLOCATE(temp(nmax,2))


	PRINT*,"Box Dimension	","   Number of boxes"
	
	!$OMP PARALLEL DEFAULT(NONE)&
	!$OMP& SHARED(totalrange,temp,list,scaling,nmax,delta)&
	!$OMP& PRIVATE(ceillist,boxnumb)
	!$OMP DO SCHEDULE(STATIC)
do1:	DO i=1,nmax
		IF(delta(i)>1) THEN
			temp(i,1) = delta(i)
			temp(i,2) = 1
			PRINT*,temp(i,1), temp(i,2)
			CYCLE
		END IF
		
		
		ALLOCATE(ceillist(SIZE(list,1),SIZE(list,2)))
		ceillist = CEILING(list/delta(i))

				
	!!!!!VERY IMPORTANT FOR FUNCTION UNIQ_ELTS.!!!!!	
	!Sorts the rest of the table.
		IF(SIZE(list,2)>1)THEN
			CALL heapsorttotal(ceillist)
		END IF

		!Count all unique elements.  
		boxnumb = uniq_elts(ceillist)

		DEALLOCATE(ceillist)

		temp(i,1) = delta(i)
		temp(i,2) = DBLE(boxnumb)

		PRINT*,temp(i,1), temp(i,2)
	END DO do1
	!$OMP END DO
	!$OMP END PARALLEL

	DEALLOCATE(delta)

		
END SUBROUTINE boxcount



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Count the number of unique elements of an ordered integer table.  

!Requires list to be sorted by HEAPSORTTOTAL from SORTERS module.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER FUNCTION uniq_elts(list)
IMPLICIT NONE

	INTEGER, DIMENSION(:,:), INTENT(INOUT) :: list
	INTEGER :: j
	LOGICAL :: same

	uniq_elts=1

	IF(SIZE(list,1)==1) RETURN

	
do2:	DO j=2,SIZE(list,1)
		same = equalcond(SIZE(list,2),j,j-1,list)
		IF(same .EQV. .FALSE.) THEN
			uniq_elts = uniq_elts+1
		END IF
	END DO do2

END FUNCTION uniq_elts


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Sets random_seed() according to computer clock to give good random numbers.


SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE init_random_seed


!********************************************************
!Function to find minimum of a real NxN array in the ith column.

DOUBLE PRECISION FUNCTION findmin(array,i)
IMPLICIT NONE

	DOUBLE PRECISION,DIMENSION(:,:) :: array
	INTEGER :: i
	DOUBLE PRECISION, DIMENSION(SIZE(array,1)) :: col
	INTEGER :: j

	DO j=1,SIZE(array,1)
		col(j)=array(j,i)
	END DO

	findmin = MINVAL(col)

END FUNCTION findmin

!********************************************************
!Function to find maximum of a real NxN array in the ith column.

DOUBLE PRECISION FUNCTION findmax(array,i)
IMPLICIT NONE

	DOUBLE PRECISION,DIMENSION(:,:) :: array
	INTEGER :: i
	DOUBLE PRECISION, DIMENSION(SIZE(array,1)) :: col
	INTEGER :: j

	DO j=1,SIZE(array,1)
		col(j)=array(j,i)
	END DO

	findmax = MAXVAL(col)

END FUNCTION findmax


!********************************************************
!Subroutine to implement the correlation dimension technique to find the fractal dimension of an object.  This involves finding the pointwise mass function and averaging to achieve the correlation integral.  Then the correlation dimension is found by taking the limit as r -> 0 in logC/logr.  

!This subroutine will output the table temp(N,2) which catalogues the logC and logr in the columns.

!This procedure is taken from Theiler public.lanl.gov/jt/Papers/est-fractal-dim.pdf





END MODULE box_count


