MODULE boundaryboxcount
USE sorters

CONTAINS

!***************************************************************************************

!Performs the boxcount on an n x m dimensional list of n points (x_1,...,x_m).  It first sorts this data, then rescales it according to the grid size.  It then finds how many boxes it takes to cover the set by doing the CEILING count below and identifying the unique elements.

!Output is the array temp(x,2) of points (delta,N) according to the input scaling.

!NOTE: the list put into boxcount requires it to be sorted from low to high in first column.

SUBROUTINE bound_boxcount(success,fail,minim,temp)
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: success, fail
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: ceilsuccess, ceilfail, joint
	DOUBLE PRECISION, INTENT(IN) :: minim
	DOUBLE PRECISION :: scaling, m, p, d, epsilon_0
	DOUBLE PRECISION, DIMENSION(SIZE(success,2)) :: maxelt, minelt, totalrange
	INTEGER :: i,i_1,j,k,l,l_1,k_1, boxnumb
	INTEGER :: nmax
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE :: temp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: delta
	LOGICAL :: same

	!Sort the table's first column.
	CALL heapsort(success)
	CALL heapsort(fail)

	m = DBLE(SIZE(success,1)+SIZE(fail,1))
	d = DBLE(SIZE(success,2)+SIZE(fail,2))
	p = 200D0	!Number of desired steps in "linear" regime.
	scaling = m**(-1D0/(p*d))

	DO j=1,SIZE(success,2)
		maxelt(j)=MAX(findmax_d(success,j),findmax_d(fail,j))
		minelt(j)=MIN(findmin_d(success,j),findmin_d(fail,j))
		totalrange(j)=maxelt(j)-minelt(j)
	END DO

	!Renormalize
	DO l=1,SIZE(success,1)
		DO k=1,SIZE(success,2)
			success(l,k)=success(l,k)-minelt(k)
			success(l,k)=success(l,k)/totalrange(k)
		END DO
	END DO
	DO l_1=1,SIZE(fail,1)
		DO k_1=1,SIZE(fail,2)
			fail(l_1,k_1)=fail(l_1,k_1)-minelt(k_1)
			fail(l_1,k_1)=fail(l_1,k_1)/totalrange(k_1)
		END DO
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
	!$OMP& SHARED(totalrange,temp,success,fail,scaling,nmax,delta)&
	!$OMP& PRIVATE(ceilsuccess,ceilfail,joint,boxnumb)
	!$OMP DO SCHEDULE(GUIDED)
do1:	DO i=1,nmax
		IF(delta(i)>1) THEN
			temp(i,1) = delta(i)
			temp(i,2) = 1
			PRINT*,temp(i,1), temp(i,2)
			CYCLE
		END IF

		
		ALLOCATE(ceilsuccess(SIZE(success,1),SIZE(success,2)))
		ALLOCATE(ceilfail(SIZE(fail,1),SIZE(fail,2)))
		ceilsuccess = CEILING(success/delta(i))
		ceilfail = CEILING(fail/delta(i))

				
	!!!!!VERY IMPORTANT FOR FUNCTION UNIQ_ELTS.!!!!!	
	!Sorts the rest of the table.
		IF(SIZE(success,2)>1)THEN
			CALL heapsorttotal(ceilsuccess)
		END IF
		IF(SIZE(fail,2)>1)THEN
			CALL heapsorttotal(ceilfail)
		END IF

		CALL INTERSECT(ceilsuccess,ceilfail,joint)

		IF(SIZE(joint,1)==0) THEN
			temp(i,1) = delta(i)
			temp(i,2) = 0D0
			DEALLOCATE(ceilsuccess,ceilfail)
			CYCLE do1
		END IF
		
		!Count all unique elements.  
		boxnumb = uniq_elts(joint)

		DEALLOCATE(ceilsuccess,ceilfail)

		temp(i,1) = delta(i)
		temp(i,2) = DBLE(boxnumb)

		PRINT*,temp(i,1), temp(i,2)
	END DO do1
	!$OMP END DO
	!$OMP END PARALLEL

	DEALLOCATE(delta)

		
END SUBROUTINE bound_boxcount



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

!*********************************************************
!SUBROUTINE to find the intersection of two integer tables, table1 and table2.  Gives the result as joint, which will be the intersection of the two tables.

!NOTE: this requires that the table1 and table2 be totally sorted, ascending.  Else, could use technique as in boundary_finder.

SUBROUTINE INTERSECT(table1,table2,joint)
IMPLICIT NONE

	INTEGER, DIMENSION(:,:), INTENT(IN) :: table1, table2
	INTEGER, DIMENSION(:,:), INTENT(OUT), ALLOCATABLE :: joint
	LOGICAL, DIMENSION(SIZE(table1,1)) :: indexer
	INTEGER :: i_1, i_2, j_1, j_2, j_3, counter, start
	LOGICAL :: same, ending

	indexer = .FALSE.
	counter = 0
	start = 1
	same = .FALSE.
	ending = .FALSE.

doi1:	DO i_1=1,SIZE(table1,1)
	doj3:	DO j_3=1,SIZE(table1,2)
			IF(i_1>1 .AND. table1(i_1,j_3)==table1(i_1-1,j_3)) THEN
				same = .TRUE.
			ELSE
				same = .FALSE.
				EXIT doj3
			END IF
		END DO doj3
	
		IF(i_1>1 .AND. same.EQV..TRUE.) THEN
			indexer(i_1) = indexer(i_1-1)
			IF(indexer(i_1).EQV..TRUE.) THEN
				counter = counter + 1
			END IF
			CYCLE
		END IF

	doj2:	DO j_2=start,SIZE(table2,1)
		doj1:	DO j_1=1,SIZE(table2,2)
				IF(table1(i_1,j_1)==table2(j_2,j_1)) THEN
					indexer(i_1) = .TRUE.
				END IF
				IF(table1(i_1,j_1)>table2(j_2,j_1)) THEN
					indexer(i_1) = .FALSE.
					start = j_2
					EXIT doj1
				END IF
				IF (table1(i_1,j_1)<table2(j_2,j_1)) THEN
					!start = j_2 -1
					indexer(i_1) = .FALSE.
					ending = .TRUE.
					EXIT doj1
				END IF
			END DO doj1
			IF(indexer(i_1).EQV..TRUE.) THEN
				counter = counter + 1
				!start = j_2
				EXIT doj2
			END IF
			IF(ending .EQV. .TRUE.) THEN
				ending = .FALSE.
				EXIT doj2
			END IF
		END DO doj2
	END DO doi1

!PRINT*,"COUNTER",counter


	ALLOCATE(joint(counter,SIZE(table1,2)))

	IF (counter == 0) RETURN

	counter = 0
doi2:	DO i_2=1,SIZE(table1,1)
		IF(indexer(i_2).EQV. .TRUE.) THEN
			counter = counter + 1
			joint(counter,:)=table1(i_2,:)
		END IF
	END DO doi2


END SUBROUTINE INTERSECT

!********************************************************
!Function to find minimum of a real NxN array in the ith column.

DOUBLE PRECISION FUNCTION findmin_d(array,i)
IMPLICIT NONE

	DOUBLE PRECISION,DIMENSION(:,:) :: array
	INTEGER :: i
	DOUBLE PRECISION, DIMENSION(SIZE(array,1)) :: col
	INTEGER :: j

	DO j=1,SIZE(array,1)
		col(j)=array(j,i)
	END DO

	findmin_d = MINVAL(col)

END FUNCTION findmin_d

!********************************************************
!Function to find maximum of a real NxN array in the ith column.

DOUBLE PRECISION FUNCTION findmax_d(array,i)
IMPLICIT NONE

	DOUBLE PRECISION,DIMENSION(:,:) :: array
	INTEGER :: i
	DOUBLE PRECISION, DIMENSION(SIZE(array,1)) :: col
	INTEGER :: j

	DO j=1,SIZE(array,1)
		col(j)=array(j,i)
	END DO

	findmax_d = MAXVAL(col)

END FUNCTION findmax_d

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


END MODULE boundaryboxcount

