module boundaryboxcount
  use types, only : dp
  use sorters

contains

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
		maxelt(j)=max(findmax_d(success,j),findmax_d(fail,j))
		minelt(j)=min(findmin_d(success,j),findmin_d(fail,j))
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

!print*,"counter",counter


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

real(dp) function findmin_d(array,i)
implicit none

	real(dp),dimension(:,:) :: array
	integer :: i
	real(dp), dimension(size(array,1)) :: col
	integer :: j

	do j=1,size(array,1)
		col(j)=array(j,i)
	end do

	findmin_d = minval(col)

end function findmin_d

!********************************************************
!Function to find maximum of a real NxN array in the ith column.

real(dp) function findmax_d(array,i)
implicit none

	real(dp),dimension(:,:) :: array
	integer :: i
	real(dp), dimension(size(array,1)) :: col
	integer :: j

	do j=1,size(array,1)
		col(j)=array(j,i)
	end do

	findmax_d = maxval(col)

end function findmax_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Sets random_seed() according to computer clock to give good random numbers.


subroutine init_random_seed()
            integer :: i, n, clock
            integer, dimension(:), allocatable :: seed

          
            call random_seed(size = n)
            allocate(seed(n))
          
            call system_clock(count=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call random_seed(put = seed)
          
            deallocate(seed)
end subroutine init_random_seed


end module boundaryboxcount

