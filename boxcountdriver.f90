program boxcountdriver
use boundaryboxcount
use box_count
use sorters
use types, only : dp
implicit none

	real(dp), dimension(:,:), allocatable :: success
	real(dp), dimension(:,:), allocatable :: fail
	real(dp), dimension(:,:), allocatable :: temp, total
	real(dp) :: epsilon_0, minim, dimn, dev, chi2
	integer :: length_s, length_f, width_s, width_f, i, k,j, check
	integer :: countboth
	namelist /tablel/ length_s, length_f, width_s, width_f, countboth

	!Reads file sizes from input file "bxfilesizes.txt".
	open(unit=1000, file="bxfilesizes.txt", status="old", delim="apostrophe")
	read(unit=1000, nml=tablel)
	allocate(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	open(unit=20,status='old',file='totalsucc.bin',form='unformatted')
	open(unit=21,status='old',file='totalfail.bin',form='unformatted')

	check = 0
doi1:	do i=1,length_s+1
		check = check + 1
		read(20,end=10) (success(i,j), j=1,width_s)
		if(check==size(success,1)) exit doi1
	end do	doi1

10	close(unit=20)

	check = 0
doi2:	do i=1,length_f+1
		check = check + 1
		read(21,end=20) (fail(i,j),j=1,width_f)
		if(check==size(fail,1)) exit doi2
	end do	doi2

	!Start boxcounting.
20	minim = .0000001_dp
	close(unit=21)


	print*,"Boxcounting success..."
	open(unit=80, file="bxsucc.80", status="new")
	call boxcount(success,minim,temp)
	do i=1,size(temp,1)
		write(unit=80,fmt=*),temp(i,1),temp(i,2)
	end do
	deallocate(temp)


	print*,"Boxcounting fail..."
	open(unit=81, file="bxfail.81", status="new")
	call boxcount(fail,minim,temp)
	do i=1,size(temp,1)
		write(unit=81,fmt=*),temp(i,1),temp(i,2)
	end do
	deallocate(temp)

	print*,"Boxcounting boundary..."
	open(unit=82, file="bxbound.82", status="new")
	call bound_boxcount(success,fail,minim,temp)
	do i=1,size(temp,1)
		write(unit=82,fmt=*),temp(i,1),temp(i,2)
	end do
	deallocate(temp)

	!boxcount both if desired.
	if (countboth==1) then
		allocate(total(length_s+length_f,width_s))

		do i=1,(length_s+length_f)
			if(k<=length_s-1) then
				do j=1,width_s
					total(i,j) = success(i,j)
				end do
			else
				do j=1, width_f
					total(i,j) = fail(i-length_s,j)
				end do
			end if
		end do
	
		deallocate(success,fail)

		print*,"Boxcounting both..."
		open(unit=83, file="bxboth.83", status="new")
		call boxcount(total,minim,temp)
		do i=1,size(temp,1)
			write(unit=83,fmt=*),temp(i,1),temp(i,2)
		end do
		deallocate(temp)
		
		deallocate(total)
	else
		deallocate(success,fail)
	end if		

end program boxcountdriver
