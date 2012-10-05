!*******************************************************************************
!Layne Price, University of Auckland, May 2012.
!*******************************************************************************

!A program that loads a success file named totalsucc.bin and a fail file named
!totalfail.bin and performs a boxcounting procedure on a) the success set b) the
!failure set c) the boundary between the two and d) the two combined together.
!
!Reads the array file sizes from a namelist stored in bxfilesizes.txt.

!Optionally takes in one set and box counts it at different resolutions.  In
!other words, it loads a set and then box counts only portions of that set.
!These portions steadily increase until the full resolution is counted.  This
!allows us to determine at which resolution the fractal dimension stabilizes.


program boxcountdriver
  use box_count
  use sorters, only : heapsort
  use types, only : dp
  use features, only : newunit
  implicit none

	real(dp), dimension(:,:), allocatable :: success, subset
	real(dp), dimension(:,:), allocatable :: fail
	real(dp), dimension(:,:), allocatable :: temp, total
	real(dp) :: epsilon_0, minim, dimn, dev, chi2
	integer :: length_s, length_f, width_s, width_f, i, k,j, check
	integer :: ier, u, uu, high, iend
  logical :: countboth, resol_count

	namelist /tablel/ length_s, length_f, width_s, width_f, countboth, &
  & resol_count, minim

	!Reads file sizes from input file "bxfilesizes.txt".
	open(unit=newunit(u), file="bxfilesizes.txt", status="old", delim="apostrophe")
	read(unit=u, nml=tablel)
  close(u)
	allocate(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	open(unit=newunit(u),status='old',file='totalsucc.bin',form='unformatted')
	check = 0
doi1:	do i=1,length_s+1
		check = check + 1
		read(u,iostat=ier) (success(i,j), j=1,width_s)
		if(check==size(success,1)) exit doi1
    if (ier<0) exit doi1
	end do	doi1
  close(u)

	open(unit=newunit(u),status='old',file='totalfail.bin',form='unformatted')
	check = 0
doi2:	do i=1,length_f+1
		check = check + 1
		read(u,iostat=ier) (fail(i,j),j=1,width_f)
		if(check==size(fail,1)) exit doi2
    if (ier<0) exit doi2
	end do	doi2
  close(u)

	!Start boxcounting.
  !-------------------------------------------------------------

  if (resol_count) then
    iend=10
    do i=1,iend
      !Count only a subset of the success set.
      high=(size(success,1)/iend)*i  !Note int div.
      print*,"Boxcounting resolution points ", high
      allocate(subset(high,size(success,2)))
      subset(:,:)=success(1:high,:)
      print*,subset(1,1)
      call boxcount(subset,minim,temp)
      !Unit to write to.
      uu=i+90
      do j=1,size(temp,1)
        write(unit=uu,fmt=*), temp(j,1), temp(j,2)
      end do
      deallocate(temp)
      deallocate(subset)
    end do
    stop
  end if


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
	call boxcount(success,fail,minim,temp)
	do i=1,size(temp,1)
		write(unit=82,fmt=*),temp(i,1),temp(i,2)
	end do
	deallocate(temp)

	!boxcount both if desired.
	if (countboth) then
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
