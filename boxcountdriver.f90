PROGRAM boxcountdriver
USE boundaryboxcount
USE box_count
USE sorters
IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: success
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fail
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: temp, total
	DOUBLE PRECISION :: epsilon_0, minim, dimn, dev, chi2
	INTEGER :: length_s, length_f, width_s, width_f, i, k,j, check
	INTEGER :: countboth
	NAMELIST /tablel/ length_s, length_f, width_s, width_f, countboth

	!Reads file sizes from input file "bxfilesizes.txt".
	OPEN(unit=1000, file="bxfilesizes.txt", status="old", delim="apostrophe")
	READ(unit=1000, nml=tablel)
	ALLOCATE(success(length_s,width_s),fail(length_f,width_f))

	!Read succ and fail sets from file.
	OPEN(unit=20,status='old',file='totalsucc.bin',form='UNFORMATTED')
	OPEN(unit=21,status='old',file='totalfail.bin',form='UNFORMATTED')

	check = 0
doi1:	DO i=1,length_s+1
		check = check + 1
		READ(20,END=10) (success(i,j), j=1,width_s)
IF (mod(check,1000000)==0) THEN
        PRINT*,	"SUCCESS",check
	PRINT*, (success(i,j),j=1,width_s)
END IF
		IF(check==SIZE(success,1)) EXIT doi1
	END DO	doi1

10	CLOSE(unit=20)

	check = 0
doi2:	DO i=1,length_f+1
		check = check + 1
		READ(21,END=20) (fail(i,j),j=1,width_f)
IF (mod(check,1000000)==0) THEN
        PRINT*, "FAIL",check
        PRINT*,	(fail(i,j),j=1,width_f)
END IF
		IF(check==SIZE(fail,1)) EXIT doi2
	END DO	doi2

	!Start boxcounting.
20	minim = .0000001D0
	CLOSE(unit=21)


	PRINT*,"Boxcounting success..."
	OPEN(unit=80, file="bxsucc.80", status="new")
	CALL boxcount(success,minim,temp)
	DO i=1,SIZE(temp,1)
		WRITE(UNIT=80,FMT=*),temp(i,1),temp(i,2)
	END DO
	DEALLOCATE(temp)


	PRINT*,"Boxcounting fail..."
	OPEN(unit=81, file="bxfail.81", status="new")
	CALL boxcount(fail,minim,temp)
	DO i=1,SIZE(temp,1)
		WRITE(UNIT=81,FMT=*),temp(i,1),temp(i,2)
	END DO
	DEALLOCATE(temp)

	PRINT*,"Boxcounting boundary..."
	OPEN(unit=82, file="bxbound.82", status="new")
	CALL bound_boxcount(success,fail,minim,temp)
	DO i=1,SIZE(temp,1)
		WRITE(UNIT=82,FMT=*),temp(i,1),temp(i,2)
	END DO
	DEALLOCATE(temp)

	!Boxcount both if desired.
	IF (countboth==1) THEN
		ALLOCATE(total(length_s+length_f,width_s))

		DO i=1,(length_s+length_f)
			IF(k<=length_s-1) THEN
				DO j=1,width_s
					total(i,j) = success(i,j)
				END DO
			ELSE
				DO j=1, width_f
					total(i,j) = fail(i-length_s,j)
				END DO
			END IF
		END DO
	
		DEALLOCATE(success,fail)

		PRINT*,"Boxcounting both..."
		OPEN(unit=83, file="bxboth.83", status="new")
		CALL boxcount(total,minim,temp)
		DO i=1,SIZE(temp,1)
			WRITE(UNIT=83,FMT=*),temp(i,1),temp(i,2)
		END DO
		DEALLOCATE(temp)
		
		DEALLOCATE(total)
	ELSE
		DEALLOCATE(success,fail)
	END IF		

END PROGRAM boxcountdriver
