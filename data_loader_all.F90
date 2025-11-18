module data_loader_all
    use mpi
    use array_dimensions
    implicit none
    private
    public :: load_checkpoint_binary, load_weightStuff_binary, load_ETDRKcoeffs_binary_named, load_mat_files

contains
!======================================================================
! 1. CHECKPOINT
!======================================================================
subroutine load_checkpoint_binary(bf, v_1,v_2,v_3,v_4,v_5,v_6,t,tStar,counter,kLength,comm)
    real(dp), intent(in) :: bf
    integer,  intent(in) :: comm
    integer,  intent(out) :: kLength, counter
    real(dp), intent(out) :: v_1(:),v_2(:),v_3(:),v_4(:),v_5(:),v_6(:)
    real(dp), intent(out) :: t,tStar
    integer :: rank,ierr,ios,kLength_file
    real(dp) :: counter_dbl
    character(len=256) :: fname,dirname
    character(len=64) :: bfstr

    call MPI_Comm_rank(comm,rank,ierr)

    !---- build folder name
    write(bfstr,'(G0.6)') bf
    bfstr  = adjustl(bfstr)
    write(dirname,'(A,A)') 'results/param_', trim(bfstr)
    fname = trim(dirname)//'/checkpoint.bin'

    if (rank==0) then
        print *, fname
        open(10,file=fname,form='unformatted',access='stream',status='old',iostat=ios)
        if (ios/=0) stop 'Error opening checkpoint file.'

        ! First thing in file = kLength - verify it matches compiled value
        read(10) kLength_file
        if (kLength_file /= KLENGTH_PARAM) then
            print *, 'ERROR: File kLength=', kLength_file, ' does not match compiled KLENGTH_PARAM=', KLENGTH_PARAM
            stop 'kLength mismatch - recompile with correct kLength'
        end if
        kLength = KLENGTH_PARAM
    end if

    call MPI_Bcast(kLength,1,MPI_INTEGER,0,comm,ierr)

    ! Arrays are fixed-size, no allocation needed

    if (rank==0) then
        ! Read exactly the expected number of elements
        read(10) v_1
        read(10) v_2
        read(10) v_3
        read(10) v_4
        read(10) v_5
        read(10) v_6
        read(10) t
        read(10) tStar
        read(10) counter_dbl
        counter = int(counter_dbl)
        close(10)
    end if

    ! Broadcast to all ranks
    call MPI_Bcast(v_1,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(v_2,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(v_3,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(v_4,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(v_5,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(v_6,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(t,       1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    print *, 't=',t
    call MPI_Bcast(tStar,   1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(counter, 1,MPI_INTEGER,        0,comm,ierr)
print *, 'load checkpoint binary v_1 = ',v_1
end subroutine load_checkpoint_binary
!======================================================================
! 2. WEIGHTSTUFF
!======================================================================
subroutine load_weightStuff_binary(bf,kLength,weight,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11,comm)
    real(dp),intent(in)::bf
    integer,intent(in)::comm
    integer,intent(out)::kLength
    real(dp),intent(out)::weight(:,:,:),CxVals(:,:,:),CyVals(:,:,:)
    integer,intent(out)::triadFlag(:,:,:),outsideCutCell(:,:,:),insideCutCell(:,:,:),Q11(:,:,:,:)
    integer :: rank,ierr,ios,kLength_file
    character(len=64)::fname,dirname
    character(len=64) :: bfstr

    call MPI_Comm_rank(comm,rank,ierr)
    if(rank==0)then
       write(bfstr,'(G0.6)') bf
bfstr  = adjustl(bfstr)
write(dirname,'(A,A)') 'results/param_', trim(bfstr)

fname = trim(dirname)//'/weightStuff.bin'
        open(10,file=fname,form='unformatted',access='stream',status='old',iostat=ios)
        if (ios /= 0) then
    print *, "Rank ", rank, " Error opening weightStuff file."
    print *, "Tried to open: ", fname
    print *, "iostat = ", ios
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    stop
end if
        read(10)kLength_file
        if (kLength_file /= KLENGTH_PARAM) then
            print *, 'ERROR: File kLength=', kLength_file, ' does not match compiled KLENGTH_PARAM=', KLENGTH_PARAM
            stop 'kLength mismatch in weightStuff - recompile with correct kLength'
        end if
        kLength = KLENGTH_PARAM
     end if
    call MPI_Bcast(kLength,1,MPI_INTEGER,0,comm,ierr)

    ! Arrays are fixed-size, no allocation needed

    if(rank==0)then
        read(10) weight
        read(10) triadFlag
        read(10) outsideCutCell
        read(10) insideCutCell
        read(10) CxVals
        read(10) CyVals
        read(10) Q11
        close(10)
    end if

    call MPI_Bcast(weight,size(weight),MPI_DOUBLE_PRECISION,0,comm,ierr)
    print *, 'weight(17,46,100) =',weight(17,46,100)
    call MPI_Bcast(triadFlag,size(triadFlag),MPI_INTEGER,0,comm,ierr)
    print *, 'triadFlag(17,46,100) =',triadFlag(17,46,100)
    call MPI_Bcast(outsideCutCell,size(outsideCutCell),MPI_INTEGER,0,comm,ierr)
    print *, 'outsideCutCell(17,46,100) =',outsideCutCell(17,46,100)
    call MPI_Bcast(insideCutCell,size(insideCutCell),MPI_INTEGER,0,comm,ierr)
    print *, 'insideCutCell(17,46,100) =',insideCutCell(17,46,100)
    call MPI_Bcast(CxVals,size(CxVals),MPI_DOUBLE_PRECISION,0,comm,ierr)
    print *, 'CxVals(17,46,100) =',CxVals(17,46,100)
    call MPI_Bcast(CyVals,size(CyVals),MPI_DOUBLE_PRECISION,0,comm,ierr)
    print *, 'CyVals(17,46,100) =',CyVals(17,46,100)
    call MPI_Bcast(Q11,size(Q11),MPI_INTEGER,0,comm,ierr)
    print *, 'Q11(1,17,46,100) =',Q11(1,17,46,100)
end subroutine load_weightStuff_binary
!======================================================================
! 3. ETDRKCOEFFS (explicit variables)
!======================================================================
subroutine load_ETDRKcoeffs_binary_named( bf, kLength,  &
    EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, &
    EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
    EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, &
    EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
    EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, &
    EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6, &
    comm )

    use mpi
    implicit none

    ! Inputs
    real(dp), intent(in) :: bf
    integer, intent(in)  :: comm

    ! Outputs
    integer, intent(out) :: kLength

    real(dp), intent(out) :: &
        EX1_1(:),EX2_1(:),Q_1(:),f1_1(:),f2_1(:),f3_1(:), &
        EX1_2(:),EX2_2(:),Q_2(:),f1_2(:),f2_2(:),f3_2(:), &
        EX1_3(:),EX2_3(:),Q_3(:),f1_3(:),f2_3(:),f3_3(:), &
        EX1_4(:),EX2_4(:),Q_4(:),f1_4(:),f2_4(:),f3_4(:), &
        EX1_5(:),EX2_5(:),Q_5(:),f1_5(:),f2_5(:),f3_5(:), &
        EX1_6(:),EX2_6(:),Q_6(:),f1_6(:),f2_6(:),f3_6(:)

    ! Locals
    integer :: rank, ierr, ios, kLength_file
    character(len=64)  :: bfstr
    character(len=256) :: dirname, fname

    call MPI_Comm_rank(comm, rank, ierr)

    !-------------------------
    ! Build file name
    !-------------------------
    write(bfstr,'(G0.6)') bf
    bfstr = adjustl(bfstr)
    write(dirname,'(A,A)') 'results/param_',trim(bfstr)
    fname = trim(dirname)//'/ETDRKcoeffs.bin'

    !-------------------------
    ! Rank 0 opens and reads
    !-------------------------
    if (rank==0) then
        open(10,file=fname,form='unformatted',access='stream', &
             status='old',iostat=ios)
        if (ios /= 0) then
            print *, "Could not open: ", fname
            stop "Error opening ETDRKcoeffs.bin"
        end if

        ! First element: kLength - verify it matches compiled value
        read(10) kLength_file
        if (kLength_file /= KLENGTH_PARAM) then
            print *, 'ERROR: File kLength=', kLength_file, ' does not match compiled KLENGTH_PARAM=', KLENGTH_PARAM
            stop 'kLength mismatch in ETDRKcoeffs - recompile with correct kLength'
        end if
        kLength = KLENGTH_PARAM
    end if

    ! Broadcast kLength
    call MPI_Bcast(kLength,1,MPI_INTEGER,0,comm,ierr)

    ! Arrays are fixed-size, no allocation needed

    !-------------------------
    ! Read all 36 arrays (rank 0)
    !-------------------------
    if (rank==0) then
        call rd(EX1_1); call rd(EX2_1); call rd(Q_1); call rd(f1_1); call rd(f2_1); call rd(f3_1)
        call rd(EX1_2); call rd(EX2_2); call rd(Q_2); call rd(f1_2); call rd(f2_2); call rd(f3_2)
        call rd(EX1_3); call rd(EX2_3); call rd(Q_3); call rd(f1_3); call rd(f2_3); call rd(f3_3)
        call rd(EX1_4); call rd(EX2_4); call rd(Q_4); call rd(f1_4); call rd(f2_4); call rd(f3_4)
        call rd(EX1_5); call rd(EX2_5); call rd(Q_5); call rd(f1_5); call rd(f2_5); call rd(f3_5)
        call rd(EX1_6); call rd(EX2_6); call rd(Q_6); call rd(f1_6); call rd(f2_6); call rd(f3_6)

        close(10)
    end if

    !-------------------------
    ! Broadcast everything
    !-------------------------
    call bcast_all()

contains

    !==============================
    subroutine rd(a)
    !==============================
        real(dp), intent(out) :: a(:)
        read(10) a
    end subroutine rd

    !==============================
    subroutine bc(a)
    !==============================
        real(dp), intent(inout) :: a(:)
        call MPI_Bcast(a, size(a), MPI_DOUBLE_PRECISION, 0, comm, ierr)
    end subroutine bc

    !==============================
    subroutine bcast_all()
    !==============================
        call bc(EX1_1); call bc(EX2_1); call bc(Q_1); call bc(f1_1); call bc(f2_1); call bc(f3_1)
        call bc(EX1_2); call bc(EX2_2); call bc(Q_2); call bc(f1_2); call bc(f2_2); call bc(f3_2)
        call bc(EX1_3); call bc(EX2_3); call bc(Q_3); call bc(f1_3); call bc(f2_3); call bc(f3_3)
        call bc(EX1_4); call bc(EX2_4); call bc(Q_4); call bc(f1_4); call bc(f2_4); call bc(f3_4)
        call bc(EX1_5); call bc(EX2_5); call bc(Q_5); call bc(f1_5); call bc(f2_5); call bc(f3_5)
        call bc(EX1_6); call bc(EX2_6); call bc(Q_6); call bc(f1_6); call bc(f2_6); call bc(f3_6)
    end subroutine bcast_all

end subroutine load_ETDRKcoeffs_binary_named

!======================================================================
! 4. WRAPPER: unify all loads into one call for main program
!     (ExtForcing will be computed separately, not loaded)
!======================================================================
subroutine load_mat_files(bf, kLength, kVals, v_1,v_2,v_3,v_4,v_5,v_6, t, tStar, counter, &
                          weight, triadFlag, outsideCutCell, insideCutCell, CxVals, CyVals, Q11, &
                          EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, &
                          EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
                          EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, &
                          EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
                          EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, &
                          EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6,ExtForcing)
  use mpi
  implicit none
  real(dp), intent(in) :: bf
  integer :: comm, ierr, rank, kLength_file
  integer, intent(out) :: kLength, counter
  real(dp), intent(out) :: t, tStar
  ! arrays
  real(dp), intent(out) :: kVals(:)
  real(dp), intent(out) :: v_1(:),v_2(:),v_3(:),v_4(:),v_5(:),v_6(:)
  real(dp), intent(out) :: weight(:,:,:),CxVals(:,:,:),CyVals(:,:,:)
  integer, intent(out) :: triadFlag(:,:,:),outsideCutCell(:,:,:),insideCutCell(:,:,:),Q11(:,:,:,:)
  real(dp), intent(out) :: &
      EX1_1(:),EX2_1(:),Q_1(:),f1_1(:),f2_1(:),f3_1(:), &
      EX1_2(:),EX2_2(:),Q_2(:),f1_2(:),f2_2(:),f3_2(:), &
      EX1_3(:),EX2_3(:),Q_3(:),f1_3(:),f2_3(:),f3_3(:), &
      EX1_4(:),EX2_4(:),Q_4(:),f1_4(:),f2_4(:),f3_4(:), &
      EX1_5(:),EX2_5(:),Q_5(:),f1_5(:),f2_5(:),f3_5(:), &
      EX1_6(:),EX2_6(:),Q_6(:),f1_6(:),f2_6(:),f3_6(:),ExtForcing(:)

  character(len=64) :: bfstr, dirname, f_kvals, f_EF
  integer :: ios

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)

  !---------------------------------------
  ! 1. checkpoint
  !---------------------------------------
  call load_checkpoint_binary(bf, v_1,v_2,v_3,v_4,v_5,v_6,t,tStar,counter,kLength,comm)

  !---------------------------------------
  ! 2. weightStuff
  !---------------------------------------
  call load_weightStuff_binary(bf,kLength,weight,triadFlag,outsideCutCell,insideCutCell,CxVals,CyVals,Q11,comm)

  !---------------------------------------
  ! 3. ETDRK coeffs
  !---------------------------------------
  call load_ETDRKcoeffs_binary_named(bf,kLength, &
       EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, &
       EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
       EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, &
       EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
       EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, &
       EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6,comm)

  !---------------------------------------
  ! 4. kVals only (ExtForcing computed separately)
  !---------------------------------------
 write(bfstr,'(G0.6)') bf
bfstr  = adjustl(bfstr)
write(dirname,'(A,A)') 'results/param_', trim(bfstr)

f_kvals = trim(dirname)//'/kVals.bin'

  if (rank==0) then
     open(11,file=f_kvals,form='unformatted',access='stream',status='old',iostat=ios)
     if(ios/=0) print *, f_kvals
     if(ios/=0) stop 'Error opening kVals.bin'
     read(11) kLength_file
     if (kLength_file /= KLENGTH_PARAM) then
         print *, 'ERROR: File kLength=', kLength_file, ' does not match compiled KLENGTH_PARAM=', KLENGTH_PARAM
         stop 'kLength mismatch in kVals - recompile with correct kLength'
     end if
     kLength = KLENGTH_PARAM
     read(11) kVals
     close(11)
  end if

  call MPI_Bcast(kLength,1,MPI_INTEGER,0,comm,ierr)
  call MPI_Bcast(kVals,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)

  if(rank==0) print *,'Loaded .bin files for bf=',bf,' kLength=',kLength,'kVals=',kVals







  write(bfstr,'(G0.6)') bf
bfstr  = adjustl(bfstr)
write(dirname,'(A,A)') 'results/param_', trim(bfstr)

f_EF = trim(dirname)//'/ExtForcing.bin'

  if (rank==0) then
     open(11,file=f_EF,form='unformatted',access='stream',status='old',iostat=ios)
     if(ios/=0) print *, f_EF
     if(ios/=0) stop 'Error opening ExtForcing.bin'
     ! Arrays are fixed-size, no allocation needed
     read(11) ExtForcing
     close(11)
  end if

  call MPI_Bcast(ExtForcing,kLength,MPI_DOUBLE_PRECISION,0,comm,ierr)

  if(rank==0) print *,'Loaded .bin files for bf=',bf,' kLength=',kLength

end subroutine load_mat_files


end module data_loader_all
