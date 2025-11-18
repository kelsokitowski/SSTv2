PROGRAM mainSST
  use mpi
  use array_dimensions
  use integrationModules
  use EDQNMstratifiedModules
  use job_parameters
  use data_loader_all
  use getForcing
  use timeRoutines
  implicit none
  !include 'mpif.h'
  !use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64

  integer :: rank, size_Of_Cluster, ierr, passiveICtrigger, SSflag
  integer :: root = 0
  !=============================================================
    !  MPI and basic parameters
    !=============================================================
    integer ::  kLength, counter
    real(dp) :: bf, t, tStar

    !=============================================================
    !  Checkpoint arrays (1D) - Fixed size using KLENGTH_PARAM
    !=============================================================
    real(dp) :: v_1(KLENGTH_PARAM), v_2(KLENGTH_PARAM), v_3(KLENGTH_PARAM)
    real(dp) :: v_4(KLENGTH_PARAM), v_5(KLENGTH_PARAM), v_6(KLENGTH_PARAM)

    !=============================================================
    !  WeightStuff arrays (3D) - Fixed size using KLENGTH_PARAM
    !=============================================================
    real(dp) :: weight(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    integer  :: triadFlag(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    integer  :: outsideCutCell(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    integer  :: insideCutCell(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    real(dp) :: CxVals(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    real(dp) :: CyVals(KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)
    integer  :: Q11(2,KLENGTH_PARAM,KLENGTH_PARAM,KLENGTH_PARAM)

    !=============================================================
    !  ETDRK coefficient arrays (1D) - Fixed size using KLENGTH_PARAM
    !=============================================================
    real(dp) :: &
        EX1_1(KLENGTH_PARAM), EX2_1(KLENGTH_PARAM), Q_1(KLENGTH_PARAM),&
        f1_1(KLENGTH_PARAM), f2_1(KLENGTH_PARAM), f3_1(KLENGTH_PARAM), &
        EX1_2(KLENGTH_PARAM), EX2_2(KLENGTH_PARAM), Q_2(KLENGTH_PARAM),&
        f1_2(KLENGTH_PARAM), f2_2(KLENGTH_PARAM), f3_2(KLENGTH_PARAM), &
        EX1_3(KLENGTH_PARAM), EX2_3(KLENGTH_PARAM), Q_3(KLENGTH_PARAM),&
        f1_3(KLENGTH_PARAM), f2_3(KLENGTH_PARAM), f3_3(KLENGTH_PARAM), &
        EX1_4(KLENGTH_PARAM), EX2_4(KLENGTH_PARAM), Q_4(KLENGTH_PARAM),&
        f1_4(KLENGTH_PARAM), f2_4(KLENGTH_PARAM), f3_4(KLENGTH_PARAM), &
        EX1_5(KLENGTH_PARAM), EX2_5(KLENGTH_PARAM), Q_5(KLENGTH_PARAM),&
        f1_5(KLENGTH_PARAM), f2_5(KLENGTH_PARAM), f3_5(KLENGTH_PARAM), &
        EX1_6(KLENGTH_PARAM), EX2_6(KLENGTH_PARAM), Q_6(KLENGTH_PARAM),&
        f1_6(KLENGTH_PARAM), f2_6(KLENGTH_PARAM), f3_6(KLENGTH_PARAM)

    !=============================================================
!  Solver/state arrays - Fixed size using KLENGTH_PARAM
!=============================================================
real(dp) :: kVals(KLENGTH_PARAM)
real(dp) :: dt(4)
real(dp) :: E(KLENGTH_PARAM), EHdir(KLENGTH_PARAM), EHdirNew(KLENGTH_PARAM), EHpol(KLENGTH_PARAM), EHpolNew(KLENGTH_PARAM)
real(dp) :: ET(KLENGTH_PARAM), Enew(KLENGTH_PARAM), ETH(KLENGTH_PARAM), ETHnew(KLENGTH_PARAM), ETnew(KLENGTH_PARAM)
real(dp) :: F(KLENGTH_PARAM), Fnew(KLENGTH_PARAM)
real(dp) :: HDIR(KLENGTH_PARAM), HPOL(KLENGTH_PARAM), HT(KLENGTH_PARAM)
real(dp) :: mu(KLENGTH_PARAM), mu1(KLENGTH_PARAM), mu3(KLENGTH_PARAM)
real(dp) :: Einitial(KLENGTH_PARAM), forcing(6), fl(KLENGTH_PARAM), fEta(KLENGTH_PARAM)
real(dp) :: RHSeqn1(KLENGTH_PARAM), RHSeqn2(KLENGTH_PARAM), RHSeqn3(KLENGTH_PARAM),&
        RHSeqn4(KLENGTH_PARAM), RHSeqn5(KLENGTH_PARAM), RHSeqn6(KLENGTH_PARAM)
real(dp) :: localRHSeqn1(KLENGTH_PARAM), localRHSeqn2(KLENGTH_PARAM), localRHSeqn3(KLENGTH_PARAM),&
        localRHSeqn4(KLENGTH_PARAM), localRHSeqn5(KLENGTH_PARAM), localRHSeqn6(KLENGTH_PARAM)
real(dp) :: duu1(KLENGTH_PARAM), duu2(KLENGTH_PARAM), duu3(KLENGTH_PARAM), duu4(KLENGTH_PARAM), duu5(KLENGTH_PARAM),&
        duu6(KLENGTH_PARAM)
real(dp) :: ExtForcing(KLENGTH_PARAM)


! Scalars (not allocatable)
real(dp) :: w(4), NtauInv
integer  :: kLengthTest, kj, v, v1, nSteps, timesteps, storeCount, y, kjl
integer  :: local_kLength, local_kStart, local_kEnd, termFlag,nPrint,qjMin,qjMax
real(dp) :: a, A1, A3, D, deltaT, deltaT_Eta, deltaTL, dt_v, epsTest, epzilon, eta
real(dp) :: Froude, NotintegerVal, k, K0, KE0, L, N, nu, Pr, Re_l, Re_taylor
real(dp) :: sigma, TauEta, TauL, TauL2, kEta,ReCorrector
real(dp) :: epsTestOnNuCubed, c0, du, kai, TauL0, gammagh,h,ySquaredAvg
integer :: comm = MPI_COMM_WORLD


 
  double precision, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)

  ! Set kLength from compile-time parameter
  kLength = KLENGTH_PARAM

  !YOU WILL NEED TO ADJUST LENGTH OF ARRAYS TO AVOID ALLOCATE FOR NOW. AFTER Y LOOP NO FURTHER EDITS NEEDED
  termFlag = 0;
  SSflag = 0;
  A1 = 0.355;
  A3 = 1.3;
  du = 10.0**(1.0/17.0);
  bf = 1;
  Pr = 1.;
   ! Initialize parameters (rank 0 reads env vars, broadcasts)
    !call initialize_parameters_mpi(bf, Pr, MPI_COMM_WORLD)
     call load_mat_files(bf, kLength, kVals, v_1,v_2,v_3,v_4,v_5,v_6, t, tStar, counter, &
                        weight, triadFlag, outsideCutCell, insideCutCell, CxVals, CyVals, Q11, &
                        EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, &
                        EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
                        EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, &
                        EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
                        EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, &
                        EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6,ExtForcing)

    if (rank == 0) then
        print *, "Initialization complete: bf =", bf, ", Pr =", Pr,' kLength =',kLength
    endif
!=============================================================
!  Arrays are now fixed-size (no allocation needed)
!=============================================================
    ! No allocation needed - arrays are declared with fixed size KLENGTH_PARAM

!print *, 'kVals=',kVals

 ! if ( rank == 0 ) then
     print *, 'bf = ',bf
     if ( (bf< 1) .and. (bf > 0.0) ) then
        Re_taylor = 400.0*(2.);
        Re_l = (Re_taylor**2.0) *3.0/20.0;
     else
        Re_taylor = 400.0*(1.+sqrt(sqrt(bf)));
        Re_l = Re_taylor**2.0 *3.0/20.0;
        ReCorrector = 0.0000047220*bf**4.0-0.0010379*bf**3.0+0.0938686817*bf**2.0 &
             +1.6641001561*bf-0.9318513730;
        Re_l = Re_l*ReCorrector
        if (bf < 7.0) then
           Re_l = Re_l*5.0;
        end if

     end if
     a= -1.0;
     L = 10.0_dp; !lengthscale of the largest eddies
   !  Re_taylor = 400.0_dp;
     print *, 'L = ',L, ' Re_taylor = ',Re_taylor
    ! Re_l = Re_taylor**2.0_dp *3.0_dp/20.0_dp;
     eta = Re_l**(-3.0_dp/4.0_dp)* L;
     nu = 0.05;
     
     print *, ' nu = ',nu
     !Re_taylor = sqrt(20*Re_l/3);
     !eta = (nu**3/epzilon)**(1/4);
     epzilon = nu**3.0/(eta**4.0);
     !nu   = ( eta**4*epzilon )**(1/3);

     sigma = 2.0;
     !kEta = (epzilon/(nu**3.0_dp))**(1.0_dp/4.0_dp);
     kEta = 1./eta;
     !K0 = 10.0_dp**(-6.0_dp)*kEta;
     K0 = 0.01; !manual choice will give kLength kVals
     !K0 = 0.0747125785428542; !manual choice will give 144 kVals,  based on below spectra where kMax is 4th value in Kvals
     !N = 10;


     D = nu/Pr;
!!$
!!$     kVals(1) = K0;
!!$
!!$     y = 2;
!!$
!!$     do while (kVals(y-1)< (10.0_dp*kEta*sqrt(Pr)) )
!!$
!!$        kVals(y) = du*kVals(y-1);
!!$        y = y+1;
!!$     end do
!!$
!!$
!!$
!!$
!!$
!!$
!!$     kLengthTest = y-1;
!!$
!!$     if (kLengthTest /= kLength) then
!!$        print *, 'kEta = ',kEta
!!$        print *, kVals
!!$        print *, 'kLength does not equal expected length, error',', kLengthTest = ',kLengthTest, ' y = ',y
!!$        print *, ' Size of kVals returns ',size(kVals,DIM=1)
!!$        print *, 'Last entry of kVals is',kVals(size(kVals,DIM=1))
!!$     end if
!!$
!!$
!!$     do kj = 1,kLength
!!$        k = kVals(kj);
!!$
!!$        fl(kj) = ( k*L/ (((k*L)**1.5+1.5-sigma/4.0)**(2.0/3.0)) )**(5.0/3.0+sigma);
!!$        fEta(kj) = exp(-5.3*(((k*eta)**4.0+0.4**4.0)**(1.0/4.0)-0.4));
!!$
!!$        !E= K0*k.**(-5/3)*epzilon**(2/3).*fl(k).*fEta(k*eta);
!!$        E(kj)=  k**(-5.0/3.0) *fl(kj)*fEta(kj);
!!$
!!$     end do
!!$
!!$     c0 = epzilon**(1./3.) / (2.*nu*trapezoidalIntegration(kVals,kVals**2.*E));
!!$     epsTest = 2*nu*trapezoidalIntegration(kVals,kVals**2.*E*epzilon**(2./3.)*c0);
!!$
!!$     print *, 'theoretical epzilon = ',epzilon,'integrated epzilon = ',epsTest
!!$     print *, 'c0 = ',c0
!!$
!!$
!!$     !print *, 'fl(1) = ',fl(1), 'fEta(1) = ', fEta(1)
!!$
!!$     !now normalize E so that epzilon is 1.
!!$
!!$     !epsTest = (2.0*nu*trapezoidalIntegration(kVals,kVals**2.0*E))**(3.0);
!!$     !epsTestOnNuCubed = (2.0*trapezoidalIntegration(kVals,kVals**2.0*E))**(3.0);
!!$     ! c0 = ( eta**4.0*epsTestOnNuCubed)**(-1./3.); !correction as eta^4 should be 1 over epsOnNuCubed
!!$     ! epsTestOnNuCubed = epsTestOnNuCubed * c0**3.0;
!!$     ! epsTest = epsTestOnNuCubed * nu**3.;
!!$     !epzilon = epsTest;
!!$
!!$     E = c0*epzilon**(2./3.)*E;
!!$     print *, 'E(100) = ',E(100),'epzilon = ',epzilon,' kvals(end)= ',k


     !close all; loglog(kVals,E);
print *, 'kVals(1) = ',kVals(1)
E = v_1
print *, 'E=',E
     KE0 = 0.5*trapezoidalIntegration(kVals,E)
print *, 'KE0',KE0

     TauL = KE0/epzilon; !eddy turnover time
     TauEta = (nu/epzilon)**0.5 ;
     TauL2 = TauEta*Re_l**0.5;
     TauL0 = TauL



     !##################################################################
     NtauInv = 0.0 ;  !this is N in units of TauEpsilon
     ! N = NtauInv/(TauL*0.63); !This is N in units of seconds.
     !Froude = 1.0/(TauL*N);
     !bf = 1.0; ! buoyancy frequency 'N' in most papers
     !Froude = 1.e13
     !##################################################################

     deltaTL = TauL/1000.0;


     !     E = 0.0; !set all IC to 0!


     Einitial = E;

     !Call getGaussExtForcing(kVals, E, eta, nu, ExtForcing)


     print *, 'External forcing at index 8 is ',ExtForcing(8),'rank= ',rank
    ! MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !ExtForcing = 0.0;
    ! E = 0.0;
    if (kLength >=130) then
    h = 1.0/(2.0*(D)*kVals(kLength)**2.0)*10.0;
    print *, ' h,kLength= ',h,kLength
    else
        h = 1.0/(2.0*(D)*kVals(kLength)**2.0)*10.0;
    end if
! end if !rank ==0 then

 qjMin = 1.0d0
    qjMax = 10.0d0
    ySquaredAvg = 1.0d0
    nPrint = 1

print *,'made it to before final bcast mainSST, rank = ',rank
    ! All arrays are now fixed-size (no allocation needed)

    ! Broadcast large arrays to all ranks
    call MPI_Bcast(kVals, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_1, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_2, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_3, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_4, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_5, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(v_6, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)
    call MPI_Bcast(ExtForcing, kLength, MPI_DOUBLE_PRECISION, root, comm, ierr)

    call MPI_Barrier(comm, ierr)
    if (rank == root) print *, 'All data broadcast complete — entering time loop.'

    !============================================================
    ! Time Loop
    !============================================================
   call run_timeloop_mpi( &
    bf, nu, D, h, &
    kVals, &
    v_1, v_2, v_3, v_4, v_5, v_6, &  ! initial state already loaded
    EX1_1,EX2_1,Q_1,f1_1,f2_1,f3_1, EX1_2,EX2_2,Q_2,f1_2,f2_2,f3_2, &
    EX1_3,EX2_3,Q_3,f1_3,f2_3,f3_3, EX1_4,EX2_4,Q_4,f1_4,f2_4,f3_4, &
    EX1_5,EX2_5,Q_5,f1_5,f2_5,f3_5, EX1_6,EX2_6,Q_6,f1_6,f2_6,f3_6, &
    weight, triadFlag, outsideCutCell, insideCutCell, CxVals, CyVals, Q11, &
    ExtForcing, du, qjMin, qjMax, ySquaredAvg, A1, A3, &
    ! restart/time bookkeeping (loaded from checkpoint)
    t, tStar, &
    ! outputs/updated scalars
    nPrint, counter, &
    comm)

    !============================================================
    ! Finalize MPI
    !============================================================
    call MPI_Barrier(comm, ierr)
    if (rank == root) then
        print *, 'Finished SLURM array task for bf=', bf, ' Pr=', Pr
    end if
    call MPI_Finalize(ierr)

end program mainSST
