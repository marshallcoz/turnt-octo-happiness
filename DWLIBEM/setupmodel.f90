      module setupmodel
      contains

      ! SET UP & READ
      subroutine setupdirectories
      use glovars, only : borrarDirectorio, rutaOut,PrintNum,saveG
      character(10) :: time
      CHARACTER(len=400) :: path
      CHARACTER(len=32) :: arg
      integer :: status
      ! preparar los directorios de trabajo y de resultados
      CALL getcwd(path)
      write(6,'(a,/,a)') 'At:',TRIM(path)
      write(path,'(a,a)') trim(path),"/WorkDir"
      CALL chdir(trim(path))
      CALL getcwd(path)
      write(6,'(a,/,a)') 'WORKDIRECTORY:',TRIM(path)
      call getMainInput
      CALL get_command_argument(1, arg)
      if ((trim(arg) .ne. '-r') .and. (trim(arg) .ne. '-f')) borrarDirectorio = .true.
      if (borrarDirectorio .eqv. .false.) then
      call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
      write(rutaOut,'(a,a)') "mkdir outs_",time
      call system(trim(adjustl(rutaOut)))
      call chdir(trim(adjustl(rutaOut)),status)
      write(rutaOut,'(a,a)') "outs_",time
      else
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then
        if (trim(arg) .eq. '-L') then
        saveG = .false.
        go to 296
        end if
      end if
      call system('rm -rf outs')
      call system('mkdir outs')
 296  CALL chdir("outs",status)
      write(rutaOut,'(a)') "outs"
      end if!
      if (status .eq. 0) call chdir("..") !workdir
      write(path,'(a,a,a)') 'cp -r ins ',trim(adjustl(rutaOut)),'/insbackup'
      call system(trim(adjustl(path)))
      call chdir(trim(adjustl(rutaOut)),status) !outs
      if (status .eq. 0) call system("rm *.*")
      if (PrintNum /= 6) open(PrintNum,FILE= "logfile.txt")
      call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
      call system('mkdir subdivs')
      call system('cp ../../DWNLIBEM.f90 DWNLIBEM.f90')
      call chdir("..")
      write(PrintNum,'(///)')
      end subroutine setupdirectories

      subroutine getMainInput
      use glovars
      use GeometryVars , only : staywiththefinersubdivision,finersubdivisionJ,&
       longitudcaracteristica_a, fraccionDeSmallestWL_segm_de_esquina
      use resultvars, only : overDeterminedSystem,OD_Jini,OD_Jend
      use waveNumVars, only : SpliK
      implicit none
      logical :: lexist
      CALL chdir("ins")
      inquire(file="maininput.txt",exist=lexist)
      if (lexist) then
        OPEN(35,FILE="maininput.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "maininput.txt" on Working directory'
      end if
      READ(35,'(I1)') verbose!; print*,"verbose =",verbose
      READ(35,'(L1)') makevideo!; print*,"make a video =",makevideo
      READ(35,'(L1)') vivaChurubusco
      READ(35,'(L1)') vivaCine
      READ(35,'(L1)') workBoundary!; print*,"boundary? ",workBoundary
      READ(35,'(L1)') flip12; if(.not. workBoundary) flip12 = .false.
      READ(35,'(L1)') plotFKS!; print*,"plotFK?",plotFKS
      READ(35,'(L1)') plotFKCK
      READ(35,'(L1)') PrintEtas
      READ(35,'(L1)') PlotFilledSabanas
      READ(35,'(L1)') saveG!; print*,"Save Green funcs?", saveG
      READ(35,*) multSubdiv!; print*,"division multiple = ", multSubdiv
      READ(35,*) staywiththefinersubdivision, finersubdivisionJ
      READ(35,*) overDeterminedSystem,OD_Jini,OD_Jend
      READ(35,*) cKbeta!; print*,"multBminIntegrando = ", cKbeta
      READ(35,*) SpliK
      read(35,*) periodicdamper!; print*,"Periodic sources damping factor = ", periodicdamper
      READ(35,*) useAzimi
      read(35,*) developerfeature
      read(35,*) borrarDirectorio
      read(35,*) PrintNum
      read(35,*) longitudcaracteristica_a
      read(35,*) fraccionDeSmallestWL_segm_de_esquina
      read(35,*) PWfrecReal
      read(35,*) comoFacDeAmpliDinamica
      close(35)
      CALL chdir("..")
      end subroutine getMainInput

      subroutine getSoilProps
      use soilVars; use waveNumVars; use waveVars, only : dt; use fitting
      use gloVars, only : verbose,PI, cKbeta,outpf => PrintNum,periodicdamper
      implicit none
      logical :: lexist
      real :: H, ALF, BET, RO
      real*8 :: maxBeta,BEALF, k1_3,kmax,frac!,m
      integer :: J,i

      CALL chdir("ins")
      inquire(file="HVDEPTH.txt",exist=lexist)
      if (lexist) then
        OPEN(7,FILE="HVDEPTH.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "HVDEPTH.txt" on Working directory'
      end if

      READ(7,*)
      READ(7,'(I2)')N   !NUMBER OF LAYERS NOT COUNTING UPPER/LOWER HALF-SPACE
      READ(7,*)
      write(outpf,'(/,A,I2,A,/)') '  ',N,' layers'
      ALLOCATE (Z(0:N+1));ALLOCATE (RHO(0:N+2)); ALLOCATE(ANU(0:N+2))
      ALLOCATE (BETA0  (0:N+2)); ALLOCATE (ALFA0(0:N+2))
      ALLOCATE (BETA   (0:N+2)); ALLOCATE (ALFA (0:N+2))
      ALLOCATE (LAMBDA0(0:N+2)); ALLOCATE (AMU0 (0:N+2))
      ALLOCATE (LAMBDA (0:N+2)); ALLOCATE (AMU  (0:N+2))
      allocate (layershadecolor(0:N+2))

      Z(0)=real(1000,8);      Z(1)=real(0,8)
      if (verbose >= 1) write(outpf,'(a)')&
      '        depth       alpha0    beta0      mu0     rho      lambda0       nu0'
      DO J=1,N
         READ(7,*) H, ALF, BET, RO
        if (H .lt. 0.0) then
          Z(0) = real(H)
          AMU0(0)=RO*BET**2.0
          BETA0(0)=BET
          ALFA0(0)=ALF
          RHO(0) = RO
          LAMBDA0(0)=RHO(0)*ALF**2.0 - real(2.)*real(AMU0(0))
          BEALF = beta0(0)/alfa0(0)
          anu(0) = (bealf**2 - 0.5)/(-1 + bealf**2)!#< b
         write(outpf,'(A,F7.1,2x, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
       ' -inf.  - ',Z(1),ALFA0(0),BETA0(0),real(AMU0(0)),  RHO(0), real(LAMBDA0(0)),real(anu(J))!#>
          READ(7,*) H, ALF, BET, RO
        end if
         Z(J+1)=Z(J)+real(H)
         AMU0(J)=RO*BET**2.0
         BETA0(J)=BET
         ALFA0(J)=ALF
         RHO(J) = RO
         LAMBDA0(J)=RHO(J)*ALF**2.0 - real(2.)*real(AMU0(J))
!        BEALF=SQRT((0.5-ANU)/(1.0-ANU)) !IF POISSON RATIO IS GIVEN
         BEALF = beta0(J)/alfa0(J)
         anu(J) = (bealf**2 - 0.5)/(-1 + bealf**2)
!        ALFA(J)=BET/BEALF !#< b
          write(outpf,&
          '(F7.1,A,F7.1,2x, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
          Z(J),' - ',Z(J+1),ALFA0(J),BETA0(J),real(AMU0(J)),&
          RHO(J), real(LAMBDA0(J)),real(anu(J)) !#>
      END DO

      READ(7,*) H, ALF, BET, RO
      if (H .lt. 0.0) then !(caso: semiespacio arriba, semiespacio abajo)
         Z(0) = real(H)
         AMU0(0)=RO*BET**2.0
         BETA0(0)=BET
         ALFA0(0)=ALF
         RHO(0) = RO
         LAMBDA0(0)=RHO(0)*ALF**2.0 - real(2.)*real(AMU0(0))
         BEALF = beta0(0)/alfa0(0)
         anu(0) = (bealf**2 - 0.5)/(-1 + bealf**2)!#< b
         write(outpf,'(A,F7.1,2x, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
         ' -inf.  - ',Z(1),ALFA0(0),BETA0(0),real(AMU0(0)),&
         RHO(0), real(LAMBDA0(0)),real(anu(0)) !#>
         READ(7,*) H, ALF, BET, RO
      end if
      AMU0(N+1)=RO*BET**2
      BETA0(N+1)=BET
      ALFA0(N+1)=ALF
      RHO(N+1) = RO
      LAMBDA0(N+1)=RHO(n+1)*ALF**2 - real(2.)*real(AMU0(n+1))
      BEALF = beta0(N+1)/alfa0(N+1)
      anu(N+1) = (bealf**2 - 0.5)/(-1 + bealf**2) !#< b
         write(outpf,'(F7.1,2x,A, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
         Z(1),' -  inf. ',ALFA0(N+1),BETA0(N+1),real(AMU0(N+1)),&
         RHO(N+1), real(LAMBDA0(N+1)),real(anu(N+1)) !#>
      i = 1;  if (Z(0) .lt. 0.0) i = 0
      minBeta = minval(beta0(i:N+1));  maxBeta = maxval(beta0(i:N+1))
      if (abs(minBeta - maxBeta) .lt. 0.01) then
       layershadecolor(i:N+1) = 0.8_4
      else
       do J=i,N+1
        layershadecolor(J)= real(0.6-(maxBeta-beta0(J))*((0.6-0.89)/(maxBeta-minBeta)),4)
!       print*,i,layershadecolor(J)
       end do
      end if
      READ(7,*)
      READ(7,*)DFREC,NFREC,NPTSTIME,nK,Qq,TW
      close(7)
      CALL chdir("..")

      !#< r ---- Calculamos NFREC resultados y hacemos zeropadding hasta NPSTIME--- !#>
         NPTSTIME = int(log10(real(NPTSTIME)) / log10(2.)+0.99 )
         NPTSTIME = 2** NPTSTIME     !adjuste to a 2**N value
         if (NPTSTIME .lt. 2*NFREC) then
           NPTSTIME=2*NFREC
           print*,"WARNING, modified NPTSTIME to" , NPTSTIME
         end if
      Dt = (1.0) / (real(NPTSTIME) * DFREC)

      !#< r ----  dk para ver estar CKbeta veces por arriba de polo Rayleigh------- !#>
      frac = 1./3.
      kmax = 0.9* (2*pi*DFREC*NFREC) / minBeta * cKbeta !1.5
      DK = kmax / nk
      k1_3 = (2*pi*DFREC*NFREC * frac) / minBeta * cKbeta

      allocate(vecNk(NFREC+1))
      ! cantidad de numeros de onda {donde se invierte la matriz} en cada frecuencia
      ! el resto hasta nmax+1(o hasta SpliK) se interpola
      vecNk = (/(i,i=1, NFREC+1)/)
      call splineIn(vecNk, & ! vector x [int]
                    vecNk, & ! vector interpolado [int]
                    int(nfrec * frac),& ! n puntos [x0 ... x1-1]
                    NFREC+1-int(nfrec * frac), & ! [x1 ...  x2 ]
                    1,                  int(0.55*nK), &  !x0,y0
                    int(nfrec * frac)+1,int(0.6*nk), &  !x1,y1
                    NFREC+1,            nk) !x2,y2
      NMAX = vecNK(NFREC)!     IF (vecNK(NFREC) .gt. nmax)
      NMAX = int(log10(real(NMAX)) / log10(2.)+0.99 )
      NMAX = 2**NMAX     !adjuste to a 2**N value

      ! complex frecuency to "smooth" the poles and be able to integrate
        ! Bouchon (2003) OMEI entre -pi/T y -2pi/T ; T= 2pi/DFREC tiempo total
        ! tengo TW la ventana de tiempo de interés.
      OMEI = - periodicdamper*PI/TW ! se cumple que exp(omei TW) << 1

      if (verbose .ge. 1) then !#< b
       write(outpf,'(a)') &
       "--- frecuency --------------------------------------------------------------------"
       write(outpf,'(a,F9.7,/,a,F15.5)') "   dt = (1.0) / (real(NPTSTIME) * DFREC) = ",Dt,"   Tmax=",Dt* NPTSTIME
       write(outpf,'(a,F8.1)') '   Atenuation Q = ',Qq
       write(outpf,'(a,I0,a,F8.4,a,F12.4,a,/)') &
           '   N. frequencies: ',NFREC,'  @',DFREC,'Hertz :. Fmax = ', &
           NFREC*DFREC,'Hertz'

       write(outpf,'(a)') &
       "--- DWN -------------------------------------------------------------------------"
       write(outpf,'(a,F15.7)') "   kmax = 0.9* (2*pi*DFREC*NFREC) / minBeta * cKbeta = ",kmax
       write(outpf,'(a,F9.7)') "   DK = kmax / nk = ",DK
       write(outpf,'(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)') "   nk a 3pt spline with (",&
          1,",",int(0.55*nK),"), (",&
          int(nfrec * frac)+1,",",int(0.6*nk),"), (",&
          NFREC+1,",",nk,")"
       write(outpf,'(a,I0)') '   nk average = ',int(sum(vecNK)/(NFREC+1))
       write(outpf,'(a,I0)') '   nmax (each sign): ',NMAX
       write(outpf,'(a,F12.7)') "   delta X = ", real(pi / (nMax * DK),4)
       write(outpf,'(a,EN14.2E2,a)') "   L = 2*pi/DK = ",2*pi/DK, "m"
       write(outpf,'(a,EN19.5E2,a)') &
       "   L/alfa = tp = ",(2*pi/DK)/maxval(abs(ALFA0)), "seconds"
       write(outpf,'(a,EN19.5E2,a)') &
       "   L/beta = ts = ",(2*pi/DK)/maxval(abs(BETA0)), "seconds"
       write(outpf,'(a,E10.2,a)') '   Frec. Imaginary part: ',OMEI,' rad/s'
       write(outpf,'(a)') &
       "---------------------------------------------------------------------------------"
      end if !#>
      end subroutine getSoilProps

      subroutine checarWisdom(a,b,c)
      use waveNumVars, only : planNmaxF,planNmaxB,&
                             planNfrecF,planNfrecB,&
                             planNtimeF,planNtimeB
      use, intrinsic :: iso_c_binding
      use glovars, only : PrintNum
      include 'fftw3.f03'
      type(C_PTR) :: plan
      integer(C_INT) :: j,a,b,c,flag
      complex*16, dimension(a) :: Ua,Va !frec
      complex*16, dimension(b) :: Ub,Vb !nmax
      complex*16, dimension(c) :: Uc,Vc !npstime
      character(len=100) :: tx
      logical :: lexist
      write(tx,'(a,I0,a,I0,a,I0,a)') "wis",a,"-",b,"-",c,".dat"
      inquire(file=trim(tx),exist=lexist)
      if (lexist .eqv. .false.) then
      flag = FFTW_PATIENT ! FFTW_MEASURE
      plan = fftw_plan_dft_1d(a, Ua,Va, FFTW_FORWARD,flag)
      plan = fftw_plan_dft_1d(a, Ua,Va, FFTW_BACKWARD,flag)
      plan = fftw_plan_dft_1d(b, Ub,Vb, FFTW_FORWARD, flag)
      plan = fftw_plan_dft_1d(b, Ub,Vb, FFTW_BACKWARD, flag)
      plan = fftw_plan_dft_1d(c, Uc,Vc, FFTW_FORWARD, flag)
      plan = fftw_plan_dft_1d(c, Uc,Vc, FFTW_BACKWARD, flag)
      j = fftw_export_wisdom_to_filename(trim(tx) // C_NULL_CHAR)
      call fftw_destroy_plan(plan)
      if (j .eq. 0) then
      print*,"i=",j," (non-zero on success)"
      stop "Error al exportar wisdom en :checarWisdom"
      else
      print*, "Plans created, at"
      print*,trim(tx)
      stop "Re run program"
      end if
      print*,"Se exportó wisdom ",trim(tx)
      else
      !importar
      j = fftw_import_wisdom_from_filename(trim(tx) // C_NULL_CHAR)
      if (j .eq. 0) then
      print*,"i=",j," (non-zero on success)"
      stop "Error al importar wisdom en :FFTW"
      else
      write(PrintNum,*)"  "
      write(PrintNum,*)"  leído wisdom ",trim(tx)
      ! crear planos
      flag = FFTW_WISDOM_ONLY
      planNfrecF = fftw_plan_dft_1d(a, Ua,Va, FFTW_FORWARD,flag)
      planNfrecB = fftw_plan_dft_1d(a, Ua,Va, FFTW_BACKWARD,flag)
      planNmaxF = fftw_plan_dft_1d(b, Ub,Vb, FFTW_FORWARD, flag)
      planNmaxB = fftw_plan_dft_1d(b, Ub,Vb, FFTW_BACKWARD, flag)
      planNtimeF = fftw_plan_dft_1d(c, Uc,Vc, FFTW_FORWARD, flag)
      planNtimeB = fftw_plan_dft_1d(c, Uc,Vc, FFTW_BACKWARD, flag)
      write(PrintNum,*)"  planos creados"
      end if
      end if
      end subroutine checarWisdom

      subroutine vaciarWisdom
      use waveNumVars, only : planNmaxF,planNmaxB,&
                             planNfrecF,planNfrecB,&
                             planNtimeF,planNtimeB
      use, intrinsic :: iso_c_binding
      include 'fftw3.f03'
      call fftw_destroy_plan(planNfrecF)
      call fftw_destroy_plan(planNfrecB)
      call fftw_destroy_plan(planNmaxF)
      call fftw_destroy_plan(planNmaxB)
      call fftw_destroy_plan(planNtimeF)
      call fftw_destroy_plan(planNtimeB)
      !call fftwf_forget_wisdom()
      !call fftwf_cleanup()
      end subroutine vaciarWisdom

      subroutine getInquirePoints
      use resultVars, only : Punto,nPtsolos,inqPoints, nIpts, iPtini,iPtfin, &
                       nSabanapts, nSecciones, Sabana,sabZeroini,sabZerofin,&
                       SabanaPlotIndividual, sabanaBajarAFrontera, &
                       n_OD,overDeterminedSystem,punEnlaFront,Punto2d,promP2d,negP2d, &
                       nBPt_topo,nBPt_cont,nBPt_vall
      use waveNumVars, only : NPTSTIME
      use GeometryVars, only: nXi,origGeom,n_topo,n_cont,n_vall,spaceBarVall
      use glovars, only : flip12
      implicit none
      interface
        function tellisoninterface(zi)
        real*8 :: zi
        end function tellisoninterface
      end interface

      integer :: thelayeris
      integer :: i,j,k,iIndex, nnsabanas,nnsecciones,nbouP,thisnsab
      logical :: lexist, tellisoninterface, ths_isoninterface,wtfk
      integer :: auxGuardarFK, ths_layer, sabanabajarmax
      real :: xini,deltax,zini,deltaz,dx,dz, SbanadeltaZ
      real*8 :: escalax,escalay,escaladx,escalady,offsetx,offsety,r,nx,nz
      logical :: cn,adentroOafuera,tellpunEnlaFront
      integer, dimension(0:2) :: reg
      type (Punto), dimension(:), allocatable :: auxInq
      CALL chdir("ins")
      ! read file
      inquire(file="interestingPoints.txt", exist=lexist)
      if(lexist)then
        open(7,FILE="interestingPoints.txt",FORM="FORMATTED")
      else
        write(6,'(a,a)') 'There is a missing input file. ',&
        'Check "interestingPoints.txt" on Working directory'
        stop 1
      end if
      READ(7,*) !Points of interest for an accelerogram
      READ(7,*) nIpts
      nPtsolos = nIpts
      READ(7,*) nnsabanas, nSabanapts, SabanaPlotIndividual, SbanadeltaZ, sabanabajarmax
      READ(7,*) nnsecciones, nsecciones
      READ(7,*) tellpunEnlaFront, spaceBarVall
      punEnlaFront=tellpunEnlaFront
      if (punEnlaFront) then; nbouP = nXI
      else; nbouP = 0
      end if
      READ(7,*) !__________________________________________
      READ(7,*) escalax,escalay
      READ(7,*) offsetx,offsety
      iPtini = 1
      if (nnsabanas .eq. 0) nSabanapts=0
      if (nnsecciones .eq. 0) nsecciones=0
      iPtfin = nIpts + nSabanapts + nsecciones + nbouP
      allocate(inqPoints(nIpts + nSabanapts + nsecciones + nbouP))
      inqPoints(:)%normal%x = 0
      inqPoints(:)%normal%z = 0
      inqPoints(:)%isOnInterface = .false.
      inqPoints(:)%isBoundary = .false.
      inqPoints(:)%guardarFK = .false.
      inqPoints(:)%guardarMovieSiblings = .false.
      inqPoints(:)%isSabana = .false.
      inqPoints(:)%isSourceSegmentForce= .false.
      inqPoints(:)%isOD= .false. ; n_OD = 0
      inqPoints(:)%region = 1
      inqPoints(:)%boundaryIndex = 0
      inqPoints(:)%atBou = .false.
      READ(7,*) !  X        Z          nx       nz     guardarFK
      do i=1, nIpts
         READ(7,*) inqPoints(i)%center%x, inqPoints(i)%center%z, &
             inqPoints(i)%normal%x, inqPoints(i)%normal%z, auxGuardarFK, &
             inqPoints(i)%isOD,inqPoints(i)%tipoFrontera
           inqPoints(i)%center%x = inqPoints(i)%center%x * escalax + offsetx
           inqPoints(i)%center%z = inqPoints(i)%center%z * escalay + offsety
           IF (inqPoints(i)%isOD) then
           n_OD = n_OD +1
           if (inqPoints(i)%tipoFrontera .eq. 1) n_OD = n_OD +1
           end if
      !encontrar el layer en el que estan o 0 si está sobre la interfaz
           inqPoints(i)%layer = thelayeris(real(inqPoints(i)%center%z,8))
           if (auxGuardarFK .ge. 1 ) inqPoints(i)%guardarFK = .true.
           inqPoints(i)%isOnInterface = &
                         tellisoninterface(real(inqPoints(i)%center%z,8))
      end do

!     addOD = 0
      iIndex = nIpts
      ! !#< r putnos sobre la geometría -------------------------- !#>
      if (punEnlaFront) then
      reg(0) = 0; reg(1)= 1; reg(2) = 2
      if (flip12) then
        reg(0) = 0; reg(1)= 2; reg(2) = 1
      end if
      nBPt_topo=0;nBPt_cont=0;nBPt_vall=0
      do i =1,nXI
       if (tellisoninterface(real(origGeom(i)%bord_A%z,8))) then
       if (origGeom(i)%tipoFrontera .eq. 1) then
       iPtfin = iPtfin -1
       nbouP = nbouP -1
       cycle
       end if;end if
       iIndex = iIndex + 1
       inqPoints(iIndex)%center = origGeom(i)%bord_A
!      print*,"------------------------";print*,inqPoints(iIndex)%center
       inqPoints(iIndex)%tipoFrontera = origGeom(i)%tipoFrontera
       IF (origGeom(i)%tipoFrontera .eq. 0) then
         inqPoints(iIndex)%region = reg(1)
       nBPt_topo = nBPt_topo+1
       if (i .eq. 1) then                            !#< r tipoFrontera 0 !#>
       inqPoints(iIndex)%normal = promP2d(origGeom(1)%normal,origGeom(n_topo)%normal)
!      print*,origGeom(1)%normal,origGeom(n_topo)%normal," ->",inqPoints(iIndex)%normal
       else
       inqPoints(iIndex)%normal = promP2d(origGeom(i)%normal,origGeom(i-1)%normal)
!      print*,origGeom(i)%normal,origGeom(i-1)%normal," ->",inqPoints(iIndex)%normal
       end if
       ELSE IF (origGeom(i)%tipoFrontera .eq. 1) then !#< r tipoFrontera 1 !#>
         inqPoints(iIndex)%region = reg(1)
       nBPt_cont=nBPt_cont+1
       if (i .eq. n_topo+1) then
       inqPoints(iIndex)%normal = negP2d(promP2d(origGeom(n_topo+1)%normal, &
                                                 origGeom(n_topo+n_cont)%normal))
!      print*,origGeom(n_topo+1)%normal,origGeom(n_topo+n_cont)%normal," ->", &
!                                                      inqPoints(iIndex)%normal
       else
       inqPoints(iIndex)%normal = negP2d(promP2d(origGeom(i)%normal,origGeom(i-1)%normal))
!      print*,origGeom(i)%normal,origGeom(i-1)%normal," ->",inqPoints(iIndex)%normal
       end if

       ELSE IF (origGeom(i)%tipoFrontera .eq. 2) then !#< r tipoFrontera 2 !#>
         inqPoints(iIndex)%region = reg(2)
       nBPt_vall=nBPt_vall+1
       if (i .eq. n_topo+n_cont+n_vall) then
       inqPoints(iIndex)%normal = promP2d(origGeom(n_topo+n_cont+1)%normal, &
                                          origGeom(n_topo+n_cont+n_vall)%normal)
!      print*,origGeom(n_topo+n_cont+1)%normal,origGeom(n_topo+n_cont+n_vall)%normal, &
!                                                     " ->",inqPoints(iIndex)%normal
       else
       inqPoints(iIndex)%normal = promP2d(origGeom(i)%normal,origGeom(i+1)%normal)
!      print*,origGeom(i)%normal,origGeom(i-1)%normal," ->",inqPoints(iIndex)%normal
       end if
       ! hacer que no esten exacto sobre la frotnera
       inqPoints(iIndex)%center%x = inqPoints(iIndex)%center%x + spaceBarVall * inqPoints(iIndex)%normal%x
       inqPoints(iIndex)%center%z = inqPoints(iIndex)%center%z + spaceBarVall * inqPoints(iIndex)%normal%z

       !#< r sobredeterminar sistema en frontera 2 !#>
!      if ( overDeterminedSystem ) then
!          n_OD = n_OD +1
!          inqPoints(iIndex)%isOD = .true.
!      end if
       END IF
       if (abs(inqPoints(iIndex)%normal%x) .lt. 0.0001) inqPoints(iIndex)%normal%x = 0
       if (abs(inqPoints(iIndex)%normal%z) .lt. 0.0001) inqPoints(iIndex)%normal%z = 0
       ! angulo para hacer la rotación a coordenadas normal y tangencial
       r = sqrt(inqPoints(iIndex)%normal%z**2 + inqPoints(iIndex)%normal%x**2)
       inqPoints(iIndex)%cosT = inqPoints(iIndex)%normal%x/r
       inqPoints(iIndex)%sinT = inqPoints(iIndex)%normal%z/r
!      print*,inqPoints(iIndex)%normal,inqPoints(iIndex)%cosT,inqPoints(iIndex)%sinT

       inqPoints(iIndex)%layer = origGeom(i)%layer
       inqPoints(iIndex)%isOnInterface = .false. !tellisoninterface(real(inqPoints(iIndex)%center%z,8))
       inqPoints(iIndex)%atBou = .true.

      end do ! iIndex
      nIpts = nIpts + nbouP
      end if ! pu en la frontera

      if (.not. overDeterminedSystem) n_OD = 0        !#< r en caso de que no hay puntos !#>
      if (n_OD .eq. 0) then
      overDeterminedSystem = .false. ! de coloc. adicionales
      print*,"Found no points to overdetermine the system"
      end if
      !
      if (nnsecciones .gt. 0) then
      read(7,*) !Secciones -------------------------
      read(7,*) escalax,escalay !; print*,escalax,escalay
      read(7,*) escaladx,escalady !; print*,escalax,escalay
      READ(7,*) offsetx,offsety !; print*, offsetx,offsety
      read(7,*) ! npuntos     Xi        deltax      Zi       deltaz       nx       nz
      do j=1,nnsecciones
        read(7,*) thisnsab,xini,deltax,zini,deltaz,nx,nz
        dx = 0.0
        dz = 0.0
        do i = 1,thisnsab
        iIndex = iIndex + 1
        inqPoints(iIndex)%isSeccion = .true.
        inqPoints(iIndex)%center%x = (xini*escalax + dx) + offsetx
        inqPoints(iIndex)%center%z = (zini*escalay + dz) + offsety
        inqPoints(iIndex)%normal%x = nx
        inqPoints(iIndex)%normal%z = nz
         if (abs(inqPoints(iIndex)%normal%x) .lt. 0.0001) inqPoints(iIndex)%normal%x = 0
         if (abs(inqPoints(iIndex)%normal%z) .lt. 0.0001) inqPoints(iIndex)%normal%z = 0
       ! angulo para hacer la rotación a coordenadas normal y tangencial
         r = sqrt(inqPoints(iIndex)%normal%z**2 + inqPoints(iIndex)%normal%x**2)
         inqPoints(iIndex)%cosT = inqPoints(iIndex)%normal%x/r
         inqPoints(iIndex)%sinT = inqPoints(iIndex)%normal%z/r
        ths_layer = thelayeris(inqPoints(iIndex)%center%z)
        ths_isoninterface = tellisoninterface(inqPoints(iIndex)%center%z)
        inqPoints(iIndex)%layer = ths_layer
        inqPoints(iIndex)%isOnInterface = ths_isoninterface
        dx = dx + deltax*escaladx
        dz = dz + deltaz*escalady
        end do
      end do
      nIpts = nIpts + nsecciones
      else
        read(7,*) !Secciones -------------------------
        read(7,*) !; print*,escalax,escalay
        READ(7,*) !; print*, offsetx,offsety
        read(7,*) ! npuntos     Xi        del
      end if !nsecciones

      !
      if (nSabanapts .gt. 0) then
      print*,"Ahora con los puntos de la Sabana"
      allocate(Sabana(nSabanapts, NPTSTIME))
      read(7,*) !Sabanapoints -------------------------
      read(7,*) escalax,escalay !; print*,escalax,escalay
      READ(7,*) offsetx,offsety
      read(7,*) !npuntos xini   deltax   zini   delta z  guardarFK
      do j=1,nnsabanas
      read(7,*) thisnsab,xini,deltax,zini,deltaz,wtfk, sabanaBajarAFrontera
      dx = 0.0
      dz = 0.0
      do i=1,thisnsab
        iIndex = iIndex + 1
        inqPoints(iIndex)%isSabana = .true.
        inqPoints(iIndex)%center%x = (xini + dx)*escalax + offsetx
        inqPoints(iIndex)%center%z = (zini + dz)*escalay

        if (sabanaBajarAFrontera) then
         do k= 1,sabanabajarmax
           cn = adentroOafuera(real(inqPoints(iIndex)%center%x,4), &
                               real(inqPoints(iIndex)%center%z,4),'void')
           if (cn .eqv. .true.) then !esta en el aire
             inqPoints(iIndex)%center%z = inqPoints(iIndex)%center%z + k* SbanadeltaZ
           else
             exit
           end if
         end do
        end if
        inqPoints(iIndex)%center%z = inqPoints(iIndex)%center%z + offsety

        ths_layer = thelayeris(inqPoints(iIndex)%center%z)
        ths_isoninterface = tellisoninterface(inqPoints(iIndex)%center%z)
        inqPoints(iIndex)%layer = ths_layer
        inqPoints(iIndex)%isOnInterface = ths_isoninterface
        dx = dx + deltax
        dz = dz + deltaz
        inqPoints(iIndex)%guardarFK = wtfk
      end do ! i
      end do ! j
      read(7,*) !hacer zeros el rango en la sabana:
      read(7,*) sabZeroini,sabZerofin
      nIpts = nIpts + nSabanapts
      else
       read(7,*) ! Sbana---
       read(7,*) !escala
       read(7,*) !offset
       read(7,*) ! encabezados
       read(7,*) ! hacer zeros
       read(7,*) ! valores
      end if !sabanas


      close(7)

      if (nIpts .lt. size(inqPoints)) then
        allocate(auxInq(nIpts))
        auxInq = inqPoints
        deallocate(inqPoints)
        allocate(inqPoints(nIpts))
        inqPoints = auxInq
        deallocate(auxInq)
      end if
      CALL chdir("..")
      print*,"chido con los puntos de interes"
      end subroutine getInquirePoints

      subroutine getPolaridad(skipdir,PSV,SH)

      logical :: skipdir(3),PSV,SH
      logical :: lexist
      integer :: input
      CALL chdir("ins")
      inquire(file="source.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="source.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "source.txt" on Working directory'
      end if
      READ(77,*);READ(77,*) input
      close(77)
      SH = .false.; PSV = .false.;skipdir = .true.

      if (input .eq. 2) then
      SH = .true.
      skipdir(2) = .false.
      return
      end if
      !
      if (input .eq. 1) then
      PSV = .true.
      skipdir(1) = .false.
      return
      end if
      !
      if (input .eq. 3) then
      PSV = .true.
      skipdir(3) = .false.
      return
      end if
      !
      if (input .eq. 4) then
      PSV = .true.
      skipdir(1) = .false.
      skipdir(3) = .false.
      return
      end if
      end subroutine getPolaridad

      subroutine getsource
      use wavevars, only: t0,maxtime! Escala,Ts,Tp, ampfunction, sigGaus,
      use sourceVars, only: Po, nFuentes!, tipoFuente, PW_pol,
      use glovars, only:pi,PrintNum
      use resultVars, only : Punto
      use Gquadrature, only: Gquad_n
      use soilVars, only : Z,N
      implicit none
      interface
        subroutine punGa (BP)
        use resultVars, only : Punto
        type (Punto), pointer :: BP
        end subroutine punGa
      end interface

      integer :: thelayeris,efsource,i,regi
      logical :: lexist, tellisoninterface, intfsource
      real    :: xfsource,zfsource,l
      real*8  :: nxfsource,nyfsource,nzfsource, PW_theta
      type (Punto), pointer :: BPi
        real :: Escala
        real :: mastimere,Ts,Tp
        integer :: ampfunction
        real :: sigGaus
        integer :: PW_pol
        integer :: tipofuente


      CALL chdir("ins")
      inquire(file="source.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="source.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "source.txt" on Working directory'
      end if
      READ(77,*);READ(77,*);READ(77,*);
      READ(77,*) nFuentes;READ(77,*);READ(77,*)

      allocate(Po(nFuentes))
      allocate(maxtime(nFuentes))
      write(Printnum,'(a)') &
       "---------------------------------------------------------------------------------"
      do i=1,nFuentes
       READ(77,*) xfsource, zfsource, nxfsource,&
                 nyfsource, nzfsource, PW_theta, l, regi,&
                 Escala, ampfunction,  mastimere, Ts, Tp, sigGaus,&
                 PW_pol, tipoFuente
       Po(i)%center%x = xfsource
       Po(i)%center%z = zfsource
       Po(i)%normal%x = nxfsource
!      Po(i)%normal%y = nyfsource
       Po(i)%normal%z = nzfsource
       Po(i)%cosT = cos((360-PW_theta)*pi/180.0)
       Po(i)%sinT = sin((360-PW_theta)*pi/180.0)
       Po(i)%gamma = PW_theta*pi/180.0 !clockwise desde el eje z (hacia abajo)
       Po(i)%length = l
       if (tipoFuente .eq. 1) then ! onda plana
         Po(i)%center%z = Z(N+1)
         efsource = N+1
         intfsource = .true.
       else
         efsource = thelayeris(real(zfsource,8))
         intfsource = tellisoninterface(real(zfsource,8))
       end if
       Po(i)%region = regi
       if (regi .eq. 2) then
         efsource = N+2
         intfsource = .false.
       end if
       Po(i)%layer = efsource
       Po(i)%isOnInterface = intfsource
       Po(i)%isBoundary = .false.
       if (tipoFuente .eq. 2) then ! fuente segmento
        Po(i)%bord_A%x = xfsource - l * 0.5 * Po(i)%cosT
        Po(i)%bord_A%z = zfsource - l * 0.5 * Po(i)%sinT
        Po(i)%bord_B%x = xfsource + l * 0.5 * Po(i)%cosT
        Po(i)%bord_B%z = zfsource + l * 0.5 * Po(i)%sinT
        Po(i)%isSourceSegmentForce = .true.
        BPi => Po(i)
        allocate(BPi%Gq_xXx_coords(Gquad_n,2))
        allocate(BPi%Gq_xXx_C(Gquad_n))
        call punGa(BPi)
       else
        Po(i)%isSourceSegmentForce = .false.
        Po(i)%length = 1.0_8
       end if
       Po(i)%Escala=Escala
       Po(i)%ampfunction=ampfunction
       maxtime(i) = mastimere
       Po(i)%Ts=Ts
       Po(i)%Tp=Tp
       Po(i)%sigGaus=sigGaus
       Po(i)%PW_pol=PW_pol
       Po(i)%tipoFuente=tipoFuente
      write(PrintNum,'(/,a,F8.2,a,F8.2,a,2x,a,F9.2,a,F9.2,a,I0,a,I0,a,L2)') &
      "   Source: (",Po(i)%center%x,",",Po(i)%center%z,")", &
      "n=[",Po(i)%normal%x,",",Po(i)%normal%z,&
      "] r= ",Po(i)%region," e=",Po(i)%layer," intf=",Po(i)%isOnInterface
      write(PrintNum,'(a,F8.2,a,F8.2)') "   cosT=",Po(i)%cosT," sinT=",Po(i)%sinT
      if (tipoFuente .eq. 2) then ! fuente segmento
        write(PrintNum,'(a,F8.2,F8.2,a,F8.2,F8.2,a)') &
        "   A=(",Po(i)%bord_A%x,Po(i)%bord_A%z,")   B=(",Po(i)%bord_B%x,Po(i)%bord_B%z,")"
        write(PrintNum,'(a,F8.2)') "   Clockwise angle from Z: ",Po(i)%gamma
      end if
      end do
      READ(77,*);READ(77,*) t0
      close(77); CALL chdir("..")

      end subroutine getsource

      subroutine getVideoPoints
      use resultVars, only : moviePoints, nMpts, &
                             iPtfin,mPtini,mPtfin
      use peli, only : coords_Z,coords_X,fotogramas,fotogramas_Region
      use meshVars
      use gloVars,only: makevideo,z0
      use waveNumVars, only : NPTSTIME
      use sourceVars, only : nFuentes
      implicit none
      integer :: iz,thelayeris
      real :: firstZ
      logical :: tellisoninterface,lexist
      CALL chdir("ins")
      inquire(file="videoParameters.txt",exist=lexist)
      if (lexist) then
        OPEN(7,FILE="videoParameters.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "videoParameters.txt" on Working directory'
      end if

      READ(7,*) !      Resolution (vertical, horizontal)
      READ(7,*) npixZ,npixX,nmarkZ,nmarkX
      READ(7,*) !      Vertical (start,end)
      READ(7,*) MeshMinZ,MeshMaxZ,firstZ
      READ(7,*) !      Horizontal (start,end)
      READ(7,*) MeshMinX,MeshMaxX
      READ(7,*) MeshVecLen
      READ(7,*) MeshVecLen2
      close(7)
      CALL chdir("..")
      if (makeVideo .eqv. .false.) return
        allocate(coords_X(npixX))
        allocate(coords_Z(npixZ))
        allocate(fotogramas(npixZ,npixX,NPTSTIME,3,nFuentes))
        allocate(fotogramas_Region(npixZ,npixX))

        fotogramas = Z0
        fotogramas_Region = 1!'estr' ! la region exterior
        nMpts = npixZ
        mPtini = iPtfin + 1
        mPtfin = mPtini + nMpts - 1
        allocate(moviePoints(nMpts))
        moviePoints(:)%isBoundary = .false.
        moviePoints(:)%isOnInterface = .false.
        moviePoints(:)%guardarFK = .false.
        moviePoints(:)%guardarMovieSiblings = .true.
        moviePoints(:)%region = 1!'estr'
        do iz = 1, npixZ
          moviePoints(iz)%center%x = MeshMinX
          moviePoints(iz)%center%z = firstZ + (MeshMaxZ - firstZ)/(npixZ-1) * (iz-1)
          coords_Z(iz) = real(moviePoints(iz)%center%z,4)
          moviePoints(iz)%layer = thelayeris(real(moviePoints(iz)%center%z,8))
          moviePoints(iz)%isOnInterface = tellisoninterface(real(moviePoints(iz)%center%z,8))
        end do!
        do iz = 1, npixX
          coords_X(iz) = real(MeshMinX + (MeshMaxX - MeshMinX)/(npixX-1) * (iz-1),4)
        end do
      end subroutine getVideoPoints

      subroutine getTopography
      !Read coordinates of collocation points and fix if there are
      !intersections with inferfaces. Also find normal vectors of segments.
      use GeometryVars, only: n_topo,n_cont,n_vall,nXI, Xcoord_ER, &
      Xcoord_Voidonly, Xcoord_Incluonly,Xcoord_flip_out,boxIncl_maxX,boxIncl_maxY, &
      boxIncl_minX,boxIncl_minY,boxVoid_maxX,boxVoid_maxY,boxVoid_minX,&
      boxVoid_minY, midPoint,normXI, origGeom, surf_poly,&
      N_de_regdionesR, N_de_segmentosR,Xcoord_Incluonly_e
      use gloVars
      use fitting
      use soilVars, only : Z,N,RHO,BETA0,ALFA0,shadecolor_inc,LAMBDA0,ANU,AMU0, layershadecolor
      use ploteo10pesos

      implicit none
      logical :: lexist, huboCambios
      real*8, dimension(:,:,:), allocatable :: auxVector
      real*8 :: l, m,BEALF
      integer :: iXI,e
      real*8 :: nuevoPx, escalax,escalay,offsetx,offsety,escala_n,minBeta,maxBeta

      real     :: errT = 0.0001
      logical, dimension(:), allocatable  :: isOnIF
      real, dimension(:), allocatable :: auxA,auxB
      CHARACTER(len=32) :: arg
      real*8, allocatable, save :: lengthXI(:),layerXI(:),cost(:),sint(:)
      logical, dimension(:), allocatable :: es_de_esquina
      integer :: nXIoriginal
      integer :: iR
      huboCambios = .false.
      allocate(auxA(1)); allocate(auxB(1))
      print*,"Reading topography file"
      CALL chdir("ins")
      !Read Surface topography
      inquire(file="boundaries.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="boundaries.txt",FORM="FORMATTED")
      else
        write(6,'(a)')'There is a missing input file. '
        stop 'Check "boundaries.txt" on Working directory'
      end if
      CALL get_command_argument(1, arg)
      IF ((LEN_TRIM(arg) .ne. 0) .and. (trim(arg) .eq. '-g')) verbose = 3
      READ(77,*)
      READ(77,*) n_topo
      READ(77,*) n_cont
      READ(77,*) n_vall
      nXI = n_topo + n_cont + n_vall !numero total de segmentos originales
      if (allocated(Xcoord_ER)) deallocate(Xcoord_ER)
      ALLOCATE (Xcoord_ER(nXI+300,2,2))
      allocate (es_de_esquina(nXI*2+300))
      allocate(auxVector(nXI+1+300,2,2))
      READ(77,*)
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
      READ(77,*) escala_n
      DO iXI = 1,n_topo
         if (escala_n .gt. 0.0) then
         READ(77,*) Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),&
         Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),es_de_esquina(iXI)
         else
         READ(77,*) Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),&
         Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),es_de_esquina(iXI)
         end if
         Xcoord_ER(iXI,1,1:2) = Xcoord_ER(iXI,1,1:2) * escalax + offsetx
         Xcoord_ER(iXI,2,1:2) = Xcoord_ER(iXI,2,1:2) * escalay + offsety
      END DO
      READ(77,*)
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
      READ(77,*) escala_n
      DO iXI = n_topo+1, n_topo+n_cont
         if (escala_n .gt. 0.0) then
         READ(77,*) Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),&
         Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),es_de_esquina(iXI)
         else
         READ(77,*) Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),&
         Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),es_de_esquina(iXI)
         end if
         Xcoord_ER(iXI,1,1:2) = Xcoord_ER(iXI,1,1:2) * escalax + offsetx
         Xcoord_ER(iXI,2,1:2) = Xcoord_ER(iXI,2,1:2) * escalay + offsety
      END DO
      READ(77,*)
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
      READ(77,*) escala_n
      DO iXI = n_topo+n_cont+1,n_topo+n_cont+n_vall !nXI
         if (escala_n .gt. 0.0) then
         READ(77,*) Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),&
         Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),es_de_esquina(iXI)
         else
         READ(77,*) Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),&
         Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1),es_de_esquina(iXI)
         end if
         Xcoord_ER(iXI,1,1:2) = Xcoord_ER(iXI,1,1:2) * escalax + offsetx
         Xcoord_ER(iXI,2,1:2) = Xcoord_ER(iXI,2,1:2) * escalay + offsety
      END DO
      READ(77,*) ! Inclusion properties
      READ(77,*) ! Alfa(m/s) Beta(m/s) Dens(T/m3)
      READ(77,*) ALFA0(N+2),BETA0(N+2),RHO(N+2)

      AMU0(N+2)=RHO(N+2)*BETA0(N+2)**2.0
      LAMBDA0(N+2)=RHO(N+2)*ALFA0(N+2)**2.0 - &
                        real(2.)*real(AMU0(N+2))
      BEALF = BETA0(N+2)/ALFA0(N+2)
      ANU(N+2) = (BEALF**2 - 0.5)/(-1 + BEALF**2)

      if (allocated(Xcoord_Incluonly)) deallocate(Xcoord_Incluonly)
      allocate (Xcoord_Incluonly(n_cont + n_vall,2,2))
      do iXI = 1,n_cont + n_vall
          Xcoord_Incluonly(iXI,:,:) = Xcoord_ER(n_topo+iXI,:,:)
          if (abs(Xcoord_Incluonly(iXI,2,1)) .le. 0.06) Xcoord_Incluonly(iXI,2,1) = 0
          if (abs(Xcoord_Incluonly(iXI,2,2)) .le. 0.06) Xcoord_Incluonly(iXI,2,2) = 0
      end do

      ! bounding box
      boxIncl_maxX = maxval(Xcoord_Incluonly(:,1,:))
      boxIncl_maxY = maxval(Xcoord_Incluonly(:,2,:))
      boxIncl_minX = minval(Xcoord_Incluonly(:,1,:))
      boxIncl_minY = minval(Xcoord_Incluonly(:,2,:))
      READ(77,*) ! and make anything inside is on the air
      READ(77,*) e
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
      if (e .ne. 0) then
      if (allocated(Xcoord_Voidonly)) deallocate(Xcoord_Voidonly)
      allocate (Xcoord_Voidonly(e,2,2))
      do iXI = 1,e
         READ(77,*) Xcoord_Voidonly(iXI,1,1),Xcoord_Voidonly(iXI,2,1),&
                    Xcoord_Voidonly(iXI,1,2),Xcoord_Voidonly(iXI,2,2)
         Xcoord_Voidonly(iXI,1,1:2) = Xcoord_Voidonly(iXI,1,1:2) * escalax + offsetx
         Xcoord_Voidonly(iXI,2,1:2) = Xcoord_Voidonly(iXI,2,1:2) * escalay + offsety
      end do
      ! bounding box
      boxVoid_maxX = maxval(Xcoord_Voidonly(:,1,:))
      boxVoid_maxY = maxval(Xcoord_Voidonly(:,2,:))
      boxVoid_minX = minval(Xcoord_Voidonly(:,1,:))
      boxVoid_minY = minval(Xcoord_Voidonly(:,2,:))
      end if

      READ(77,*) !Cualquier receptor dentro de estas superficies está en la region homogenea
      READ(77,*) N_de_regdionesR !N de regiones
      READ(77,*) N_de_segmentosR !Total de segmentos
      if (N_de_regdionesR .ne. 0) then
      if (allocated(Xcoord_Incluonly)) deallocate(Xcoord_Incluonly)
      allocate (Xcoord_Incluonly(N_de_segmentosR,2,2))
      if (allocated(Xcoord_Incluonly_e)) deallocate(Xcoord_Incluonly_e)
      allocate (Xcoord_Incluonly_e(N_de_regdionesR+1))
      Xcoord_Incluonly_e(1) = 1
      do iR = 1,N_de_regdionesR
      READ(77,*) Xcoord_Incluonly_e(iR+1)
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
!     print*,"from ",sum(Xcoord_Incluonly_e(1:iR))," to ",sum(Xcoord_Incluonly_e(2:iR+1))
      do iXI = sum(Xcoord_Incluonly_e(1:iR)),sum(Xcoord_Incluonly_e(2:iR+1))
         READ(77,*) Xcoord_Incluonly(iXI,1,1), Xcoord_Incluonly(iXI,2,1),&
                    Xcoord_Incluonly(iXI,1,2), Xcoord_Incluonly(iXI,2,2)

         Xcoord_Incluonly(iXI,1,1:2) = Xcoord_Incluonly(iXI,1,1:2) * escalax + offsetx
         Xcoord_Incluonly(iXI,2,1:2) = Xcoord_Incluonly(iXI,2,1:2) * escalay + offsety
      end do
      end do
!     do ixi = 1, N_de_segmentosR; print*,ixi,Xcoord_Incluonly(iXI,1,1), Xcoord_Incluonly(iXI,2,1),&
!        Xcoord_Incluonly(iXI,1,2), Xcoord_Incluonly(iXI,2,2); end do
!     do ixi = 1, N_de_regdionesR+1; print*,ixi,Xcoord_Incluonly_e(ixi)
!     end do
!     stop
      ! bounding box
      boxIncl_maxX = maxval(Xcoord_Incluonly(:,1,:))
      boxIncl_maxY = maxval(Xcoord_Incluonly(:,2,:))
      boxIncl_minX = minval(Xcoord_Incluonly(:,1,:))
      boxIncl_minY = minval(Xcoord_Incluonly(:,2,:))
      else
       READ(77,*)
       READ(77,*)
       READ(77,*)
      end if
      !
      if (flip12) then
      READ(77,*) ! Borde de región E cuando flip12 .true.
      READ(77,*) e
      READ(77,*) escalax,escalay
      READ(77,*) offsetx,offsety
      if (e .ne. 0) then
      if (allocated(Xcoord_flip_out)) deallocate(Xcoord_flip_out)
      allocate (Xcoord_flip_out(e,2,2))
      do iXI = 1,e
         READ(77,*) Xcoord_flip_out(iXI,1,1),Xcoord_flip_out(iXI,2,1),&
                    Xcoord_flip_out(iXI,1,2),Xcoord_flip_out(iXI,2,2)
         Xcoord_flip_out(iXI,1,1:2) = Xcoord_flip_out(iXI,1,1:2) * escalax + offsetx
         Xcoord_flip_out(iXI,2,1:2) = Xcoord_flip_out(iXI,2,1:2) * escalay + offsety
      end do
      end if
      end if
      close(77)
      CALL chdir("..")

      minBeta = minval(BETA0(1:N+2))
      maxBeta = maxval(BETA0(1:N+2))

      if (abs(minBeta - maxBeta) .lt. 0.01) then
        shadecolor_inc = layershadecolor(1)
      else
        e=N+2
!         layershadecolor(e)= real(0.4-(maxBeta-beta0(e))*((0.4-0.89)/(maxBeta-minBeta)),4)
          layershadecolor(e)= real(0.6-(maxBeta-beta0(e))*((0.6-0.89)/(maxBeta-minBeta)),4)
          !print*,e,layershadecolor(e)
!       end do
          shadecolor_inc = layershadecolor(N+2)
      end if

!     go to 384
      nXIoriginal = nXI
      ! cortar en intersección con estratos y determinar estrato de cada segemento
      if (verbose >= 4) write(PrintNum,*)"begin slice with layers"
      iXI = 1
      DO !para cada segmento

      !#< r we check if there is a change of medium along segment iXI  !#>
      DO e = 2,N+1 !en cada interfaz (excepto la superficie)
      if (verbose >= 4) write(PrintNum,*) "e= ",e
        if ((abs(Xcoord_ER(iXI,2,1) - Z(e)) < errT) .or. &
            (abs(Xcoord_ER(iXI,2,2) - Z(e)) < errT)) then
!       if ((abs(anint(Xcoord_ER(iXI,2,1) * 1000) - anint(Z(e) * 1000)) < errT) .or. &
!           (abs(anint(Xcoord_ER(iXI,2,2) * 1000) - anint(Z(e) * 1000)) < errT)) then
            if (verbose >= 4) write(PrintNum,*)"already sliced here"
        else ! not already
          if (verbose >= 4) write(PrintNum,*) "divi"
            ! si es un segmento que se forma recorriendo los puntos hacia Z+
            if (Xcoord_ER(iXI,2,1) < Z(e)  .AND. Xcoord_ER(iXI,2,2) > Z(e)) then
               if (verbose >= 4) write(PrintNum,*)"hay un cambio de estrato hacia z+"
               if (abs(Xcoord_ER(iXI,1,2) - Xcoord_ER(iXI,1,1)) .le. 0.001) then ! es vertical
                 nuevoPx = Xcoord_ER(iXI,1,1)
               else ! es inclinado
                 m = (Xcoord_ER(iXI,2,2) - Xcoord_ER(iXI,2,1))/(Xcoord_ER(iXI,1,2) - Xcoord_ER(iXI,1,1))
                 surf_poly = real((/Xcoord_ER(iXI,2,1) - Xcoord_ER(iXI,1,1)*m , m /),4)
                 !are the polynomial coefficients from lower to top. A0 A1 A2 ... An
                 nuevoPx = splitatY(surf_poly,1,real(Z(e),4),&
                 real(Xcoord_ER(iXI,1,1),4),real(Xcoord_ER(iXI,1,2),4))
               end if !
               if ((sqrt((Xcoord_ER(iXI,1,1)-nuevoPx)**2 + (Xcoord_ER(iXI,2,1) - Z(e))**2) .lt. 0.001) .or. &
                   (sqrt((Xcoord_ER(iXI,1,2)-nuevoPx)**2 + (Xcoord_ER(iXI,2,2) - Z(e))**2) .lt. 0.001)) then
                   if (verbose >= 4) write(PrintNum,*) "muy corto"
                   cycle ! si es muy corto no lo queremos
               end if!
             if (verbose >= 4) then
              write(Printnum,'(a,F10.4,F10.4,a,F10.4,F10.4,a,F10.4,a,F10.4,a)') &
              "(dn)Segment (x,z):[",Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1)," to" &
              ,Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),"] crosses interface at (x,z)=(",nuevoPx,",",Z(e),")"
             end if

               ! insertamos el nuevo punto en el vector de puntos
!              if (allocated(auxVector)) deallocate(auxVector)
!              allocate(auxVector(nXI+1,2,2)) ! una arista nueva
               auxVector(1:iXI,1:2,1:2) = Xcoord_ER(1:iXI,1:2,1:2)
               auxVector(iXI,1,2) = nuevoPx
               auxVector(iXI,2,2) = Z(e)
               auxVector(iXI+1,1,1) = nuevoPx
               auxVector(iXI+1,2,1) = Z(e)
               auxVector(iXI+1,1,2) = Xcoord_ER(iXI,1,2)
               auxVector(iXI+1,2,2) = Xcoord_ER(iXI,2,2)
               auxVector(iXI+2:nXI+1,1:2,1:2) = Xcoord_ER(iXI+1:nXI,1:2,1:2)
!              es_de_esquina(1:iXI) !se queda igual
               es_de_esquina(iXI+2:nXI+1) = es_de_esquina(iXI+1:nXI) !recorre
               es_de_esquina(iXI+1) = es_de_esquina(iXI) ! se duplica

               nXI = nXI + 1
               Xcoord_ER(1:nXI,1:2,1:2) = auxVector(1:nXI,1:2,1:2)

!              if (allocated(Xcoord_ER))  deallocate(Xcoord_ER)
!              nXI = nXI + 1
!              allocate(Xcoord_ER(nXI,2,2))
!              Xcoord_ER(1:nXI,1:2,1:2) = auxVector(1:nXI,1:2,1:2)

               if(n_topo+n_cont+1 .le. iXI) n_vall = n_vall + 1
               if(n_topo+1 .le. iXI .and. iXI .le. n_cont+n_topo) n_cont = n_cont + 1
               if(iXI .le. n_topo) n_topo = n_topo + 1
               if (verbose >= 4) write(PrintNum,*) "fin hacia z+"
               cycle !estrato
            end if !(Xcoord_ER(iXI,2,1) < Z(e)  .AND. Xcoord_ER(iXI,2,2) > Z(e))
      ! si es un segmento que se forma recorriendo los puntos hacia Z-
            if (Xcoord_ER(iXI,2,1) > Z(e)  .AND. Xcoord_ER(iXI,2,2) < Z(e)) then
               if (verbose >= 4) write(PrintNum,*)"hay un cambio de estrato hacia z-"
               if (abs(Xcoord_ER(iXI,1,2) - Xcoord_ER(iXI,1,1)) .le. 0.001) then ! es vertical
                 nuevoPx = Xcoord_ER(iXI,1,1)
               else ! es inclinado
                 m = (Xcoord_ER(iXI,2,2) - Xcoord_ER(iXI,2,1))/(Xcoord_ER(iXI,1,2) - Xcoord_ER(iXI,1,1))
                 surf_poly = real((/Xcoord_ER(iXI,2,1) - Xcoord_ER(iXI,1,1)*m , m /),4)
                 nuevoPx = splitatY(surf_poly,1,real(Z(e),4),&
                 real(Xcoord_ER(iXI,1,1),4),real(Xcoord_ER(iXI,1,2),4))
               end if!
               if ((sqrt((Xcoord_ER(iXI,1,1)-nuevoPx)**2 + (Xcoord_ER(iXI,2,1) - Z(e))**2) .lt. 0.001) .or. &
                   (sqrt((Xcoord_ER(iXI,1,2)-nuevoPx)**2 + (Xcoord_ER(iXI,2,2) - Z(e))**2) .lt. 0.001)) then
                   if (verbose >= 4) write(PrintNum,*) "muy corto"
                   cycle ! si es muy corto no lo queremos
               end if!
             if (verbose >= 4) then
               write(Printnum,'(a,F10.4,F10.4,a,F10.4,F10.4,a,F10.4,a,F10.4,a)') &
              "(dn)Segment (x,z):[",Xcoord_ER(iXI,1,1),Xcoord_ER(iXI,2,1)," to" &
              ,Xcoord_ER(iXI,1,2),Xcoord_ER(iXI,2,2),"] crosses interface at (x,z)=(",nuevoPx,",",Z(e),")"
             end if

               ! insertamos el nuevo punto en el vector de puntos
!              if (allocated(auxVector)) deallocate(auxVector)
!              allocate(auxVector(nXI+1,2,2)) ! un segmento más
               auxVector(1:iXI,1:2,1:2) = Xcoord_ER(1:iXI,1:2,1:2)
               auxVector(iXI,1,2) = nuevoPx
               auxVector(iXI,2,2) = Z(e)
               auxVector(iXI+1,1,1) = nuevoPx
               auxVector(iXI+1,2,1) = Z(e)
               auxVector(iXI+1,1,2) = Xcoord_ER(iXI,1,2)
               auxVector(iXI+1,2,2) = Xcoord_ER(iXI,2,2)
               auxVector(iXI+2:nXI+1,1:2,1:2) = Xcoord_ER(iXI+1:nXI,1:2,1:2)
!              es_de_esquina(1:iXI) !se queda igual
               es_de_esquina(iXI+2:nXI+1) = es_de_esquina(iXI+1:nXI) ! se recorre
               es_de_esquina(iXI+1) = es_de_esquina(iXI) ! se duplica

               nXI = nXI + 1
               Xcoord_ER(1:nXI,1:2,1:2) = auxVector(1:nXI,1:2,1:2)


!              if (allocated(Xcoord_ER)) deallocate(Xcoord_ER)
!              nXI = nXI + 1
!              allocate(Xcoord_ER(nXI,2,2))
!              Xcoord_ER(1:nXI,1:2,1:2) = auxVector(1:nXI,1:2,1:2)

               if(n_topo+n_cont+1 .le. iXI) n_vall = n_vall + 1
               if(n_topo+1 .le. iXI .and. iXI .le. n_cont+n_topo) n_cont = n_cont + 1
               if(iXI .le. n_topo) n_topo = n_topo + 1
               if (verbose >= 4) write(PrintNum,*) "fin hacia z-"
               cycle !estrato
            end if !(Xcoord_ER(iXI,2,1) > Z(e)  .AND. Xcoord_ER(iXI,2,2) < Z(e))
        end if ! already
        if (verbose >= 4) write(PrintNum,*) "fin e=",e
      end do ! e
      !n_topo+n_cont+1,n_topo+n_cont+n_vall
      if (iXI .eq. n_topo+n_cont) exit ! salidr de Do de iXI
!     if (iXI .eq. nXI) exit ! salidr de Do de iXI
      ixI = iXI+1
      end do !iXI

      deallocate(auxVector)
      if (verbose >= 4) write(PrintNum,*)"out of slice"
      if (nXI .le. 0) stop "errror nXI <= 0"

      ! Right hand normals
      ALLOCATE (normXI(nXI,2))
      do iXI=1,nxi
        normXI(ixi,1) = Xcoord_ER(iXI,2,1) - Xcoord_ER(iXI,2,2)
        normXI(ixi,2) = Xcoord_ER(iXI,1,2) - Xcoord_ER(iXI,1,1)
        l = sqrt(normxi(ixi,1)**2. + normxi(ixi,2)**2.)
        normxi(ixi,1) = normxi(ixi,1)/ l
        normxi(ixi,2) = normxi(ixi,2)/ l
      end do

      ! midpoints and lenght
      allocate (midPoint(nXI,2))
      allocate (lengthXI(nXI))
      allocate (cost(nXI))
      allocate (sint(nXI))
      do iXI=1,nxi
        midPoint(ixi,1) = (Xcoord_ER(iXI,1,2) + Xcoord_ER(iXI,1,1))/2.0
        midPoint(ixi,2) = (Xcoord_ER(iXI,2,2) + Xcoord_ER(iXI,2,1))/2.0
        lengthXI(ixi) = sqrt((Xcoord_ER(iXI,1,2)-Xcoord_ER(iXI,1,1))**2.0 &
                             +(Xcoord_ER(iXI,2,2)-Xcoord_ER(iXI,2,1))**2.0)
      end do
      cost = (Xcoord_ER(1:nXi,1,2)- midPoint(1:nXi,1))/(lengthXI(1:nXi)/2.)
      sint = (Xcoord_ER(1:nXi,2,2)- midPoint(1:nXi,2))/(lengthXI(1:nXi)/2.)

      ! estrato
      allocate (layerXI(nXI))
      layerxi(:) = N+1
      do e=1,N
        forall(ixi=1:nxi, &
          (z(e) < midpoint(ixi,2) .and. midpoint(ixi,2) <= z(e+1))) &
          layerxi(ixi) = e
      end do

      allocate (isOnIF(nXI)) ; isOnIF = .false.
      do e=1,N+1
        forall(ixi=1:nxi, &
          ((Z(e)-errT .lt. midPoint(iXI,2)) &
            .and. (midPoint(iXI,2) .lt. Z(e)+errT))) &
          isOnIF(iXI) = .true.
      end do

      ! save original Geometry for the posterity
      allocate(origGeom(nXI))
      ! central point
      origGeom(1:nXi)%center%x = midPoint(1:nXi,1)
      origGeom(1:nXi)%center%z = midPoint(1:nXi,2)
      !borders of segment
      origGeom(1:nXi)%bord_A%x = Xcoord_ER(1:nXi,1,1)
      origGeom(1:nXi)%bord_A%z = Xcoord_ER(1:nXi,2,1)
      origGeom(1:nXi)%bord_B%x = Xcoord_ER(1:nXi,1,2)
      origGeom(1:nXi)%bord_B%z = Xcoord_ER(1:nXi,2,2)
      !add normal
      origGeom(1:nXi)%normal%x = normXI(1:nXi,1)
      origGeom(1:nXi)%normal%z = normXI(1:nXi,2)
      !add length
      origGeom(1:nXi)%length = lengthXI(1:nXi)
      origGeom(1:nXi)%segmentoDeEsquina = es_de_esquina(1:nXI)
      origGeom(1:nXi)%cosT = cost(1:nXi)
      origGeom(1:nXi)%sinT = sint(1:nXi)
      !add layer
      origGeom(1:nXi)%layer = int(layerXI(1:nXi))
      origGeom(1:nXi)%isBoundary = .true.
      origGeom(1:nXi)%isOnInterface = isOnIF(1:nXi)
      origGeom(1:nXi)%guardarFK = .false.
      origGeom(1:nXi)%guardarMovieSiblings = .false.

      !tipo de frontera
      !  TE^0 + TE^d = 0
      origGeom(1:n_topo)%tipoFrontera = 0
      !  TE^0 + TE^d = TR^r; uE^0 + uE^d = uR^r
      origGeom(n_topo+1:n_cont+n_topo)%tipoFrontera = 1
      !  TR^r = 0
      origGeom(n_cont+n_topo+1:nXI)%tipoFrontera = 2

      write(Printnum,'(A, F7.1,2x, F7.1,2x, E8.2,2x, E8.2,2x, E8.2,2x, E8.2)') &
         'inclusion  ',ALFA0(N+2),BETA0(N+2),real(AMU0(N+2)),&
         RHO(N+2), real(LAMBDA0(N+2)),real(anu(N+2))

        if (verbose .ge. 3) then
         DO iXI = 1,nXI
        write(Printnum,'(I0,a,F5.2,a,F5.2,a,F5.2,a,F5.2,a,F5.2)', ADVANCE = "NO") iXI,&
         " c[",origGeom(iXI)%center%x,",",origGeom(iXI)%center%z,"] n[",&
         origGeom(iXI)%normal%x,",",origGeom(iXI)%normal%z,"] l:",origGeom(iXI)%length
        write(Printnum,'(a,F5.2,a,F5.2,a,i0,a,i0)') &
        "  co",origGeom(iXI)%cosT," si",origGeom(iXI)%sinT," e",origGeom(iXI)%layer, &
        " t",origGeom(iXI)%tipoFrontera
         end do
        end if
      deallocate(auxA);deallocate(auxB);deallocate(es_de_esquina);deallocate(isonif)
      deallocate(lengthXI);deallocate(layerXI);deallocate(cost);deallocate(sint)
      end subroutine getTopography

      subroutine setVideoPointsRegions
      use peli, only : coords_Z,coords_X,fotogramas_Region
      use meshVars, only : npixX,npixZ
      use glovars, only : flip12,PrintNum,verbose
      use debugstuff
      implicit none
      integer :: i,j
      integer, dimension(0:2) :: reg
      logical :: cn, adentroOafuera
      ! !1 estr, 2 incl, 0 void
      reg(0) = 0; reg(1)= 1; reg(2) = 2
      if (flip12) then
      reg(0) = 0; reg(1)= 2; reg(2) = 1
      end if
      ! asigmar la region a cada pixel de la pelicula
       do i=1,npixZ
        do j=1,npixX
          fotogramas_Region(i,j) = reg(1)!'estr'
          cn = adentroOafuera(coords_X(j), coords_Z(i),'void')
          if (cn .eqv. .true.) then
            fotogramas_Region(i,j) = reg(0)!'void'
            cycle
          end if
          cn = adentroOafuera(coords_X(j), coords_Z(i),'incl')
          if (cn .eqv. .true.) then
            fotogramas_Region(i,j) = reg(2)!'incl'
          end if
        end do
       end do
      if (verbose .ge. 2) then
      print*,"";print*,"fotogramas_Region"
      call showMNmatrixI(npixZ,npixX,fotogramas_Region,"fotRe",PrintNum)
      end if

      end subroutine setVideoPointsRegions

      subroutine setInqPointsRegions
      use resultVars, only : allpoints,nPts,n_OD,nIpts
      use soilVars, only : N
      use glovars, only : verbose, Printnum, flip12
      implicit none
      integer :: i
      integer, dimension(0:2) :: reg
      logical :: cn,adentroOafuera
      !! 0 void, 1 estr, 2 incl
      reg(0) = 0; reg(1)= 1; reg(2) = 2
      if (flip12) then
        reg(0) = 0; reg(1)= 2; reg(2) = 1
      end if
!     allpoints(1:nPts)%region = reg(1)!'estr'
      if (verbose .ge. 1) then
      write(PrintNum,'(a)') "------------------------------------------------"
      write(PrintNum,'(a,I0)') "nPts=",nPts
      write(PrintNum,'(a,I0)') "nIpts=",nIpts
      write(PrintNum,'(a,I0,a)') "There are ",n_OD," points to overdetermine the system"
      write(PrintNum,*) "center,region,layer,isOD"
      end if
      do i = 1, nPts
        if (allpoints(i)%atBou .eqv. .false.) then
        allpoints(i)%region = reg(1)!'estr'
        cn = adentroOafuera(real(allpoints(i)%center%x,4), &
                            real(allpoints(i)%center%z,4),'void')
!       print*,i,cn;cycle
        if (cn .eqv. .true.) then
            allpoints(i)%region = reg(0)!'void'
            cycle
        end if
        cn = adentroOafuera(real(allpoints(i)%center%x,4), &
                            real(allpoints(i)%center%z,4),'incl')
        if (cn .eqv. .true.) then
            allpoints(i)%region = reg(2)!'incl'
        end if
!       if (abs(z(0)) .gt. 0.0001) then ! si no hay un semiespacio arriba
!         if (real(allpoints(i)%center%z,4) .lt. 0.0) then
!           allpoints(i)%region = reg(0)!'void'
!         end if
!       end if!
        end if!
        !if (allpoints(i)%region .eq. 2) allpoints(i)%layer = N+2
        if (verbose .ge. 1) then
        write(PrintNum,*)i,"[",allpoints(i)%center%x, &
        ",",allpoints(i)%center%z,"] is ",allpoints(i)%region,allpoints(i)%layer,allpoints(i)%isOD
        end if
      end do!;stop "setInqPointsRegions"
      end subroutine setInqPointsRegions
      end module setupmodel

      function thelayeris(zi)
      use soilVars, only : Z,N
      implicit none
      integer :: i,e,thelayeris
      real*8 :: zi
      i = 1
      if (Z(0) .lt. 0.0) i = 0
      do e=i,N
         if(real(zi) .lt. Z(e+1)) then
            exit
         end if
      end do
      thelayeris = e
      end function thelayeris

      function tellisoninterface(zi)
      use soilVars, only : Z,N
      implicit none
      real*8 ::  errT = 0.0001 !0.01_8
      integer :: e
      real*8 :: zi
      logical :: tellisoninterface
      tellisoninterface = .false.
      do e=1,N+1
        if((Z(e)-errT .lt. real(zi)) &
             .and. (real(zi) .lt. Z(e)+errT)) then
           tellisoninterface = .true.
        end if
      end do
      end function tellisoninterface



      function adentroOafuera(x,y,region)
       use geometryvars, only : Xcoord_Voidonly,Xcoord_Incluonly, &
            boxIncl_maxX,boxIncl_maxY,boxIncl_minX,boxIncl_minY, &
            boxVoid_maxX,boxVoid_maxY,boxVoid_minX,boxVoid_minY, &
            N_de_regdionesR,Xcoord_Incluonly_e

       implicit none
       logical :: adentroOafuera, crossingNumber
       character(LEN=4) :: region
       real*4,intent(in) :: x,y
!      real*16 :: vt
       integer :: nXI,iR!cn,k,
       real*8 :: XcooBox_maxX,XcooBox_minX,XcooBox_maxY,XcooBox_minY
       real*8, dimension(:,:,:), pointer :: Xcoord_ER
       nullify(Xcoord_ER)
       !print*,x,y,region
       adentroOafuera = .false.
       if (region .eq. 'void') then
       XcooBox_maxX = boxVoid_maxX
       XcooBox_maxY = boxVoid_maxY
       XcooBox_minX = boxVoid_minX
       XcooBox_minY = boxVoid_minY
       if (allocated( Xcoord_Voidonly)) then
          Xcoord_ER => Xcoord_Voidonly
          nXI = size(Xcoord_ER(:,1,1))
          if (XcooBox_minX .le. x) then
          if (x .le. XcooBox_maxX) then
          if (XcooBox_minY .le. y) then
          if (y .le. XcooBox_maxY) then
              adentroOafuera  = crossingNumber(x,y,nXI,Xcoord_ER)
          end if
          end if
          end if
          end if
       else
       return
       end if
       else if (region .eq. 'incl') then
       XcooBox_maxX = boxIncl_maxX
       XcooBox_maxY = boxIncl_maxY
       XcooBox_minX = boxIncl_minX
       XcooBox_minY = boxIncl_minY
       if (allocated( Xcoord_Incluonly)) then
       if (N_de_regdionesR .ne. 0) then
       do iR = 1, N_de_regdionesR
          Xcoord_ER => &
          Xcoord_Incluonly(sum(Xcoord_Incluonly_e(1:iR)): sum(Xcoord_Incluonly_e(2:iR+1)),1:2,1:2)
!         print*,sum(Xcoord_Incluonly_e(1:iR)),sum(Xcoord_Incluonly_e(2:iR+1)),x,y
          XcooBox_maxX = maxval(Xcoord_ER(:,1,:))
          XcooBox_maxY = maxval(Xcoord_ER(:,2,:))
          XcooBox_minX = minval(Xcoord_ER(:,1,:))
          XcooBox_minY = minval(Xcoord_ER(:,2,:))
          nXI = size(Xcoord_ER(:,1,1))

          if (XcooBox_minX .le. x) then
          if (x .le. XcooBox_maxX) then
          if (XcooBox_minY .le. y) then
          if (y .le. XcooBox_maxY) then
              adentroOafuera  = crossingNumber(x,y,nXI,Xcoord_ER)
              if (adentroOafuera) return
          end if
          end if
          end if
          end if

       end do
       else
       return
!      Xcoord_ER => Xcoord_Incluonly
       end if
       else
       return
       end if
       else
       return
       end if
!      print*,XcooBox_maxX,XcooBox_maxY,XcooBox_minX,XcooBox_minY
      end function adentroOafuera

       function crossingNumber(x,y,nXI,Xcoord_ER)
       implicit none
       logical :: crossingNumber
       real*4,intent(in) :: x,y
       real*16 :: vt
       integer :: cn,k,nXI
       real*8, dimension(nXI,2,2) :: Xcoord_ER
       crossingNumber = .false.
         cn = 0
          do k=1,nXI
            if (Xcoord_ER(k,2,1) .le. y) then
            if (Xcoord_ER(k,2,2) .gt. y) then
            vt = ( y -  Xcoord_ER(k,2,1)) / &
                 ( Xcoord_ER(k,2,2) - Xcoord_ER(k,2,1))
            if (x .lt. Xcoord_ER(k,1,1) + vt * (Xcoord_ER(k,1,2) -  Xcoord_ER(k,1,1))) then
              cn = cn + 1
              cycle
            end if
            end if
            end if
            !
            if (Xcoord_ER(k,2,1) .gt. y) then
            if (Xcoord_ER(k,2,2) .le. y) then
            vt = ( y -  Xcoord_ER(k,2,1)) / &
                 ( Xcoord_ER(k,2,2) - Xcoord_ER(k,2,1))
            if (x .lt. Xcoord_ER(k,1,1) + vt * (Xcoord_ER(k,1,2) -  Xcoord_ER(k,1,1))) then
              cn = cn + 1
              cycle
            end if
            end if
            end if

            ! crossing number:
            !  http://geomalgorithms.com/a03-_inclusion.html
          end do !k
          if (cn .eq. 1 .or. cn .eq. 3 .or. cn .eq. 5 .or. cn .eq. 7) then
            crossingNumber = .true. !el punto pertenece a la ragion de adentro
          end if

      end function crossingNumber
