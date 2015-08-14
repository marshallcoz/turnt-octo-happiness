      PROGRAM DWNLIBEM 
      use gloVars
      use refSolMatrixVars, only : Ak,B,BparaGa,BparaNu,CoefparGa, CoefparNu,subMatD0,subMatS0
      use waveNumVars
      use soilVars, only : alfa,beta, alfa0,beta0,minbeta,n,z,amu,lambda,rho,qq
      use debugStuff
      use resultVars
      use ploteo10pesos
      use sourceVars, only : nFuentes,currentiFte,Po,psv,sh
      use MeshVars, only : npixX
      use dislin
      use waveVars, only : t0,Uo
      use peli, only : fotogramas,fotogramas_region
      use GeometryVars, only : longitudcaracteristica_a
      use cine
      use setupmodel
      implicit none
      interface
        include 'interfaz.f'
        subroutine diffField_at_iz(i_zF,dir_j,J,cOME)
          integer, intent(in) :: i_zF,dir_j,J
          complex*16, intent(in),target  :: cOME
        end subroutine diffField_at_iz
        
        subroutine intrplr_gloMat(k0,n,pt_cOME_i,pt_ipivA,pt_workA)
          use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
          use waveNumVars, only : k_vec, NMAX
          use fitting
          implicit none
          integer :: k0,n
          complex*16, pointer :: pt_cOME_i
          integer, dimension(:), pointer :: pt_ipivA
          complex*16, dimension(:),pointer :: pt_workA
        end subroutine intrplr_gloMat
        
        subroutine drawBoundary(BP, nbpts, titleN, extension, & 
         zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente)
          use resultVars, only : Punto
          type (Punto), dimension(:), pointer :: BP
          integer, intent(in) :: nbpts
          character(LEN=100) :: titleN
          character(LEN=3) :: extension
          logical, intent(in) :: zoomGeom, plotReceptoresA, & 
          plotReceptoresB,plotFuente
        end subroutine drawBoundary
        
        subroutine G0estr(MecElem,p_x,J,cOME_in,dir)
          use resultVars, only : Punto,FFres
          complex*16, dimension(1:5) :: MecElem
          type(Punto), pointer :: p_X
          integer, intent(in) :: J,dir
          complex*16, intent(in),target  :: cOME_in
        end subroutine G0estr
      end interface
      integer :: J,l,iP,ik,iz,i,ip_x,ipxi,m
      integer :: dir,iPhi,iPhi_I,iPhi_F,ipxi_I,ipxi_F
      character(LEN=100) :: titleN,yax,CTIT,tt
      character(LEN=9)  :: xAx  
      CHARACTER(len=32) :: arg
      character(LEN=3) :: extension
      logical :: skipdir(3)
      logical :: thereisavirtualsourceat,onlythisJ
      real, dimension(2) :: tarray
      real :: result, lastresult
      complex*16, dimension(:),allocatable :: auxGvector
      type (Punto), dimension(:), pointer :: BP
      integer, dimension(:),allocatable,target :: ipivA
      integer, dimension(:),pointer :: pt_ipivA
      complex*16, dimension(:),allocatable,target :: workA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16, dimension(:,:), pointer :: pointAp
      real*8, pointer :: pt_k
      complex*16, pointer :: pt_cOME_i
      integer :: frecIni, frecEnd,tam
      integer :: info,k0
      integer :: Mi,Ni
      complex*16,dimension(:,:),allocatable :: AF
      complex*16,dimension(:,:),allocatable :: Xbem
      complex*16,dimension(:),allocatable :: WORKbem
      real*8,dimension(:),allocatable :: Rbem,Cbem,RWORK
      real*8 :: RCOND,FERR,BERR
      character*1 :: EQUED
      complex*16,dimension(400,5) :: OUTVAR
      complex*16, dimension(1:5) :: MecElem
      type(Punto), pointer :: p_X
!       integer, dimension(275) :: pasaNopasa 
!       pasaNopasa(1:275) = 1 !1 pasa, 0 no pasa
      
      outvar = 0
      !#< blue
      call system('clear')
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .eq. 0) then 
        print*,"usage DWNIBEM.out [option]"
        print*,""
        print*,"     -r        : Correr programa normal."
        print*,"     -r FF       : Correr programa normal sólo en campo libre"
        print*,"     -g        : Solo imprimir geometría."
        print*,"     -a        : Solo Imprimir funcion de amplitud."
        print*,"     -f [frec] : Solo Ejecutar la frecuencia [frec]."
        print*,"     -f [frec] FF : Solo Ejecutar la frecuencia [frec] en campo libre."
        print*,"     -L        : Cargar fotogramas de video/G.bin"
        print*,""
        stop 
      else
        if (trim(arg) .eq. '-r' .or. &
            trim(arg) .eq. '-g' .or. &
            trim(arg) .eq. '-a' .or. &
            trim(arg) .eq. '-f' .or. &
            trim(arg) .eq. '-L') then 
           print*,"Correr programa"
        else
           print*,"unknown input argument"
           stop
        end if
      end if !#>
      call ETIME(tarray, result)
      call setupdirectories ! llama getmaininput
      call getSoilProps 
      frecIni = NFREC+1
      frecEnd = 1
      onlythisJ = .false. !#< blue
      IF (LEN_TRIM(arg) .ne. 0) then
        if (trim(arg) .eq. '-f') then
          CALL get_command_argument(2, arg)
          IF (LEN_TRIM(arg) .ne. 0) then
            read(arg,*) frecIni
            if (frecIni .lt. 1) stop "-f lower than 1"
            if (frecIni .gt. NFREC+1) stop "-f higher than NFREC+1 "
            frecEnd = frecIni
            onlythisJ = .true.
            CALL get_command_argument(3, arg)
            IF (LEN_TRIM(arg) .ne. 0) then
              if(trim(arg) .eq. 'FF') then
                print*,"Override: Only Free Field" 
                write(PrintNum,'(a)') "Override: Only Free Field"
                workBoundary = .false.
              end if
            end if
          else
            print*," "; print*," "
            print*,"no single frequency specified on argument -f"
            stop
          end if
        end if!
        CALL get_command_argument(2, arg)
            IF (LEN_TRIM(arg) .ne. 0) then
              if(trim(arg) .eq. 'FF') then
                print*,"Override: Only Free Field" 
                write(PrintNum,'(a)') "Override: Only Free Field"
                workBoundary = .false.
              end if
            end if
      end if !#>
      call checarWisdom(2*nfrec,2*nmax,NPTSTIME) ! FFTw
        nIpts=0; nMpts=0; nBpts = 0; iPtfin = 0; mPtfin = 0
        
      call getsource
      call getPolaridad(skipdir,PSV,SH)
      if(.not. onlythisJ) allocate(Uo(NPTSTIME,nFuentes)) 
      
      if (PSV) write(PrintNum,*) "  this a P-SV case"
      if (SH)  write(PrintNum,*) "  this a SH case"
      if (workBoundary .eqv. .true.) call getTopography
      call getInquirePoints
      allocate(XF(nSabanapts,Nfrec+1,2))
      
      call getVideoPoints
      
        Npts = nIpts + nMpts
        write(PrintNum,'(a,I0)') '   Number of fixed receptors: ',Npts
        allocate (allpoints(Npts))
        allpoints(iPtini:iPtfin)= inqPoints(iPtini:iPtfin); deallocate(inqPoints)
      if (makeVideo) then
        allpoints(mPtini:mPtfin) = moviePoints; deallocate(moviePoints); end if
        
      call setInqPointsRegions
      do iP_x = 1,nIpts
       if (allpoints(iP_x)%isOD) then 
       if (allpoints(iP_x)%tipoFrontera .eq. 0) allpoints(iP_x)%region = 1
       if (allpoints(iP_x)%tipoFrontera .eq. 1) allpoints(iP_x)%region = 1
       if (allpoints(iP_x)%tipoFrontera .eq. 2) allpoints(iP_x)%region = 2
       end if
      end do
      
      if (makeVideo) then
      call setVideoPointsRegions
      call chdir(trim(adjustl(rutaOut)))
      call Churubusco(.true.)
      call chdir("..")
      end if !#< blue
       call chdir(trim(adjustl(rutaOut)))
       call system("mkdir matrices")
       do currentiFte = 1,nFuentes 
       write(arg,'(a,I0)') 'mkdir phi',currentiFte
       call system(trim(adjustl(arg)))
       
       write(arg,'(a,I0)') 'mkdir traces',currentiFte
       call system(trim(adjustl(arg)))
       
       write(arg,'(a,I0)') 'mkdir video',currentiFte
       call system(trim(adjustl(arg))) 
       end do !
       do currentiFte = 1,nfuentes
         if(.not. onlythisJ) call sourceAmplitudeFunction
         write(titleN,'(a,I0,a)') '0___Geometry',currentiFte,'.pdf'
         write(extension,'(a)') 'PDF'
         BP => BouPoints
         call drawBoundary(BP,nbpts,titleN,extension,.false.,.false.,.false.,.true.)
         write(titleN,'(a,I0,a)') '0___ModeloCompleto',currentiFte,'.pdf'
         call drawBoundary(BP,nbpts,titleN,extension,.false.,.true.,.true.,.true.)
       end do !
       currentiFte = 0
       
       write(extension,'(a)') 'PDF'
       BP => BouPoints
       write(titleN,'(a)') '0___SensorsA.pdf'
       call drawBoundary(BP,nbpts,titleN,extension,.false.,.true.,.false.,.false.)
       if (nsecciones .gt. 0) then
       write(titleN,'(a)') '0___SensorsB.pdf'
       call drawBoundary(BP,nbpts,titleN,extension,.false.,.false.,.true.,.false.)
       end if!
       
       if (workBoundary .eqv. .true.) then 
         write(titleN,'(a)') '0___Inclusion.pdf'
         write(extension,'(a)') 'PDF'
         BP => BouPoints
         call drawBoundary(BP,nbpts,titleN,extension,.true.,.false.,.false.,.false.)
       end if
       call chdir("..")
              
       call get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then 
        if (trim(arg) .eq. '-g') stop "argumento -g : Geometría"
      end if !#>
! predimensionamientos
      write(PrintNum,'(a)') &
       "---------------------------------------------------------------------------------"
      if (SH) then 
        i = 2*N+1
        if (Z(0) .lt. 0.0) i = i + 1 !HS arriba
      end if !
      if (PSV) then 
        i = 4*N+2
        if (Z(0) .lt. 0.0) i = i + 2 !HS arriba
      end if
      ik = 0 !indice para onda plana
        allocate (Ak(i,i,0:2*nmax)); allocate (B(i,0:2*nmax)) 
        allocate(ipivA(i)); allocate(workA((i)*(i)))
   !   to store the strata diffracted displacement: W,U,V,...
      do iP=iPtini,iPtfin
        allocate(allpoints(ip)%resp(NFREC+1,nFuentes)) ! para prescindir de W
        if (comoFacDeAmpliDinamica ) & 
        allocate(allpoints(ip)%facAmpli(NFREC+1,nFuentes)) ! para prescindir de W
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%U = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%V = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%W = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Tx = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Ty = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%Tz = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%szz = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%szx = z0
        allpoints(ip)%resp(1:NFREC+1,1:nFuentes)%sxx = z0
         ! Plot integrand FK of Green(l) function
        if(plotFKS) then; if(allpoints(iP)%guardarFK)then
         allocate(allpoints(iP)%FK(NFREC+1,NMAX,3)) !W,U,V
         allpoints(iP)%FK = Z0
        end if;end if
        ! para los puntos en la frontera en los que se quiere graficar
        if (allpoints(iP)%atBou) allocate(allpoints(iP)%S(NPTSTIME,5)) !(traza,componente)
      end do
      
       if (makeVideo) then
        do iP=mPtini,mPtfin
         allocate(allpoints(iP)%Wmov(NFREC+1,3,npixX,nFuentes)) !W,U,V
         allpoints(iP)%Wmov = z0
        end do
       end if
      ! un indice de punto para todos
      do i=1,Npts
        allpoints(i)%pointIndex = i
      end do ! i
      allocate( gamma_E(0:2*nmax,N+1))
      allocate( nu_E(0:2*nmax,N+1))
      allocate( eta_E(0:2*nmax,N+1))
      allocate( subMatD0(2,4,N+1,0:2*nmax))
      allocate( subMatS0(3,4,N+1,0:2*nmax))      
      allocate( k_vec(2*nmax))
      k_vec(1) = real(dk * 0.01,8)
      do ik = 2,NMAX+1
        k_vec(ik) = real(ik-1,8) * dk
      end do!
      do ik = nmax+2,2*NMAX
        k_vec(ik) = (ik - 2*NMAX - 1) * dk
      end do
      allocate(t_vec(NPTSTIME))
      t_vec(1:NPTSTIME) = z0
      t_vec(1) = exp(cmplx(0.0,-0.01* dfrec * t0*(2*pi),8))
      do i = 2,nfrec+1
        t_vec(i) = exp(cmplx(0.0,-(i-1)* dfrec * t0*(2*pi),8))
      end do!
      do i = NPTSTIME-NFREC+2,NPTSTIME
        t_vec(i) = exp(cmplx(0.0,-(i-NPTSTIME-1)* dfrec * t0*(2*pi),8))
      end do
      
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then 
        if (trim(arg) .eq. '-L') then
        ! cargar G para misma fuente,geometria,receptores
        ! Cuando cambiamos la función de amplitud
          print*,""
          print*,"********************"
          print*,"*** loading data ***"
          print*,"********************"
          print*,""
          write(PrintNum,'(a)') "*** loading data ***"
         !#< r  TODO: recalcular campo incidente !#>
         call loadW(PSV,SH)
!        nPtsolos = 3
!        ! hacer ceros
!        do i=1,275
!          if (pasaNopasa(i) .eq. 0) then
!            print*," zeros to the ",i
!            allpoints(i)%resp(:,:)%U = 0
!            allpoints(i)%resp(:,:)%V = 0
!            allpoints(i)%resp(:,:)%W = 0
!            allpoints(i)%resp(:,:)%Tx = 0
!            allpoints(i)%resp(:,:)%Ty = 0
!            allpoints(i)%resp(:,:)%Tz = 0
!            allpoints(i)%resp(:,:)%szz = 0
!            allpoints(i)%resp(:,:)%szx = 0
!            allpoints(i)%resp(:,:)%sxx = 0
!          end if
!        end do
         ! borr´e facamp
         
         call chdir(trim(adjustl(rutaOut)))
         if (makeVideo) call loadG_fotogramas
      do currentiFte = 1,nfuentes
           print*,"fuente=,",currentiFte
           write(arg,'(a,I0)') 'traces',currentiFte
           CALL chdir(trim(arg))
           

           call plotSisGram(PSV,SH,.false.)
           if (plotFKCK) call F_K_exp(XF)
           CALL chdir("..")
           
           if (makeVideo) then 
             !call loadG_fotogramas
             write(arg,'(a,I0)') 'video',currentiFte
             CALL chdir(trim(arg))
             if (PSV .and. vivaChurubusco) call Churubusco(.false.)
             if (SH) call Hollywood(3)
             CALL chdir("..")
           end if
           !
           if (workboundary .and. punEnlaFront) then
             write(arg,'(a,I0)') 'video',currentiFte
             CALL chdir(trim(arg))
             if (PSV) call CINETECA
             CALL chdir("..")
           end if!
      end do
         stop "done"
        end if ! load 
      end if! argument
      
      !#< g
!     if (developerfeature .ne. 0) then
!     allocate(developerAUXvec(NFREC,2))
!     end if !#>
      call preparePointerTable(pota,.true.,0.1_8) !#< b
      call ETIME(tarray, result)
      call system('clear')
      write(PrintNum,*)"Pre-prosses took ", result,"seconds"
      write(PrintNum,'(/,/)')
      lastresult = result
      if (Verbose .le. 1) then
      write(6,'(A)', ADVANCE = "NO") " "
      write(6,'(A)', ADVANCE = "NO") repeat("_",58)
      write(6,'(A)', ADVANCE = "NO") " "
      print*,""
      write(6,'(A)', ADVANCE = "NO") repeat(" ",60)
      end if !#>
      DO J = frecIni,frecEnd,-1 !NFREC+1,1,-1
        FREC=DFREC*real(J-1); if (J .eq. 1)  FREC = 0.5_8 * DFREC  ! Hz
        OME=2.0*PI*FREC !rad/s
        COME = CMPLX(OME, OMEI,8)!periodic sources damping
        COME = COME * cmplx(1.0, -1.0/2.0/Qq,8) !histeretic damping
        pt_come_i => COME
        ik = N+1; i = 1
        if (workBoundary) ik = N+2 !peros solo si hay inclusion (con topografia no)
        if (Z(0) .lt. 0.0) i = 0
        if (useAzimi) then
        ! Azimi attenuation: see Kennett, The seismic wavefield, vol1 pag 148-150
          do l=i,ik !N+1 estratos y una inclusión
          ! "Normally we expect that loss in dilatation is very small compared 
          ! with that in shear so that Qp^-1 << Qs^-1 , and then "
!         Qs(l) = Qq !* 4./3. * (beta0(l)/alfa0(l))**2. !; print*,"Qs ",Qs(l)
!         Qp(l) = Qq !; print*,"Qp ",Qp(l)
          ! alfa0 y beta0 son velocidades de referencia a T=1sec
           beta(l) = cmplx((beta0(l))*(1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
           alfa(l) = cmplx((alfa0(l))*(1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
           aMU(l) = RHO(l) * beta(l)**2.
           Lambda(l) = RHO(l)*alfa(l)**2. - real(2.)*aMU(l)
          end do 
        else
          do l=i,ik !N+1 estratos y una inclusión
           beta(l) = beta0(l); alfa(l) = alfa0(l)
           aMU(l) = RHO(l) * beta(l)**2.
           Lambda(l) = RHO(l)*alfa(l)**2. - real(2.)*aMU(l)
          end do 
        end if
      !#< b
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),60)
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),17) !eta
        if (workBoundary) write(6,'(A)', ADVANCE = "NO") repeat(char(8),10) !nbp
        write(6,'(A)', ADVANCE = "NO") "["
        write(6,'(A,A)', ADVANCE = "NO") & 
        repeat("X",int((58.0/NFREC)*(NFREC+1-J))),&
        repeat("_",58-int((58.0/NFREC)*(NFREC+1-J)))
        write(6,'(A)', ADVANCE = "NO") "]" 
      if(PrintNum .ne. 6)then
        call ETIME(tarray, result)
        write(PrintNum,'(A,I0,A,EN18.2,A,EN18.2,A)', ADVANCE = "YES") &
        'w(',J,') | ',FREC," | ",result,"sec"
      end if
      if (onlythisJ) write(PrintNum,*) "Frec = ",COME, " eta = ",&
      (ome*longitudcaracteristica_a)/(pi*beta0(N+1))!#>
      
      if (workBoundary) then ! Subsegment the topography if neccesssary
         call subdivideTopo(J,FREC, onlythisJ,minBeta,beta0(i:N+2),nbpts,BouPoints)
          write(6,'(A,I5)', ADVANCE = "NO") "nbp= ", nbpts
         call preparePointerTable(pota,.false.,smallestWL) !(solo DWNs)
    
       l = n_top_sub + 2* n_con_sub + n_val_sub
       if (PSV) ik = 2 ;  if (SH) ik = 1
       i = 0
       if (.not. allocated(ibemMat))then
         i = 1
       else
         if (overDeterminedSystem) then
           m = ik * (l+n_OD)
         else
           m = ik * l
         end if!
         if (size(ibemMat,1) .ne. m) then 
         i = 1
         deallocate(ibemMat);deallocate(trac0vec)
         deallocate(IPIVbem);deallocate(auxGvector)
         deallocate(ibemMatS)
         end if
       end if!
       if (i .eq. 1) then
        if ((overDeterminedSystem) .and. & 
            (OD_Jini .le. J) .and. &
            (J .le. OD_Jend)) then
         allocate(ibemMat(ik * (l+n_OD),ik * l));allocate(trac0vec  (ik * (l+n_OD)))
         !allocate(IPIVbem(1));                   allocate(auxGvector(ik * (l+n_OD)))
         allocate(IPIVbem(ik * l)); allocate(auxGvector(ik * l))
         allocate(ibemMatS(ik * (l+n_OD),ik * l))
        else
         allocate(ibemMat(ik * l,ik * l));allocate(trac0vec  (ik * l))
         allocate(IPIVbem(ik * l));       allocate(auxGvector(ik * l))
         allocate(ibemMatS(ik * l,ik * l))
        end if
       end if
       ibemMat = z0 ; trac0vec = z0 ; auxGvector = z0
      end if!workbou
      
         pt_ipivA => ipivA
         pt_workA => workA
         
         call makeGANU (J) !los numeros de onda horizontales para P y S 
      if (PSV) then
         tam = size(Ak,1)
!      print*,"vecNK(J)=",vecNK(J)
!      print*,"   pos",min(int(vecNK(J)*SpliK),nmax)+1," neg", 2*nmax-((min(int(vecNK(J)*SpliK),nmax))-2) 
!      k0 = min(int(vecNK(J)*SpliK),nmax)+1 !vecNK(J) 
!      k0 = vecNK(J); 
!      print*,"ik= 1:",k0
       do ik = 1,vecNK(J) ! k positivo (Aorig -> nmax+1)
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)            
!        print*,ik, 2*nmax - (ik-2)   
         call gloMat_PSV(pointAp,pt_k,ik)
!        print*,ik,"A   =",sum(pointAp)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
!        print*,ik,"A-1 =",sum(pointAp)
       end do ! ik
       k0 = vecNK(J); 
       call intrplr_gloMat(k0,15,pt_cOME_i,pt_ipivA,pt_workA)         
!      do ik = 1,min(int(vecNK(J)*SpliK),nmax)+1
!        print*,ik,sum(Ak(1:tam,1:tam,ik))
!      end do
!      print*,"sum=",sum(Ak(1:tam,1:tam,1:min(int(vecNK(J)*SpliK),nmax)+1))
!      stop
       call parImpar_gloMat ! k negativo
       
       
       
       ! Funciones de Green para propagar a la frontera de cada estrato
       ! (sin fase vertical, la fase vertical se agrega para cada fuente)
       if(.not. allocated(BparaGa)) allocate(BparaGa(tam,N+1,2*nmax,2));BparaGa=0
       if(.not. allocated(BparaNu)) allocate(BparaNu(tam,N+1,2*nmax,2));BparaNu=0
       call PSVpaGaNU(J)
       ! Multiplicar partes de la m
       if(.not. allocated(CoefparGa)) allocate(CoefparGa(tam,N+1,2*nmax,2,2))
       if(.not. allocated(CoefparNu)) allocate(CoefparNu(tam,N+1,2*nmax,2,2))
       call PSVMatAporGaNU(J)
       
      end if!psv ............................................
      if (SH) then
         !Ak = Z0
         tam = size(Ak,1)
         !print*,"N=",N," tam=",tam
!        k0 = min(int(vecNK(J)*SpliK),nmax)+1
      Do ik = 1,vecNK(J)!k0!nmax+1
      ! k positivo
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)  
!        print*,ik, 2*nmax - (ik-2) 
         call globalmatrix_SH(pointAp,pt_k,ik)
!        print*,ik,"A   =",sum(pointAp)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
!        print*,ik,"A-1 =",sum(pointAp)
      end do ! ik
       k0 = vecNK(J); 
      call intrplr_gloMat(k0,15,pt_cOME_i,pt_ipivA,pt_workA) 
      ! k negativo (es simétrico)
         Ak(1:tam,1:tam,Nmax+2:2*nmax) = &
         Ak(1:tam,1:tam,Nmax:2:-1)
!     k0 = min(int(vecNK(J)*SpliK),nmax); k0 = 2*nmax-(k0-2)
!        Ak(1:tam,1:tam,k0:2*nmax) = &
!        Ak(1:tam,1:tam,k0-2:2:-1)
      end if!sh
     
      ! matriz IBEM 
      if (workboundary) then
      currentiFte = 0
      do dir= 1,3 !x,y,z direction of force application
        if(dir .eq. 2) then
         if(skipdir(dir)) cycle
        else ! 1 o 3
         if(skipdir(1) .and. skipdir(3)) cycle
        end if! dir
      if (verbose .ge. 2) print*,""
      if (verbose .ge. 2) print*,'********(campo difractda en medio estratificado)*********************'
         do iz = 1,nZs !por cada Z de fuentes en tabla de segmentos tipo 0 y 1
           if (thereisavirtualsourceat(iz)) then ! (ipxi = 1,n_top_sub)
             call diffField_at_iz(iz,dir,J,cOME)
           end if 
         end do !iz
      if (verbose .ge. 2) print*,'********(campo refractado en inclusión y FF de inclusion (columnas xi=d2))***' 
         if (n_con_sub .gt. 0) then 
           do iPxi = n_top_sub +1, n_top_sub + n_con_sub 
             call reffField_by_(iPxi,dir,cOME)
           end do ! iPxi
         end if !n_cont
      if (verbose .ge. 2) print*,'********(frontera libre en inclusión (columnas de xi=d1))*****************'
         if (n_val_sub .gt. 0) then
           do iPxi = n_top_sub + n_con_sub +1, n_top_sub + n_con_sub + n_val_sub
             call reffField_by_(iPxi,dir,cOME)
           end do ! iPxi
         end if !n_vall
      if (verbose .ge. 2) print*,'*****(funcions de Green de desplazamiento en la región R)*****************'
         call GreenReg_R(J,dir,cOME)
      end do !dir
      
      ! guardar la matriz ibem para que no se estropee al invertir distintas fuentes
      ibemMatS = ibemMat
      end if !workboundary
!     campo incidente 
       do currentiFte = 1,nFuentes !para cada una de las incidencias
       if (verbose .ge. 2) print*,'*****( Incidencia )*************************************'
       call makeGANU0 !gamma y nu con el número de onda de ésta incidencia plana
       if (workboundary) then
         ibemMat = ibemMatS !guardada por si hay varias incidencias
         trac0vec = 0
       end if
       
       do dir= 1,3 !x,y,z direction of force application
!       print*,"dir=",dir
        if(dir .eq. 2) then
           if(skipdir(dir)) cycle
        else ! 1 o 3
           if(skipdir(1) .and. skipdir(3)) cycle
        end if! dir
        if(.not. skipdir(dir)) then 
           if (Po(currentiFte)%region .eq. 1) then
              call diffField_at_iz(0,dir,J,cOME)
           else ! Po en región 2
!          if (Po(currentiFte)%region .eq. 2) then
              if (workboundary) then
                 call termIndepR(dir,cOME)
              else
                  stop "indicó región 2 en un problema de 1 región"
              end if
           end if
        end if
       end do !dir 
       if (workboundary) then
      if (PSV) then 
      ik = 2
      l = n_top_sub + 2* n_con_sub + n_val_sub  
      if ((overDeterminedSystem) .and. & 
          (OD_Jini .le. J) .and. &
          (J .le. OD_Jend)) then !#< r -.-.-.-.-.-.-.-.-.-.-.-.-.- !#>
      Mi = ik*(l + n_OD)
      Ni = ik*l
      call overDetermineSystem(J,ik*l+1,PSV)
      call solveOverDetermineSystem(Mi,Ni)
      Mi = Ni
      else !#< r -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- !#>
      iPIVbem = ik*l
      Mi = ik*l
      Ni = ik*l
       
!#< b
!      call chdir(trim(adjustl(rutaOut)))
!      call chdir('matrices')
!     print*,n_top_sub,"n_top_sub"
!     print*,n_con_sub,"n_con_sub"
!     print*,n_val_sub,"n_val_sub"
!      write(arg,'(a,I0,a)') "outA",J,".m"
!      open(421,FILE= trim(arg),action="write",status="replace")
!      write(arg,'(a)') "BiSIN"
!      call scripToMatlabMNmatrixZ(size(trac0vec,1),1,trac0vec,arg,421)
!      write(arg,'(a,I0)') "Mi",J
!      call scripToMatlabMNmatrixZ(size(ibemMat,1),size(ibemMat,2),ibemMat,arg,421)
!      close(421)
!      CALL chdir(".."); CALL chdir("..")
!      stop 456
!      
!     if (verbose .ge. 1) call showMNmatrixZ(size(ibemMat,1),size(ibemMat,2), ibemMat ," mat ",6);stop
!     if (verbose .ge. 1) call showMNmatrixZ(size(trac0vec,1),1 , trac0vec,"  b  ",6) ;stop
!      call chdir(trim(adjustl(rutaOut))) 
!      open(421,FILE= "outA.m",action="write",status="replace")
!      write(arg,'(a)') "Bf"
!      call scripToMatlabMNmatrixZ(Mi,1,trac0vec(1:Mi),arg,421)
!      write(arg,'(a)') "Mf"
!      call scripToMatlabMNmatrixZ(Mi,Ni,ibemMat(1:Mi,1:Ni),arg,421)
!      close(421)
!      CALL chdir("..")
!      stop 456 
!#>
            
      ! driver simple:
!     call zgesv(Mi,1,ibemMat,Mi,IPIVbem,trac0vec,Mi,info)
      
      ! driver experto:
      allocate(AF(Mi,Mi))
      allocate(Xbem(Mi,1))
      allocate(Rbem(Mi))
      allocate(Cbem(Mi))
      allocate(WORKbem(2*Mi))
      allocate(RWORK(2*Mi))
      call ZGESVX('E','N',Mi,1,ibemMat,Mi,AF,Mi,IPIVbem,EQUED,Rbem, &
      Cbem,trac0vec,Mi,Xbem,Mi,RCOND,FERR,BERR,WORKbem,RWORK,INFO)
      
      if(info .ne. 0) stop "problem with ibem system"
      if (any(isnan(real(trac0vec)))) stop "751 valió madres el ibem"
      trac0vec = Xbem(:,1)
      deallocate(Af);deallocate(Xbem);deallocate(Rbem);deallocate(Cbem)
      deallocate(WORKbem);deallocate(RWORK)
      end if !#< r -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- !#>
      
      
      !#< b
!     call showMNmatrixZ(Mi,1, trac0vec,"phi  ",6)
      if (verbose .ge. 2) then
         call chdir(trim(adjustl(rutaOut)))
         write(arg,'(a,I0)') 'phi',currentiFte
         CALL chdir(trim(arg))
         write(titleN,'(a,I0,a)')'h_phi_',J,'_E.pdf'
         call drawPHI(titleN,1)
         write(titleN,'(a,I0,a)')'v_phi_',J,'_E.pdf'
         call drawPHI(titleN,2)
         if (n_con_sub .gt. 0) then
            write(titleN,'(a,I0,a)')'h_phi_',J,'_R.pdf'
            call drawPHI(titleN,-1)
            write(titleN,'(a,I0,a)')'v_phi_',J,'_R.pdf'
            call drawPHI(titleN,-2)
         end if
         if(verbose .ge. 3)  call showMNmatrixZabs(Mi,1, trac0vec," phi ",6)
         CALL chdir(".."); CALL chdir("..")
      end if!
      if(verbose .ge. 4) write(PrintNum,'(a)') "add diffracted field by topography"
      !#>
            
      do iP_x = 1,nIpts !cada receptor X 
!         print*,""
!         print*,allpoints(iP_x)%pointIndex," iP_x=(",allpoints(iP_x)%center,")"
!         print*,allpoints(iP_x)%resp(J,currentiFte)%W," W"
!         print*,allpoints(iP_x)%resp(J,currentiFte)%U," U"
!         print*,allpoints(iP_x)%resp(J,currentiFte)%Tx," Tx"
!         print*,allpoints(iP_x)%resp(J,currentiFte)%Tz," Tz"
          if (allpoints(iP_x)%region .eq. 1) then !'estr'
            iPhi_I = 1
            iPhi_F = (n_top_sub + n_con_sub) * 2
            ipxi_I = 1
            ipxi_F =  n_top_sub + n_con_sub
          else if (allpoints(iP_x)%region .eq. 2) then !'incl'
            iPhi_I = (n_top_sub + n_con_sub) * 2 + 1
            iPhi_F = (n_top_sub + 2* n_con_sub + n_val_sub) * 2
            ipxi_I =  n_top_sub + 1
            ipxi_F =  n_top_sub + n_con_sub + n_val_sub
          else ! 'void'
            cycle
          end if
        do i=1,2 !dirección de desplazamiento {W,U} en ip_X !PSV
          auxGvector = z0
          iPhi = iPhi_I 
          do iPxi = ipxi_I,ipxi_F ! recopilamos  G_ij
            auxGvector(iPhi)   = boupoints(iPxi)%G(iP_X,i,1) ! por fzas horizontales:
            auxGvector(iPhi+1) = boupoints(iPxi)%G(iP_X,i,3) ! por fzas verticales:
            iPhi = iPhi + 2
          end do         
          if(i .eq. 1) allpoints(iP_x)%resp(J,currentiFte)%W = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &!sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%W
          if(i .eq. 2) allpoints(iP_x)%resp(J,currentiFte)%U = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &!sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%U
        end do !i
        
        do i=4,5 !tracciones
          auxGvector = z0
          iPhi = iPhi_I 
          do iPxi = ipxi_I,ipxi_F ! recopilamos  G_ij
            auxGvector(iPhi)   = boupoints(iPxi)%G(iP_X,i,1) ! por fzas horizontales:
            auxGvector(iPhi+1) = boupoints(iPxi)%G(iP_X,i,3) ! por fzas verticales:
            iPhi = iPhi + 2
          end do 
          if(i .eq. 4) allpoints(iP_x)%resp(J,currentiFte)%Tz = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%Tz
          if(i .eq. 5) allpoints(iP_x)%resp(J,currentiFte)%Tx = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%Tx
        end do !i 
        
        do i=7,9 !esfuerzos
          auxGvector = z0
          iPhi = iPhi_I 
          do iPxi = ipxi_I,ipxi_F ! recopilamos  G_ij
            auxGvector(iPhi)   = boupoints(iPxi)%G(iP_X,i,1) ! por fzas horizontales:
            auxGvector(iPhi+1) = boupoints(iPxi)%G(iP_X,i,3) ! por fzas verticales:
            iPhi = iPhi + 2
          end do 
          if(i .eq. 7) allpoints(iP_x)%resp(J,currentiFte)%sxx = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%sxx
          if(i .eq. 8) allpoints(iP_x)%resp(J,currentiFte)%szx = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%szx
          if(i .eq. 9) allpoints(iP_x)%resp(J,currentiFte)%szz = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J,currentiFte)%szz
        end do !i
      
      end do !iP_x
      !
      if (makeVideo) then 
      do m = 1,npixX
      do iP_x = mPtini,mPtfin !cada receptor X 
      do i=1,2 !dirección de desplazamient {W,U} en ip_X !PSV
          if (fotogramas_Region(iP_x-nIpts,m) .eq. 1) then !'estr'
!         print*,m,iP_x-nIpts,"estr"
            iPhi_I = 1
!           iPhi_F = (n_top_sub + n_con_sub) * 2
            ipxi_I = 1
            ipxi_F = n_top_sub + n_con_sub
          else if (fotogramas_Region(iP_x-nIpts,m) .eq. 2) then !'incl'
!         print*,m,iP_x-nIpts,"incl"
            iPhi_I = (n_top_sub + n_con_sub) * 2 + 1
!           iPhi_F = (n_top_sub + 2* n_con_sub + n_val_sub) * 2
            ipxi_I = n_top_sub + 1
            ipxi_F = n_top_sub + n_con_sub + n_val_sub
          else ! 'void'
!         print*,m,iP_x-nIpts,"void" 
            cycle
          end if
          auxGvector = z0
          iPhi = iPhi_I 
          do iPxi = ipxi_I,ipxi_F ! recopilamos  G_ij
            auxGvector(iPhi)   = boupoints(iPxi)%Gmov(iP_X-nIpts,i,1,m) ! por fzas horizontales:
            auxGvector(iPhi+1) = boupoints(iPxi)%Gmov(iP_X-nIpts,i,3,m) ! por fzas verticales:
            iPhi = iPhi + 2
          end do 
          fotogramas(iP_x-nIpts,m,J,i,currentiFte) = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                                         allpoints(iP_x)%Wmov(J,i,m,currentiFte)
          
      end do !i
      end do !iP_x
      end do !m
      end if !makevideo
      end if !PSV
      
      
!     stop 2125
      if (SH) then
      ik = 1
      l = n_top_sub + 2* n_con_sub + n_val_sub !cantidad de segmentos
      iPIVbem = ik * l
      Mi = ik * l
      Ni = ik * l
!     !#< b
!     if (verbose .ge. 3) call showMNmatrixZ(Mi,Ni, ibemMat," mat ",6)
!     if (verbose .ge. 3) call showMNmatrixZ(Mi,1, trac0vec,"  b  ",6)
!      call chdir(trim(adjustl(rutaOut)))
!      call chdir('matrices')
!     print*,n_top_sub,"n_top_sub"
!     print*,n_con_sub,"n_con_sub"
!     print*,n_val_sub,"n_val_sub"
!      write(arg,'(a,I0,a)') "outA",J,".m"
!      open(421,FILE= trim(arg),action="write",status="replace")
!      write(arg,'(a)') "Bincidente"
!      call scripToMatlabMNmatrixZ(size(trac0vec,1),1,trac0vec,arg,421)
!      write(arg,'(a,I0)') "Mibem",J
!      call scripToMatlabMNmatrixZ(size(ibemMat,1),size(ibemMat,2),ibemMat,arg,421)
!      close(421)
!      CALL chdir(".."); CALL chdir("..")
!     !#>
      call zgesv(Mi,1,ibemMat,Mi,IPIVbem,trac0vec,Mi,info)
      if(info .ne. 0) stop "problem with ibem system"
      
      if (any(isnan(real(trac0vec)))) then ! NAN is not equal even to itself
      stop "891 valio madres el ibem"; end if!
      !#< b
      !call showMNmatrixZ(Mi,1, trac0vec," phi ",6)
      if (verbose .ge. 2) then  
         call chdir(trim(adjustl(rutaOut))) 
         write(arg,'(a,I0)') 'phi',currentiFte
         CALL chdir(trim(arg))
         write(titleN,'(a,I0,a)')'n_phi_',J,'_E.pdf'
         call drawPHI(titleN,3)
         if (n_con_sub .gt. 0) then
            write(titleN,'(a,I0,a)')'n_phi_',J,'_R.pdf'
            call drawPHI(titleN,-3)
         end if!
         if (verbose .ge. 3) then 
           call showMNmatrixZ(Mi,1, trac0vec," phi ",6)
           write(titleN,'(a,I0,a)')'n_phiVal_',J,'.pdf'
           write(CTIT,'(I0)') Mi
           write(yax,'(a)') 'abs(phi) '
           call plotXYcomp(trac0vec(1:Mi),1.0,Mi,titleN, &
           '  phi n  ',yax, CTIT,1200,800,0.0)
         end if
        CALL chdir(".."); CALL chdir("..")
      end if! 
      if(verbose .ge. 4) write(PrintNum,'(a)') "add diffracted field by topography"
      !#>
      do iP_x = 1,nIpts !cada receptor X 
          if (allpoints(iP_x)%region .eq. 1) then !'estr'
            iPhi_I = 1
            iPhi_F = n_top_sub + n_con_sub
            ipxi_I = 1
            ipxi_F = n_top_sub + n_con_sub
          else if (allpoints(iP_x)%region .eq. 2) then !'incl'
            iPhi_I = n_top_sub + n_con_sub + 1
            iPhi_F = n_top_sub + 2* n_con_sub + n_val_sub
            ipxi_I = n_top_sub + 1
            ipxi_F = n_top_sub + n_con_sub + n_val_sub
          else ! 'void'
            cycle
          end if
      do i= 3,6,3!dirección de desplazamientraccion V,Ty !SH
         auxGvector = z0 !(1:ik)
!        print*, iPhi_I,iPhi_F, " -- ", ipxi_I,ipxi_F
         iPxi = ipxi_I
         do iPhi = iPhi_I,iPhi_F
           auxGvector(iPhi) = boupoints(iPxi)%G(iP_X,i,2)
!          if(i.eq.3) write(6,'(a,i0,a,i0,a,F10.6,F10.6)') & 
!          "G(",iP_X,",",iPxi,")22=",& 
!          real(boupoints(iPxi)%G(iP_X,i,2)),aimag(boupoints(iPxi)%G(iP_X,i,2))
           iPxi = iPxi + 1
         end do !iPhi         
          if (verbose .ge. 4) call showMNmatrixZ(Mi,1,auxGvector," auxG",6)
          
          
          if (i .eq. 3) allpoints(iP_x)%resp(J,currentiFte)%V = & 
                        allpoints(iP_x)%resp(J,currentiFte)%V + &
                        sum(trac0vec * auxGvector)
          if (i .eq. 6) allpoints(iP_x)%resp(J,currentiFte)%Ty = & 
                        allpoints(iP_x)%resp(J,currentiFte)%Ty + &
                        sum(trac0vec * auxGvector)
      end do !i
      end do !iP_x
      
      if (makeVideo) then
      do m = 1,npixX
      do iP_x = mPtini,mPtfin !cada receptor X 
!     do i=3,3 !dirección de desplazamient V !SH
          i = 3
          if (fotogramas_Region(iP_x-nIpts,m) .eq. 1) then !'estr'
!         print*,m,iP_x-nIpts,"estr"
            iPhi_I = 1
            iPhi_F = n_top_sub + n_con_sub
            ipxi_I = 1
            ipxi_F = n_top_sub + n_con_sub
          else if (fotogramas_Region(iP_x-nIpts,m) .eq. 2) then !'incl'
!         print*,m,iP_x-nIpts,"incl"
            iPhi_I = n_top_sub + n_con_sub + 1
            iPhi_F = n_top_sub + 2* n_con_sub + n_val_sub
            ipxi_I = n_top_sub + 1
            ipxi_F = n_top_sub + n_con_sub + n_val_sub
          else ! 'void'
!         print*,m,iP_x-nIpts,"void"
            cycle
          end if
          auxGvector = z0 !(1:ik)
!        print*, iPhi_I,iPhi_F, " -- ", ipxi_I,ipxi_F
         iPxi = ipxi_I
         do iPhi = iPhi_I,iPhi_F
!          print*, iPxi ,iP_X-nIpts, m, boupoints(iPxi)%Gmov(iP_X-nIpts,i,2,m) 
           auxGvector(iPhi) = boupoints(iPxi)%Gmov(iP_X-nIpts,i,2,m)
           iPxi = iPxi + 1
         end do !iPhi 
         if (verbose .ge. 4) call showMNmatrixZ(Mi,1,auxGvector," auxG",6)
          fotogramas(iP_x-nIpts,m,J,i,currentiFte) = sum(trac0vec * auxGvector) + &
          allpoints(iP_x)%Wmov(J,i,m,currentiFte) 
!     end do !i
      end do !iP_x
      end do !m
      end if !makevideo
      end if !SH
!     stop 2214
      else !not workboundary
      if (makeVideo) then
      if (PSV) then
      do m = 1,npixX
      do iP_x = mPtini,mPtfin  
      do i=1,2 !dirección de desplazamient {W,U} en ip_X !PSV
          fotogramas(iP_x-nIpts,m,J,i,currentiFte) = allpoints(iP_x)%Wmov(J,i,m,currentiFte) 
      end do !i
      end do !iP_x
      end do !m
      end if !PSV
      if (SH) then
      do m = 1,npixX
      do iP_x = mPtini,mPtfin
          fotogramas(iP_x-nIpts,m,J,3,currentiFte) = allpoints(iP_x)%Wmov(J,3,m,currentiFte) 
      end do !iP_x
      end do !m
      end if !SH
      end if !makevideo
      end if !workboundary
      
      if (comoFacDeAmpliDinamica) then
      do dir= 1,3 !x,y,z direction of force application
        if(dir .eq. 2) then
           if(skipdir(dir)) cycle
        else ! 1 o 3
           if(skipdir(1) .and. skipdir(3)) cycle
        end if! dir
        if(.not. skipdir(dir)) then 
      
      if (PSV) then
      do iP_x = iPtini,iPtfin !cada receptor X  
          p_x => allpoints(iP_x)
          call G0estr(MecElem,p_x,J,cOME,dir)
          allpoints(iP_x)%facAmpli(J,currentiFte)%W = allpoints(iP_x)%resp(J,currentiFte)%W / MecElem(1)
          allpoints(iP_x)%facAmpli(J,currentiFte)%U = allpoints(iP_x)%resp(J,currentiFte)%U / MecElem(2)
          allpoints(iP_x)%facAmpli(J,currentiFte)%sxx = allpoints(iP_x)%resp(J,currentiFte)%sxx / MecElem(5)
          allpoints(iP_x)%facAmpli(J,currentiFte)%szx = allpoints(iP_x)%resp(J,currentiFte)%szx / MecElem(4)
          allpoints(iP_x)%facAmpli(J,currentiFte)%szz = allpoints(iP_x)%resp(J,currentiFte)%szz / MecElem(3)
      end do ! iP_x
      end if ! PSV 
      
         end if !skipdir
      end do !dir
      end if ! comoFacDeAmpliDinamica
     
      if (onlythisJ) then
          write(PrintNum,'(a)') ""
          write(PrintNum,'(a)') "Rultados en puntos receptores: ------"
          do iP_x = 1,nPts
           if ((.not. allpoints(iP_x)%isSabana) .and. &
               (.not. allpoints(iP_x)%guardarMovieSiblings)) then
            write(PrintNum,'(a,i0,a,F5.1,a,F5.1,a,a,F3.1,a,F3.1,a)')&
            "ip",ip_x,"[",allpoints(iP_x)%center%x,",",allpoints(iP_x)%center%z,"] ", &
            "n = [",allpoints(iP_x)%normal%x,",",allpoints(iP_x)%normal%z,"]"
            write(PrintNum,'(a,i0,a,i0)')&
            "reg",allpoints(iP_x)%region," layer",allpoints(iP_x)%layer
            if (allpoints(iP_x)% guardarMovieSiblings) then
            if (PSV) then
            write(PrintNum,*)&
            "   W=",allpoints(iP_x)%Wmov(J,1,1:npixX,currentiFte)
            write(PrintNum,*)&
            "   U=",allpoints(iP_x)%Wmov(J,2,1:npixX,currentiFte)
            else !SH
            write(PrintNum,*)&
            "   V=",allpoints(iP_x)%Wmov(J,3,1:npixX,currentiFte)
            end if
            else
            if (PSV) then
            write(PrintNum,*)&
            "   W=",allpoints(iP_x)%resp(J,currentiFte)%W
            write(PrintNum,*)&
            "   U=",allpoints(iP_x)%resp(J,currentiFte)%U
            write(PrintNum,*)&
            "   Tz=",allpoints(iP_x)%resp(J,currentiFte)%Tz
            write(PrintNum,*)&
            "   Tx=",allpoints(iP_x)%resp(J,currentiFte)%Tx
            write(PrintNum,*)&
            "   sxx=",allpoints(iP_x)%resp(J,currentiFte)%sxx
            write(PrintNum,*)&
            "   szx=",allpoints(iP_x)%resp(J,currentiFte)%szx
            write(PrintNum,*)&
            "   szz=",allpoints(iP_x)%resp(J,currentiFte)%szz
            
            ! a polares locales y guardar:
!           print*," Aguas con la transf de coordenadas lin929"
!           allpoints(iP_x)%sinT = (980 - allpoints(iP_x)%center%z)/5.5
!           allpoints(iP_x)%cosT = (10 - allpoints(iP_x)%center%x)/5.5
!           
!           outvar(iP_x,1) = UR*10. - allpoints(iP_x)%center%x
!           outvar(iP_x,2) = UR*980. - allpoints(iP_x)%center%z
!           
!           outvar(iP_x,3) = allpoints(iP_x)%resp(J,currentiFte)%sxx
!           outvar(iP_x,4) = allpoints(iP_x)%resp(J,currentiFte)%szx
!           outvar(iP_x,5) = allpoints(iP_x)%resp(J,currentiFte)%szz
            else !SH
            write(PrintNum,*)&
            "   V=",allpoints(iP_x)%resp(J,currentiFte)%V
            write(PrintNum,*)&
            "   Ty=",allpoints(iP_x)%resp(J,currentiFte)%Ty
            end if
            
!           outvar(iP_x,3) = abs( &
!            allpoints(iP_x)%sinT **2 * allpoints(iP_x)%resp(J,currentiFte)%sxx + &
!            allpoints(iP_x)%cosT **2 * allpoints(iP_x)%resp(J,currentiFte)%szz - &
!            2*allpoints(iP_x)%sinT*allpoints(iP_x)%cosT*allpoints(iP_x)%resp(J,currentiFte)%szx &
!            ) !stt
            end if ! guardarMovieSiblings
           end if ! isSabana
          end do !iP_x
!     if (PSV) then    
!     open(169,FILE= 'outvar.m',action="write",status="replace")
!     write(arg,'(a)') "out"
!     call scripToMatlabMNmatrixZ(nPts-1,5,outvar(2:nPts,1:5),arg,169)
!     close(169)
!     end if
      end if! onlythisJ
      
      
      
      ! sabana en frecuencia eta
      if (PrintEtas .or. onlythisJ) then 
        call chdir(trim(adjustl(rutaOut)))
        write(arg,'(a,I0)') 'phi',currentiFte
        CALL chdir(trim(arg))
        if (PSV) then
        write(tt,'(a,I0,a)') "W_eta",J,".pdf"
        call plot_at_eta(frecIni,tt)
        write(tt,'(a,I0,a)') "U_eta",J,".pdf"
        call plot_at_eta(frecIni,tt)
        else
!       print*,"print at eta SH not implemented line966"
        write(tt,'(a,I0,a)') "V_eta",J,".pdf"
        call plot_at_eta(frecIni,tt)
        end if
        CALL chdir("..");CALL chdir("..")
        if ((onlythisJ .eqv. .true.) .and. (currentiFte .eq. nFuentes))then 
          write(PrintNum,'(/,a,F7.3)')&
         "at normalized frequency eta=",abs(come*1.0/(pi*beta(N+1)))
          Stop "-f argum finish"
        end if!
        if (onlythisJ) cycle
      end if
      
      ! w->time & plot
      if (J .eq. frecEnd) then 
      print*,"w->time & plot"
      if (plotFKS) then  
           write(xAx,'(a)')"frec [Hz]"
           write(yAx,'(a)')" K [1/m] "
           call chdir(trim(adjustl(rutaOut)))
           write(arg,'(a,I0)') 'traces',currentiFte
           CALL chdir(trim(arg))
         do iP = iPtini,iPtfin
           if (allpoints(iP)% guardarFK) then
             if (SH) then
             write(tt,'(a,I0)')"0_FK_",iP
             call plotFK(allpoints(iP)%FK(1:NFREC+1,1:NMAX,1:3), &
                         real(allpoints(iP)%center%x,4), & 
                         real(allpoints(iP)%center%z,4), & 
                         tt,xAx,yAx,PrintNum,3,3, onlythisJ, frecEnd)
             end if !
             if (PSV) then
             write(tt,'(a,I0)')"0_FK_",iP
             call plotFK(allpoints(iP)%FK(1:NFREC+1,1:NMAX,1:3), &
                         real(allpoints(iP)%center%x,4), & 
                         real(allpoints(iP)%center%z,4), & 
                         tt,xAx,yAx,PrintNum,1,2, onlythisJ, frecEnd)
             
             open(421,FILE= "outA_w.m",action="write",status="replace")
             write(arg,'(a,I0,a)')"FK_",iP,"_w"
             call scripToMatlabMNmatrixZ(NFREC+1,NMAX,& 
                                         allpoints(iP)%FK(1:NFREC+1,1:NMAX,2),& 
                                         arg,421)
             close(421)
             end if
           end if
         end do
           CALL chdir("..");CALL chdir("..")
      end if
      
      write(6,'(a)')"Printing seismograms, etc..."
      
      ! mostrar sismogramas en los puntos de interés
           call chdir(trim(adjustl(rutaOut)))
           write(arg,'(a,I0)') 'traces',currentiFte
           CALL chdir(trim(arg))
       
      call plotSisGram(PSV,SH,.true.)    
      if (plotFKCK) call F_K_exp(XF)
      CALL chdir("..")
         
      if (makeVideo) then 
        write(arg,'(a,I0)') 'video',currentiFte
        CALL chdir(trim(arg))
        call crepa_four_fotogramas
        if (PSV .and. vivaChurubusco) call Churubusco(.false.)
        if (SH .and. vivaChurubusco) call Hollywood(3)
        call chdir("..")
      end if
      !
      if (workboundary .and. punEnlaFront .and. vivaCine) then
        write(arg,'(a,I0)') 'video',currentiFte
        CALL chdir(trim(arg))
        if (PSV) call CINETECA
        CALL chdir("..")
      end if!
           
      end if !J = frecEn
      end do !currentiFte
      END DO ! J: frequency loop
      
!     deallocate(B);deallocate(IPIV)
      if(verbose >= 1) write(PrintNum,'(a)')" done"
      ! finish
      call ETIME(tarray, result)
      if (result .ge. 60) then
      write(PrintNum,'(a,f10.3,a)') "Elapsed ",result/60,"minutes"
      write(6,'(a,f10.3,a)') "Elapsed ",result/60,"minutes"
      else
      write(PrintNum,'(a,f10.3,a)') "Elapsed ",result,"seconds"
      write(6,'(a,f10.3,a)') "Elapsed ",result,"seconds"
      end if
      
      call vaciarWisdom
      Write(PrintNum,'(a)') ' done '
      Write(6,'(a)') ' done '
      close(PrintNum)
      
      END program
      
      
! pointer table 
      subroutine preparePointerTable(pota,firstTime,smallestWL)
      use resultVars, only : nPts,allpoints,nBpts,BouPoints, & 
                             fixedpota,nZs,Punto,n_top_sub,n_con_sub
      use debugstuff
      use Gquadrature, only : Gquad_n
      use glovars, only : verbose,workBoundary, rutaOut
      use waveNumVars, only : frec
      ! tabla con las coordenadas Z (sin repetir).
      implicit none
      logical,intent(in) :: firstTime
      real*8,intent(in) :: smallestWL
      integer :: i,j,iP!,e
      logical :: nearby,esnueva
      type(Punto), pointer :: PX
      integer, allocatable, dimension(:,:) :: pota,auxpota
      real*8 :: bola
      character(LEN=100) :: titleN
      bola = min(0.01_8*smallestWL,0.01_8)
      ! si es la primera vez que corre sólo agregamos los allpoints 
      if (firstTime) then
       if (verbose .ge. 4) print*, "generating master table"
       allocate(auxpota(npts,npts+2))
       auxpota = 0
       ! siempre agregamos el primer punto [ allpoints(1)% ]
       nZs = 1
       auxpota(1,1) = 1 !nXs allpoints
       auxpota(1,3) = 1 !4,5,6 ... allpoints -7,-8,-9 ... boupoints
       !         ^--- el (3) siempre está
       ! agregamos los demás sin repetir
       if (nPts>=2) then
       do iP = 2,nPts
         if (allpoints(ip)%region .ne. 1) then 
          if (.not. allpoints(ip)%isOD) &
           cycle !nada mas agregar receptores en medio estratificado
         end if
         esnueva = .true.
         do i = 1,nZs
           if (nearby(allpoints(iP)%center%z, & 
               allpoints(auxpota(i,3))%center%z,0.01_8)) then
             !agregar a coordenada Z existente
             auxpota(i,1) = 1 + auxpota(i,1)
             auxpota(i,auxpota(i,1) + 2) = iP
             esnueva = .false.
             exit
           end if
         end do
         if (esnueva) then
           !nueva coordenada Z
           nZs = nZs + 1
           auxpota(nZs,1) = 1
           auxpota(nZs,3) = iP
         end if
       end do
       end if
       j = maxval(auxpota(:,1))
       
       allocate(PoTa(nZs,2+j))
       Pota = auxpota(1:nZs,1:2+j)
       deallocate(auxpota)
       if (workBoundary) then
         allocate(fixedPoTa(nZs,2+j))
         fixedPota = Pota
       end if
!      call showMNmatrixI(nZs,2+j,fixedPota,"po_ta",6);stop
      else !--first time----------------------------------------------------------------------------
       ! Dada la tabla de los puntos fijos.
       ! Agrega los centros de los segmentos que se resuelven con DWN. 
       if (verbose .ge. 4) print*, "updating master table"
       if (allocated(pota)) deallocate(pota)
       if (allocated(auxpota)) deallocate(auxpota)
       allocate(auxpota(nZs + nBpts * Gquad_n, 2 + & 
                        maxval(fixedPota(:,1)) + nBpts * Gquad_n))
       auxpota = 0
       nZs = size(fixedPota,1)
!      call showMNmatrixI(size(fixedPota,1),size(fixedPota,2),fixedPota,"fixed",6)
       auxpota(1:size(fixedPota,1),1:size(fixedPota,2)) = fixedPota
       ! *****************************************************************************
       do iP = 1,n_top_sub + n_con_sub ! (Sólo lo que se resuelve con DWN) antes:nBpts
         esnueva = .true.
         do i = 1,nZs !para cada profundidad ya identificada en la tabla
          ! diferenciar de que grupo es la coordenda
           if (associated(PX)) nullify(PX)
           if (auxpota(i,3) .gt. 0) then
             PX => allpoints(auxpota(i,3))!; print*,"allp "
           else
             PX => boupoints(abs(auxpota(i,3)))!; print*,"boup "
           end if
          ! revisar si en la tabla ya hay una profundidad cercana 
           if (nearby(boupoints(iP)%center%z,PX%center%z,bola)) then
!            print*,"estan cerca"
             if (auxpota(i,3) .lt. 0) then !negativo => boundary
               if (.not.(nearby(PX%length,boupoints(iP)%length,bola)) .or. &
              (.not.(nearby(abs(PX%normal%x), & 
              abs(boupoints(iP)%normal%x),bola)))) then
                 esnueva = .true. !;print*," they are different"
                 exit
               end if
             end if
             !inscribir a coordenada Z existente
             auxpota(i,2) = 1 + auxpota(i,2)
             auxpota(i,auxpota(i,1) + auxpota(i,2) + 2) = - iP
             esnueva = .false.
!            print*,"inscrito. renglon:",i,"(",auxpota(i,1)," ",auxpota(i,2),")"
             exit
           end if !nearby
         end do ! i,nZs
         if (esnueva) then
           !nueva coordenada Z
           nZs = nZs + 1
           auxpota(nZs,1) = 0
           auxpota(nZs,2) = 1
           auxpota(nZs,3) = - iP
!          print*,"nuevo Z ahora son ",nZs
         end if
       end do ! iP
       j = maxval(auxpota(:,1) + auxpota(:,2))
       allocate(PoTa(nZs,2+j))
       Pota = auxpota(1:nZs,1:2+j)
       deallocate(auxpota)
       if(verbose .ge. 3) then
         call chdir(trim(adjustl(rutaOut)))
         CALL chdir("insbackup")
         call system("mkdir pota")
         call chdir("pota")
         write(titleN,'(a,I0,a)') "pota_",int(frec*100),"E-100HZ.txt"
         open(35303,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         call showMNmatrixI(nZs,2+j,pota,"po_ta",35303)
         close(35303)
         call chdir("..");call chdir("..");call chdir("..")
       end if
      end if !firsttime
      if (verbose .ge. 4) call showMNmatrixI(nZs,2+j,pota,"po_ta",6)
      if (verbose .ge. 4) print*, "out of preparePointerTable"
      end subroutine preparePointerTable
      
      function nearby(a,b,bola)
      implicit none
      real*8, intent(in) :: a,b,bola
      logical :: nearby
      nearby = .false.
      if ((b-bola .le. a) .and. (a .le. b+bola)) then 
        nearby = .true.
      end if
!     print*, b-bola, "<=", a, "<=",b+bola ," :: ", nearby
      end function nearby 
      
      subroutine asociar(PX, i_Fuente ,itabla_z, itabla_x)
      use resultVars, only : pota,Punto,allpoints,boupoints
      use sourcevars, only : Po
      implicit none
      type(Punto), pointer :: PX
      integer, intent(in) :: itabla_x, itabla_z, i_Fuente
        nullify(PX)
        if (itabla_z .ne. 0) then
          if (pota(itabla_z, itabla_x) .gt. 0) then
             PX => allpoints(pota(itabla_z, itabla_x))
          else
             PX => boupoints(abs(pota(itabla_z, itabla_x)))
          end if
        else
          if (i_Fuente .eq. 0) stop "i_fuente = 0 en asociar"
          PX => Po(i_Fuente)
        end if
      end subroutine asociar
      
      function porLoMenosUnoEsEstr(itabla_z)
      use resultVars, only : pota,Punto,allpoints
      implicit none
      type(Punto), pointer :: PX
      logical :: porLoMenosUnoEsEstr
      integer :: itabla_z ,itabla_x
      
      porLoMenosUnoEsEstr = .true.
      if (pota(itabla_z,2) .gt. 0) return !si hay una frontera
      
      porLoMenosUnoEsEstr = .false.
      do itabla_x = 3,2+ pota(itabla_z,1) ! si hay un receptor en medio estratificado
       nullify(PX)
       PX => allpoints(pota(itabla_z,itabla_x))
       if (PX%region .eq. 1) porLoMenosUnoEsEstr = .true. !'estr'
      end do
      
      end function porLoMenosUnoEsEstr
      subroutine subdivideTopo(iJ,frec, onlythisJ ,minBeta,bet,nbpts,BouPoints)
      !optimized segementation of geometry.
      use gloVars, only : multSubdiv,V => verbose,makevideo,z0, rutaOut
      use GeometryVars, only : origGeom,nXI, segemntedcoords, & 
       n_topo,n_cont,staywiththefinersubdivision, finersubdivisionJ,fraccionDeSmallestWL_segm_de_esquina
      use waveNumVars, only : dfrec,nfrec, smallestWL
      use resultVars, only : Punto,nIpts,nMpts,n_top_sub,n_con_sub,n_val_sub
      use meshVars, only : npixX
      use ploteo10pesos
      use Gquadrature, only : Gquad_n
      use soilvars, only : N
      implicit none
      interface
         subroutine punGa (BP) 
          use resultVars, only : Punto 
          type (Punto), pointer :: BP
         end subroutine punGa
         
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente
      end subroutine drawBoundary
         
      end interface 
      
      integer, intent(in) :: iJ
      real*8, intent(in) :: frec
      logical :: onlythisJ
      integer :: nbpts
      real*8  :: minBeta
      real*8,dimension(:) :: BET
      type (Punto), dimension(:), allocatable, target ::  BouPoints
      type (Punto), dimension(:), pointer :: BP
      type (Punto), pointer :: BPi
      
      integer :: iXI,nMxDeDivi,iDi,nsubdivs,ndivsiguales,nSegmeTotal
      integer :: bou_conta_ini, bou_conta_fin,nsubSegme
      real*8 :: sixthofWL,f_frec,resi,deltax,deltaz
      character(LEN=100) :: txt
      character(LEN=3) :: extension
      real*8     :: errT = 0.01_8
      type (segemntedcoords), dimension(:), allocatable :: subdiv
      
      f_frec = frec
      smallestWL = minBeta / f_frec
      ! que la subdivision minima sea la correspondiente a 1/3 de fmax
      if ((onlythisJ .eqv. .false.) .and. (iJ * 1.0 < 0.33 * nfrec)) then
        f_frec = 0.33 * nfrec * dfrec
        smallestWL = minBeta / f_frec
      end if
      
      ! para forzar que el refinamiento sea al menos tan fino como en la frecuencia finersubdivisionJ
      if (staywiththefinersubdivision) then 
        f_frec = min(finersubdivisionJ , nfrec)*dfrec
        if ((allocated(boupoints)) .and. (iJ .le. finersubdivisionJ) .and. (iJ .lt. nfrec)) go to 938
      end if
      ! if we keep segmenting
      if (allocated(subdiv)) then
        nsubdivs = size(subdiv)
        do ixi = 1,nsubdivs
          if(allocated(subdiv(ixi)%x)) deallocate(subdiv(ixi)%x)
          if(allocated(subdiv(ixi)%z)) deallocate(subdiv(ixi)%z)
        end do
      else
        allocate(subdiv(nxi))
      end if!
      if(V .ge. 4) print*,"nXI=",nXI
      nSegmeTotal = 0
      ! dividimos cada elemento original
      do iXI = 1,nXI
      if(V .ge. 4) print*,"dividing big segment",iXI
      if (iXI .le. n_topo) then 
        sixthofWL = real(BET(origGeom(iXI)%layer),8) / f_frec / multSubdiv 
      else if ((n_topo+1 .le. iXI) .and. (iXI .le. n_topo+ n_cont)) then
        sixthofWL = min(real(BET(origGeom(iXI)%layer),8),real(BET(N+2),8)) / f_frec / multSubdiv
      else ! n_topo+n_cont+1 .le. iXI 
        sixthofWL = real(BET(N+2),8) / f_frec / multSubdiv
      end if
      ! pero si el elemento es de esquina, que sea más manchado
      if (origGeom(iXI)%segmentoDeEsquina) sixthofWL = sixthofWL / fraccionDeSmallestWL_segm_de_esquina
      
      if(V .ge. 4) print*,"sixthofWL =", sixthofWL
      
      ! delta normalizado y luego del tamaño sixthofWL
      deltaX = (origGeom(iXI)%bord_B%x - origGeom(iXI)%bord_A%x)
      deltaX = deltaX / origGeom(iXI)%length * real(sixthofWL,4)
      deltaZ = (origGeom(iXI)%bord_B%z - origGeom(iXI)%bord_A%z)
      deltaZ = deltaZ / origGeom(iXI)%length * real(sixthofWL,4)
      nMxDeDivi = ceiling(origGeom(iXI)%length/sixthofWL)
      if (V .ge. 4) then
      print*,""
      print*,frec,'Hz'
      print*,"on layer=", origGeom(iXI)%layer
      print*,real(BET(origGeom(iXI)%layer),8),'m/s'
      print*,"forsegment ",iXI
      print*,"of lenght ",origGeom(iXI)%length
      print*, sixthofWL*multSubdiv,'m is the wavelenght'
      print*, "fract. of wavelen = ",sixthofWL,'m'
      print*,"the it will be diveded nmxdedivi =", nMxDeDivi
      end if!
      
      if(V .ge. 4) print*,"segm,",iXI,"of L =",origGeom(iXI)%length, & 
      "  wavelen/6=",sixthofWL,  " nMxDeDivi=", nMxDeDivi
      
      if (origGeom(iXI)%length .le. sixthofWL ) then! no hace falta segm
        if(V .ge. 4) print*,"no hace falta segmentar"
        allocate(subdiv(ixi)%x(2));allocate(subdiv(ixi)%z(2))
        subdiv(ixi)%x(1) = origGeom(iXI)%bord_A%x !x
        subdiv(ixi)%z(1) = origGeom(iXI)%bord_A%z !z
        subdiv(ixi)%x(2) = origGeom(iXI)%bord_B%x !x
        subdiv(ixi)%z(2) = origGeom(iXI)%bord_B%z !z
        nsubdivs = 2
        nSegmeTotal = nSegmeTotal + 1
        if(V .ge. 4) print*,"nSegmeTotal=", nSegmeTotal
      else ! origGeom(iXI)%length .gt. sixthofWL       ! si hace falta segmentar 
        if(V .ge. 4) print*,"si hace falta segmentar"
        if (nMxDeDivi .lt. 3) nMxDeDivi = 3
        ! ¿Sólo dividir por la mitad u optimizando?
        if (origGeom(iXI)%length .le. 2.* sixthofWL) then
        if(V .ge. 4) print*,"     sólo por la mitad"
                                                             ! sólo por la mitad
        allocate(subdiv(ixi)%x(3));allocate(subdiv(ixi)%z(3))
        subdiv(ixi)%x(1) = origGeom(iXI)%bord_A%x !x
        subdiv(ixi)%z(1) = origGeom(iXI)%bord_A%z !z
        subdiv(ixi)%x(2) = (origGeom(iXI)%bord_A%x + origGeom(iXI)%bord_B%x)/2.
        subdiv(ixi)%z(2) = (origGeom(iXI)%bord_A%z + origGeom(iXI)%bord_B%z)/2.
        subdiv(ixi)%x(3) = origGeom(iXI)%bord_B%x !x
        subdiv(ixi)%z(3) = origGeom(iXI)%bord_B%z !z
        nsubdivs = 3
        nSegmeTotal = nSegmeTotal + 2
        else ! origGeom(iXI)%length .gt. 2* sixthofWL
        if(V .ge. 4) print*,"     segmentar chido"
                                                             ! segmentar chido
        ndivsiguales = floor(origGeom(iXI)%length/2./sixthofWL)
        resi = origGeom(iXI)%length - 2. * ndivsiguales * sixthofWL
!       print*,"ndivsiguales", ndivsiguales
!       print*,"resi", resi
        if(resi .gt. sixthofWL) then
          if(V .ge. 4) print*,"     dividir segmento central por la mitad"
          ! dividir segmento central por la mitad
          nsubdivs = 2*ndivsiguales + 2 + 1
          allocate(subdiv(ixi)%x(nsubdivs));allocate(subdiv(ixi)%z(nsubdivs))
          subdiv(ixi)%x(ndivsiguales+1) = origGeom(iXI)%bord_A%x + deltaX*ndivsiguales
          subdiv(ixi)%z(ndivsiguales+1) = origGeom(iXI)%bord_A%z + deltaZ*ndivsiguales
          subdiv(ixi)%x(ndivsiguales+3) = origGeom(iXI)%bord_B%x - deltaX*ndivsiguales
          subdiv(ixi)%z(ndivsiguales+3) = origGeom(iXI)%bord_B%z - deltaZ*ndivsiguales 
          subdiv(ixi)%x(ndivsiguales+2) = (subdiv(ixi)%x(ndivsiguales+1) + subdiv(ixi)%x(ndivsiguales+3))/2.
          subdiv(ixi)%z(ndivsiguales+2) = (subdiv(ixi)%z(ndivsiguales+1) + subdiv(ixi)%z(ndivsiguales+3))/2.
        else ! es más grande de 2.*sixthofWL
          if(V .ge. 4) print*,"     odd number of segments"
          ! odd number
          if (resi .gt. errT) then
           nsubdivs = 2*ndivsiguales + 1 + 1
           allocate(subdiv(ixi)%x(nsubdivs));allocate(subdiv(ixi)%z(nsubdivs))
           subdiv(ixi)%x(ndivsiguales+1) = origGeom(iXI)%bord_A%x + deltaX*ndivsiguales
           subdiv(ixi)%z(ndivsiguales+1) = origGeom(iXI)%bord_A%z + deltaZ*ndivsiguales
           subdiv(ixi)%x(ndivsiguales+2) = origGeom(iXI)%bord_B%x - deltaX*ndivsiguales
           subdiv(ixi)%z(ndivsiguales+2) = origGeom(iXI)%bord_B%z - deltaZ*ndivsiguales 
          else !por fin si solo por la mitad
           nsubdivs = 2*ndivsiguales + 1
           allocate(subdiv(ixi)%x(nsubdivs));allocate(subdiv(ixi)%z(nsubdivs))
           subdiv(ixi)%x(ndivsiguales+1) = (origGeom(iXI)%bord_A%x + deltaX*ndivsiguales + & 
                                            origGeom(iXI)%bord_B%x - deltaX*ndivsiguales)/2.
           subdiv(ixi)%z(ndivsiguales+1) = (origGeom(iXI)%bord_A%z + deltaZ*ndivsiguales + &
                                            origGeom(iXI)%bord_B%z - deltaZ*ndivsiguales)/2.
          end if
        end if ! resi .le. sixthofWL
        nSegmeTotal = nSegmeTotal + nsubdivs -1
        if(V .ge. 4) print*,"     los que son iguales a cada lado"
        do idi = 1, ndivsiguales ! los que son iguales a cada lado
            subdiv(ixi)%x(idi) = origGeom(iXI)%bord_A%x + deltaX *(idi-1)
            subdiv(ixi)%z(idi) = origGeom(iXI)%bord_A%z + deltaZ *(idi-1)
            subdiv(ixi)%x(nsubdivs +1 - idi) = origGeom(iXI)%bord_B%x - deltaX *(idi-1)
            subdiv(ixi)%z(nsubdivs +1 - idi) = origGeom(iXI)%bord_B%z - deltaZ *(idi-1)
        end do !idi
        end if ! (origGeom(iXI)%length .le. 2.* sixthofWL)
      end if ! origGeom(iXI)%length .le. sixthofWL
!     if (ij .lt. 161) stop "3895"
      
      if(V .ge. 4) then; DO idi = 1, nsubdivs
       print*,"(",subdiv(ixi)%x(idi),",",subdiv(ixi)%z(idi),")"
      end do; end if
      end do !iXI
      ! para no chorrear memoria
      if (allocated(boupoints)) then
        do idi = 1,size(boupoints)
          if (allocated(boupoints(idi)%Gq_xXx_coords)) deallocate(BouPoints(idi)%Gq_xXx_coords)
          if (allocated(boupoints(idi)%Gq_xXx_C)) deallocate(BouPoints(idi)%Gq_xXx_C)
        end do 
        deallocate(boupoints)
      end if
      
      allocate(boupoints(nSegmeTotal)) !numero de segmentos total
        do idi = 1,size(boupoints)
          allocate(BouPoints(idi)%Gq_xXx_coords(Gquad_n,2))
          allocate(BouPoints(idi)%Gq_xXx_C(Gquad_n))
        end do 
      nBpts = nSegmeTotal
      
      ! Heredar
      bou_conta_ini = 0
      bou_conta_fin = 0
      n_top_sub=0; n_con_sub=0; n_val_sub=0
      do ixi = 1,nXI !nuevamente para cada segmento original
      nsubSegme = size(subdiv(ixi)%x)-1
      if(V .ge. 4) print*,"Gau pst for segms in Big segm",ixi,"with",nsubSegme,"subs"
      bou_conta_fin = bou_conta_fin + nsubSegme
      
      BouPoints(1+bou_conta_ini:bou_conta_fin)%bord_A%x = subdiv(ixi)%x(1:nsubSegme)
      BouPoints(1+bou_conta_ini:bou_conta_fin)%bord_A%z = subdiv(ixi)%z(1:nsubSegme)
      BouPoints(1+bou_conta_ini:bou_conta_fin)%bord_B%x = subdiv(ixi)%x(2:nsubSegme+1)
      BouPoints(1+bou_conta_ini:bou_conta_fin)%bord_B%z = subdiv(ixi)%z(2:nsubSegme+1)
      BouPoints(1+bou_conta_ini:bou_conta_fin)%center%x = &
         (subdiv(ixi)%x(2:nsubSegme+1) + subdiv(ixi)%x(1:nsubSegme))/2.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%center%z = &
         (subdiv(ixi)%z(2:nsubSegme+1) + subdiv(ixi)%z(1:nsubSegme))/2.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%normal%x = origGeom(iXI)%normal%x
!     BouPoints(1+bou_conta_ini:bou_conta_fin)%normal%y = 0.0_8
      BouPoints(1+bou_conta_ini:bou_conta_fin)%normal%z = origGeom(iXI)%normal%z
      BouPoints(1+bou_conta_ini:bou_conta_fin)%length = sqrt(&
           (subdiv(ixi)%x(2:nsubSegme+1) - subdiv(ixi)%x(1:nsubSegme))**2 + &
           (subdiv(ixi)%z(2:nsubSegme+1) - subdiv(ixi)%z(1:nsubSegme))**2)
      
      BouPoints(1+bou_conta_ini:bou_conta_fin)%cosT = origGeom(iXI)%cosT
      BouPoints(1+bou_conta_ini:bou_conta_fin)%sinT = origGeom(iXI)%sinT
      BouPoints(1+bou_conta_ini:bou_conta_fin)%layer = origGeom(iXI)%layer
      BouPoints(1+bou_conta_ini:bou_conta_fin)%region = 1!'estr'
      BouPoints(1+bou_conta_ini:bou_conta_fin)%tipoFrontera = origGeom(iXI)%tipoFrontera
      if (origGeom(iXI)%tipoFrontera .eq. 0 ) n_top_sub= n_top_sub + nsubSegme
      if (origGeom(iXI)%tipoFrontera .eq. 1 ) n_con_sub = n_con_sub + nsubSegme
      if (origGeom(iXI)%tipoFrontera .eq. 2 ) n_val_sub = n_val_sub + nsubSegme
      BouPoints(1+bou_conta_ini:bou_conta_fin)%isBoundary = .true.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%isOnInterface = .false.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%guardarFK = .false.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%guardarMovieSiblings = .false.
      BouPoints(1+bou_conta_ini:bou_conta_fin)%isSourceSegmentForce = .false.
      
      BouPoints(1+bou_conta_ini:bou_conta_fin)%length = ( & 
        (subdiv(ixi)%x(2:nsubSegme+1) - subdiv(ixi)%x(1:nsubSegme))**2. + &
        (subdiv(ixi)%z(2:nsubSegme+1) - subdiv(ixi)%z(1:nsubSegme))**2. )** 0.5
          
      bou_conta_ini = bou_conta_ini + nsubSegme
      end do !iXI
      
      if (V .ge. 4) then;do idi = 1, nSegmeTotal
        print*,"-------division de puntos------------"
        print*,idi,":(", BouPoints(idi)%bord_A,") - (", & 
                         BouPoints(idi)%bord_B,") : ",BouPoints(idi)%length
        print*,"              (",BouPoints(idi)%center,")"
        print*,"       n = ", BouPoints(idi)%normal
        print*,"       e = ", BouPoints(idi)%layer
      end do;end if

 938  do ixi = 1,nBpts
        if(allocated(boupoints(ixi)%G)) deallocate(boupoints(ixi)%G)
        !                                 ,--- 1 a 3 direccion de la fuerza
        allocate(boupoints(ixi)%G(nIpts,9,3)); boupoints(ixi)%G = z0
        !                               `--- 1 a 6 elemento mecánico
        !                                          1 : W      4 : Tz    7 : sxx
        !                                          2 : U      5 : Tx    8 : szx
        !                                          3 : V      6 : Ty    9 : szz
      end do!
!     print*,n_top_sub,n_con_sub,n_val_sub,nbpts;stop 4126
      do ixi = 1,nBpts
        BouPoints(iXi)%boundaryIndex = iXi ! para indice de la matriz 
        BPi => BouPoints(iXi)
        call punGa(BPi) ! the Gaussian integration points
      end do
      
      if (makeVideo) then
      do ixi = 1,nBpts
        if(allocated(boupoints(ixi)%Gmov)) deallocate(boupoints(ixi)%Gmov)
        allocate(boupoints(ixi)%Gmov(nMpts,3,3,npixX)); boupoints(ixi)%Gmov = z0
      end do
      end if
      
      ! draw the subdivision
      if (V .ge. 2) then
       call chdir(trim(adjustl(rutaOut)))
       CALL chdir("subdivs")
       write(txt,'(a,I0,a)') 'Division_at[J=',iJ,'].pdf'
       write(extension,'(a)') 'PDF'
       BP => BouPoints
       call drawBoundary(BP,nBpts,txt, extension,.true.,.false.,.false.,.false.)
       CALL chdir("..") !out of subdivs
       CALL chdir("..") !out of outs
      end if!
      if (V .ge. 4) print*,"ended subdivideTopo"
      end subroutine subdivideTopo
      
      subroutine punGa (BP) 
      ! Coordenadas de los puntos de integración Gaussiana.
      use resultVars, only : Punto 
      use glovars, only : verbose
      use Gquadrature, only: Gqu_n => Gquad_n, & 
                             Gqu_t => Gqu_t_30, & 
                             Gqu_A => Gqu_A_30
      implicit none
      real*8, dimension(2) :: A,B
      real*8, dimension(:), pointer :: G_c
      real*8, pointer :: L
      real*8, dimension(:,:), pointer :: Gq
      real*8, dimension(2) :: norm_comp
      real*8, dimension(2) :: ABp !x or y coords of points A and B
      integer :: i, xory, xoryOtro
      real*8 :: xA,yA,xB,yB
      real*8 :: interpol
      type (Punto), pointer :: BP
      
      A = (/ BP%bord_A%x, BP%bord_A%z /)
      B = (/ BP%bord_B%x, BP%bord_B%z /)
      L => BP%length
      Gq => BP%Gq_xXx_coords(1:Gqu_n,1:2)
      G_c => BP%Gq_xXx_C(1:Gqu_n)
      
      norm_comp(1)=abs(B(1)-A(1)) / L
      norm_comp(2)=abs(B(2)-A(2)) / L
      
      if (norm_comp(2) > norm_comp(1)) then
            !print*," la pendiente es mayormente vertical "
         xory = 2 
         xoryOtro = 1
      else 
            !print*," la pendiente es mayormente horizontal "
         xory = 1 
         xoryOtro = 2
      end if
      
      ABp(1) = A(xory)
      ABp(2) = B(xory)
      
      do i = 1,Gqu_n !ceros de Legendre (una coordenada):      
        Gq(i,xory) = (ABp(2)+ABp(1))/2. + (ABp(2)-ABp(1))/2. * Gqu_t(i)
        G_c(i) = abs(L)/2. * Gqu_A(i)
      end do
      ! la otra coordenada:
        xA = ABp(1)
        yA = A(xoryOtro)
        xB = ABp(2)
        yB = B(xoryOtro)
        do i = 1,Gqu_n
          Gq(i,xoryOtro) = interpol(xA,yA,xB,yB,Gq(i,xory))
        end do
      
      if (verbose .ge. 4) then
        write(6,'(a,F12.5,a,F12.5,a,F12.5,a,F12.5,a,F12.8,a,F12.6,a,F12.6,a)') & 
               "[",A(1),",",A(2),"]-[", & 
                 B(1),",",B(2), & 
                 "] L:", L, &
                 " n:[", BP%normal%x,",", BP%normal%z,"]"
        if (xory .eq. 1) print*,"mayormente horizontal"
        if (xory .eq. 2) print*,"mayormente vertical"          
        !print*,"{",xA,",",yA,"}-{",xB,",",yB,"} Gquad points:"
        do i = 1,Gqu_n
          print*,"Gq",i,"[", BP%Gq_xXx_coords(i,1), " , ", &
          BP%Gq_xXx_coords(i,2), "] :: ", BP%Gq_xXx_C(i)
        end do
        print*,""
      end if
      end subroutine punGa
      
      function interpol(xA,yA,xB,yB,x)
      implicit none
      real*8 :: interpol
      real*8, intent(in) :: xA,yA,xB,yB,x
      real*8 :: m
      !interpolación lineal de la ordenada de x que está entre A y B
      m = (yB-yA) / (xB-xA)
      interpol = yA + m * (x-xA)
      end  function interpol
      
      function thereisavirtualsourceat(iz)
      use resultVars, only : pota
      integer :: iz
      logical :: thereisavirtualsourceat
      thereisavirtualsourceat = .false.
      if (pota(iz,2) .gt. 0) then
         thereisavirtualsourceat = .true.
      end if
      end function thereisavirtualsourceat
      
      subroutine sourceAmplitudeFunction
      use wavelets !las funciones: ricke
      use waveVars, only : Dt,Uo,dt,maxtime!, ampfunction,Escala
      use waveNumVars, only : DFREC,nfrec, NPTSTIME,OMEI
      use gloVars, only : ve => verbose,Ur,Ui,z0,Printnum,rutaOut
      use sourceVars, only:Po,iFte=>currentiFte, nFuentes
      use ploteo10pesos
      implicit none
      integer  :: i,nval
      character(LEN=9)          :: logflag
      character(LEN=100)        :: titleN,xAx,yAx,CTIT
      CHARACTER(len=32) :: arg
      logical :: inquire, lexist, argumA
      real*8 :: val1,val2
      complex*16, dimension(NPTSTIME) :: S
      !complex*16 :: FFTW!(NPTSTIME)
      integer :: n_maxtime
      argumA = .false.
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then 
        if (trim(arg) .eq. '-a') then 
        iFte = 1
        ve = 6
        argumA = .true.
        end if
      end if
      
       n_maxtime = int(maxtime(iFte)/dt)
       if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
       if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
      
      ! grafica del factor de correción
      S = exp(- OMEI * Dt*((/(i,i=0, NPTSTIME-1)/)))
      write(titleN,'(a)') 'x-LevantonDWN.pdf'
      write(CTIT,'(a)') 'exp(-omei t)'
      write(xAx,'(a)') 't[sec]'
      write(yAx,'(a)') ''
        call plotXYcomp(S(1:n_maxtime),real(Dt,4), & 
         int(n_maxtime),titleN,xAx,yAx, CTIT ,1200,800,0.0)
         
         
!     real :: factor
      ! Amplitude function of incident wave
      ! prepare the signal we will use as incident wave amplitude.
      ! el tamaño del ricker es 2*NFREC porque incluirá el conjugado
      
 153  Uo(:,iFte)=z0
      if(Po(iFte)%ampfunction .eq. 1) then
        call ricker(Uo(:,iFte),Po(iFte)%Ts,Po(iFte)%Tp) ! Ricker wavelet saved on Uo
          if (ve .ge. 1) then
            write(Printnum,'(a)')'   Incident wave amplitude function: Ricker'
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-t.pdf'
            write(CTIT,'(a)') 'WaveAmplitude of Ricker wavelet'
            xAx = 'time[sec]'
            write(yAx,'(a)') 'amplitude'
            call plotXYcomp(Uo(1:n_maxtime,iFte),real(Dt,4), & 
                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
          end if
        ! forward
        Uo(:,iFte) = FFTW(NPTSTIME,Uo(:,iFte),-1,Dt)     
        Uo(1,iFte) = 0   
          if (ve .ge. 1) then
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
            xAx = 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'      
!           logflag = 'logx     '      
!           logflag = 'none     '
            call plotXYcomp(Uo(:,iFte),real(DFREC,4), & 
                 size(Uo(:,iFte)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), & 
!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
          end if
     !********************************************************************************
      elseif(Po(iFte)%ampfunction .eq. 2) then ! Gaussian
        call gaussian(Uo(:,iFte),Po(iFte)%sigGaus)
          if (ve .ge. 1) then
           write(Printnum,'(a)')'   Incident wave amplitude function: Gaussian'
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
            write(CTIT,'(a)') 'WaveAmplitude of Gaussian wavelet'
            write(xAx,'(a)') 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'
!           call plotXYcomp(Uo(:,iFte),real(DFREC,4), & 
!                n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!           print*,Uo(:,iFte)
            call plotXYcomp(Uo(:,iFte),real(DFREC,4), & 
                 size(Uo(:,iFte)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), & 
!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
          end if
      elseif(Po(iFte)%ampfunction .gt. 20) then ! GaussianCúbico
        call gaussian(Uo(:,iFte),Po(iFte)%sigGaus)
        Uo(:,iFte) = Uo(:,iFte) ** (Po(iFte)%ampfunction-20)
          if (ve .ge. 1) then
           write(Printnum,'(a,I0)')'   Incident wave amplitude function: Gaussian ',int((Po(iFte)%ampfunction-20))
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
            write(CTIT,'(a,I0,a)') 'WaveAmplitude of Gaussian',int((Po(iFte)%ampfunction-20)),' wavelet'
            write(xAx,'(a)') 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'
!           call plotXYcomp(Uo(:,iFte),real(DFREC,4), & 
!                n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
!           print*,Uo(:,iFte)
            call plotXYcomp(Uo(:,iFte),real(DFREC,4), & 
                 size(Uo(:,iFte)),titleN,xAx,yAx,CTIT,1200,800,0.0)
!           call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), & 
!           titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
          end if
      elseif(Po(iFte)%ampfunction .eq. 3) then ! inDispl.txt
     !********************************************************************************
        CALL chdir("..")
        CALL chdir("ins")
        inquire(file="inAmplitude.txt",exist=lexist)
        if (lexist) then
        OPEN(77,FILE="inAmplitude.txt",FORM="FORMATTED")
        else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "inAmplitude.txt" on input directory'
        end if
        READ(77,*)
        READ(77,*) inquire
        READ(77,*) nval
        if (nval .gt. NPTSTIME) then 
        print*, "**** warning inAmplitude trimmed to ", NPTSTIME
        nval = NPTSTIME
        end if
        Uo(:,iFte) = 0
        if ( inquire) then
        do i = 1,nval
        READ(77,*) val1
        Uo(i,iFte) = UR * val1
        end do
        close(77)
          if (ve .ge. 1) then
            CALL chdir("..")
            call chdir(trim(adjustl(rutaOut)))
            write(Printnum,'(a)')'   Incident wave amplitude function from file'
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-t.pdf'
            write(CTIT,'(a)') 'WaveAmplitude'
            xAx = 'time[sec]'
            write(yAx,'(a)') 'amplitude'
            call plotXYcomp(Uo(1:NPTSTIME,iFte),real(Dt,4), & 
                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
          end if
        
        ! forward
        Uo(:,iFte) = FFTW(NPTSTIME,Uo(:,iFte),-1,Dt) 
        
          if (ve .ge. 1) then
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
            xAx = 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'      
            logflag = 'logx     '      
!           logflag = 'none     '
            call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), & 
            titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
            CALL chdir("..")
            CALL chdir("ins")
          end if
        else
        do i = 1,nval
        READ(77,*) val1, val2
        Uo(i,iFte) = UR * val1 + UI * val2
        end do
        if (ve .ge. 1) then
            CALL chdir("..")
            call chdir(trim(adjustl(rutaOut)))
            write(titleN,'(a,I0,a)') 'x-amp',iFte,'-f.pdf'
            xAx = 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'      
            logflag = 'logx     '      
!           logflag = 'none     '
            call plotSpectrum(Uo(:,iFte),real(DFREC,4), size(Uo(:,iFte)),int(size(Uo(:,iFte))/2), & 
            titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
            CALL chdir("..")
            CALL chdir("ins")
          end if
        close(77)
        end if
        
        CALL chdir("..")
        call chdir(trim(adjustl(rutaOut)))
        
      else  ! DIRAC -----------------------------------------
       if (ve .ge. 1) write(Printnum,'(a)')'   Incident wave amplitude function: Dirac delta'
        Uo(:,iFte)=UR
      end if !tipoPulso
      Uo(:,iFte) = Uo(:,iFte)*Po(iFte)%Escala !escala de la señal
      write(Printnum,'(a)') &
       "---------------------------------------------------------------------------------"
      if (argumA .eqv. .true.) then 
      if (iFte .eq. nfuentes) then
      stop "argumento -a"
      else
      iFte = iFte + 1
      go to 153
      end if
      end if
      end subroutine sourceAmplitudeFunction
      subroutine diffField_at_iz(i_zF,dir_j,J,cOME_in)
!#define ver 1
      ! esta función es llamada con cada profundidad donde hay
      ! por lo menos una fuente.
      use gloVars, only: z0, plotFKS,UI,UR,PWfrecReal
      use resultVars, only : pota,Punto,nZs,MecaElem,FFres
      use refSolMatrixVars, only : B,Ak
      use waveNumVars, only : NMAX,k_vec,dk,vecNK,SpliK,OME!,DFREC
      use wavelets 
      use dislin
      use sourceVars, only: Po,iFte=>currentiFte!nFuentes,tipofuente, , PW_pol
      use soilvars, only:N,Z,alfa0,beta0,alfa,beta
      use, intrinsic :: iso_c_binding!, only : C_INT
      use debugStuff
      implicit none
      interface
        include 'interfaz.f'
        function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k,ik)
          use soilvars, only:N
          complex*16, dimension(1:5) :: PSVdiffByStrata
          real*8, intent(in)           :: z_i,k
          complex*16, intent(in)       :: cOME_i  
          integer, intent(in)          :: e,ik
          complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
        end function PSVdiffByStrata
        
        function SHdiffByStrata(coefOndas_SH,z_i,e,cOME_i,k,ik)
          use soilvars, only:N
          complex*16, dimension(1:3) :: SHdiffByStrata
          real*8, intent(in)           :: z_i,k
          complex*16, intent(in)       :: cOME_i  
          integer, intent(in)          :: e,ik
          complex*16, dimension(2*N+1),intent(in) :: coefOndas_SH
        end function SHdiffByStrata
        
        subroutine  eGAeNU(i_zF,ik,pXI,dj)
        use resultVars, only : Punto
          integer :: ik,i_zF,dj
          type(Punto), pointer :: pXi
        end subroutine  eGAeNU
      end interface
      integer, intent(in) :: i_zF,dir_j,J
      complex*16, intent(in),target  :: cOME_in
      logical,pointer :: intf
      integer, pointer :: ef
      real*8,target :: k
      real*8, pointer :: zf,xf,pt_k
      complex*16, dimension(:,:), allocatable, target :: auxK,savedAuxK
      complex*16, target  :: cOME,alf,bet
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:),allocatable,target :: ipivA
      integer, dimension(:),pointer :: pt_ipivA
      complex*16, dimension(:),allocatable,target :: workA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16,dimension(:,:),pointer :: pointA
      type(Punto), pointer :: pXi,p_X
      logical :: auxLogic ,porLoMenosUnoEsEstr,isPW
      integer :: ik,tam,itabla_z,itabla_x,iMec,mecS,mecE,&
                 nXis,n_Xs,iXi,dj,pos,ne!,i_Fuente,i_FuenteFinal
!     type(MecaElem)  :: Meca_diff, SHdiffByStrata
#ifdef ver
      character(LEN=32) :: arg
      real :: result,lastresult
      real, dimension(2) :: tarray
#endif
      allocate(auxK(2*nmax,5)); allocate(savedAuxK(2*nmax,5))
      isPW = .false.
      if (i_zF .eq. 0) then
      if (Po(iFte)%tipofuente .eq. 1) then
      isPW = .true.
      end if
      end if
      
      cOME = cOME_in 
      dj = dir_j; if(dj .eq. 3) dj = 2 
      if (i_zF .eq. 0) then
         itabla_x = 3 ! En la tabla (pota) de índices: la fuente real -> (0,3)
!        i_FuenteFinal = 1 ! Cantidad de fuentes reales
         if (isPW) then ! onda plana incidente·p
            ! con incidencia de onda plana no usamos atenuación
           if (PWfrecReal) then 
             cOME = OME * UR  !real(cOME_in) * UR!
             alf = alfa0(N+1)
             bet = beta0(N+1)
           else
             cOME = cOME_in  !real(cOME_in) * UR!
             alf = alfa(N+1)
             bet = beta(N+1)
           end if
                                      k = real(cOME/bet)*sin(Po(iFte)%gamma) !SV,SH
           if(Po(iFte)%PW_pol .eq. 2) k = real(cOME/alf)*sin(Po(iFte)%gamma) ! P
!          if(Po(iFte)%PW_pol .eq. 1) k = real(cOME/bet)*sin(Po(iFte)%gamma)
!          if(Po(iFte)%PW_pol .eq. 2) k = real(cOME/alf)*sin(Po(iFte)%gamma)
!          if(Po(iFte)%PW_pol .eq. 3) k = real(cOME/bet)*sin(Po(iFte)%gamma)
         end if! ·································································n
      else; itabla_x = 2 + pota(i_zF,1) + 1 !    la primer fuente virtual
!           i_FuenteFinal = 1; 
      end if       !        a esa profundidad
      ! ------- para cada una de las fuentes en la tabla -------------------------
!     do i_Fuente = 1,i_FuenteFinal
       call asociar(pXi,iFte,i_zF,itabla_x) ! asociar apuntador a fuente [pXi]
!      if (J .eq. 1 .and. i_zF .eq. 0) print*,"zf=pXi%center%z",pXi%center%z
       xf=>pXi%center%x;zf=>pXi%center%z;ef=>pXi%layer;intf=>pXi%isOnInterface
!      print*,"pXi: x",pXi%center,pXi%normal,pXi%layer,pXi%isoninterface,pxi%tipofuente,pxi%region
#ifdef ver
       call ETIME(tarray, result)
       print*,"pXi: x",pXi%center%x,"z",pXi%center%z,& 
       "isboundary",pXi%isboundary,"isOnInterface",intf, result
       lastresult = result
#endif
      ! Si es la fuente real y es una onda plana no se usa el DWN. Se calcula para
      ! el número de onda horizontal asociado al ángulo de incidencia ············
      if (isPW) then ! onda plana incidente   ·
          if (dir_j .eq. 2) then! SH                                             ·o
            tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1!                       ·n
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·d
            pt_k => k; pt_come_i => cOME!                                        ·a
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call globalmatrix_SH(pointA,pt_k,0)!                                 ·p
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·l
            call SHvectorB_ondaplana(B(:,0),pxi%gamma)!                          ·a
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·n
          else!  P-SV                                                            ·a
            tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2!                       ·
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·
            pt_k => k; pt_come_i => cOME!                                        ·
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call gloMat_PSV(pointA,pt_k,0)!                                      ·
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·
            call PSVvectorB_ondaplana(B(:,0),pxi%gamma)!                         ·
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·
          end if!                                                                ·
          pos = 0; ne = 2*nmax+1
      else ! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · 
      ! La fuente es cilíndrica, puede tratarse de la fuente real o una virtual ·· 
        pos = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(pos-2)                    !
        B(:,pos:ne) = 0                                                           !
        if (dir_j .eq. 2) then                                                   !
         tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1                           !o
         if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) return
         
         ! si la fuente está sobre una interfaz ...............
!        if (pXi%isOnInterface) then
         do ik = 1,pos+1                                                            !n
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik),ik)             !d
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !a
         end do                                                                  !                                     
         do ik = ne,2*NMAX                                                       !c
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik),ik)             !i
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !l
         end do                                                                  !i
!        else ! la fuente está entre interfaces ..............
!        stop "fuerza entre interfaces no implementado para SH"
!          do ik = 1,pos+1
!            stop 1990
!            call eGAeNU(i_zF,ik,pXI,dj)
!          end do!                                                                 
!          do ik = ne,2*NMAX
!            call eGAeNU(i_zF,ik,pXI,dj)
!          end do
!        end if
        else!  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .     !n        
         tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2                           !d
         if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) return         
         
         ! si la fuente está sobre una interfaz ...............
         if (pXi%isOnInterface) then
           do ik = 1,pos+1                                                         !r
             call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !i
!            stop 2046
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !c
           end do                                                                  !a
           do ik = ne,2*NMAX                                                       !
             call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !
           end do
         else ! la fuente está entre interfaces ..............
           do ik = 1,pos+1
             call eGAeNU(i_zF,ik,pXI,dj)
!            print*,"Bik1=",B(:,ik)
!            if (ik .eq. 443) stop
           end do!                                                                 
           do ik = ne,2*NMAX
             call eGAeNU(i_zF,ik,pXI,dj)
           end do
!            print*,sum(B(:,1:pos+1)),sum(B(:,ne:2*nmax))
         end if
        end if                                                                   !
      end if! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·
#ifdef ver
       call ETIME(tarray, result)
       print*,"vector B,  x=inv(A)B (todas k) ",result-lastresult
       lastresult = result
       print*,"resultados para cada Z donde hay receptores:"
#endif
      ! resultados para cada Z donde hay receptores en el medio estratificado ...
      do itabla_z = 1,nZs
      !(en la tabla todos los puntos son receptores)
       auxLogic = porLoMenosUnoEsEstr(itabla_z)
#ifdef ver 
       print*,"  itabla_z",itabla_z,"de",nZs," -----------------------------"
       print*,"   pota(itabla_z,:)=",pota(itabla_z,:),auxLogic
#endif
       !
       if (auxLogic .eqv. .false.) cycle
        ! (el primer receptor de la tabla --------| a esa profundidad) 
            call asociar(p_X, 0, itabla_z, 3)
!           print*,"p_X: x",p_x%center%x,"z",p_x%center%z,& 
!      "isboundary",p_x%isboundary,"isOnInterface",p_x%isoninterface;stop
#ifdef ver 
       call ETIME(tarray, result)
       print*,"  p_x%center%z",p_x%center%z,"[m]  e=",p_X%layer
       lastresult = result
#endif
            if (p_X%layer .eq. N+2) cycle
            if (dir_j .eq. 2) then; mecS = 1; mecE = 3 !V,s32,s12
            else;                   mecS = 1; mecE = 5 !W,U,s33,s31,s11
            end if
            
      ! ... elementos mecánicos a la profundidad p_X%center%z ........
      ! .... usando los coeficientes de las ondas en el estrato ......
!     savedauxk = z0
!     savedauxk(po+1:ne-1,:) = 0
      if (isPW) then ! onda plana·············
          if (dir_j .eq. 2) then!                                               ·
!           if(Po(iFte)%PW_pol .eq. 3) k = OME/beta0(N+1)*sin(pXi%gamma)!       ·
!                Meca_diff = SHdiffByStrata(B(:,0), &!                          ·o
!                            p_X%center%z, p_X%layer, & !                       ·n
!                            cOME,k,mecS,mecE,outpf) !                          ·d
!                savedauxK(1,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE) !          ·a
                 savedauxK(1,1:3) = SHdiffByStrata(B(:,0), & 
                              p_X%center%z, p_X%layer,cOME,k,0)!
!               print*,"savedauxk(1,1:3)=",savedauxk(1,1:3)
          else !                                                                ·
!           if(Po(iFte)%PW_pol .eq. 1) k = OME/beta0(N+1)*sin(pXi%gamma)!       ·p
!           if(Po(iFte)%PW_pol .eq. 2) k = OME/alfa0(N+1)*sin(pXi%gamma)!       ·l
                 savedauxk(1,1:5) = PSVdiffByStrata(B(:,0), &!                  ·a
                              p_X%center%z, p_X%layer,cOME,k,0)!                ·a
!               print*,"savedauxk(1,1:5)=",savedauxk(1,1:5)
          end if!                                                               ·
      else ! onda plana incidente / onda cilíndrica circular ····················
        do ik = 1,pos+1                                                         !
          if (dir_j .eq. 2) then !SH                                            !
!             Meca_diff = SHdiffByStrata(B(:,ik), &                             !
!                            p_X%center%z, p_X%layer, &                         !o
!                            cOME,k_vec(ik),mecS,mecE,outpf)                    !n
!             savedauxK(ik,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE)              !d
              savedauxK(ik,1:3) = SHdiffByStrata(B(:,ik), & 
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)
          else !PSV                                                             !a
              savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                    !c
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)        !i
!             print*,ik,"savedauxk(ik,1.+5)=",sum(savedauxk(ik,1:5))
          end if                                                                !l
        end do ! ik                                                             !i
        do ik = ne,2*Nmax                                                       !n
          if (dir_j .eq. 2) then !SH                                            !d
!             Meca_diff = SHdiffByStrata(B(:,ik), &                             !r
!                            p_X%center%z, p_X%layer, &                         !i
!                            cOME,k_vec(ik),mecS,mecE,outpf)                    !c
!            savedauxK(ik,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE)               !a
             savedauxK(ik,1:3) = SHdiffByStrata(B(:,ik), & 
                                 p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)
          else !PSV                                                             !
             savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                     !
                                 p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)     !
!             print*,ik,"savedauxk(ik,1.+5)=",sum(savedauxk(ik,1:5))
          end if                                                                !
        end do ! ik                                                             !
!       print*,"sum ik1=",sum(savedauxk(1:pos+1,1)) + sum(savedauxk(ne:2*nmax,1))
!       print*,"sum ik2=",sum(savedauxk(1:pos+1,2)) + sum(savedauxk(ne:2*nmax,2))
!       print*,"sum ik3=",sum(savedauxk(1:pos+1,3)) + sum(savedauxk(ne:2*nmax,3))
!       print*,"sum ik4=",sum(savedauxk(1:pos+1,4)) + sum(savedauxk(ne:2*nmax,4))
!       print*,"sum ik5=",sum(savedauxk(1:pos+1,5)) + sum(savedauxk(ne:2*nmax,5))
!       print*,"pos=",pos,"  ne=",ne
!       stop
      end if! onda cilíndrica circular ··········································
      ! fase horizontal (fuente..receptor)
      ! para cada fuente a la profundidad iz ................... 
        if (i_zF .eq. 0) then
           nXis = 1 ! la fuente real es sólo una
        else
           nXis = pota(i_zF,2) ! fuente virtual (pueden ser varias 
        end if                 !                 a la misa profundidad)
#ifdef ver
       call ETIME(tarray, result)
       print*,"  obtenidos el mecanicos (todos los k)",result-lastresult
       lastresult = result
       print*,"   ahora la fase horizontal:"
       print*,"   pota(i_zF,:)=",pota(i_zF,:)
#endif
        do iXi = 1,nXis
          if (i_zF .ne. 0) then ! fuentes virtuales con igual Z
            itabla_x = 2 + pota(i_zF,1) + iXi
            call asociar(pXi, 0, i_zF, itabla_x)
            xf => pXi%center%x
          end if! i_zF ne 0
#ifdef ver
      print*,"     faseHorzDeFuente iXi",ixi,"de",nXis
      print*,"     fuente virtual pXi%center%x  %z", & 
                   xf,pXi%center%z," l=",pxi%length
      print*,"     cada horz de receptores"
#endif
      ! para cada X receptor a la profundidad itabla_z ......................
        n_Xs = pota(itabla_z,1) + pota(itabla_z,2)
        do itabla_x =1,n_Xs 
          call asociar(p_x, 0, itabla_z, 2+ itabla_x)
          if (p_x%isboundary .eqv. .false.) then 
              if (p_x%region .ne. 1) cycle; end if!'estr'
#ifdef ver
       print*,"        faseHorzDeReceptor itabla_x",itabla_x,"de",n_Xs
       print*,"        region: ",p_x%region
       print*,"        p_x%center%x",p_x%center%x
       print*,"        is boundary ",p_x%isBoundary
       print*,"        is movie ",p_x%guardarMovieSiblings
#endif
      ! fase horizontal
            ! reponer auxK original sin fase horizontal
!         auxK(1:2*nmax,mecS:mecE) = savedAuxK(1:2*nmax,mecS:mecE)
          auxK(1:pos+1,mecS:mecE)      = savedAuxK(1:pos+1,mecS:mecE)
          auxK(ne:2*nmax,mecS:mecE) = savedAuxK(ne:2*nmax,mecS:mecE)
            ! agregar información fase horizontal de fuente y receptor 
!           print*,"xf=",xf,"  x=",p_x%center%x," k=",k
          do imec = mecS,mecE !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
            if (isPW) then 
                if (p_x%guardarMovieSiblings .eqv. .false.) then
                auxk(1,imec) = auxk(1,imec) * &           ! onda plana      ·
                exp(-UI*k*(p_x%center%x - xf))            !                 ·
!               print*,"auxk(1,",imec,")=",auxk(1,imec)
                CYCLE ! imec                              !                 ·
                end if
            end if ! ························································
!           print*,k_vec(:)
            do ik = 1,pos+1                                                 !
                auxk(ik,imec) = auxk(ik,imec) * &                           !
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!c
!               print*,auxk(ik,imec);stop
            end do !  ik                                                    !i
            do ik = ne,2*Nmax                                               !l
                auxk(ik,imec) = auxk(ik,imec) * &                           !i
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!n
            end do !  ik                                                    !d
#ifdef ver
       call ETIME(tarray, result)                                           !i
       print*,"        agregado fase horizontal",result-lastresult          !c
       lastresult = result                                                  !a
#endif
                ! guardar el FK por fuente real ---------------------       !
                 if(plotFKS) then
                 if (i_zF .eq. 0) then 
                 if (p_x%guardarFK .eqv. .true.) then                       !
                 if (dir_j .eq. 2) then 
                    if ( imec .eq. 1) then
                      p_x%FK(J,1:nmax,3) = &                                !
                      p_x%FK(J,1:nmax,3) + auxK(1:nmax,imec) 
                    end if
                 else
                    if( imec .le. 2) then      
                      p_x%FK(J,1:pos,imec) = auxK(1:pos,imec) 
                      p_x%FK(J,pos+1:nmax,imec) = 0
                    end if
                end if;end if;end if;end if!-------------------------       !
#ifdef ver
       call ETIME(tarray, result)
       print*,"        guardado el integrando",result-lastresult            !
       lastresult = result
#endif
      ! K -> X  .........................................................   !
         if (p_x%guardarMovieSiblings .eqv. .false.) then
#ifdef ver         
        write(arg,'(a,I0,I0,a)') "k_at",p_x%pointIndex,imec,".m"
        OPEN(739,FILE=trim(arg),FORM="FORMATTED",ACTION='WRITE')
        write(arg,'(a,I0,a,I0)') "kat",p_x%pointIndex,"_",imec
        CALL scripToMatlabMNmatrixZ(2*nmax,1,auxK(:,imec),arg,739)
        close(739)
#endif 
             
             
!            auxK(:,iMec) = FFTW(2*nmax,auxK(:,iMec),+1,dk) !backward       !
!            print*,""
!            print*,j,imec
!            print*,auxK(1,iMec)
!            print*, 'p=',sum(auxK(2:pos+1,iMec))
!            print*, 'n=',sum(auxK(ne:2*nmax,iMec))
             auxK(1,iMec) = sum(auxK(1:pos+1,iMec))+sum(auxK(ne:2*nmax,iMec))
             auxK(1,iMec) = auxK(1,iMec)*dk
!            print*,pos,ne
!            print*,auxK(1,iMec);stop
             
#ifdef ver
       call ETIME(tarray, result)                                           !
       print*,"        fork non movie", result-lastresult                   !
       lastresult = result                                                  !
#endif
         end if                                                             !
      end do !imec !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
         ! print*,"auxK",auxK(1,1:5)
              !
              !                    |        RECEPTOR
              !                    | allpoint  |   boupoint
              !   ----------------------------------------------
              !   fuente real      |     W,U   | G,T: term. indep.
              !   fuente virtual   |      G    | G,T: ibem.mat.
          ! Se ha calculado el campo difractado por estratificación, 
          ! falta el campo directo (analítico) 
              if(p_x%isBoundary) then
                if(i_zF .eq. 0) then ! pXi es la fuente real
                  ! term. indep.
!                 print*,pXi%center
                  call fill_termindep(auxK(1,mecS:mecE),& 
                                      come,mecS,mecE,p_x,pXi,dir_j)
                else ! pXi es una fuente virtual
                  ! func. Green tracciones para matriz IBEM
                  call fill_ibemMat(i_zF,auxK(1,mecS:mecE),& 
                                    come,mecS,mecE,p_x,pXi,dir_j)
                end if
              else ! p_x no es un punto de colocación, es un receptor
                call fill_diffbyStrata(i_zF,J,auxK(1:2*nmax,mecS:mecE),& 
                                       come,mecS,mecE,p_x,pXi,dir_j)
              end if ! isBoundary
#ifdef ver 
      print*,"La función de Green resultado se distribuye";print*,"" 
#endif
          end do !itabla_x
        end do !iXi
      end do !itabla_z
!     end do !i_Fuente
#ifdef ver
       stop "2921 end of diffField_at_iz"
#endif
      end subroutine diffField_at_iz
      
      subroutine G0estr(MecElem,p_x,J,cOME_in,dir_j)
!#define ver 1
      ! funcion de Green en medio estratificado. 
      ! cálculo lento uno por uno.
      ! Para onda planas y para fuerza en una dirección x o z
      use gloVars, only : z0,UI,UR,PWfrecReal
      use resultVars, only : Punto,MecaElem,FFres
      use refSolMatrixVars, only : B,Ak
      use waveNumVars, only : NMAX,k_vec,dk,vecNK,SpliK,OME
      use wavelets 
      use dislin
      use sourceVars, only: Po,iFte=>currentiFte
      use soilvars, only:N,Z,alfa0,beta0,alfa,beta
      use, intrinsic :: iso_c_binding
      use debugStuff
      implicit none
      interface
        include 'interfaz.f'
         
        function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k,ik)
          use soilvars, only:N
          complex*16, dimension(1:5) :: PSVdiffByStrata
          real*8, intent(in)           :: z_i,k
          complex*16, intent(in)       :: cOME_i  
          integer, intent(in)          :: e,ik
          complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
        end function PSVdiffByStrata
        
        subroutine  eGAeNU(i_zF,ik,pXI,dj)
          use resultVars, only : Punto
          integer :: ik,i_zF,dj
          type(Punto), pointer :: pXi
        end subroutine  eGAeNU
      end interface
      type(FFres) :: FFdirecto
      complex*16, dimension(1:5) :: MecElem
      type(Punto), pointer :: p_X
      integer, intent(in) :: J,dir_j
      complex*16, intent(in),target  :: cOME_in
      
      real*8,target :: k
      real*8, pointer :: zf,xf,pt_k
      complex*16, dimension(:,:), allocatable, target :: auxK!,savedAuxK
      complex*16, target  :: cOME,alf,bet
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:),allocatable,target :: ipivA
      integer, dimension(:),pointer :: pt_ipivA
      complex*16, dimension(:),allocatable,target :: workA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16,dimension(:,:),pointer :: pointA
      type(Punto), pointer :: pXi
      logical :: isPW
      integer :: ik,tam,iMec,dj,pos,ne
            
      allocate(auxK(2*nmax,5))
      isPW = .false.
      if (Po(iFte)%tipofuente .eq. 1) isPW = .true.
            
      cOME = cOME_in 
      dj = dir_j; if(dj .eq. 3) dj = 2 
         if (isPW) then ! onda plana incidente
           if (PWfrecReal) then !  no usamos atenuación
             cOME = OME * UR  !real(cOME_in) * UR!
             alf = alfa0(N+1)
             bet = beta0(N+1)
           else
             cOME = cOME_in 
             alf = alfa(N+1)
             bet = beta(N+1)
           end if
           ! numeros de onda horizontales 
           if(Po(iFte)%PW_pol .eq. 1) k = real(cOME/bet)*sin(Po(iFte)%gamma)
           if(Po(iFte)%PW_pol .eq. 2) k = real(cOME/alf)*sin(Po(iFte)%gamma)
         end if! ································································
      call asociar(pXi,iFte,0,3) ! asociar apuntador a fuente [pXi]
#ifdef ver 
      print*,"En G0estr" 
      print*,"p_x:",p_x%center
      print*,"pXi:",pXi%center ,"ispw", isPW
#endif 
       xf=>pXi%center%x;zf=>pXi%center%z!;ef=>pXi%layer;intf=>pXi%isOnInterface
      ! Si es la fuente real y es una onda plana no se usa el DWN. Se calcula para
      ! el número de onda horizontal asociado al ángulo de incidencia ············
      if (isPW) then ! onda plana incidente   ·
            tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2!                       ·
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·
            pt_k => k; pt_come_i => cOME!                                        ·
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                · 
            call gloMat_PSV(pointA,pt_k,0)!                                      ·
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·
            call PSVvectorB_ondaplana(B(:,0),pxi%gamma)!                         ·
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·
            pos = 0; ne = 2*nmax+1
      else ! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · 
      ! La fuente es cilíndrica 
        pos = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(pos-2)                !
        B(:,pos:ne) = 0                                                         !
        tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2                           !
         ! si la fuente está sobre una interfaz ...............                 !
         if (pXi%isOnInterface) then                                            !
           do ik = 1,pos+1                                                      !
             call PSVvectorB_force(0,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                               !
           end do                                                               !
           do ik = ne,2*NMAX                                                    !
             call PSVvectorB_force(0,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik),ik)   !
             B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                               !
           end do                                                               !
         else ! la fuente está entre interfaces ..............                  !
           do ik = 1,pos+1                                                      !
             call eGAeNU(0,ik,pXI,dj)                                           !
           end do!                                                              !             
           do ik = ne,2*NMAX                                                    !
             call eGAeNU(0,ik,pXI,dj)                                           !
           end do                                                               !
         end if                                                                 !
      end if! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · ·      
      ! resultado en el receptor
      if (isPW) then ! onda plana·············                                  .
              auxK(1,1:5) = PSVdiffByStrata(B(:,0), &!                          ·
                              p_X%center%z, p_X%layer,cOME,k,0)!                ·
      else ! onda plana incidente / onda cilíndrica circular ····················
        do ik = 1,pos+1                                                         !
              auxK(ik,1:5) = PSVdiffByStrata(B(:,ik), &                         !
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)        !
        end do ! ik                                                             !
        do ik = ne,2*Nmax                                                       !
              auxK(ik,1:5) = PSVdiffByStrata(B(:,ik), &                         !
                              p_X%center%z, p_X%layer,cOME,k_vec(ik),ik)        !
        end do ! ik                                                             !
      end if! onda cilíndrica circular ··········································     
      ! agregar información fase horizontal de fuente y receptor 
       do imec = 1,5 !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
            if (isPW) then 
                auxk(1,imec) = auxk(1,imec) * &           ! onda plana      ·
                exp(-UI*k*(p_x%center%x - xf))            !                 ·
                CYCLE ! imec                              !                 ·
            end if ! ························································
            do ik = 1,pos+1                                                 !
                auxk(ik,imec) = auxk(ik,imec) * &                           !
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!c
            end do !  ik                                                    !i
            do ik = ne,2*Nmax                                               !l
                auxk(ik,imec) = auxk(ik,imec) * &                           !i
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!n
            end do !  ik                                                    !d
      ! K -> X  .........................................................   !
             auxK(1,iMec) = sum(auxK(1:pos+1,iMec))+sum(auxK(ne:2*nmax,iMec))
             auxK(1,iMec) = auxK(1,iMec)*dk             
       end do !imec !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.     
      ! campo directo si se está en el mismo estrato que la fuente
      call FFpsv(0,FFdirecto,dir_j,p_X,pXi,cOME,1,5)
#ifdef ver 
      print*,"directo_compW  =", FFdirecto%W
      print*,"directo_compU  =", FFdirecto%U
      print*,"directo_compsxx  =", FFdirecto%sxx
      print*,"directo_compszx  =", FFdirecto%szx
      print*,"directo_compszz  =", FFdirecto%szz
#endif       
      ! resultados
      MecElem(1) = FFdirecto%W + auxk(1,1)
      MecElem(2) = FFdirecto%U + auxk(1,2)
      MecElem(3) = FFdirecto%sxx + auxk(1,5)
      MecElem(4) = FFdirecto%szx + auxk(1,4)
      MecElem(5) = FFdirecto%szz + auxk(1,3)
      end subroutine G0estr
! G_stra - matrix pointAp,pt_k,pt_cOME_i
      subroutine makeGANU (J)
      use waveNumVars, only : vecNK,SpliK,nmax,cOME,k_vec, & 
         gamma=>gamma_E,nu=>nu_E,eta=>eta_E
      use soilVars, only : ALFA,BETA,N
!     use dislin
      implicit none
      integer :: J,po,ne,e,ii,ik,ikI(2),ikF(2)
      real*8  :: k
      complex*16 :: omeAlf,omeBet
      
      ! gamma y nu en esta frecuencia
       po = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(po-2)
       ikI(1) = 1  ! positivos
       ikF(1) = po+1 ! 
       ikI(2) = ne     ! negativos
       ikF(2) = 2*nmax !
       do e = 1,N+1
          omeAlf = cOME**2.0/ALFA(e)**2.0
          omeBet = cOME**2.0/BETA(e)**2.0
       
       ! de ondas cilíndricas
       do ii = 1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
          k = k_vec(ik)
          ! algunas valores constantes para todo el estrato          
          gamma(ik,e) = sqrt(omeAlf - k**2.0)
          nu(ik,e) = sqrt(omeBet - k**2.0)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z crece.
          if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
          if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
          eta(ik,e) = 2.0*gamma(ik,e)**2.0 - cOME**2.0 / BETA(e)**2.0
       end do ! ik
       end do ! ii
!      call qplot(real(k_vec,4),real(aimag(gamma(:,e)),4), 2*nmax)
!      stop
       end do ! e 
      end subroutine makeGANU
      
      subroutine makeGANU0
      use waveNumVars, only : come,ome, & 
         gamma=>gamma_E,nu=>nu_E,eta=>eta_E
      use soilVars, only : N,alfa0,beta0,alfa,beta
      use sourceVars, only : PoFte=>Po,iFte=>currentiFte
      use glovars, only: Ur,PWfrecReal
      implicit none
      integer :: e,ik
      real*8  :: k
       do e = 1,N+1
       ! de la onda plana
       ik = 0
       if (PoFte(iFte)%PW_pol .eq. 1) then
        if (PWfrecReal) then
        k = real(OME/beta0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(COME/beta(e)*sin(PoFte(iFte)%gamma))
        end if
      elseif (PoFte(iFte)%PW_pol .eq. 2) then 
        if (PWfrecReal) then
        k = real(OME/alfa0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(cOME/alfa(e)*sin(PoFte(iFte)%gamma))
        end if
      elseif (PoFte(iFte)%PW_pol .eq. 3) then 
        if (PWfrecReal) then
        k = real(OME/beta0(e)*sin(PoFte(iFte)%gamma))
        else
        k = real(COME/beta(e)*sin(PoFte(iFte)%gamma))
        end if
      end if
       if (PWfrecReal) then
       gamma(ik,e) = sqrt(UR*OME**2.0/ALFA0(e)**2.0 - k**2.0)
!      print*,gamma(ik,e),"gamma(ik,e)"
       nu(ik,e) = sqrt(UR*OME**2.0/BETA0(e)**2.0 - k**2.0)
       if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
!      print*,gamma(ik,e),"gamma(ik,e)"
       if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
       eta(ik,e) = 2.0*gamma(ik,e)**2.0 - UR*OME**2.0 / BETA0(e)**2.0
       else
       gamma(ik,e) = sqrt(cOME**2.0/ALFA(e)**2.0 - k**2.0)
       nu(ik,e) = sqrt(cOME**2.0/BETA(e)**2.0 - k**2.0)
       if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
       if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
       eta(ik,e) = 2.0*gamma(ik,e)**2.0 - cOME**2.0 / BETA(e)**2.0
       end if
       end do
       end subroutine makeGANU0
       

      subroutine gloMat_PSV(this_A,k,ik)
      ! Calcular para +k. Una vez invertida la matrix, hay paridades para -k.
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
      use debugStuff
      use waveNumVars, only : gamma_E,nu_E,eta_E,nmax
      use refSolMatrixVars, only : subMatD0,subMatS0
      implicit none
      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      integer :: ik,neik
      real*8     :: k2
      complex*16, dimension(2,4) :: subMatD,subMatS
      complex*16, dimension(4) :: diagMat
      complex*16 :: gamma,nu,xi,eta,ega,enu,ck
      integer    :: i,iR,iC,e,bord
      this_A = 0 
      ck = -k*UR
      iR= 0;iC= 0;i=1
      if (Z(0) .lt. 0.0) then !Half-Space por arriba de z=0 
        i = 0;iR= -2; iC= -2
      end if
      
      DO e = i,N+1
!         if (ik .ne. 0) then
          gamma = gamma_E(ik,e)
          nu = nu_E(ik,e)
!         else
!         gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
!         nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
!         if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
!         if(aimag(nu).gt.0.0)nu= conjg(nu)
!         end if
          xi = k**2.0 - nu**2.0 
!         eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
          eta = eta_E(ik,e)
!         print*,ik,e,k,gamma,nu,eta," L2184"
          ! en fortran los elementos se indican por columnas:
          subMatD0(1:2,1:4,e,ik) = RESHAPE((/ -gamma, ck, ck,nu,& 
                                               gamma, ck, ck,-nu /), &
                           (/ 2,4 /))
          subMatD0(1:2,1:4,e,ik) = subMatD0(1:2,1:4,e,ik) * UI
          
          k2 = 2.0*k
!         subMatS0 = RESHAPE((/ xi,-k2*gamma,-k2*nu,-xi,& 
!                              xi,k2*gamma,k2*nu,-xi /),&
!                          (/2,4/)) 
          subMatS0(1:3,1:4,e,ik) = RESHAPE(& 
                        (/ xi,      -k2*gamma,     eta,     &
                          -k2*nu,     -xi,        k2*nu,   &
                           xi,       k2*gamma,     eta,       &
                           k2*nu,     -xi,       -k2*nu /),&
                           (/3,4/)) 
          subMatS0(1:3,1:4,e,ik) = amu(e) * subMatS0(1:3,1:4,e,ik)
          ! y para k negativo aprovechando la simetria de gamma y nu
          neik = 2*nmax - (ik-2)
          if (neik .le. 2*nmax) then
          subMatD0(1:2,1:4,e,neik) = RESHAPE((/ -gamma, -ck, -ck,nu,& 
                                               gamma, -ck, -ck,-nu /), &
                           (/ 2,4 /))
          subMatD0(1:2,1:4,e,neik) = subMatD0(1:2,1:4,e,neik) * UI
          ! y para el k negativo aprovechando simetria de gamma,nu,xi,eta,
          subMatS0(1:3,1:4,e,neik) = RESHAPE(& 
                        (/ xi,      k2*gamma,     eta,     &
                           k2*nu,     -xi,       -k2*nu,   &
                           xi,      -k2*gamma,     eta,       &
                          -k2*nu,     -xi,        k2*nu /),&
                           (/3,4/)) 
          subMatS0(1:3,1:4,e,neik) = amu(e) * subMatS0(1:3,1:4,e,neik)
          end if
          
          
          
          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1 
        do bord = 0,1
          if ((e .eq. 0) .and. (bord .eq. 0)) cycle
          if (e+bord > N+1) exit 
          ! si 1+0;1+1;2+0;[2+1] > 2
        if (bord .eq. 0) then !---->---->---->---->---->---->
          if (e /= N+1) then !(radiation condition lower HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e)))
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))
          else
            ega = Z0
            enu = Z0
          end if
!         diagMat = RESHAPE((/ UR, Z0, Z0, Z0, & 
!                              Z0, UR, Z0, Z0, & 
!                              Z0, Z0, ega, Z0, & 
!                              Z0, Z0, Z0, enu /), &
!                          (/ 4,4 /))
          diagMat = (/ UR, UR, ega, enu /)
        else !bord .eq. 1 -----------------------------------
          if (e /= 0) then !(radiation condition upper HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e))) 
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))       
          else
            ega = Z0
            enu = Z0
          end if !<----<----<----<----<----<----<----<----<--
!         diagMat = RESHAPE((/ ega, Z0, Z0, Z0, & 
!                              Z0, enu, Z0, Z0, & 
!                              Z0, Z0, UR, Z0, & 
!                              Z0, Z0, Z0, UR /), &
!                          (/ 4,4 /))
          diagMat = (/ ega, enu, UR, UR /)
        end if
          ! desplazamientos
!         subMatD = matmul(subMatD0(1:2,1:4,e,ik),diagMat)
          subMatD(1:2,1) = subMatD0(1:2,1,e,ik)*diagMat(1)
          subMatD(1:2,2) = subMatD0(1:2,2,e,ik)*diagMat(2)
          subMatD(1:2,3) = subMatD0(1:2,3,e,ik)*diagMat(3)
          subMatD(1:2,4) = subMatD0(1:2,4,e,ik)*diagMat(4)
          ! esfuerzos
!         subMatS = matmul(subMatS0(1:2,1:4,e,ik),diagMat)
          subMatS(1:2,1) = subMatS0(1:2,1,e,ik)*diagMat(1)
          subMatS(1:2,2) = subMatS0(1:2,2,e,ik)*diagMat(2)
          subMatS(1:2,3) = subMatS0(1:2,3,e,ik)*diagMat(3)
          subMatS(1:2,4) = subMatS0(1:2,4,e,ik)*diagMat(4)
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0 .AND. e /= 0 ) then 
           if (e /= N+1) then !(radiation condition)
            if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+4 ) = -subMatD
            end if
             this_A( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatS
           else
           if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+2 ) = -subMatD(:,1:2)
           end if
             this_A( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatS(:,1:2)
            ! exit
           end if
          end if
          
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
           if (e /= 0) then !(radiation condition upward)
            ! ambas ondas
            this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatD
            this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatS
           else
            ! solo onda hacia arriba
            this_A( iR+3 : iR+4 , iC+3 : iC+4 ) = subMatD(:,3:4)
            this_A( iR+5 : iR+6 , iC+3 : iC+4 ) = subMatS(:,3:4)
           end if
          end if
        end do !bord loop del borde i superior o nferior
          iR= iR+4 
          iC= iC+4
      END DO !{e} loop de las macro columnas para cada estrato
!            call showMNmatrixZ(size(this_A,1),& 
!            size(this_A,2),this_A,"A    ",6) 
!            stop "gloMat_PSV"
      end subroutine gloMat_PSV
      
      
      subroutine globalmatrix_SH(this_A,k,ik)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
      use debugStuff  
      use waveNumVars, only : nu_E
      implicit none
      
      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      integer :: ik
!     complex*16, intent(in),pointer     :: cOME_i  
      
      real*8     :: z_i
      complex*16, dimension(1,2) :: subMatD, subMatS, subMatD0, subMatS0
      complex*16, dimension(2,2) :: diagMat
      complex*16 :: nu,enuN,enuP
      integer    :: i,iR,iC,e,bord
      this_A = 0
      z_i = k ! (nada mas para que no chiste)
      iR= 0;iC= 0;i=1
      if (Z(0) .lt. 0.0) then !Half-Space por arriba de z=0 
        i = 0;iR= -1; iC= -1
      end if
      
      DO e = i,N+1!cada estrato
          nu = nu_E(ik,e)!; print*,e,N,nu
          ! algunas valores constantes para todo el estrato
!         nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z es mayor que cero.
!         if(aimag(nu).gt.0.0)nu= conjg(nu)

          ! en fortran los elementos se indican por columnas:
          subMatD0 = RESHAPE((/ UR,UR /), (/ 1,2 /))
          subMatS0 = RESHAPE((/ -UI*nu,UI*nu /),(/1,2/)) 
          subMatS0 = amu(e) * subMatS0                  
          
          ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1 
        do bord = 0,1
          if ((e .eq. 0) .and. (bord .eq. 0)) cycle
          if (e+bord > N+1) then ! si 1+0;1+1;2+0;[2+1] > 2
            exit
          end if                           
          ! la profundidad z de la frontera superior del estrato
                z_i = Z(e+bord)   ! e=1 , bord=0  ->  z = z0 = 0
                                  ! e=1 , bord=1  ->  z = Z1 = h1
          !downward waves
          if (e /= 0) then !(radiation condition upper HS)
            enuN = exp(-UI * nu * (z_i-Z(e)))
          else
            enuN = Z0
          end if
          !upward waves    
          if (e /= N+1) then !(radiation condition lower HS)
            enuP = exp(UI * nu * (z_i-Z(e+1)))
          else
            enuP = Z0
          end if
          diagMat = RESHAPE((/ enuN, Z0,   & 
                               Z0,   enuP /), &
                              (/ 2,2 /))
          subMatD = matmul(subMatD0,diagMat)
          subMatS = matmul(subMatS0,diagMat)
!         call showMNmatrixZ(2,2,diagmat,"diag ",6)
!         call showMNmatrixZ(1,2,submatd,"  D  ",6)
!         call showMNmatrixZ(1,2,submats,"  S  ",6)
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0 .AND. e /= 0 ) then 
           if (e /= N+1) then !(radiation condition downward)
            if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR   : iR   , iC+1 : iC+2 ) = -subMatD!; print*,"b1"
            end if !(e==2,3,...)
             this_A( iR+1 : iR+1 , iC+1 : iC+2 ) = -subMatS!; print*,"b2"
           else
           if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR   : iR   , iC+1 : iC+1 ) = -subMatD(:,1:1)!; print*,"b3"
           end if
             this_A( iR+1 : iR+1 , iC+1 : iC+1 ) = -subMatS(:,1:1)!; print*,"b4"
!            exit
           end if
          end if
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
           if (e /= 0) then !(radiation condition upward)
            ! ambas ondas
            this_A( iR+2 : iR+2 , iC+1 : iC+2 ) = subMatD!; print*,"b5",iR+2, iC+1
            this_A( iR+3 : iR+3 , iC+1 : iC+2 ) = subMatS
           else
            ! solo onda hacia arriba
            this_A( iR+2 : iR+2 , iC+2 : iC+2 ) = subMatD(:,2:2)!; print*,"b5",iR+2, iC+1
            this_A( iR+3 : iR+3 , iC+2 : iC+2 ) = subMatS(:,2:2)
           end if
!           this_A( iR+2 : iR+2 , iC+1 : iC+2 ) = subMatD!; print*,"b5",iR+2, iC+1
!           this_A( iR+3 : iR+3 , iC+1 : iC+2 ) = subMatS
          end if
!         print*,"bord=",bord,"e=",e,"IR=",IR,"IC=",IC
!         call showMNmatrixZ(2*N+2,2*N+2,this_A,"Ash  ",6) 
        end do !bord loop del borde i superior o nferior
          iR= iR+2 
          iC= iC+2
      END DO !{e} loop de las macro columnas para cada estrato
!         call showMNmatrixZ(size(this_A,1),size(this_A,2),this_A,"Ash  ",6) 
!         stop 5248
!     print*,"glomat",ik,"line2789"
      end subroutine globalmatrix_SH
      
      subroutine inverseA(A,ipiv,work,n)
!     use debugstuff
      integer, intent(in) :: n
      complex*16, dimension(:,:), intent(inout),pointer :: A
      integer, dimension(:), intent(inout),pointer :: ipiv
      complex*16, dimension(:),intent(inout),pointer :: work
      integer :: info
      integer :: lwork
      lwork = n*n
      call zgetrf(n,n,A,n,ipiv,info)
      if(info .ne. 0) stop "inverseA :Problem at LU factorization of matrix Try more Q"
      ! Para más de 3 estratos !--------
      ! 1) convertir a almacenamiento bandeado
!       2) call ZGBTRF( n, n, 5, 5, A, n, IPIV, INFO )
!     if(info .ne. 0) stop "Problem at LU factorization of matrix "
      ! luego con cada B
!       3) call ZGBTRS( 'N', N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
      call zgetri(n,a,n,ipiv,work,lwork,info)
      if(info .ne. 0) stop "inverseA :Problem at inverse of matrix "
      end subroutine inverseA
      
      subroutine intrplr_gloMat(k0,n,pt_cOME_i,pt_ipivA,pt_workA)
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use waveNumVars, only : k_vec, NMAX
      use fitting
      use soilvars, only:alfa,beta,Nestr => N
      use waveNumVars, only : gamma=>gamma_E,nu=>nu_E
      use sourceVars, only  : PSV
      implicit none
      interface
      include 'interfaz.f'
      end interface
      integer :: i,k0,k1,k2,n,tam,r,c,ik,e
!     integer :: d != 4
      real*8 :: ratio,step
      real*8, dimension(n)     :: xdat
      integer,dimension(n)     :: idat
      complex*16, dimension(n) :: ydat
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:), pointer :: pt_ipivA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16, dimension(:,:), pointer :: pointAp
      real*8, pointer :: pt_k
      complex*16 :: m
      ! armar los vectores de datos
      ratio = 1.0 ! mayor o igual a 1
      tam = size(Ak,1)
!     step = ((nmax+1) - (k0))/(n-1) !iguales
      if (int(n/2) .lt. real(n)/2.) then
      step = (int((n-1)/2)) * n  !impar
      else
      step = (int((n-1)/2) + 1) * n !par
      end if
      step = ((nmax+1) - (k0))/(ratio*step) !crecim geometrico
      idat(1) = k0
      do i=2,n-1
         idat(i) =  idat(i-1) + int((i-1)*step)
         xdat(i) = k_vec(idat(i))
      end do
      idat(n) = nmax+1
      xdat(n) = k_vec(idat(n))
!     do i=1,n;print*,i,idat(i);end do;stop
      do i=2,n
         ik = idat(i)
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)
            do e=1,Nestr+1
              gamma(ik,e) = sqrt(pt_cOME_i**2.0/ALFA(e)**2.0 - pt_k**2.0)
              nu(ik,e) = sqrt(pt_cOME_i**2.0/BETA(e)**2.0 - pt_k**2.0)
              if(aimag(gamma(ik,e)).gt.0.0_8)gamma(ik,e) = conjg(gamma(ik,e))
              if(aimag(nu(ik,e)).gt.0.0_8)nu(ik,e)=conjg(nu(ik,e))
            end do ! e
         if (PSV) then
         call gloMat_PSV(pointAp,pt_k,ik)
         else
         call globalmatrix_SH(pointAp,pt_k,ik)
         end if
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
      end do
      
      ! interpolar cada elemento (rectas)
      do r=1,tam
      do c=1,tam
        ydat(1) = Ak(r,c,idat(1))
!       print*,ydat(1)
        k1 = idat(1)+1
        do i=2,n
           k2 = idat(i) !k0 + int((i-1)*step)
           ydat(i) = Ak(r,c,k2)
           m = (ydat(i)-ydat(i-1))/(xdat(i)-xdat(i-1))
           do ik=k1,k2
             Ak(r,c,ik) = ydat(i-1) + m * (k_vec(ik) - xdat(i-1))
           end do
           k1 = k2+1
        end do
!       print*,""
!       do i=k0-2,nmax+1
!         print*,i,k_vec(i),Ak(r,c,i)
!       end do       
!       stop
      end do !c
      end do !r
      end subroutine intrplr_gloMat
      
      
      subroutine intrplr_zpoly_gloMat(k0,n,pt_ipivA,pt_workA)
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use waveNumVars, only : k_vec, NMAX
      use fitting
      implicit none
      interface
      include 'interfaz.f'
      end interface
      integer :: i,k0,k1,k2,n,tam,r,c
      integer :: d != 4
      real*8 :: step
      real*8, dimension(n)     :: xdat
      complex*16, dimension(n) :: ydat
!     complex*16, pointer :: pt_cOME_i
      integer, dimension(:), pointer :: pt_ipivA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16, dimension(:,:), pointer :: pointAp
      real*8, pointer :: pt_k
      complex*16, dimension(10+1) :: po
      d = 10
      ! armar los vectores de datos
      tam = size(Ak,1)
      step = ((nmax+1) - (k0))/(n-1)
      xdat(1) = k_vec(k0)
      xdat(n) = k_vec(nmax+1)
      do i=2,n-1
         k1 = k0 + int((i-1)*step)
         xdat(i) = k_vec(k1)
      end do!
      do i=2,n
         k1 = k0 + int((i-1)*step)
         pointAp => Ak(1:tam,1:tam,k1)
         pt_k => k_vec(k1)
         call gloMat_PSV(pointAp,pt_k,k1)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
      end do
      ! interpolar cada elemento
      do r=1,tam
      do c=1,tam
        ydat(1) = Ak(r,c,k0)
        do i=2,n
           k1 = k0 + int((i-1)*step)
           ydat(i) = Ak(r,c,k1)
         if (abs(ydat(i)) .lt. 1E-5) then
           k2 = k1
           ydat(i:n) = 0
           exit
         else
           k2 = nmax+1
         end if
        end do
        
        po = Zpolyfit(xdat(1:n),ydat(1:n),n,d)
        do i=k0+1,k2
           Ak(r,c,i) = Zevapol(po,d,k_vec(i))
        end do
        do i=k2,nmax+1
           Ak(r,c,i) = 0
        end do
      end do !c
      end do !r
      end subroutine intrplr_zpoly_gloMat
      
      
      ! de los elementos a_ij
      ! pares en k:  i={1,3,5...} con j={1,3,4...} ; i={2,4,6...} con j={2,4,6...}
      ! impares  k:  i={1,3,5...} con j={2,4,6...} ; i={2,4,6...} con j={1,3,5...}
      subroutine parImpar_gloMat
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use waveNumVars, only : NMAX
      use debugStuff
      implicit none
      integer :: r,c,ik,ikp,tam
      integer,dimension(:,:),allocatable :: Comal
      logical :: s !empieza en priemr columna
      tam = size(Ak,1)
      allocate(Comal(tam,tam))
      Comal =1
      ! los odd
      s = .false.
      do r = 1,tam
      if (s) then
      Comal(r,1:tam:2) = -1
      else
      Comal(r,2:tam:2) = -1
      end if
      s = (.not. s)
      end do
!     call showMNmatrixI(tam,tam, Comal," Que ",6)
      ! negativos
      do ik=nmax+2,2*NMAX
        ikp = (2*NMAX +2-ik)
!     print*,ik,ikp!;cycle
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ikp),"lado positivo",6)
        do r=1,tam
        do c=1,tam
        Ak(r,c,ik) = Ak(r,c,ikp) * Comal(r,c)
        end do
        end do
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),"lado negativo",6)
!     stop
      end do
!     stop "parImpar_gloMat"
      end subroutine parImpar_gloMat
      
! G_stra - termIndPSV
      subroutine PSVvectorB_ondaplana(this_B,gamma)
      use soilvars, only : n,lambda0,amu0,lambda,amu,alfa0,beta0,Z,beta,alfa
      use glovars, only:UI,z0,PWfrecReal
      use sourceVars, only: Po,iFte=>currentiFte
      use waveNumVars, only : cOME,ome
      use debugStuff
      implicit none
      complex*16, intent(inout), dimension(1:4*N+2) :: this_B
!     complex*16, intent(in)    :: come ! no trae amortiguamiento
      real*8, intent(in) :: gamma
      integer :: i,e
      real*8,dimension(1:2) :: theta
      complex*16 :: kx,kz,U,W,c,la,am
      real*8 :: z_loc
      !     Colocamos la onda incindente en la interfaz
      !     con el semiespacio de abajo.
      e = Po(iFte)%layer! N+1 
      z_loc = Z(N+1)
      
      if (Po(iFte)%PW_pol .eq. 1) then
        if (PWfrecReal) then
        c = beta0(e) !SV
        else
        c = beta(e) !SV
        end if
        theta(1) = cos(gamma)
        theta(2) = sin(gamma)
      elseif (Po(iFte)%PW_pol .eq. 2) then 
        if (PWfrecReal) then
        c = alfa0(e) !SV
        else
        c = alfa(e) !SV
        end if
        theta(1) = sin(gamma)
        theta(2) = -cos(gamma)
      end if
      !
      if (PWfrecReal) then
      kx = ome/c * sin(gamma)
      kz = ome/c * cos(gamma)
      la = LAMBDA0(e)
      am = AMU0(e)
      else
      kx = come/c * sin(gamma)
      kz = come/c * cos(gamma)
      la = LAMBDA(e)
      am = AMU(e)
      end if
      U = (theta(1))* exp(UI * kz * (z_loc))
      W = (theta(2))* exp(UI * kz * (z_loc))
      
      i=0
      this_B(1:4*N+2) = Z0 
      if (Z(0) .lt. 0.0) then ! Semiespacio en z<0 ···········
        i = 2                                                !
        this_B(1+4*(e-1)-2 + i) = W !  w                     !
        this_B(1+4*(e-1)-1 + i) = U !  u                     !
      end if                                                 !  
      ! ······················································
      if (e .ne. 1) then                                     ! 
        this_B(1+4*(e-1)-2 + i) = W !  w                     !
        this_B(1+4*(e-1)-1 + i) = U !  u                     !
      end if                                                 !
      !.......................................................
      
      ! Tracciones en la frontera de la región de la fuente.........
      this_B(1+4*(e-1)   + i) = UI * ( &                           !
                            ( W * kz * (la + 2.0 * am)) & !
                          - ( U * kx * la)) ! szz           !
      this_B(1+4*(e-1)+1 + i) = UI * am &                      !
                          * ( kz * U - kx * W ) ! szx              !
      !                   sxx = UI * ( &                           !
      !                   - ( U * kx * (LAMBDA(e) + 2.0*AMU(e))) & !
      !                   + ( W * kz * LAMBDA(e)))                 !
      !............................................................!
!     call showMNmatrixZ(4*N+2,1, this_B ,"  B  ",6)
      end subroutine PSVvectorB_ondaplana
      
      
      subroutine PSVvectorB_force(i_zF,this_B,tam,pXi,direction,cOME,k,ik)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO,NPAR
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only: Po,iFte=>currentiFte!use sourceVars, only: tipofuente
      use waveNumVars, only : gamma_E,nu_E
      implicit none
      
      integer, intent(in) :: i_zF,tam,ik
      complex*16, intent(inout), dimension(tam) :: this_B
      integer,    intent(in)    :: direction
      real*8,     intent(in)    :: k
      complex*16, intent(in)    :: cOME
      type(Punto),intent(in),target    :: pXi
      
      integer, pointer :: e_f
      real*8, pointer :: z_f
      logical, pointer :: fisInterf
      
      integer :: iIf,nInterf
      real    :: SGNz
      complex*16 :: gamma,nu,DEN,L2M, argum,sincmod
      real*8     :: errT = 0.0001_8
!     complex*16 :: omeAlf,omeBet
      
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: egamz,enuz,gamz,nuz
      complex*16, dimension(2) :: sincGamma, sincNu
      
                                  !  1 para Greeni1 (horizontal),
                                  !  3 para Greeni3 (vertical)
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s331,s131,s333,s313
      
      real*8 :: a
      real*8, pointer :: cose,seno
      integer :: el_tipo_de_fuente
      
      this_B = Z0
      sincGamma = Z0
      sincNu = Z0
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
         go to 349
      end if   
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1)  z_loc(2) = Z(e_f+1) - z_f !upward (+)
      
      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if  ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "PSVvectorB_force iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente 
      
      DEN = 4.0*PI*RHO(e_f)*cOME**2.0
      L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
      gamma = gamma_E(ik,e_f)
      nu = nu_E(ik,e_f)
          
      do iIf = 1,2
          egamz(iIf) = exp(-UI*gamma*ABS(z_loc(iIf)))
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))
          
          gamz(iIf) = gamma*ABS(z_loc(iIf))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do 
      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincGamma(1) = UR * egamz(1)
        sincGamma(2) = UR * egamz(2)
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
      ! ahora estamos involucrando información del receptor en el vector de fuente
      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
      
        A = pxi%length * 0.5_8  ! 2a=lenght
        cose => pxi%cosT
        seno => pxi%sinT
      
      ! en cada interfaz 
        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
        argum = gamma*A*seno + UR * (k*A*cose)
        sincGamma(2) = sincmod(argum,gamz(2))*(2*A)
        
        argum = nu*A*seno + UR * (k*A*cose)
        sincNu(2) = sincmod(argum,nuz(2))*(2*A)
        
        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
        argum = -gamma*A*seno + UR * (k*A*cose)
        sincGamma(1) = sincmod(argum,gamz(1))*(2*A)
        
        argum = -nu*A*seno + UR * (k*A*cose)
        sincNu(1) = sincmod(argum,nuz(1))*(2*A)
        
      end if !fuente tipo 0 o 1
      
      do iIf = 1,2    
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if

      G31(iIf) = -UI/DEN * SGNz*k*(sincGamma(iIf) & 
                           - sincNu(iIf))
                           
      if (direction .eq. 1) then !  G31, G11, S331, S131
          
      G11(iIf) = -UI/DEN * (k**2.0/gamma*sincGamma(iIf) &
                           +nu*sincNu(iIf)) 
      
      s331(iIf) = -UR/DEN * ( &
                    (k*gamma*L2M + lambda(e_f)*k**3.0/gamma)* sincGamma(iIf)&
                  + (-2.0*amu(e_f)*k*nu)* sincNu(iIf) &
                    ) !
                    
      s131(iIf) = -UR/DEN * amu(e_f)*SGNz * ( &             
                    (2.0*k**2.0)* sincGamma(iIf) &
                    + (nu**2.0-k**2.0)* sincNu(iIf) &
                    ) !
      end if 
      !
      if (direction .eq. 3) then ! G33,G31,S333,S313
      
      G33(iIf) = -UI/DEN * (gamma*sincGamma(iIf) & 
                           +k**2.0/nu*sincNu(iIf))
                    
      s333(iIf) = -UR/DEN * SGNz * ( &              
                    (gamma**2.0*L2M + k**2.0*lambda(e_f))* sincGamma(iIf) &
                  + (2.0*amu(e_f)*k**2.0)* sincNu(iIf) &
                    ) !              
                    
      s313(iIf) = -UR/DEN * amu(e_f) * ( &              
                    (2.0*k*gamma)* sincGamma(iIf) &
                  - (k/nu*(nu**2.0-k**2.0))* sincNu(iIf) &
                    ) !
      end if
      end do !iIf interface
      
 349  if (fisInterf) then
      if(direction .eq. 1) then
        s331(1) = 0
        S131(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
      elseif (direction .eq. 3) then
        s333(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
        S313(1) = 0
      end if
        G33(1) = 0
        G31(1) = 0
        G11(1) = 0
      end if
      
      ! El vector de términos independientes genera el campo difractado
      if (z(0) .gt. 0.0) then !--------------------------
      if (direction .eq. 1) then ! fuerza HORIZONTAL
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+4*(e_f-1)-2) = G31(1)!  w
        this_B(1+4*(e_f-1)-1) = G11(1)!  u
       end if 
        this_B(1+4*(e_f-1)  ) = S331(1)! szz
        this_B(1+4*(e_f-1)+1) = S131(1)! szx   ! delta
      
      if (.not. fisInterf) then ! la fuerza no en la interfaz
      !                     =      (2) interfaz de abajo
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+2) = - G31(2)!  w
        this_B(1+4*(e_f-1)+3) = - G11(2)!  u
        this_B(1+4*(e_f-1)+4) = - S331(2)! szz
        this_B(1+4*(e_f-1)+5) = - S131(2)! szx
       end if
      end if
      
      elseif (direction .eq. 3) then ! fuerza VERTICAL
      !                     =     (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+4*(e_f-1)-2) = G33(1)!  w 
        this_B(1+4*(e_f-1)-1) = G31(1)!  u 
       end if 
        this_B(1+4*(e_f-1)  ) = S333(1)! szz   ! delta
        this_B(1+4*(e_f-1)+1) = S313(1)! szx 
    
      if (.not. fisInterf) then
      !                     =    (2) interfaz de abajo (con Fza en el HS no entra aqui)
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+2) = - G33(2)!  w 
        this_B(1+4*(e_f-1)+3) = - G31(2)!  u
        this_B(1+4*(e_f-1)+4) = - S333(2)! szz 
        this_B(1+4*(e_f-1)+5) = - S313(2)! szx 
       end if
      end if
      end if ! direction
      else !Z(0) < 0 un semiespacio arriba --------------------------
      if (direction .eq. 1) then ! fuerza HORIZONTAL
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+4*(e_f-1)  ) = G31(1)!  w
        this_B(1+4*(e_f-1)+1) = G11(1)!  u
        this_B(1+4*(e_f-1)+2) = S331(1)! szz
        this_B(1+4*(e_f-1)+3) = S131(1)! szx   ! delta
       end if 
      !
      if (.not. fisInterf) then ! la fuerza no en la interfaz
      !                     =      (2) interfaz de abajo
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+4) = - G31(2)!  w
        this_B(1+4*(e_f-1)+5) = - G11(2)!  u
        this_B(1+4*(e_f-1)+6) = - S331(2)! szz
        this_B(1+4*(e_f-1)+7) = - S131(2)! szx
       end if
      end if
      
      elseif (direction .eq. 3) then ! fuerza VERTICAL
      !                     =     (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+4*(e_f-1)  ) = G33(1)!  w 
        this_B(1+4*(e_f-1)+1) = G31(1)!  u 
        this_B(1+4*(e_f-1)+2) = S333(1)! szz   ! delta
        this_B(1+4*(e_f-1)+3) = S313(1)! szx 
       end if 
      !
      if (.not. fisInterf) then
      !                     =    (2) interfaz de abajo (con Fza en el HS no entra aqui)
       if (e_f .ne. N+1) then
        this_B(1+4*(e_f-1)+4) = - G33(2)!  w 
        this_B(1+4*(e_f-1)+5) = - G31(2)!  u
        this_B(1+4*(e_f-1)+6) = - S333(2)! szz 
        this_B(1+4*(e_f-1)+7) = - S313(2)! szx 
       end if
      end if
      end if ! direction
      end if 
      ! cuand la fuente es una onda plana, sólo un número de onda
      ! es distinto de cero, k = omega/c y no se hace fft en k
!     print*,"e_f=",e_f
!     print*,this_B
      end subroutine PSVvectorB_force
      
      
      subroutine PSVpaGaNU(J) 
      use waveNumVars, only : cOME,vecNK,k_vec,nmax,SpliK,gamma_E,nu_E
      use soilVars
      use gloVars, only : UR,UI,PI
      use refSolMatrixVars, only : BparaGa,BparaNu
      implicit none
      integer, intent(in) :: J
      complex*16, dimension(:,:,:,:),pointer :: this_B
      complex*16 :: gamma,nu,DEN,L2M!,omeAlf,omeBet
      integer :: ii,ik,e,Ga_o_Nu,dir,pos,neg,ikI(2),ikF(2)
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s331,s131,s333,s313
      real*8 :: k
       pos = min(int(vecNK(J)*SpliK),nmax); neg = 2*nmax-(pos-2)
       ikI(1) = 1
       ikF(1) = pos+1
       ikI(2) = neg
       ikF(2) = 2*nmax
       do ii=1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
       k = k_vec(ik)
       do e=1,N+1 ! estrato que contiene la fuerza
         
         DEN = 4.0*PI*RHO(e)*cOME**2.0
!        omeAlf = cOME**2.0/ALFA(e)**2.0
!        omeBet = cOME**2.0/BETA(e)**2.0
         L2M = LAMBDA(e) + 2.0*AMU(e)
      
!         ! algunas valores constantes para todo el estrato          
!         gamma = sqrt(omeAlf - k**2.0)
!         nu = sqrt(omeBet - k**2.0)
!         ! Se debe cumplir que la parte imaginaria del número de onda 
!         ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
!         ! nencial a medida que z crece.
!         if(aimag(gamma).gt.0.0_8)gamma = conjg(gamma)
!         if(aimag(nu).gt.0.0_8)nu=conjg(nu)
          gamma = gamma_E(ik,e)
          nu = nu_E(ik,e)
         
         ! the green function indexes go for the part
         ! that is multiplied by exp(gamma) : 1
         ! and                by exp(nu) : 2
         G31(1) = -UI/DEN * k 
         G31(2) =  UI/DEN * k 
         
         G11(1) = -UI/DEN * k**2.0/gamma
         G11(2) = -UI/DEN * nu
         
         S331(1) = -UR/DEN * (k*gamma*L2M + lambda(e)*k**3.0/gamma)
         S331(2) = -UR/DEN * (-2.0*amu(e)*k*nu)
         
         S131(1) = -UR/DEN * amu(e)* (2.0*k**2.0)
         S131(2) = -UR/DEN * amu(e)* (nu**2.0-k**2.0)
         
         G33(1) = -UI/DEN * gamma
         G33(2) = -UI/DEN * (k**2.0/nu)
         
         s333(1) = -UR/DEN * (gamma**2.0*L2M + k**2.0*lambda(e))
         s333(2) = -UR/DEN * (2.0*amu(e)*k**2.0)
         
         S313(1) = -UR/DEN * amu(e) * (2.0*k*gamma)
         S313(2) =  UR/DEN * amu(e) * (k/nu*(nu**2.0-k**2.0))
         
       do Ga_o_Nu = 1,2
         if (Ga_o_Nu .eq. 1) this_B=>BparaGa
         if (Ga_o_Nu .eq. 2) this_B=>BparaNu
         ! para dir 1
         dir = 1
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2,e,ik,dir) = G31(Ga_o_Nu) * (-1)!  w
        this_B(1+4*(e-1)-1,e,ik,dir) = G11(Ga_o_Nu)!  u
       end if 
        this_B(1+4*(e-1)  ,e,ik,dir) = S331(Ga_o_Nu)! szz
        this_B(1+4*(e-1)+1,e,ik,dir) = S131(Ga_o_Nu) * (-1)! szx   ! delta
      !                     =      (2) interfaz de abajo
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2,e,ik,dir) = - G31(Ga_o_Nu)* (1)!  w
        this_B(1+4*(e-1)+3,e,ik,dir) = - G11(Ga_o_Nu)!  u
        this_B(1+4*(e-1)+4,e,ik,dir) = - S331(Ga_o_Nu)! szz
        this_B(1+4*(e-1)+5,e,ik,dir) = - S131(Ga_o_Nu) * (1)! szx
       end if
         
         ! para dir 3
         dir = 2
       if (e .ne. 1) then
        this_B(1+4*(e-1)-2,e,ik,dir) = G33(Ga_o_Nu)!  w 
        this_B(1+4*(e-1)-1,e,ik,dir) = G31(Ga_o_Nu)* (-1)!  u 
       end if 
        this_B(1+4*(e-1)  ,e,ik,dir) = S333(Ga_o_Nu)* (-1)! szz   ! delta
        this_B(1+4*(e-1)+1,e,ik,dir) = S313(Ga_o_Nu)! szx 
      !                     =    (2) interfaz de abajo 
       if (e .ne. N+1) then
        this_B(1+4*(e-1)+2,e,ik,dir) = - G33(Ga_o_Nu)!  w 
        this_B(1+4*(e-1)+3,e,ik,dir) = - G31(Ga_o_Nu)* (1)!  u
        this_B(1+4*(e-1)+4,e,ik,dir) = - S333(Ga_o_Nu)* (1)! szz 
        this_B(1+4*(e-1)+5,e,ik,dir) = - S313(Ga_o_Nu)! szx 
       end if
       
       end do !Ga_o_Nu
       end do ! e
       end do ! ik
       end do !ii
      end subroutine PSVpaGANU
      
      subroutine PSVMatAporGaNU(J) ! multiplicar marices A y (casi) B
      use waveNumVars, only : vecNK, nmax,SpliK
      use soilVars, only:N
      use refSolMatrixVars, only : BparaGa,BparaNu,CoefparGa,CoefparNu,Ak
      implicit none
      integer, intent(in) :: J
      complex*16, dimension(:,:,:,:),pointer :: this_B
      complex*16, dimension(:,:,:,:,:),pointer ::this_coef
      integer :: ii,ik,e,Ga_o_Nu,dir,pos,neg,ikI(2),ikF(2),tam,& 
                coIu,coFu,coId,coFd
       tam = 4*N+2
       pos = min(int(vecNK(J)*SpliK),nmax); neg = 2*nmax-(pos-2)
       ikI(1) = 1
       ikF(1) = pos+1
       ikI(2) = neg
       ikF(2) = 2*nmax
       do ii=1,2 ! los num de onda positivos, negativos
       do ik = ikI(ii),ikF(ii) !  todos los num de onda
       do Ga_o_Nu = 1,2
         if (Ga_o_Nu .eq. 1) then 
           this_B=>BparaGa ;  this_coef=>CoefparGa ; end if!
         if (Ga_o_Nu .eq. 2) then 
           this_B=>BparaNu ;  this_coef=>CoefparNu ; end if
         do dir = 1,2 ! para cada dirección de la fuerza
           do e=1,N+1 ! estrato que contiene la fuerza
               coIu = 1+4*(e-1)   ! sz
             if (e .ne. 1) then
               coIu = 1+4*(e-1)-2 ! w
               coFu = 1+4*(e-1)-1 ! u
             end if ! e!= 1
               coFu = 1+4*(e-1)+1 ! sx
             if (e .ne. N+1) then
               coId = 1+4*(e-1)+2 ! w
!              1+4*(e-1)+3        ! u
!              1+4*(e-1)+4        ! sz
               coFd = 1+4*(e-1)+5 ! sx
             end if ! e!= HS
!            print*,ii,ik,Ga_o_nu,dir,e;print*,"   ",coIu,coFu,coId,coFd
               this_coef(1:tam,e,ik,dir,1) = matmul(Ak(1:tam,coIu:coFu,ik),&
                  this_B(coIu:coFu,e,ik,dir))
             if (e .ne. N+1) then
               this_coef(1:tam,e,ik,dir,2) = matmul(Ak(1:tam,coId:coFd,ik),&
                  this_B(coId:coFd,e,ik,dir))
             end if
           end do ! e
         end do ! dir
       end do ! Ga_o_Nu
       end do ! ik
       end do ! ii
!      print*,CoefparGa(:,1,10,1,1);stop "PSVMatAporGaNU"
      end subroutine PSVMatAporGaNU
      
      ! ****  la fuente cilíndrica está entre interfaces de los estratos ****
      subroutine  eGAeNU(i_zF,ik,pXI,dj) 
      !hacer sincGamma y sincNu para cada interfaz y devolver B final
        use waveNumVars, only : k_vec,gamma_E,nu_E!cOME,
        use resultVars, only : Punto
        use soilVars
        use sourceVars, only: Po,iFte=>currentiFte
        use glovars, only : UR,UI
        use refSolMatrixVars, only : B, CoefparGa, CoefparNu
        implicit none
        integer :: ik,i_zF,dj
        type(Punto), pointer :: pXi
        integer, pointer :: e_f
        real*8, pointer :: z_f
        logical, pointer :: fisInterf
        integer :: iIf
        complex*16 :: gamma,nu,argum,sincmod
        
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
        real*8,     dimension(2) :: z_loc
        complex*16, dimension(2) :: egamz,enuz,gamz,nuz
        complex*16, dimension(2) :: sincGamma, sincNu
        real*8 :: a,k
        real*8, pointer :: cose,seno
        integer :: el_tipo_de_fuente,tam
      
      tam = 4*N+2
      k = k_vec(ik)
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      if (fisInterf) stop "eGAeNU: force on interface"
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)    (interfaz de arriba)
      if (e_f .ne. N+1)  z_loc(2) = Z(e_f+1) - z_f !upward (+) (interfaz de abajo )
      
      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "eGAeNU iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente 
      gamma = gamma_E(ik,e_f)
      nu = nu_E(ik,e_f)
      
      do iIf = 1,2
          egamz(iIf) = exp(-UI*gamma*ABS(z_loc(iIf)))
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))
          
          gamz(iIf) = gamma*ABS(z_loc(iIf))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do 
      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincGamma(1) = UR * egamz(1)
        sincGamma(2) = UR * egamz(2)
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
      ! ahora estamos involucrando información del receptor en el vector de fuente
      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
        A = pxi%length * 0.5_8  ! 2a=lenght
        cose => pxi%cosT
        seno => pxi%sinT
      
      ! en cada interfaz 
        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
        argum = gamma*A*seno + UR * (k*A*cose)
        sincGamma(2) = sincmod(argum,gamz(2))*(2*A)
        
        argum = nu*A*seno + UR * (k*A*cose)
        sincNu(2) = sincmod(argum,nuz(2))*(2*A)
        
        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
        argum = -gamma*A*seno + UR * (k*A*cose)
        sincGamma(1) = sincmod(argum,gamz(1))*(2*A)
        
        argum = -nu*A*seno + UR * (k*A*cose)
        sincNu(1) = sincmod(argum,nuz(1))*(2*A)
        
      end if !fuente tipo 0 o 1
      
      !#< r obtener vector de Coeficientes de ondas final !#>
      ! multiplicar por los coeficientes de cada interfaz para Ga y Nu
      B(1:tam,ik) = CoefparGa(1:tam,e_f,ik,dj,1)* sincGamma(1) + &
                    CoefparNu(1:tam,e_f,ik,dj,1)* sincNu(1) 
      if (e_f .ne. N+1) then
      B(1:tam,ik) = B(1:tam,ik) + & 
                    CoefparGa(1:tam,e_f,ik,dj,2)* sincGamma(2) + &
                    CoefparNu(1:tam,e_f,ik,dj,2)* sincNu(2)
      end if
!     print*,"eGAeNU",ik,sum(B(1:tam,ik))!;stop
      end subroutine  eGAeNU
      
      

! G_stra - termIndSH
      subroutine SHvectorB_ondaplana(this_B,gamma)
      use soilvars, only : n,amu0,amu,beta0,Z,beta
      use glovars, only:UR,UI,z0,PWfrecReal
      use waveNumVars, only : cOME,ome
      use debugStuff
      use sourceVars, only: Po,iFte=>currentiFte
      implicit none
      complex*16, intent(inout), dimension(1:4*N+2) :: this_B
      real*8, intent(in) :: gamma
      integer :: e
      complex*16 :: kx,kz,V,am,c
      real*8 :: z_loc
            
      e = Po(iFte)%layer ! puede ser 1 o N+1
      
        if (PWfrecReal) then
        c = beta0(e) !SV
        else
        c = beta(e) !SV
        end if
      
      if (PWfrecReal) then
      kx = ome/c * sin(gamma)
      kz = ome/c * cos(gamma)
      am = AMU0(e)
      else
      kx = come/c * sin(gamma)
      kz = come/c * cos(gamma)
      am = AMU(e)
      end if
      
      this_B(1:2*N+1) = Z0 
      if (Z(0) .lt. 0.0) then ! Semiespacio en z<0 ···········
!       i = 1                                                !
!       this_B(1+2*(e-1)-1 + i) = V !  v                     !
        stop "SHvectorB_ondaplana Semiespacio arriba no implementado"
      end if                                                 !  
      ! ······················································
     
!     if (e .eq. 1) then ! Colocamos la onda incindente en la frontera libre.
!         z_loc = 0
!     else
          z_loc = Z(N+1) ! La incidencia en la última interfaz
!     end if
      
       ! Tracciones libres en z=0
       V = UR * exp(UI * kz * (z_loc))
       this_B(1+2*(e-1)  ) = UI * am * V * (- kx * 0 + kz * 1)
       
      ! Continuidad en z = z1
      if (N .ge. 1) then
       this_B(1+2*(e-1)-1) = V
      end if
      
      end subroutine SHvectorB_ondaplana
      
!     subroutine SHvectorB_ondaplana(this_B,gamma)
!     use soilvars, only : n,amu0,amu,beta0,Z,beta
!     use glovars, only:UR,UI,z0,PWfrecReal
!     use waveNumVars, only : cOME,ome
!     use debugStuff
!     implicit none
!     complex*16, intent(inout), dimension(1:4*N+2) :: this_B
!     real*8, intent(in) :: gamma
!     integer :: i,e
!     complex*16 :: kx,kz,V,am,c
!     real*8 :: z_loc
!     !     Colocamos la onda incindente en la interfaz
!     !     con el semiespacio de abajo.
!     e = N+1 
!     z_loc = 0! (Z(N+1)- Z(N+1)) 
!     
!       if (PWfrecReal) then
!       c = beta0(N+1) !SV
!       else
!       c = beta(N+1) !SV
!       end if
!     
!     if (PWfrecReal) then
!     kx = ome/c * sin(gamma)
!     kz = ome/c * cos(gamma)
!     am = AMU0(N+1)
!     else
!     kx = come/c * sin(gamma)
!     kz = come/c * cos(gamma)
!     am = AMU(N+1)
!     end if
!     V = UR * exp(UI * kz * (z_loc))
!     
!     i=0
!     this_B(1:2*N+1) = Z0 
!     if (Z(0) .lt. 0.0) then ! Semiespacio en z<0 ···········
!       i = 1                                                !
!       this_B(1+2*(e-1)-1 + i) = V !  v                     !
!     end if                                                 !  
!     ! ······················································
!     ! Desplazamientos en la frontera de la región de la fuente....
!     if (e .ne. 1) then                                     ! 
!       this_B(1+2*(e-1)-1 + i) = V !  v                     !
!     end if                                                 !
!     !.......................................................
!     ! Tracciones en la frontera de la región de la fuente.........
!     this_B(1+2*(e-1)   + i) = UI * am * V * (- kx * 0 + kz * 1)
!     !............................................................!
!     end subroutine SHvectorB_ondaplana
      
      
      subroutine SHvectorB_force(i_zF,this_B,tam,pXi,cOME,k,ik)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO,NPAR
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only: Po,iFte=>currentiFte! use sourceVars, only : tipofuente!,PW_theta
      use waveNumVars, only : nu_E
      implicit none
      
      integer, intent(in) :: i_zF,tam,ik
      complex*16,    intent(inout), dimension(tam) :: this_B
      real*8,     intent(in)    :: k
      complex*16, intent(in)    :: cOME
      type(Punto),intent(in),target    :: pXi
      
      integer, pointer :: e_f
      real*8, pointer :: z_f,cose,seno
      logical, pointer :: fisInterf
      integer :: iIf,nInterf, el_tipo_de_fuente
      real    :: SGNz
      complex*16 :: nu,DEN, argum,sincmod
      real*8     :: errT = 0.0001_8
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: G22,S22,enuz,nuz,sincNu
      real*8 :: a
      argum = cOME ! (ana mas para que no chiste)
      this_B = Z0
      sincNu = Z0
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
         go to 459
      end if 
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1) z_loc(2) = Z(e_f+1) - z_f !upward (+)
      
      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "SHvectorForce iFte=0"
      if (i_zF .eq. 0) el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual u onda plana)
      !  0  puntual 
      !  2  segmento
      
      G22=Z0;S22=z0 
      DEN = (4.0*AMU(e_f)*UI*PI)
!     omeBet = cOME**2.0/BETA(e_f)**2.0
!     nu = sqrt(omeBet - k**2.0)
!     if(aimag(nu).gt.0.0_8)nu= conjg(nu)
      nu = nu_E(ik,e_f)
      
      do iIf = 1,2
          enuz(iIf) = exp(-UI*nu*ABS(z_loc(iIf)))
          nuz(iIf) = nu*ABS(z_loc(iIf))
      end do
       
      if (el_tipo_de_fuente .le. 1) then ! fuente puntual 0 u onda plana 1
        sincNu(1) = UR * enuz(1)
        sincNu(2) = UR * enuz(2)
      elseif (el_tipo_de_fuente .eq. 2) then ! la fuente es un segmento
      ! ahora estamos involucrando información del receptor en el vector de fuente
      ! en particular, si el receptor (la interfaz) está arriba o abajo de la fuente.
      
        A = pxi%length * 0.5_8 ! 2a=lenght
        cose => pxi%cosT
        seno => pxi%sinT
      
      ! en cada interfaz 
        ! caso (receptor abajo de la fuente) : interfaz de abajo (2)
        argum = nu*a*seno + UR * (k*a*cose)
        sincNu(2) = sincmod(argum,nuz(2))*A*2.0_8
!       !sincNu(2) = enuz(2) * (sin(k*A)/(k*A)*(2*A)) !no se inclina
        
        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
        argum = -nu*a*seno + UR * (k*a*cose)
        sincNu(1) = sincmod(argum,nuz(1))*A*2.0_8
!       !sincNu(1) = enuz(1) * (sin(k*A)/(k*A)*(2*A)) !no se inclina
      end if !fuente tipo 2
      
      ! en cada interfaz (1 arriba) y (2 abajo)
      do iIf = 1,2! nInterf
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if
          
          G22(iIf) = UR/DEN * sincNu(iIf) / nu
          S22(iIf) = -UR / (4.0_8*pi) * sincNu(iIf) * SGNz
          
      end do !iIf interf    
  459 if (fisInterf) then
          S22(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
          G22(1) = 0
      end if
      
      ! El vector de términos independientes genera el campo difractado
      if (z(0) .gt. 0.0) then !--------------------------
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 1) then
        this_B(1+2*(e_f-1)-1) = G22(1)
       end if 
        this_B(1+2*(e_f-1)  ) = S22(1)
      
       if (.not. fisInterf) then ! la fuerza no calló en la interfaz
      !                     =      (2) interfaz de abajo
        if (e_f .ne. N+1) then
         this_B(1+2*(e_f-1)+1) = - G22(2)
         this_B(1+2*(e_f-1)+2) = - S22(2)
        end if
       end if
      else !Z(0) < 0 un semiespacio arriba
      !                     =      (1) interfaz de arriba
       if (e_f .ne. 0) then
        this_B(1+2*(e_f-1)  ) = G22(1)
        this_B(1+2*(e_f-1)+1) = S22(1)
       end if !
      !                     =      (2) interfaz de abajo
       if (.not. fisInterf) then ! la fuerza no calló en la interfaz
        if (e_f .ne. N+1) then
         this_B(1+2*(e_f-1)+2) = - G22(2)
         this_B(1+2*(e_f-1)+3) = - S22(2)
        end if
       end if
      end if
!     print*,"e_f=",e_f
!     print*,this_B;stop 3759
      end subroutine SHvectorB_force
      
      FUNCTION SINCMOD(ARG,ARGZ)  
      use glovars, only : UI    
      COMPLEX*16 :: SINCMOD,ARG, ARGZ
      SINCMOD=EXP(-UI*ARGZ)
      IF(ABS(ARG).LE.0.0001_8)RETURN
      SINCMOD=(EXP(UI*(ARG-ARGZ))-EXP(-UI*(ARG+ARGZ)) )/2.0_8/UI/ARG
      END function SINCMOD
! G_stra - results 
      ! coeficientes de las ondas planas emitidias en cada interfaz 
      ! para representar el campo difractado por estratigrafía
      function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k,ik)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0!,PrintNum,verbose
      use resultVars, only : MecaElem
      use waveNumVars, only : gamma_E,nu_E,ome
      use refSolMatrixVars, only : subMatD0,subMatS0
      implicit none
      
      complex*16, dimension(1:5) :: PSVdiffByStrata
      real*8, intent(in)           ::  z_i,k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e,ik
      complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
      complex*16 :: gamma,nu,xi,eta
      complex*16 :: egammaN,enuN,egammaP,enuP
      complex*16, dimension(2,4) :: subMatD
      complex*16, dimension(3,4) :: subMatS
      complex*16, dimension(1:4) :: coeffsPSV
      integer :: i 
      
      PSVdiffByStrata = 0
      if (ik .ne. 0) then
       gamma = gamma_E(ik,e)
       nu = nu_E(ik,e)
       xi = k**2.0 - nu**2.0
       eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
!     print*,ik,gamma,nu
      else !ik = 0
       gamma = sqrt(UR*OME**2.0_8/ALFA0(e)**2.0_8 - k**2.0_8)
       nu = sqrt(UR*OME**2.0_8/BETA0(e)**2.0_8 - k**2.0_8)
       if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
       if(aimag(nu).gt.0.0)nu= conjg(nu)
       xi = k**2.0 - nu**2.0
       eta = 2.0*gamma**2.0 - OME**2.0 / BETA0(e)**2.0
      end if
      !downward waves
        if (e /= 0) then !(radiation condition upper HS)
          egammaN = exp(-UI * gamma * (z_i-Z(e))) 
          enuN = exp(-UI * nu * (z_i-Z(e)))
        else
          egammaN = Z0
          enuN = Z0
        end if
      !upward waves 
        if (e /= N+1) then !(radiation condition)
          egammaP = exp(UI * gamma * (z_i-Z(e+1)))
          enuP = exp(UI * nu * (z_i-Z(e+1)))
        else
          egammaP = Z0
          enuP = Z0
        end if
          
      !la matrix diagonal
      !coeficientes de las ondas en el estrato
      i = 0
      if (z(0) .lt. 0.0) i = 2
      
        if (e .eq. 0) then! semiespacio de arriba
          coeffsPSV(1:2) = (/z0,z0/)
          coeffsPSV(3:4) = coefOndas_PSV(1 : 2)
        elseif (e .eq. N+1) then ! semiespacio de abajo
          coeffsPSV(1:2) = coefOndas_PSV(4*(e-1)+1+i : 4*(e-1)+2+i)
          coeffsPSV(3:4) = (/z0,z0/)
        else! estrato
          coeffsPSV(1:4) = coefOndas_PSV(4*(e-1)+1+i : 4*(e-1)+4+i)
        end if      
  
      ! desplazamientos
        ! {W}
        ! {U}
        subMatD = subMatD0(1:2,1:4,e,ik)
!       print*,ik,e,sum(subMatD0(1:2,1:4,e,ik))
        subMatD(:,1) = subMatD(:,1)*egammaN*coeffsPSV(1)
        subMatD(:,2) = subMatD(:,2)*enuN*coeffsPSV(2)
        subMatD(:,3) = subMatD(:,3)*egammaP*coeffsPSV(3)
        subMatD(:,4) = subMatD(:,4)*enuP*coeffsPSV(4)
        
        PSVdiffByStrata(1) = sum(subMatD(1,:)) !W
        PSVdiffByStrata(2) = sum(subMatD(2,:)) !U
!       print*,ik,"e",egammaN,enuN
!       print*,ik,"e",egammaP,enuP
!       print*,ik,coeffsPSV(1)
!       print*,ik,coeffsPSV(2)
!       print*,ik,coeffsPSV(3)
!       print*,ik,coeffsPSV(4)
!       print*,ik,"(1)_",PSVdiffByStrata(1)
!       print*,ik,"(2)_",PSVdiffByStrata(2)
        
       
      ! esfuerzos
        subMatS = subMatS0(1:3,1:4,e,ik)
        subMatS(:,1) = subMatS(:,1)*egammaN*coeffsPSV(1)
        subMatS(:,2) = subMatS(:,2)*enuN*coeffsPSV(2)
        subMatS(:,3) = subMatS(:,3)*egammaP*coeffsPSV(3)
        subMatS(:,4) = subMatS(:,4)*enuP*coeffsPSV(4)     
        
        PSVdiffByStrata(3) = sum(subMatS(1,:)) !s33
        PSVdiffByStrata(4) = sum(subMatS(2,:)) !s31
        PSVdiffByStrata(5) = sum(subMatS(3,:)) !s11 
        
!       if (sum(abs(PSVdiffByStrata)) .gt. 10000) then
!       print*,ik,e
!       print*,"coeffsPSV ",coeffsPSV(1:4) !ok
!       print*,"subMatD0(1:2,1:4,e,ik)",subMatD0(1:2,1:4,e,ik)
!       print*,"subMatS0(1:3,1:4,e,ik)",subMatS0(1:3,1:4,e,ik)
!       print*,"egammaN",egammaN
!       print*,"enuN",enuN
!       print*,"egammaP",egammaP
!       print*,"enuP",enuP
!       print*,"   W,U",PSVdiffByStrata(1:2)
!       print*,"   str",PSVdiffByStrata(3:5)
!       
!       call system("echo -e ""\033[31m problem \033[0m here""")
!       end if
      end function PSVdiffByStrata
      
      function SHdiffByStrata(coefOndas_SH,z_i,e,cOME_i,k,ik)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
      use waveNumVars, only : nu_E,ome
!     use resultVars, only : MecaElem
      
      implicit none
      complex*16, dimension(1:3)   :: SHdiffByStrata
      real*8, intent(in)           :: z_i,k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e,ik
      complex*16, dimension(2*N+1),intent(in) :: coefOndas_SH
      complex*16 :: nu
      complex*16 :: enuN,enuP
      complex*16, dimension(1,2) :: subMatDsh
      complex*16, dimension(2,2) :: subMatSsh
      complex*16, dimension(2,2) :: diagMatSH
      complex*16, dimension(2,1) :: coeffsSH
      complex*16, dimension(1,1) :: resDsh
      complex*16, dimension(2,1) :: resSsh
      integer :: i 
      SHdiffByStrata = z0
      ! algunas valores constantes para todo el estrato
      !nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
      !if(aimag(nu).gt.0.0)nu= conjg(nu)
      if (ik .ne. 0) then
       nu = come_i ! (nada mas para que no chiste)
       nu = nu_E(ik,e)
      else ! ik = 0
       nu = sqrt(UR*OME**2.0_8/BETA0(e)**2.0_8 - k**2.0_8)
       if(aimag(nu).gt.0.0)nu= conjg(nu)
      end if
          !downward waves
          if (e /= 0) then !(radiation condition upper HS)
             enuN = exp(-UI * nu * (z_i-Z(e)))! enuN = exp(-UI * nu * abs(z_i-Z(e)))
          else
            enuN = Z0
          end if
          !upward waves 
          if (e /= N+1) then !(radiation condition)
            enuP = exp(UI * nu * (z_i-Z(e+1)))! enuP = exp(-UI * nu * abs(z_i-Z(e+1)))
          else
            enuP = Z0
          end if
          
          !la matrix diagonal
          diagMatSH = RESHAPE((/ enuN, Z0,   & 
                                 Z0,   enuP /), &
                              (/ 2,2 /))
      !coeficientes de las ondas en el estrato
      i = 0
      if (z(0) .lt. 0.0) i = 1
        if (e .eq. 0) then! semiespacio de arriba
          coeffsSH(1:1,1) = z0
          coeffsSH(2:2,1) = coefOndas_SH(1:1)
        elseif (e .eq. N+1) then! semiespacio de abajo
          coeffsSH(1:1,1) = coefOndas_SH(2*(e-1)+1+i : 2*(e-1)+1+i)
          coeffsSH(2:2,1) = z0
        else! estrato
          coeffsSH(1:2,1) = coefOndas_SH(2*(e-1)+1+i : 2*(e-1)+2+i)
        end if
        
      ! desplazamientos
!     if (mecStart .eq. 1)then
        subMatDsh = RESHAPE((/ UR,UR /), (/ 1,2 /))
        subMatDsh = matmul(subMatDsh,diagMatSH)
        resDsh = matmul(subMatDsh, coeffsSH)
        SHdiffByStrata(1) = resDsh(1,1) !V
!     end if ! desplazamientos
      
      ! esfuerzos
!     if (mecEnd .eq. 3) then
        subMatSsh = RESHAPE((/ -UI*nu,-UI*k,UI*nu,-UI*k/),(/2,2/)) 
        subMatSsh = amu(e) * subMatSsh
        subMatSsh = matmul(subMatSsh,diagMatSH)
        resSsh = matmul(subMatSsh, coeffsSH)
        SHdiffByStrata(2) = resSsh(1,1) !s32
        SHdiffByStrata(3) = resSsh(2,1) !s12
!     end if ! esfuerzos
          
      end function SHdiffByStrata
      
! G_full PSV   
      Subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
!#define ver 1
      ! Sanchez-Sesma y Campillo, 1991.  Mismo resultado que:
      ! Kaussel, Fundamental solutions in elastodynamics... pag 38
      use soilVars ,only : alfa0,beta0,Lambda0,AMU0,alfa,beta,amu,lambda,rho,N!,Z
      use gloVars, only:UI,UR,one,z0,PWfrecReal
      use hank !     use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: Po,iFte=>currentiFte!tipofuente, PW_pol
      use waveNumVars, only : OME
      use Gquadrature, only : Gquad_n      
      implicit none
      interface
        subroutine greenexPSV(greenex,dir_j,p_X,pXi,estrato,cOME)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: greenex
         integer,intent(in) :: dir_j, estrato
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16,intent(in) :: cOME
        end subroutine greenexPSV
      end interface
      type (FFres), intent(out) :: FF
      integer,    intent(in)  :: mecS,mecE, i_zF
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16, intent(in) :: cOME
      integer, intent(in) :: dir_j
      real*8 :: r,gamma(2)
      complex*16 :: A,B,C,Dqr,Dkr,kx,kz
      complex*16 :: omeP,omeS
      complex*16 :: H0s,H1s,H2s,H0p,H1p,H2p !Hankel 
      complex*16 :: szz,szx,sxx
      complex*16 :: la,am
      integer :: i,j
      integer, pointer :: e
      real*8 :: nX(2)
      integer :: iGq,nGq
      real*8, pointer :: xf,zf,GqC
      real*8 :: deltaij
      real*8,dimension(2) :: theta
      integer :: el_tipo_de_fuente
      integer, target :: estrato
      logical :: shouldI,XinoEstaEnInterfaz,usarGreenex,estratosIguales
#ifdef ver 
      print*,""     
      print*,"-------------------------------------------"
      print*,"i_zF=",i_zF
      print*,"px",p_x%center,"e=",p_x%layer
      print*,"pxi",pXi%center,"e=",pXi%layer
      print*,"dir_j=",dir_j
#endif     
      FF%W=z0; FF%U=z0;FF%Tz=z0;FF%Tx=z0
      FF%sxx = z0;FF%szx = z0;FF%szz = z0 
      estratosIguales = .false.     
      XinoEstaEnInterfaz = .false.
      usarGreenex = .false.
      shouldI = .false.
      if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
      if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
      j = dir_j ! <***********  dir_j = 3  (vertical)
      if(j .eq. 3) j = 2 ! para coincidir con los indicies
      el_tipo_de_fuente = 2 !(fuente segmento)
      nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
      xf => one; zf => one
      
      if (i_zF .eq. 0) then ! es la fuente real
         if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "FFPSV iFte=0" !failsafe
         el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual:0, onda plana:1, segmento:2)
         if (pXi%region .eq. 2) then !en la inclusión R
            shouldI = .true.
            XinoEstaEnInterfaz = .true.
            estrato = N+2
            e => estrato
         else !pXi%region = 1 E
            if (p_x%layer .eq. pXi%layer) shouldI = .true.
            if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
            estrato = p_x%layer
            e => p_x%layer
         end if
      else ! i_zf .ne. 0  una fuente virtual de ibem
         el_tipo_de_fuente = 2 !(fuente segmento)
         XinoEstaEnInterfaz = .true. !(se ha encargado la geometría)
         if (i_zF .eq. -1) then !(campo refractado en inclusion R)
            shouldI = .true.
            estrato = N+2 !(para tomar las propiedades de la inclusión)
            e => estrato
         else !(normal: campo refractado en medio estratificado E)
            if (p_x%layer .eq. pXi%layer) then 
               shouldI = .true. ! si están en el mismo estrato
               estrato = p_x%layer
               e => p_x%layer
            end if
         end if
      end if
#ifdef ver
      print*,"shouldI=",ShouldI," estrato=",estrato
      print*,"el_tipo_de_fuente=",el_tipo_de_fuente
#endif
      if (shouldI) then  
      if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then ! onda plana !
#ifdef ver
      print*,"onda plana"
#endif
!        if (p_x%isOnInterface .eqv. .true.) return  ! creo
         if (Po(iFte)%PW_pol .eq. 1) then !SV
            if (PWfrecReal) then
               c = UR*beta0(e)
            else
               c = beta(e)
            end if
            theta(1) = cos(Po(iFte)%gamma)
            theta(2) = sin(Po(iFte)%gamma)
!           print*,theta(1),cos(pxi%gamma)
!           print*,theta(2),sin(pxi%gamma);stop
         elseif (Po(iFte)%PW_pol .eq. 2) then
            if (PWfrecReal) then !P
                c = UR*alfa0(e)
            else
                c = alfa(e) !SV
            end if 
            theta(1) = sin(Po(iFte)%gamma)
            theta(2) = -cos(Po(iFte)%gamma)
!           theta(1) = sin(pxi%gamma)
!           theta(2) = -cos(pxi%gamma); stop
         end if!
         if (PWfrecReal) then
            kx = UR*real(ome/c * sin(Po(iFte)%gamma))
            kz = UR*real(ome/c * cos(Po(iFte)%gamma))
            la = UR*real(LAMBDA0(e))
            am = UR*real(AMU0(e))
         else
            kx = UR*real(come/c * sin(Po(iFte)%gamma))
            kz = UR*real(come/c * cos(Po(iFte)%gamma))
            la = UR*real(LAMBDA(e))
            am = UR*real(AMU(e))
         end if
!       print*,"kxkzlaam",kx,kz,la,am
        FF%U = (theta(1))* exp(UI * kz * (p_x%center%z - 0)) &         !
              * exp(-UI * kx * (p_x%center%x - Po(iFte)%center%x))          !
        FF%W = (theta(2))* exp(UI * kz * (p_x%center%z - 0)) &         !
              * exp(-UI * kx * (p_x%center%x - Po(iFte)%center%x))          !
        szz = UI * ( &                                                      !
                    ( FF%W * kz * (la + 2.0* am)) &                         !
                  - ( FF%U * kx * la))                                      !
        szx = UI * am * ( kz * FF%U - kx * FF%W )                           !
        sxx = UI * ( &                                                      !
                  - ( FF%U * kx * (la + 2.0*am)) &                          !
                  + ( FF%W * kz * la))                                      !
        FF%Tx = sxx * nx(1) + szx * nx(2)                                   !
        FF%Tz = szx * nx(1) + szz * nx(2)                                   !
        FF%sxx = sxx                                                        !
        FF%szx = szx                                                        !
        FF%szz = szz                                                        !
#ifdef ver
      print*,FF%W,"FF%W"
      print*,FF%U,"FF%U"
      print*,FF%Tx,"FF%Tx"
      print*,FF%Tz,"FF%Tz" 
#endif
        return                                                              !
      end if! fin onda plana -´-´-´-´-´-´-´-´-´-´-´-´-´-´--´-´ fin onda plana          
      if (XinoEstaEnInterfaz .eqv. .true.) then 
         xf => pXi%center%x
         zf => pXi%center%z
         r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.) 
         
        ! para funciones de Greeen en frontera con cont. de desplazamiento
        if ((p_x%isboundary .eqv. .true.) .and. &
            (pXi%isboundary .eqv. .true.) .and. & 
            (p_x%boundaryIndex .eq. pXi%boundaryIndex)) usarGreenex = .true.
!       print*,usarGreenex
        ! para campo cercano por fuente segmento
        if ((i_zF .eq. 0) .and. &
            (el_tipo_de_fuente .eq. 2) .and. &      !fuente real
            (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
!       print*,usarGreenex    
!       if ((p_X%isboundary .eqv. .false.) .and. &
!           (pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
!           (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
        if ((pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
            (abs(r) .lt. pXi%length/2)) usarGreenex = .true.

!       print*,usarGreenex    
        ! si fuente real y cilindrica entonces no
        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false. 
!       print*,usarGreenex
        ! si es un puntno para sobredeterminar el sistema entonces no
        if (p_X%isOD) usarGreenex = .false.
!       print*,usarGreenex
!       !#< r
!       usarGreenex = .false.
!       !#>
#ifdef ver
      print*,"onda cilíndrica r=",r, "usarGreenex=",usarGreenex
#endif
      !Con integración de Gauss
      if (usarGreenex .eqv. .false.) then
         if (el_tipo_de_fuente .eq. 0) then
            nGq = 1
            GqC => ONE
         else !.eq. 2 (fuente segmento suficientemente alejado)
            if (pXi%isBoundary) nGq = Gquad_n ! IBEM
            if (pXi%isSourceSegmentForce) nGq = Gquad_n 
            if (nGq .ne. Gquad_n) stop "chin 6600"
         end if
         do iGq = 1,nGq
            if (nGq .gt. 1) then
                xf => pXi%Gq_xXx_coords(iGq,1)
                zf => pXi%Gq_xXx_coords(iGq,2)
                GqC => pXi%Gq_xXx_C(iGq) 
            end if
            r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
#ifdef ver
      print*,"iGq=",iGq,"xf,zf,r,Gqc=",xf,zf,r,Gqc
#endif
      gamma(1) = (p_X%center%x - xf) / r ! gamma x
      gamma(2) = (p_X%center%z - zf) / r ! gamma z
      omeP = cOME * r / alfa(e)
      omeS = cOME * r / beta(e)
      ! funcs de Hankel de segunda especie
      call hankels(omeP,H0p,H1p)
      H2p = -H0p + 2./omeP * H1p
      call hankels(omeS,H0s,H1s)
      H2s = -H0s + 2./omeS * H1s
      
      A = H0p/alfa(e)**2. + H0s/beta(e)**2.
      B = H2p/alfa(e)**2. - H2s/beta(e)**2.
#ifdef ver
      print*,"gamma,omeP,omeS,A,B",gamma(1),gamma(2)
      print*,omeP,omeS
      print*,A,B
#endif      
      ! desplazamientos 
      if (mecS .eq. 1) then
      ! W
      i = 2
      FF%W = FF%W + (-UI/8.0/rho(e)*& 
      (A*deltaij(i,j)-(2.*gamma(i)*gamma(j)-deltaij(i,j))*B)) * GqC
      ! U
      i = 1
      FF%U = FF%U + (-UI/8.0/rho(e)*& 
      (A*deltaij(i,j)-(2.*gamma(i)*gamma(j)-deltaij(i,j))*B)) * GqC
      end if!mecs1
      
      if (mecE .eq. 5) then
      nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
      ! tracciones
      Dqr = omeP*H1p
      Dkr = omeS*H1s
      C = Dqr/alfa(e)**2. - Dkr/beta(e)**2.
      
      ! TZ
      i = 2      
      FF%Tz = FF%Tz + (& 
      amu(e)*UI /(2.*rho(e)*r)*(& 
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* & 
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))& 
      ) * GqC
      
      ! TX
      i = 1
      FF%Tx = FF%Tx + ( & 
      amu(e)*UI /(2.*rho(e)*r)*(& 
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* & 
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))& 
      ) * GqC
      
      ! sxx
      i = 1
      nX(1) = 1; nX(2) = 0
      FF%sxx = FF%sxx + ( & 
      amu(e)*UI /(2.*rho(e)*r)*(& 
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* & 
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))& 
      ) * GqC
      ! szx
      i = 1
      nX(1) = 0; nX(2) = 1
      FF%szx = FF%szx + ( & 
      amu(e)*UI /(2.*rho(e)*r)*(& 
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* & 
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))& 
      ) * GqC
      ! szz
      i = 2
      nX(1) = 0; nX(2) = 1
      FF%szz = FF%szz + (& 
      amu(e)*UI /(2.*rho(e)*r)*(& 
      (B+(lambda(e)*Dqr)/(2.*amu(e)*alfa(e)**2.))*gamma(j)*nX(i) &
      + (B + Dkr/(2.*beta(e)**2.))* & 
      (gamma(i)*nX(j) + (gamma(1)*nX(1) + gamma(2)*nX(2))*deltaij(i,j)) + &
      (C-4.*B)*gamma(i)*gamma(j)*(gamma(1)*nX(1) + gamma(2)*nX(2)))& 
      ) * GqC
      
      end if !mece5
      end do ! iGq
      else !. usarGreenex ..................................................
#ifdef ver
      print*,"Greenex"
#endif
        call greenexPSV(FF,dir_j,p_X,pXi,estrato,cOME) ! W,U
!     print*,FF%W,"FF%W"
!     print*,FF%U,"FF%U"
!     print*,FF%Tx,"FF%Tx"
!     print*,FF%Tz,"FF%Tz"  
!     print*,"Greenex";stop
!        if (p_X%isOD) then 
!          print*,p_x%center,"  ->  ",pXi%center
!          print*,"tx=",FF%Tx
!          print*,"tz=",FF%Tz
!          stop "OD in FFpsv"
!          if (dir_j .eq. 1) then; FF%Tx=0.5*UR; FF%Tz=z0;end if!
!          if (dir_j .eq. 3) then; FF%Tx=z0;     FF%Tz=0.5*UR; end if
!        end if
      end if ! greenex o gaussiana
      end if ! XinoEstaEnInterfaz: on the interface
#ifdef ver
      print*,FF%W,"FF%W"
      print*,FF%U,"FF%U"
      print*,FF%Tx,"FF%Tx"
      print*,FF%Tz,"FF%Tz"
      print*,""
!     stop "FFpsv"
#endif
      end if !should I?
      end subroutine FFpsv
           
      subroutine greenexPSV(FF,dir_j,p_X,pXi, e,cOME)
      use glovars, only : UR,UI,pi,z0
      use resultvars, only : Punto,FFres
      use soilVars ,only : alfa,beta,rho,amu,lambda
      implicit none
      type (FFres), intent(out) :: FF !U y W dados dir_j 1 o 3
      integer,intent(in) :: dir_j,e
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16,intent(in) :: cOME !frecuencia compleja (rad/s)
      complex*16 :: AQA,AKA,BEALF,BEALF2,ALFBE
      complex*16 :: argk,ark2,h0kr,argq,arq2,h0qr,AA,BB,h2kr,h2qr
      complex*16 :: GN11,GN31,GN13,GN33
      complex*16 :: Dqr,Dkr,CC,tt1,tt2,tt3,tt4
      real*8 :: G1,G3,depi,RNM,DSM,c1,c2,c13,c23,geu
!     stop "Greeneer"
      FF%W=z0;FF%U=z0;FF%Tx=z0;FF%Tz=z0
!     if (dir_j .eq. 1) then; FF%Tx=0.5*UR; FF%Tz=z0;end if
!     if (dir_j .eq. 3) then; FF%Tx=z0;     FF%Tz=0.5*UR; end if
!     FF%sxx = 0;FF%szx = 0;FF%szz = 0 
      GEU=0.5772156649_8
      
!     e = p_X%layer !N+2
      BEALF = beta(e)/alfa(e)
      ALFBE = 1.0/BEALF
      RNM = sqrt((p_x%center%x-pXi%center%x)**2. + & 
                 (p_x%center%z-pXi%center%z)**2.) ! distancia x a xi
      AKA = cOME/beta(e)! cOME !cOME/beta(e) !
      AQA = cOME/alfa(e)! AKA*BEALF !cOME/alfa(e) !
      DSM = pXI%length 
      IF(RNM .lt. 0.001_8)THEN
         ! vector tangente a la normal del segmento
         G1 = - pXi%normal%z 
         G3 = + pXi%normal%x
      ELSE 
         ! cosenos directores 
         G1=(p_X%center%x - pXi%center%x)/ RNM
         G3=(p_X%center%z - pXi%center%z)/ RNM   
      END IF
      
      DEPI=2.0/PI
      BEALF2=BEALF*BEALF
         !new green function (with analytic integration)                                               
                                         
         C1=1.0+RNM/DSM*2.0
         C2=1.0-RNM/DSM*2.0
         C13=C1*C1*C1
         C23=C2*C2*C2
         GEU=0.5772156649_8
         ARGK=AKA*DSM 
         ARK2=ARGK*ARGK
         H0KR=UR*(1.0-ARK2*(C13+C23)/96.0) &
          -UI*DEPI*( GEU-1.0 + 0.5*C1*LOG(ARGK*C1/4.0) &
                             + 0.5*C2*LOG(ARGK*C2/4.0) & 
                    +(4./3.-GEU)*ARK2*(C13+C23)/96.0 &
                    -ARK2*C13*LOG(ARGK*C1/4.0)/96.0 & 
                    -ARK2*C23*LOG(ARGK*C2/4.0)/96.0 )  
         ARGQ=AQA*DSM
         ARQ2=ARGQ*ARGQ
         H0QR=UR*(1.0-ARQ2*(C13+C23)/96.0) &
          -UI*DEPI*( GEU-1.0 + 0.5*C1*LOG(ARGQ*C1/4.0) &
                             + 0.5*C2*LOG(ARGQ*C2/4.0) &
                    +(4./3.-GEU)*ARQ2*(C13+C23)/96.0 &
                    -ARQ2*C13*LOG(ARGQ*C1/4.0)/96.0 &
                    -ARQ2*C23*LOG(ARGQ*C2/4.0)/96.0 )  
            
         AA=H0QR /(alfa(e)**2)+ H0KR /(beta(e)**2)
         AA = AA * pXI%length ! ----------------------------------
         
         h2qr = -UI/PI   
         h2qr = h2qr+ARQ2/192.*(C13+C23)*(UR-UI*DEPI*(GEU-13./12.))
         h2qr = h2qr-ARQ2/192.*UI*DEPI*C13*LOG(ARGQ*C1/4.0)
         h2qr = h2qr-ARQ2/192.*UI*DEPI*C23*LOG(ARGQ*C2/4.0)
         
         h2kr = -UI/PI
         h2kr = h2kr+ARK2/192.*(C13+C23)*(UR-UI*DEPI*(GEU-13./12.))
         h2kr = h2kr-ARK2/192.*UI*DEPI*C13*LOG(ARGK*C1/4.0)
         h2kr = h2kr-ARK2/192.*UI*DEPI*C23*LOG(ARGK*C2/4.0)
         
         BB = h2qr /(alfa(e)**2) - h2kr /(beta(e)**2)
         BB = BB * pXI%length ! ----------------------------------
         
         GN11=-UI/(8.*rho(e))*(AA-(2.0*G1*G1-1.0)*BB)         
         GN33=-UI/(8.*rho(e))*(AA-(2.0*G3*G3-1.0)*BB)         
         GN13=+UI/(4.*rho(e))*G1*G3*BB
         GN31=+UI/(4.*rho(e))*G1*G3*BB
         
         !para esfuerzos
         Dqr = (ARQ2/2)*( H0QR + h2qr )
         Dkr = (ARK2/2)*( H0KR + h2kr )
         CC = Dqr/alfa(e)**2 - Dkr/beta(e)**2
         
         tt1 = (UI*amu(e))/(2*rho(e)*DSM)
         tt2 = BB + (lambda(e)*Dqr)/(2*amu(e)*alfa(e)**2)
         tt3 = BB + (Dkr)/(2*beta(e)**2)
         tt4 = CC - 4*BB
         
         if (dir_j .eq. 1) then
           FF%U = GN11
           FF%W = GN13
           ! i=1 j=1
           FF%Tx = - tt1 * (&
           (tt2)* G1* p_X%normal%x + & 
           (tt3)* (G1* p_X%normal%x + G1* p_X%normal%x + G3* p_X%normal%z) + &
           (tt4)* (G1* G1* (G1* p_X%normal%x + G3* p_X%normal%z)) &
           )
           ! i=3 j=1
           FF%Tz = - tt1 * (&
           (tt2)* G1* p_X%normal%z + & 
           (tt3)* (G3* p_X%normal%x) + &
           (tt4)* (G3* G1* (G1* p_X%normal%x + G3* p_X%normal%z)) &
           )
         else if (dir_j .eq. 3) then
           FF%U = GN31
           FF%W = GN33
           ! i=1 j=3
           FF%Tx = - tt1 * (&
           (tt2)* G3* p_X%normal%x + & 
           (tt3)* (G1* p_X%normal%z) + &
           (tt4)* (G1* G3* (G1* p_X%normal%x + G3* p_X%normal%z)) &
           )
           ! i=3 j=3
           FF%Tz = - tt1 * (&
           (tt2)* G3* p_X%normal%z + & 
           (tt3)* (G3* p_X%normal%z + G1* p_X%normal%x + G3* p_X%normal%z) + &
           (tt4)* (G3* G3* (G1* p_X%normal%x + G3* p_X%normal%z)) &
           )
         end if!
         if (.not. p_x%isOD) then
           FF%Tx = z0;   FF%Tz = z0
!        else
!          stop "greenexPSV"
         end if
!     print*,dir_j,RNM,p_x%center,pXI%center,FF%U,FF%W
      end subroutine greenexPSV
      
! G_full SH 
      subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
!#define ver 1
      use soilVars ,only : amu,amu0,beta,beta0,N!,z
      use gloVars, only:UR,UI,z0,ONE,PWfrecReal
      use hank !use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: Po,iFte=>currentiFte!use sourceVars, only: tipofuente
      use waveNumVars, only : OME
      use Gquadrature, only : Gquad_n
      implicit none
      interface
        subroutine greenexSH(FF,p_X,pXi,e,cOME)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF 
         integer,intent(in) :: e
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16,intent(in) :: cOME
        end subroutine greenexSH
      end interface
      type (FFres), intent(out) :: FF !V,Ty
      integer,    intent(in)     :: mecS,mecE, i_zF
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16, intent(in)     :: cOME
      
      real*8 :: r,gamma(2)!,longonda
      complex*16 :: omeS,c,am
      complex*16 :: H0s,H1s,kx,kz,szy,sxy
      integer, pointer :: e
      integer :: iGq,nGq
      real*8 ::  nX(2)
      real*8, pointer :: xf,zf,GqC
      integer :: el_tipo_de_fuente
      integer, target :: estrato
      logical :: shouldI,estratosIguales,XinoEstaEnInterfaz,usarGreenex
#ifdef ver 
      print*,""     
      print*,"-------------------------------------------"
      print*,"i_zF=",i_zF
      print*,"px",p_x%center,"e=",p_x%layer
      print*,"pxi",pXi%center,"e=",pXi%layer
#endif     
      FF%V=z0;FF%Ty=z0
      estratosIguales = .false.
      XinoEstaEnInterfaz = .false.
      usarGreenex = .false.
      shouldI = .false.
      if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
      if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
      el_tipo_de_fuente = 2 ! fuente segmento para IBEM
      nx(1) = p_X%normal%x;nx(2) = p_X%normal%z
      xf => one;zf => one
      
      if (i_zF .eq. 0) then
         if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "FFSH iFte=0"
         el_tipo_de_fuente = Po(iFte)%tipofuente !(puntual u onda plana) 
            if (pXi%region .eq. 2) then !en la inclusión R
               shouldI = .true.
               XinoEstaEnInterfaz = .true.
               estrato = N+2
               e => estrato
            else ! en el exterior E
               if (p_x%layer .eq. pXi%layer) shouldI = .true.
               if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
               estrato = p_x%layer
               e => p_x%layer
            end if 
      else !i_zf .ne. 0  una fuente virtual de ibem
         el_tipo_de_fuente = 2 !(fuente segmento)
         XinoEstaEnInterfaz = .true. !(se ha encargado la geometría)
         if (i_zF .eq. -1) then !(Fza distr en segmento en región R) 
            shouldI = .true.
            estrato = N+2 !(para tomar las propiedades de la inclusión)
            e => estrato
         else !(normal: campo refractado en medio estratificado E)
            if (p_x%layer .eq. pXi%layer) then 
               shouldI = .true. ! si están en el mismo estrato
               estrato = p_x%layer
               e => p_x%layer
            end if
         end if
      end if! 
#ifdef ver
      print*,"shouldI=",ShouldI," estrato=",estrato
      print*,"el_tipo_de_fuente=",el_tipo_de_fuente
#endif
      if (shouldI) then  
      if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then ! onda plana !
#ifdef ver
      print*,"onda plana"
#endif     
!        if (Po(iFte)%PW_pol .eq. 3) then ! SH
            if (PWfrecReal) then
               c = UR*beta0(e)
               kx = UR*real(ome/c * sin(Po(iFte)%gamma))
               kz = UR*real(ome/c * cos(Po(iFte)%gamma))
               am = UR*real(AMU0(e))
            else
               c = beta(e)
               kx = UR*real(come/c * sin(Po(iFte)%gamma))
               kz = UR*real(come/c * cos(Po(iFte)%gamma))
               am = UR*real(AMU(e))
            end if
!           if (PWfrecReal) then
!              c = UR*beta0(N+1)
!              kx = UR*real(ome/c * sin(Po(iFte)%gamma))
!              kz = UR*real(ome/c * cos(Po(iFte)%gamma))
!              am = UR*real(AMU0(N+1))
!           else
!              c = beta(N+1)
!              kx = UR*real(come/c * sin(Po(iFte)%gamma))
!              kz = UR*real(come/c * cos(Po(iFte)%gamma))
!              am = UR*real(AMU(N+1))
!           end if
!        end if
        ! las expresiones se encuentran en todos lados, por ejemplo:
        ! Gil-Zepeda et al 2003
        ! SV:
        FF%V = (UR)* exp(-UI*(kx * (p_x%center%x - Po(iFte)%center%x) - &
                              kz * (p_x%center%z - 0)))
!       szy, sxy
        sxy = am * FF%V * (-UI * kx)
        szy = am * FF%V * (UI * kz)
        FF%Ty = sxy * nx(1) + szy * nx(2)
!       !#< r  
!       FF%V = 0
!       FF%Ty = 0
!       print*,"lin 4713"
!       !#>
#ifdef ver
      print*,FF%V,"FF%V"
      print*,FF%Ty,"FF%Ty"
#endif 
        return
      end if! onda plana
      if (XinoEstaEnInterfaz .eqv. .true.) then !..............................
        xf => pXi%center%x
        zf => pXi%center%z
        r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
        
        ! para funciones de Greeen en frontera con cont. de desplazamiento
        if ((p_x%isboundary .eqv. .true.) .and. &
            (pXi%isboundary .eqv. .true.) .and. & 
            (p_x%boundaryIndex .eq. pXi%boundaryIndex)) usarGreenex = .true.
!       print*,usarGreenex
        ! para campo cercano por fuente segmento
        if ((i_zF .eq. 0) .and. &
            (el_tipo_de_fuente .eq. 2) .and. &      !fuente real
            (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
!       print*,usarGreenex    
!       if ((p_X%isboundary .eqv. .false.) .and. &
!           (pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
!           (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
        if ((pXi%isboundary .eqv. .true.) .and. &   !fuente virtual
            (abs(r) .lt. pXi%length/2)) usarGreenex = .true.

!       print*,usarGreenex    
        ! si fuente real y cilindrica entonces no
        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false. 
!       print*,usarGreenex
        ! si es un puntno para sobredeterminar el sistema entonces no
        if (p_X%isOD) usarGreenex = .false.
!       print*,usarGreenex
!       !#< r
!       usarGreenex = .false.
!       !#>
#ifdef ver
      print*,"onda cilíndrica r=",r, "usarGreenex=",usarGreenex
#endif
      !Con integración de Gauss
      if (usarGreenex .eqv. .false.) then
         if (el_tipo_de_fuente .eq. 0) then                             !
           nGq = 1                                                      !
           GqC => ONE                                                   ! 
         else ! eq. 2; print*,"fuente segmento "                        !
           if (pXi%isBoundary) nGq = Gquad_n ! IBEM                     !
           if (pXi%isSourceSegmentForce) nGq = Gquad_n                  !
           if (nGq .ne. Gquad_n) stop "chin 6600"                       !
         end if                                                         !
         do iGq = 1,nGq !..........................................     !
            if (nGq .gt. 1) then                                  !     !
               xf => pXi%Gq_xXx_coords(iGq,1)                     !     !
               zf => pXi%Gq_xXx_coords(iGq,2)                     !     !
               GqC => pXi%Gq_xXx_C(iGq)                           !     !
            end if                                                !     !
      r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)     !     !
#ifdef ver
      print*,"iGq=",iGq,"xf,zf,r,Gqc=",xf,zf,r,Gqc
#endif
      gamma(1) = (p_X%center%x - xf) / r ! gamma x                !     !
      gamma(2) = (p_X%center%z - zf) / r ! gamma z                !     !
      omeS = cOME * r / beta(e)                                   !     !
      call hankels(omeS,H0s,H1s) ! Hankel de segunda especie      !     !  
#ifdef ver
      print*,"omeS=",omeS
      print*,"H0s,H1s=",H0s,H1s
#endif 
      if(mecS .eq. 1) FF%V = FF%V + (-UI/(4.*amu(e))*H0s) * GqC   !     !
      if(mecE .eq. 3) then                                        !     !
!     FF%s32 = FF%s32 + (UI*cOME*(p_x%center%z-zf)) &             !     !
!                        /(4.0_8*beta(e)*r)*H1s * C               !     !
!     FF%s12 = FF%s12 + (UI*cOME*(p_x%center%x-xf))/ &            !     !
!                         (4.0_8*beta(e)*r)*H1s * C               !     !
      FF%Ty = FF%Ty + (UI*cOME/beta(e)/4.0*H1s* &                 !     !
               (gamma(1)*nx(1) + gamma(2)*nx(2))) * GqC           !     !
      end if !mec3                                                !     !
      end do ! iGq ................................................     !
      else !. usarGreenex ...............................
#ifdef ver
      print*,"Greenex"
#endif
        call greenexSH(FF,p_X,pXi,e,cOME)
      end if ! greenex o gaussiana .....................................!
      end if ! XinoEstaEnInterfaz: on the interface
#ifdef ver
      print*,FF%V,"FF%V"
      print*,FF%Ty,"FF%Ty"
      print*,""
!     stop "FFsh"
#endif
      end if ! should I? ....................................................
      end subroutine FFsh
      
      subroutine greenexSH(FF,p_X,pXi,e,cOME)
             
      !funcion de Green G22 analitica en la fuente
      ! se usan las series ascendentes de las funciones de Bessel
      use glovars, only : UR,UI,pi,z0
      use resultvars, only : Punto,FFres
      use soilVars ,only : beta,amu
      implicit none
      type (FFres), intent(out) :: FF 
      integer,intent(in) :: e
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16,intent(in) :: cOME !frecuencia compleja (rad/s)
      
      real*8 :: RNM ! distancia (muy pequeña o 0)
      real*8 :: DSM ! longitud del segmento de integración
      real*8 :: c1,c2,c13,c23,geu,depi
      complex*16 :: aka,argk,ark2,h0kr
      
      FF%V=z0;FF%Ty=z0
      geu=0.5772156649_8
      RNM = sqrt((p_x%center%x-pXi%center%x)**2. + & 
                 (p_x%center%z-pXi%center%z)**2.) ! distancia x a xi
      AKA = cOME/beta(e)! cOME !cOME/beta(e) !
      DSM = pXI%length 
      
      DEPI=2.0/PI
      c1=1.0+RNM/DSM*2.0
      c2=1.0-RNM/DSM*2.0
      c13= c1**3.0
      c23= c2**3.0
      argk= AKA * DSM
      ark2= argk**2.0
      
      h0kr = ur*(1.0-ark2*(c13+c23)/96.0) &
          - ui*depi*( geu-1.0 + 0.5*c1*log(argk*c1/4.0) &
                              + 0.5*c2*log(argk*c2/4.0) &
          +(4.0/3.0-geu)*ark2*(c13+c23)/96.0 &
            - ark2*c13*log(argk*c1/4.0)/96.0 &
            - ark2*c23*log(argk*c2/4.0)/96.0)
      
      FF%V = - ui/(4.0*amu(e)) * h0kr * DSM 
      
      !#< r  
      !Falta Ty  
      !#>
      if (.not. p_x%isOD) then
           FF%Ty = z0
      end if
!     print*,RNM,p_x%center,pXI%center,FF%V
      end subroutine greenexSH
                  
      function deltaij(i,j)
      integer :: i,j
      real*8 :: deltaij
      ! i : dirección del receptor {1;2} x,z
      ! j : dirección de la fuente {1;2} x,z
      deltaij = real(0.,8)
      if (i .eq. j) then
        deltaij = real(1.,8)
      end if
      end function deltaij
!#< r IBEM - termIndep !#>
      subroutine fill_termindep(auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres,trac0vec,n_con_sub,n_top_sub
      use sourceVars, only: Po,iFte=>currentiFte!use sourceVars, only: tipofuente
      implicit none
      interface
         include 'interfazFF.f'
      end interface
      
      type(Punto), pointer :: p_X,pXi
      integer :: mecS,mecE
      integer, intent(in) :: dir_j
      complex*16, dimension(1,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      real*8 :: nf(3)
      complex*16 :: TractionPSV, TractionSH
      
      ! aquí siempre la fuente es la fuente real.
      ! esta función se ejecuta se ejecuta tantas veces como las necesarias para llenar trac0vec
      
      !print*,p_x%boundaryIndex," # ",p_x%center,p_x%normal,p_x%length
      
      nf(1) = pXi%normal%x
      nf(2) = 1!pXi%normal%y
      nf(3) = pXi%normal%z
      ! si es onda plana no importa nf conque sea unitario.
      if (iFte .eq. 0) stop "fill_termindep iFte=0"
      if (Po(iFte)%tipofuente .eq. 1) nf(1:3) = 1 ! fuente real y onda plana
      
      if (pXi%region .eq. 2) then 
        auxK(1,mecS:mecE) = 0
        if (p_x%boundaryIndex .le. n_top_sub) return
        nf(1:3) = -nf(1:3)
      end if !
      
      if (dir_j .eq. 2) then ! SH
        call FFsh(0,FF,p_X,pXi,cOME,1,3)
!     print*,"------------- fill_termindep"
!     print*,p_X%boundaryindex,p_X%center
!     print*,pXi%center
!     print*,auxk(1,1),FF%V
!     print*,""
      !  | Ty |
        trac0vec(p_x%boundaryIndex) = & 
        trac0vec(p_x%boundaryIndex) - & 
          (TractionSH(auxk(1,2),auxk(1,3),p_x%normal) + FF%Ty) * nf(dir_j)
           
        if (p_X%tipoFrontera .eq.1) then !los desplazamientos
      !  |   V   |
        trac0vec(p_x%boundaryIndex+ n_con_sub) = & 
        trac0vec(p_x%boundaryIndex+ n_con_sub) - & 
          (auxk(1,1) + FF%V)* nf(dir_j)
        end if
      else !PSV 
      call FFpsv(0,FF,dir_j,p_X,pXi,cOME,1,5)
!     print*,"------------- fill_termindep"
!     print*,p_X%boundaryindex,p_X%center,p_X%normal
!     print*,p_x%length,p_x%cost,p_x%sint
!     print*,TractionPSV(auxk(1,3:5), p_x%normal,0),FF%Tx
!     print*,TractionPSV(auxk(1,3:5), p_x%normal,1),FF%Tz
!     print*,auxk(1,1),FF%W
!     print*,auxk(1,1),FF%U
!     print*,""
      !  | Tx |
      !  | Tz |
        trac0vec(p_x%boundaryIndex *2 - (1 - 0)) = &
        trac0vec(p_x%boundaryIndex *2 - (1 - 0)) - (&
          (TractionPSV(auxk(1,3:5), p_x%normal,0) + FF%Tx) * nf(dir_j)) !ok
          
        trac0vec(p_x%boundaryIndex *2 - (1 - 1)) = &
        trac0vec(p_x%boundaryIndex *2 - (1 - 1)) - (&
          (TractionPSV(auxk(1,3:5), p_x%normal,1) + FF%Tz) * nf(dir_j)) !ok
          
       if (p_X%tipoFrontera .eq.1) then !los desplazamientos
      !  |   W   |
      !  |   U   |
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) - &
          (auxK(1,1) + FF%W)* nf(dir_j)
                   
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) - &
          (auxK(1,2) + FF%U)* nf(dir_j)
       end if
      end if !dir_j
      end subroutine fill_termindep
      
      subroutine termIndepR(dir_j,cOME)
      use sourceVars, only : Po,iFte=>currentiFte
      use resultvars, only : trac0vec,Punto,FFres,& 
            boupoints,n_con_sub,n_top_sub,n_val_sub
      implicit none
      interface
        include 'interfazFF.f'
      end interface
      integer, intent(in) :: dir_j
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      type(Punto), pointer :: p_X
      real*8 :: nf(3)
      integer :: ip_X
      
!     do iPo =1,nFuentes
!     print*,",ipo",ipo
      nf(1) = Po(iFte)%normal%x
      nf(2) = 1!Po(iFte)%normal%y
      nf(3) = Po(iFte)%normal%z
!     print*,"nf=",nf
      if (iFte .eq. 0)  stop "termIndR iFte=0"
      if (Po(iFte)%tipofuente .eq. 1) nf(1:3) = 1 ! fuente real y onda plana
      
      if (dir_j .eq. 2) then ! SH
        stop "termIndepR SH _ falta line4637"
      else !PSV
        do ip_X = n_top_sub + 1, n_top_sub + n_con_sub + n_val_sub
          nullify(p_X)
          p_X => boupoints(ip_X)
!         print*,ip_x,p_x%center,p_x%tipoFrontera
          call FFpsv(0,FF,dir_j,p_X,Po(iFte),cOME,3,5)
!         print*,FF
      !  | Tx |
      !  | Tz |
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) - (FF%Tx * nf(dir_j)) 
          
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) - (FF%Tz * nf(dir_j))
        end do
      end if !dir_j
!     PRINT*,trac0vec
      end subroutine termIndepR
!#< r IBEM - reg E        !#>
      subroutine fill_ibemMat(i_zF,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres,ibemMat,n_con_sub
      use glovars, only : z0,UR
!     use debugstuff
      implicit none
      interface
        include 'interfazFF.f'
      end interface
      
      type(Punto), pointer :: pXi,p_X
      integer :: mecS,mecE,i_zF,dj
      integer, intent(in) :: dir_j
      complex*16, dimension(1,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      complex*16 :: TractionPSV, TractionSH
      if (p_X%tipoFrontera .gt. 1) return !o sea =2: (tracciones libres desde inclusión)
!     print*,"ibemmat",auxK(1,1:5)
!     call showMNmatrixZ(size(ibemMat,1),size(ibemMat,2), ibemMat ," mat ",6)
      if (dir_j .eq. 2) then ! SH
        call FFsh(i_zF,FF,p_X,pXi,cOME,1,3)
      !  |  Tyy  |
        if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
           ibemMat(p_x%boundaryIndex,pXi%boundaryIndex) = 0.5_8*UR
        else
           ibemMat(p_x%boundaryIndex,pXi%boundaryIndex) = &
              (TractionSH(auxk(1,2),auxk(1,3),p_x%normal) + FF%Ty)
        end if
           ! y si es un p. de coloc. en frontera con continuidad
        if (p_X%tipoFrontera .eq.1) then 
      !  |   V   |
           ibemMat(p_x%boundaryIndex + n_con_sub,pXi%boundaryIndex) = &
              (auxk(1,1) + FF%V)
        end if
      else !PSV -------------------------------------------------------------
      !  |  Txx Txz  |   
      !  |  Tzx Tzz  |
!       call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,1,5)
        dj = dir_j
        if (dj .eq. 3) dj = 2
        
        if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
          if (dir_j .eq. 1) then
            ibemMat(p_x%boundaryIndex *2 - 1, & 
                    pXi%boundaryIndex *2 - 1) = 0.5_8*UR
            ibemMat(p_x%boundaryIndex *2    , & 
                    pXi%boundaryIndex *2 - 1) = z0
          else if(dir_j .eq. 3) then
            ibemMat(p_x%boundaryIndex *2 - 1, & 
                    pXi%boundaryIndex *2 ) = z0
            ibemMat(p_x%boundaryIndex *2 , & 
                    pXi%boundaryIndex *2 ) = 0.5_8*UR
          end if!
        else
            if (i_zF .le. 0) then 
            print*,i_zF;print*,p_X%center;print*,pXi%center
            stop "fill_ibemMat: (i_zF .le. 0)"; end if
          call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,3,5)
!         !#< r 
!         auxk(1,1:5) = 0
!         print*,"lin 5143"
!         FF%Tx = 0;   FF%Tz = 0
!         stop
!         !#>
          
          ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                  pXi%boundaryIndex *2 -(2 - dj)) = &
                  (TractionPSV(auxK(1,3:5),p_x%normal,0) + FF%Tx)
          ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                  pXi%boundaryIndex *2 -(2 - dj)) = &
                  (TractionPSV(auxK(1,3:5),p_x%normal,1) + FF%Tz)
        end if
           ! y si es un p. de coloc. en frontera con continuidad
        if (p_X% tipoFrontera .eq.1) then 
      !  |  Wx   Wz  |
      !  |  Ux   Uz  |
          call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,1,2)
          ibemMat((p_x%boundaryIndex *2 -(1 - 0)) + 2* n_con_sub, & 
                   pXi%boundaryIndex *2 -(2 - dj)) = (auxK(1,1) + FF%W)
                       
          ibemMat((p_x%boundaryIndex *2 -(1 - 1)) + 2* n_con_sub, & 
                   pXi%boundaryIndex *2 -(2 - dj)) = (auxK(1,2) + FF%U)
        end if
      end if !dir_j
!     if (maxval(abs(ibemMat)) .gt. 100000) then
!     print*,"fill_ibemMat ######################"
!     print*,dir_j,i_zF,p_x%boundaryIndex,pXi%boundaryIndex
!     print*,FF%Tx,FF%Tz,FF%W,FF%U
!     print*,auxk(1,1)
!     print*,auxk(1,2)
!     print*,auxk(1,3)
!     print*,auxk(1,4)
!     print*,auxk(1,5)
!     call showMNmatrixZ(size(ibemMat,1),size(ibemMat,2), ibemMat ," mat ",6)
!     stop "fill_ibemMat"
!     end if
      end subroutine fill_ibemMat
            
      subroutine fill_diffbyStrata(i_zf,J,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres, nIpts
      use waveNumVars, only : NMAX,k_vec,dk,vecNK,SpliK,ome
      use meshvars, only: npixX,MeshMaxX,MeshMinX 
      use wavelets !fork
      use soilvars, only:alfa0,beta0,alfa,beta,N
      use sourceVars, only: PoFte =>  Po, iFte => currentiFte!tipofuente, PW_pol
      use glovars, only : UR,UI,PWfrecReal
      use peli, only : fotogramas_Region
      implicit none
      interface   
       include 'interfazFF.f'
      end interface 
      type(Punto), pointer :: p_X,pXi,p_Xmov
      type(Punto),target :: p_xaux
      integer :: i,ik,i_zf,iMec,mecS,mecE,J,mecaElemEnd,po,ne
      integer, intent(in) :: dir_j
      complex*16, dimension(2*nmax,mecS:mecE), target :: auxK
      complex*16, dimension(2*nmax,mecS:mecE) :: auxKmo
      real*8 :: mov_x
      complex*16 :: c,kx!,omega,bet,alf
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      real*8 :: nf(3)
      logical :: PW
      nf(1) = pXi%normal%x
      nf(2) = 1!pXi%normal%y
      nf(3) = pXi%normal%z
!     print*,i_zf,dir_j
      mecaElemEnd = 2 !PSV
      if (dir_j .eq. 2) mecaElemEnd = 1 !SH
      if (i_zF .eq. 0) then
        if ((i_zF .eq. 0) .and. (iFte .eq. 0)) stop "fill_diffbyStra iFte=0"
        if (PoFte(iFte)%tipofuente .eq. 1) nf(1:3) = 1.0 !fuente real y onda plana
        if (abs(nf(dir_j)) .lt. 0.0001) return !fuente real fuerza; componente nulo
      end if
      ! diffracted field due to the stratification U,V,W or G
      if (p_x%guardarMovieSiblings) then
               p_xaux%center%x = p_x%center%x
               p_xaux%center%z = p_x%center%z
               p_xaux%normal = p_x%normal
               p_xaux%isOnInterface = p_x%isOnInterface
               p_xaux%layer = p_x%layer
               p_Xmov => p_xaux
        po = min(int(vecNK(J)*SpliK),nmax); ne = 2*nmax-(po-2)
        do i = 1,npixX ! para cada hermanito
          mov_x = MeshMinX + (MeshMaxX - MeshMinX)/(npixX-1) * (i-1)! la coordenada x
          p_Xmov%center%x = mov_x
          
          ! si no es de la region 0 no vale la pena calcularlo
          if(fotogramas_Region(p_x%pointIndex-nIpts,i) .ne. 1) then !'estr'
!           print*,"cycled [",p_Xmov%center%x,",",p_xaux%center%z,"]"
            cycle
          end if
         
         ! agregar fase horizontal
      PW = .false.
      if (i_zF .eq. 0) then
       if (PoFte(iFte)%tipofuente .eq. 1) then ! onda plana !
         PW = .true.
       end if
      end if!
!      if ((i_zF .eq. 0) .and. (PoFte(iFte)%tipofuente .eq. 1)) then ! onda plana !
      if (PW) then
!      if (p_x%isOnInterface .eqv. .true.) return  ! creo                   !
       if ((PoFte(iFte)%PW_pol .eq. 1) .or. (PoFte(iFte)%PW_pol .eq. 3)) then
         if (PWfrecReal) then
           c = UR*beta0(N+1) !SV
         else
           c = beta(N+1) !SV
         end if
       elseif (PoFte(iFte)%PW_pol .eq. 2) then
         if (PWfrecReal) then
           c = UR*alfa0(N+1) !SV
         else
           c = alfa(N+1) !SV
         end if
       end if!
       if (PWfrecReal) then
         kx = UR*real(ome/c * sin(PoFte(iFte)%gamma))
       else
         kx = UR*real(come/c * sin(PoFte(iFte)%gamma))
       end if 
        
        ! aplicar la fase
        do imec = 1,mecaElemEnd
          auxkmo(1,imec) = auxk(1,imec) * & 
             exp(-UI*kx*(p_Xmov%center%x - PoFte(iFte)%center%x))
        end do !imec
      else !fuente cilndrica
        do imec = 1,mecaElemEnd
           do ik = 1,po+1!2*Nmax
            auxKmo(ik,imec) = auxk(ik,imec) * &
            exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_Xmov%center%x-pXi%center%x), 8))
           end do!  ik
           do ik = ne,2*Nmax
            auxKmo(ik,imec) = auxk(ik,imec) * &
            exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_Xmov%center%x-pXi%center%x), 8))
           end do!  ik
           
           !K->X
           !           auxKmo(:,imec) = FFTW(2*nmax,auxKmo(:,iMec),+1,dk)
            auxKmo(1,imec) = sum(auxKmo(1:po+1,imec)) + sum(auxKmo(ne:2*nmax,imec))
            auxKmo(1,imec) = auxKmo(1,imec) * dk
        end do !imec
      end if
      
      
      ! almacenar
          if (dir_j .eq. 2) then ! SH
          call FFsh(i_zf,FF,p_Xmov,pXi,cOME,1,1)
            if(i_zf .eq.0) then
              p_x%Wmov(J,3,i,iFte) = & 
              p_x%Wmov(J,3,i,iFte) + (auxKmo(1,1) + FF%V) * nf(dir_j) ! V
            else
              pXi%Gmov(p_x%pointIndex-nIpts,3,dir_j,i) = & 
              pXi%Gmov(p_x%pointIndex-nIpts,3,dir_j,i) + auxKmo(1,1) + FF%V ! V
            end if
          else !PSV
          call FFpsv(i_zF,FF,dir_j,p_Xmov,pXi,cOME,1,1)
            if(i_zf .eq.0) then
              p_x%Wmov(J,1,i,iFte) = p_x%Wmov(J,1,i,iFte) + (auxKmo(1,1) + FF%W) * nf(dir_j) !W
              p_x%Wmov(J,2,i,iFte) = p_x%Wmov(J,2,i,iFte) + (auxKmo(1,2) + FF%U) * nf(dir_j) !U
            else
              pXi%Gmov(p_x%pointIndex-nIpts,1,dir_j,i) = & 
              pXi%Gmov(p_x%pointIndex-nIpts,1,dir_j,i) + auxKmo(1,1) + FF%W ! W
              pXi%Gmov(p_x%pointIndex-nIpts,2,dir_j,i) = & 
              pXi%Gmov(p_x%pointIndex-nIpts,2,dir_j,i) + auxKmo(1,2) + FF%U ! U
            end if
          end if !dir_j
        end do ! i
      else !not a movie point .......................................................................
      if (dir_j .eq. 2) then ! SH
      if (pXi%region .ne. p_x%region) then 
      print*,"no on the same region", pXi%region,pXi%center," - ",p_X%region,p_X%center
      return
      end if
       call FFsh(i_zf,FF,p_X,pXi,cOME,1,3) !incidencia directa
       
!      !#< r
!      print*,p_X%center
!      print*,auxk(1,1) ,abs(auxk(1,1)),FF%V,abs(FF%V),abs(auxk(1,1)+FF%V)
!      FF%V = 0
!      auxk(1,1) = 0
!      !#>
       
       if(i_zf .eq.0) then
!         print*,p_x%resp(J,iFte)%V,auxk(1,1),FF%V,nf(dir_j);stop 4898
!         if (pXi%region .ne. p_x%region) return
          p_x%resp(J,iFte)%V = & 
          p_x%resp(J,iFte)%V + (auxk(1,1) + FF%V) !* nf(dir_j) ! V
          p_x%resp(J,iFte)%Ty = & 
          p_x%resp(J,iFte)%Ty + & 
          ((auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty) !* nf(dir_j)
          !                     s12                       s32
       else
          pXi%G(p_x%pointIndex,3,dir_j) = & 
          pXi%G(p_x%pointIndex,3,dir_j) + auxK(1,1) + FF%V ! V
          pXi%G(p_x%pointIndex,6,dir_j) = & 
          pXi%G(p_x%pointIndex,6,dir_j) + & 
          (auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty ! Ty
          !                     s12                       s32
       end if
      else !PSV
!     if (pXi%region .ne. p_x%region) return
      if (pXi%region .ne. p_x%region) then 
      print*,"no on the same region", pXi%region,pXi%center," - ",p_X%region,p_X%center
      return
      end if
       call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,1,5)  !incidencia directa       
!      print*,"x= ",p_x%center
!      print*,p_x%resp(J,iFte)%W
!      print*,auxK(1,1)
!      print*,FF%W
!      print*," "
       if(i_zf .eq. 0) then
!        if (pXi%region .ne. p_x%region) return
         p_x%resp(J,iFte)%W = p_x%resp(J,iFte)%W + (auxK(1,1) + FF%W) * nf(dir_j) !W
         p_x%resp(J,iFte)%U = p_x%resp(J,iFte)%U + (auxK(1,2) + FF%U) * nf(dir_j) !U
         p_x%resp(J,iFte)%Tz = p_x%resp(J,iFte)%Tz + &
         ((auxk(1,4)* p_x%normal%x + auxk(1,3)* p_x%normal%z) + FF%Tz) * nf(dir_j)
         !      s31                       s33
         p_x%resp(J,iFte)%Tx = p_x%resp(J,iFte)%Tx + &
         ((auxk(1,5)* p_x%normal%x + auxk(1,4)* p_x%normal%z) + FF%Tx) * nf(dir_j)
         !      s11                       s31
         p_x%resp(J,iFte)%sxx = p_x%resp(J,iFte)%sxx + & 
         (auxk(1,5) + FF%sxx) * nf(dir_j) !sxx
         p_x%resp(J,iFte)%szx = p_x%resp(J,iFte)%szx + & 
         (auxk(1,4) + FF%szx) * nf(dir_j) !szx
         p_x%resp(J,iFte)%szz = p_x%resp(J,iFte)%szz + & 
         (auxk(1,3) + FF%szz) * nf(dir_j) !szz
       else
         pXi%G(p_X%pointIndex,1,dir_j) = & 
         pXi%G(p_X%pointIndex,1,dir_j) + (auxK(1,1) + FF%W) !W
         pXi%G(p_X%pointIndex,2,dir_j) = & 
         pXi%G(p_X%pointIndex,2,dir_j) + (auxK(1,2) + FF%U) !U
         pXi%G(p_X%pointIndex,4,dir_j) = & 
         pXi%G(p_X%pointIndex,4,dir_j) + &
         ((auxk(1,4)* p_x%normal%x + auxk(1,3)* p_x%normal%z) + FF%Tz) !Tz
         pXi%G(p_X%pointIndex,5,dir_j) = & 
         pXi%G(p_X%pointIndex,5,dir_j) + &
         ((auxk(1,5)* p_x%normal%x + auxk(1,4)* p_x%normal%z) + FF%Tx) !Tx
         
         pXi%G(p_X%pointIndex,7,dir_j) = pXi%G(p_X%pointIndex,7,dir_j) + & 
         auxk(1,5) + FF%sxx !sxx
         pXi%G(p_X%pointIndex,8,dir_j) = pXi%G(p_X%pointIndex,8,dir_j) + & 
         auxk(1,4) + FF%szx !szx
         pXi%G(p_X%pointIndex,9,dir_j) = pXi%G(p_X%pointIndex,9,dir_j) + & 
         auxk(1,3) + FF%szz !szz
       end if
      end if !dir_j
      end if ! guardarMovieSiblings
      end subroutine fill_diffbyStrata
                

      function TractionPSV(RW,normal,l)
      use resultVars, only: Punto2d
      implicit none
      complex*16 :: TractionPSV
      complex*16, intent(in) :: RW(3:5)
      type(Punto2d), intent(in) :: normal
      integer, intent(in) :: l ! l=0-> Tx ; l=1-> Tz
      TractionPSV = RW(5-l) * normal%x + &
                    RW(4-l) * normal%z           
      end function TractionPSV        
      ! las tracciones:
      !   ,--- componente de la tracción : m
      !   |     ,--- cara
      !   |     |
      !  Tx = Sxx nx + Szx nz  |___ (0) fza real (fuente)  .
      !  Tz = Szx nx + Szz nz  |                           .  
      !
      !   ,--- componente de la tracción : l
      !   |,--- (1),(3) dirección de la fuerza : m
      !   ||    ,--- cara
      !   ||    |,--- fza
      !  Txx = Sxx1 nx + Szx1 nz  |___ (1) fza horizontal  .
      !  Tzx = Szx1 nx + Szz1 nz  |                        . 
      !  Txz = Sxx3 nx + Szx3 nz |___ (3) fza vertical     .
      !  Tzz = Szx3 nx + Szz3 nz |                         .
      
      !  T_lm = S_lkm * n_k
      !       = S_l1m * n_1 + S_l3m * n3
      !         __|__         __|__
      !      s11     s31   s13    s33
      !      s11     s31   s31    s33  (son equivalentes)
      !       5       4     4      3   (indice en RW(_,i) )
      !       0       1     0      1   (indice de submatriz: l )
      
      function TractionSH(s32,s12,normal)
      use resultVars, only: Punto2d
      implicit none
      complex*16 :: TractionSH
      complex*16, intent(in) :: s32,s12 
      type(Punto2d), intent(in) :: normal
      TractionSH = s12 * normal%x + &
                   s32 * normal%z
      end function TractionSH
                        
            
!#< r IBEM - reg R        !#>
      subroutine reffField_by_(ipXi,dir_j,cOME)
      ! llenamos una columna de la matriz del ibem en la región
      ! de condiciones de continuidad entre E y R y front libre en R (misma columna)
      use glovars, only : z0,UR
      use resultVars, only : Punto,ibemMat,FFres,n_top_sub,n_con_sub,n_val_sub, boupoints
      implicit none
      interface !porque recibe un puntero como argumento
         include 'interfazFF.f'
      end interface
      integer, intent(in) :: ipXi,dir_j ! ip_X : n_topo+1,n_topo+n_cont
      complex*16, intent(in),target  :: cOME
      type(Punto), pointer :: pXi,p_X
      integer :: ip_X,dj
      type(FFres),target :: FF
      
      nullify(pXi)
      pXi => boupoints(ipXi)
      ! la frontera de la regi´on R
      do ip_X = n_top_sub +1, n_top_sub + n_con_sub + n_val_sub
        nullify(p_X)
        p_X => boupoints(ip_X)
        
      if (dir_j .eq. 2) then ! SH
        call FFsh(-1,FF,p_X,pXi,cOME,1,3)
       if(p_X%tipoFrontera .lt. 2) then
      !  |  Tyy  |
        if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
           ibemMat(p_x%boundaryIndex, &
                   pXi%boundaryIndex + n_con_sub) = 0.5_8*UR
        else
           ibemMat(p_x%boundaryIndex, &
                   pXi%boundaryIndex + n_con_sub) = - FF%Ty
        end if
           ! y si es un p. de coloc. en frontera con continuidad
        if (p_X%tipoFrontera .eq.1) then 
      !  |   V   |
!       r = sqrt((p_x%center%x-pXi%center%x)**2. + (p_x%center%z-pXi%center%z)**2.)
!       FF%V = greenexG22(cOME/beta(N+2),r,pXi%length*2,N+2)
      
           ibemMat(p_x%boundaryIndex + n_con_sub, & 
                   pXi%boundaryIndex + n_con_sub) = - FF%V
        end if
           ! y si es un p. de coloc. en frontera libre
       else !p_X%tipoFrontera .eq. 2
      !  |  Tyy  |
           if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
               ibemMat(p_x%boundaryIndex + n_con_sub, &
                   pXi%boundaryIndex + n_con_sub) = - 0.5_8*UR
           else
               ibemMat(p_x%boundaryIndex + n_con_sub, & 
                   pXi%boundaryIndex + n_con_sub) = FF%Ty
           end if
       end if
      else !PSV **********************************************************************
        dj = dir_j
        call FFpsv(-1,FF,dir_j,p_X,pXi,cOME,1,5)
        if (dj .eq. 3) dj = 2
      if(p_X%tipoFrontera .lt. 2) then
      !  |  Txx Txz  |   
      !  |  Tzx Tzz  |
        if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
          if (dir_j .eq. 1) then
            ibemMat(p_x%boundaryIndex *2 - 1, & 
                   (pXi%boundaryIndex *2 - 1) + 2* n_con_sub) = 0.5_8*UR
            ibemMat(p_x%boundaryIndex *2    , & 
                   (pXi%boundaryIndex *2 - 1) + 2* n_con_sub) = z0
          else if(dir_j .eq. 3) then
            ibemMat(p_x%boundaryIndex *2 - 1, & 
                   (pXi%boundaryIndex *2 - 0) + 2* n_con_sub) = z0
            ibemMat(p_x%boundaryIndex *2 , & 
                   (pXi%boundaryIndex *2 - 0) + 2* n_con_sub) = 0.5_8*UR
           end if
        else
          ibemMat(p_x%boundaryIndex *2 -(1 - 0), & 
                 (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = - FF%Tx
          ibemMat(p_x%boundaryIndex *2 -(1 - 1), & 
                 (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = - FF%Tz
        end if
           ! y si es un p. de coloc. en frontera con continuidad
        if (p_X% tipoFrontera .eq.1) then
      !  |   W   |
      !  |   U   |
!      if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
!         call greenexPSV(FF,dir_j,p_X,pXi,cOME)
!      end if
          ibemMat((p_x%boundaryIndex *2 -(1 - 0)) + 2* n_con_sub, & 
                  (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = - FF%W
                       
          ibemMat((p_x%boundaryIndex *2 -(1 - 1)) + 2* n_con_sub, & 
                  (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = - FF%U
        end if
          ! y si es un p. de coloc. en frontera libre
      else !p_X%tipoFrontera .eq. 2
      !  |  Txx Txz  |   
      !  |  Tzx Tzz  |
           if (pXi%boundaryIndex .eq. p_x%boundaryIndex) then
               if (dir_j .eq. 1) then
            ibemMat((p_x%boundaryIndex *2 - 1)+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 - 1)+ 2* n_con_sub) = -0.5_8*UR
            ibemMat((p_x%boundaryIndex *2    )+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 - 1)+ 2* n_con_sub) = z0
               else if(dir_j .eq. 3) then
            ibemMat((p_x%boundaryIndex *2 - 1)+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 - 0)+ 2* n_con_sub) = z0
            ibemMat((p_x%boundaryIndex *2    )+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 - 0)+ 2* n_con_sub) = -0.5_8*UR
                end if
           else
            ibemMat((p_x%boundaryIndex *2 -(1 - 0) )+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = FF%Tx
            ibemMat((p_x%boundaryIndex *2 -(1 - 1) )+ 2* n_con_sub, & 
                    (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub) = FF%Tz
           end if
      end if
      end if !dir_j
      end do !ip_X
      
      end subroutine reffField_by_
      
      subroutine GreenReg_R(Jfrec,dir_j,cOME) !Región R
      use resultvars, only:Punto,FFres, nIpts, nBpts,allpoints, & 
                           boupoints, n_top_sub
      use glovars, only: makeVideo
      use peli, only : coords_Z,coords_X,fotogramas_Region
      use meshVars, only : npixX,npixZ
      use sourceVars, only : Po,nFuentes,currentiFte
      use soilvars, only : N
      implicit none
      interface
         include 'interfazFF.f'
      end interface
      integer, intent(in) :: Jfrec,dir_j
      complex*16, intent(in)  :: cOME
      integer :: iP_x,iPXi,iFte
      type(Punto), pointer :: p_X,pXi,p_Xmov
      integer :: i,j,iPo
      type(Punto),target :: p_xaux
      type(FFres),target :: FF
      real*8 :: nf(3)
      ! ciclar en todos los puntos y averiguar se le corresponde
      do iP_x = 1, nIpts ! (todos los puntos receptores menos las fronteras de integración)
        p_X => allpoints(iP_x)
        if (p_x%isboundary) stop "GreenReg_R: p_x is boundary "
        if (p_x%region .ne. 2) cycle !'incl'
        
        ! para todas las fuentes en la región R
      do iPXi = n_top_sub+1,nBpts ! las fuentes virtual en la región R
        pXi => boupoints(iPXi)
        if (dir_j .eq. 2) then ! SH
         call FFsh(-1,FF,p_X,pXi,cOME,1,3)
          pXi%G(p_x%pointIndex,3,dir_j) = FF%V ! V
          pXi%G(p_x%pointIndex,6,dir_j) = FF%Ty ! Ty
        else !PSV
         call FFpsv(-1,FF,dir_j,p_X,pXi,cOME,1,5) 
          pXi%G(p_X%pointIndex,1,dir_j) = FF%W !W 
          pXi%G(p_X%pointIndex,2,dir_j) = FF%U !U 
          pXi%G(p_X%pointIndex,4,dir_j) = FF%Tz !Tz
          pXi%G(p_X%pointIndex,5,dir_j) = FF%Tx !Tx
          pXi%G(p_X%pointIndex,7,dir_j) = FF%sxx !sxx
          pXi%G(p_X%pointIndex,8,dir_j) = FF%szx !szx
          pXi%G(p_X%pointIndex,9,dir_j) = FF%szz !szz
        end if !dir_j
      end do !iPXi
      
      ! y la fuente real (para cada una de ellas)
      do iPo = 1,nfuentes
      if (Po(iPo)%region .eq. 2) then ! la fuente real en la región R
      pXi => Po(iPo)
      nf(1) = pXi%normal%x
      nf(2) = 1!pXi%normal%y
      nf(3) = pXi%normal%z
      iFte = iPo; currentiFte = iFte
!     print*,"iFte=",iFte
      if (dir_j .eq. 2) then ! SH
       call FFsh(0,FF,p_X,pXi,cOME,1,3)
          p_x%resp(Jfrec,iFte)%V = p_x%resp(Jfrec,iFte)%V + (FF%V) * nf(dir_j) ! V
          p_x%resp(Jfrec,iFte)%Ty = p_x%resp(Jfrec,iFte)%Ty + (FF%Ty) * nf(dir_j)
      else !PSV
       call FFpsv(0,FF,dir_j,p_X,pXi,cOME,1,5) 
!        if (pXi%region .ne. p_x%region) cycle
         p_x%resp(Jfrec,iFte)%W = p_x%resp(Jfrec,iFte)%W + (FF%W) * nf(dir_j) !W
         p_x%resp(Jfrec,iFte)%U = p_x%resp(Jfrec,iFte)%U + (FF%U) * nf(dir_j) !U
         p_x%resp(Jfrec,iFte)%Tz = p_x%resp(Jfrec,iFte)%Tz + (FF%Tz) * nf(dir_j) !Tx
         p_x%resp(Jfrec,iFte)%Tx = p_x%resp(Jfrec,iFte)%Tx + (FF%Tx) * nf(dir_j) !Tx
         p_x%resp(Jfrec,iFte)%sxx = p_x%resp(Jfrec,iFte)%sxx + (FF%sxx) * nf(dir_j) !sxx
         p_x%resp(Jfrec,iFte)%szx = p_x%resp(Jfrec,iFte)%szx + (FF%szx) * nf(dir_j) !szx
         p_x%resp(Jfrec,iFte)%szz = p_x%resp(Jfrec,iFte)%szz + (FF%szz) * nf(dir_j) !szz
      end if !dir_j
      end if !Po 2
      end do !iPo
      currentiFte = 0
      end do !iP_x
      
      if (makeVideo) then
      do i=1,npixZ
        p_X => allpoints(nIpts+i)
        do j=1,npixX ! para cada hermanito
          if (fotogramas_Region(i,j) .eq. 2) then !'incl'
            p_xaux%center%x = coords_X(j)
            p_xaux%center%z = coords_Z(i)
!     print*,""
!     print*,"this is mov incl",j,i," (",p_xaux%center%x,",",p_xaux%center%z,")"
!           p_xaux%normal%x = 1.0_8 !
!           p_xaux%normal%y = 1.0_8 !  no se usan 
!           p_xaux%normal%z = 1.0_8 !
            p_xaux%region = 2; p_xaux%layer = N+2 !'incl'
            p_xaux%boundaryindex = -10000
            p_xaux%isOD = .false.
            p_Xmov => p_xaux
               
            ! para todas las fuentes virtuales en la frontera de R
            do iPXi = n_top_sub +1,nBpts
              pXi => boupoints(iPXi)
              if (dir_j .eq. 2) then ! SH
                call FFsh(-1,FF,p_Xmov,pXi,cOME,1,1)
                pXi%Gmov(i,3,dir_j,j) = FF%V ! V
!               print*,ipXi,i+nIpts,j,"=",FF%V
              else !PSV
                call FFpsv(-1,FF,dir_j,p_Xmov,pXi,cOME,1,1)
!               print*,ipXi,j,i,"=",FF%W,FF%U
                pXi%Gmov(i,1,dir_j,j) = FF%W ! W
                pXi%Gmov(i,2,dir_j,j) = FF%U ! U
              end if !dir_j
            end do !ipXi   
            
            ! y la fuente real
            do iPo = 1,nfuentes
              if (Po(iPo)%region .eq. 2) then
                pXi => Po(iPo)
                nf(1) = pXi%normal%x
                nf(2) = 1!pXi%normal%y
                nf(3) = pXi%normal%z
                iFte = iPo; currentiFte = iFte
                if (abs(nf(dir_j)) .lt. 0.0001) cycle
                if (dir_j .eq. 2) then ! SH
                  call FFsh(0,FF,p_Xmov,pXi,cOME,1,1)
         p_Xmov%Wmov(Jfrec,3,j,iFte) = p_Xmov%Wmov(Jfrec,3,j,iFte) + FF%V * nf(dir_j) !V
                else !PSV
                  call FFpsv(0,FF,dir_j,p_Xmov,pXi,cOME,1,1) 
         p_x%Wmov(Jfrec,1,j,iFte) = p_x%Wmov(Jfrec,1,j,iFte) + FF%W * nf(dir_j) !W
         p_x%Wmov(Jfrec,2,j,iFte) = p_x%Wmov(Jfrec,2,j,iFte) + FF%U * nf(dir_j) !U
                end if !dir_j
              end if ! Po 2
            end do ! iPo
            currentiFte = 0
          end if ! si 
        end do !j
      end do !i
      end if ! makeVideo
      end subroutine GreenReg_R
      
!#< r IBEM - overDet        !#> 
      subroutine overDetermineSystem(J,ren,PSV)
      use resultvars, only: nIpts,& 
      ibemMat,trac0vec,allpoints,boupoints,Punto,FFres,& 
      n_top_sub,n_con_sub,n_val_sub
      use waveNumVars, only : cOME
      use sourceVars, only : iFte => currentiFte
!     use debugStuff
      use glovars, only : z0!,rutaOut
      implicit none
      interface
      include 'interfazFF.f'
      end interface
      integer, intent(in) :: J
      integer :: ren,renStep,col
      logical, intent(in) :: PSV
      type(Punto), pointer :: pXi,p_X
      integer :: iP_x,iPxi,ipxi_I,ipxi_F,dj,dir_j!,iPhi_I,iPhi_F
      type(FFres),target :: FF
!     CHARACTER(len=32) :: arg
      
!     print*,"ren=",ren !primer renglon del primer elemento de sobredeterminad
!     call showMNmatrixZ(size(ibemMat,1),size(ibemMat,2), ibemMat ," mat ",6)
!     call showMNmatrixZ(size(trac0vec,1),1 , trac0vec,"  b  ",6)
!     call chdir(trim(adjustl(rutaOut))) 
!      open(421,FILE= "outA.m",action="write",status="replace")
!      write(arg,'(a)') "Bi"
!      call scripToMatlabMNmatrixZ(size(trac0vec,1),1,trac0vec,arg,421)
!      write(arg,'(a)') "Mi"
!      call scripToMatlabMNmatrixZ(size(ibemMat,1),size(ibemMat,2),ibemMat,arg,421)
!      close(421)
!      CALL chdir("..")
      if (iFte .eq. 0) stop "overDertSys iFte=0"
      
      do iP_x = 1,nIpts  !cada receptor X
      p_X => allpoints(iP_x) 
!     col = 1
      if (.not. p_X%isOD) cycle !; print*,"iP_X = ",iP_x,"  iP_XCen=",p_X%center
      if (PSV) then
       if (p_X% tipoFrontera .eq.0) then!#< r frontera tipo trontera libre en medio estrat !#>
          ipxi_I = 1
          ipxi_F = n_top_sub
          do iPxi = ipxi_I,ipxi_F!; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi)!; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2!; print*,"dir_j ",dir_j," dj ",dj
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
            col = pXi%boundaryIndex *2 -(2 - dj)
!           print*,"ibemMat",pXi%G(p_X%pointIndex,5,dir_j),pXi%G(p_X%pointIndex,4,dir_j)
              ibemMat(ren,   col) = pXi%G(p_X%pointIndex,5,dir_j) !Tx
              ibemMat(ren+1, col) = pXi%G(p_X%pointIndex,4,dir_j) !Tz
            end do !dj
          end do !iPxi
            !  | Tx |
            !  | Tz |
!           print*,"tra0",p_x%resp(J,iFte)%Tx,p_x%resp(J,iFte)%Tz
              trac0vec(ren  ) = p_x%resp(J,iFte)%Tx  !0
              trac0vec(ren+1) = p_x%resp(J,iFte)%Tz  !0
          renStep = 2
       else if (p_X%tipoFrontera .eq.1) then !#< r  frontera tipo continuidad !#>
          ipxi_I = 1
          ipxi_F = n_top_sub + n_con_sub
          do iPxi = ipxi_I,ipxi_F!; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi)!; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2!; print*,"dir_j ",dir_j," dj ",dj
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
            col = pXi%boundaryIndex *2 -(2 - dj)
              ibemMat(ren,   col) = pXi%G(p_X%pointIndex,5,dir_j) !Tx
              ibemMat(ren+1, col) = pXi%G(p_X%pointIndex,4,dir_j) !Tz
            !  |  Wx   Wz  |
            !  |  Ux   Uz  |
              ibemMat(ren+2, col) = pXi%G(p_X%pointIndex,1,dir_j) !W
              ibemMat(ren+3, col) = pXi%G(p_X%pointIndex,2,dir_j) !U
            end do !dj
          end do !iPxi
          ipxi_I = n_top_sub + 1
          ipxi_F = n_top_sub + n_con_sub + n_val_sub
          do iPxi = ipxi_I,ipxi_F!; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi)!; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2!; print*,"dir_j ",dir_j," dj ",dj
            call FFpsv(-1,FF,dir_j,p_X,pXi,cOME,1,5)
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
            col = (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub
              ibemMat(ren,   col) = -FF%Tx !Tx
              ibemMat(ren+1, col) = -FF%Tz !Tz
            !  |  Wx   Wz  |
            !  |  Ux   Uz  |
              ibemMat(ren+2, col) = -FF%W !W
              ibemMat(ren+3, col) = -FF%U !U
            end do !dj
          end do !iPxi
            !  | Tx |
            !  | Tz |
            !  | W  |
            !  | U  |
              trac0vec(ren  ) = p_x%resp(J,iFte)%Tx
              trac0vec(ren+1) = p_x%resp(J,iFte)%Tz
              trac0vec(ren+2) = p_x%resp(J,iFte)%W
              trac0vec(ren+3) = p_x%resp(J,iFte)%U
          renStep = 4
       else if (p_X%tipoFrontera .eq.2) then !#< r  frontera tipo frontera libre en inclusión elástica !#>
          ! los segmentos en la frontera de continuidad y en
          ! la frontera libre en la inclusión
          ipxi_I = n_top_sub + 1
          ipxi_F = n_top_sub + n_con_sub + n_val_sub
!         print*,"ipxi_I=",ipxi_I,"ipxi_F=",ipxi_F
          do iPxi = ipxi_I,ipxi_F !; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi) !; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2 !; print*,"dir_j ",dir_j," dj ",dj
            call FFpsv(-1,FF,dir_j,p_X,pXi,cOME,3,5)
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
!           col = pXi%boundaryIndex *2 -(2 - dj)+ 2* n_con_sub + n_top_sub 
            col = (pXi%boundaryIndex *2 -(2 - dj))+ 2* n_con_sub
!           print*,"col=",col,"  ren=",ren
!           print*,"| Tx |",FF%Tx
!           print*,"| Tz |",FF%Tz
            
              ibemMat(ren,   col) =  FF%Tx !Tx
              ibemMat(ren+1, col) =  FF%Tz !Tz
            end do !dj
          end do !iPxi
!             print*,"trac0vec(",ren,ren+1,")"
!             print*,"| Tx |",p_x%resp(J,iFte)%Tx
!             print*,"| Tz |",p_x%resp(J,iFte)%Tz
              trac0vec(ren  ) = z0!p_x%resp(J,iFte)%Tx
              trac0vec(ren+1) = z0!p_x%resp(J,iFte)%Tz
          renStep = 2
       end if !%tipoFrontera
      else !PSV / SH
         stop "overDetermineSystem: Falta SH sobreDet"
      end if !PSV/SH
      ren = ren + renStep
      end do !iP_x
      
!      write(arg,'(a)') "BfGx"
!      call scripToMatlabMNmatrixZ(size(trac0vec,1),1,trac0vec,arg,421)
!      write(arg,'(a)') "MfGx"
!      call scripToMatlabMNmatrixZ(size(ibemMat,1),size(ibemMat,2),ibemMat,arg,421)
!      close(421)
!      CALL chdir("..")
!     call showMNmatrixZ(size(ibemMat,1),size(ibemMat,2), ibemMat ," mat ",6)
!     call showMNmatrixZ(size(trac0vec,1),1 , trac0vec,"  b  ",6)
!     stop "overDetermineSystem"
      end subroutine overDetermineSystem
      
      subroutine solveOverDetermineSystem(M,N)
      use resultvars, only: A => ibemMat, B => trac0vec
      use debugStuff
      use glovars, only : v => verbose
      implicit none
      integer :: M,N,NRHS,LDA,LDB,LWORK,LWMAX,INFO
      complex*16, dimension(:), allocatable :: WORK
      complex*16, dimension(:,:), allocatable :: copyA
      CHARACTER(len=32) :: arg
      EXTERNAL  :: ZGELS
      INTRINSIC :: INT,MIN
!     real*8,dimension(:),allocatable :: S,Rwork
!     real*8 :: Rcond
      NRHS = 1
      LWMAX = 500
      LDA = M
      LDB = M
      allocate(WORK(LWMAX))
      
      if (v .gt. 2) then !#< b
      write(arg,*) "  A  "
      call showMNmatrixZ(M,N   , A,arg,6)
!     call scripToMatlabMNmatrixZ(M,N,A(1:M,1:N),arg,6)
      write(arg,*) "  B  "
      call showMNmatrixZ(M,NRHS, B,"  b  ",6) 
!     call scripToMatlabMNmatrixZ(M,1,B(1:M),arg,6)
      allocate(copyA(size(A,1),size(A,2)))
      copyA = A; end if!#>
      ! Using QR ----------------------------------------
      
      ! Query the optimal workspace.............
      LWORK = -1                               !
      call ZGELS('No transpose', M, N, NRHS, & !
                 A, LDA,                     & !
                 B, M, WORK,                 & !
                 LWORK,INFO )                  !
      LWORK = int(WORK(1))                     !
      deallocate(WORK); allocate(WORK(LWORK))  !
      !········································· 
      WORK = 0 !#< b
!     print*,"M",M
!     print*,"M",N
!     print*,"NRHS",NRHS
!     print*,"LDA",LDA
!     print*,"LDB",LDB
!     print*,'LWORK=',LWORK !#>
      
      call ZGELS('No transpose', M, N, NRHS, & 
                 A, LDA,                     & 
                 B, M, WORK,                 & 
                 LWORK,INFO ) 
      
      !Check for the full rank.
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
         WRITE(*,*)'factor of A is zero, so that A does not have full '
         WRITE(*,*)'rank; the minimum norm solution could not be '
         WRITE(*,*)'computed.'
         STOP "solveOverDetermineSystem"
      END IF
      
      ! Using single value decomposition
!     allocate(Rwork(5*MN))
!     allocate(S(MN))
!     Rcond = dble(-1) !for machine precision
!     Rcond = 1E-8
!     call ZGELSS(Mi,Ni,1,ibemMat,Mi,trac0vec,Mi,& 
!                 S, RCOND, RANK,WORK, LWORK, RWORK, INFO )
!     print*,"Rank",Rank;print*,"S",S
      
      if (v .gt. 2) then  !#< b
      write(arg,*) "phi  "
      call showMNmatrixZ(N,NRHS,B(1:N),arg,6)  
      call scripToMatlabMNmatrixZ(N,NRHS,B(1:N),arg,6)
      call showMNmatrixZ(M,NRHS, matmul(copyA(1:M,1:N),B(1:N)),"  bBA",6)
      end if!#>
      end subroutine solveOverDetermineSystem

      subroutine crepa_four_fotogramas
      use glovars, only : saveG,verbose
      use peli, only : fotogramas 
      use sourceVars, only : SH,PSV,ifTe => currentiFte
      use waveNumVars, only : t_vec,omei, nF => NFREC, nT => NPTSTIME
      use meshVars, only : npixX,npixZ
      use waveVars, only : Uo,Dt
      use wavelets
      use dislin
      implicit none
      integer :: i,ix,iz,imec,mecS,mecE
      complex*16, dimension(:), pointer :: p_fot
      character(LEN=100) :: nam
!     real*8 :: factor
      if (verbose .ge. 1) print*,"frec -> time"
      if (iFte .eq. 0) stop "crepa_four_fotogr iFte=0"
!     factor = sqrt(1.0*nT)
      mecS = 3; mecE = 2; if (PSV) mecS =1; if (SH) mecE =3
      if (saveG .eqv. .true.) then
       write(nam,'(a,I0,a)') "G",iFte,".bin"
       OPEN(6374,FILE=trim(nam),STATUS='UNKNOWN', ACCESS='STREAM',ACTION='WRITE')
       write(6374) mecS
       write(6374) mecE
       write(6374) npixZ
       write(6374) npixX
       write(6374) nT
      end if
      do imec = mecS,mecE
      do iz=1,npixZ
        do ix=1,npixX
          fotogramas(iz,ix,nT-nF+2:nT,imec,iFte) = conjg(fotogramas(iz,ix,nF:2:-1,imec,iFte))
          p_fot => fotogramas(iz,ix,1:nT,imec,iFte)
          p_fot = p_fot * t_vec
          if (saveG .eqv. .true.) then 
           do i=1,nT
            write(6374) p_fot(i)
           end do
          end if
          p_fot = p_fot * Uo(:,iFte)
        ! al tiempo:
          p_fot = FFTW(nT,p_fot,+1,1/(nT*dt))
        ! remover efecto de la frecuencia imaginaria
          p_fot = p_fot * & 
          exp(-OMEI * Dt*((/(i,i=0,nT-1)/)))
!         if (verbose .ge. 3) then
!           CALL chdir("perPixelTraces")
!           write(titleN,"(i0,a,i0,a)") ix,"_",iz,".pdf"
!           CALL SETFIL(trim(titleN))
!           call qplot(real((/((i-1)*Dt,i=1,800)/),4),real(p_fot(1:800),4),800)
!           CALL chdir("..")
!         end if
        end do !ix
      end do !iz
      end do !imec
      if (saveG .eqv. .true.) THEN 
      close(6374)
      print*,"saved G function"
!     CALL chdir("..")
      end if
      end subroutine crepa_four_fotogramas
      
      subroutine loadG_fotogramas
      use glovars, only:verbose,rutaOut
      use peli, only : fotogramas 
      use waveNumVars, only : omei, nT => NPTSTIME!,dfrec
      use meshVars, only : npixX,npixZ
      use waveVars, only : Uo,Dt
!     use soilVars, only : Qq
      use wavelets
      use sourceVars, only : nFuentes 
      use dislin
      implicit none
      integer :: i,ix,iz,imec,mecS,mecE,iFte
      complex*16, dimension(:), pointer :: p_fot
      character(LEN=100) :: titleN,nam
      call chdir(trim(adjustl(rutaOut)))
      do iFte = 1,nFuentes
      write(nam,'(a,I0)') "video",iFte
      CALL chdir(trim(nam))
       write(nam,'(a,I0,a)') "G",iFte,".bin"
       OPEN(6375,FILE=trim(nam),STATUS='UNKNOWN', ACCESS='STREAM',ACTION='READ')
       read(6375) mecS
       read(6375) mecE
       read(6375) npixZ
       read(6375) npixX
       read(6375) nT
      do imec = mecS,mecE
      do iz=1,npixZ
      do ix=1,npixX
      do i=1,nT      
         read(6375) fotogramas(iz,ix,i,imec,iFte)
      end do
          p_fot => fotogramas(iz,ix,1:nT,imec,iFte)
          p_fot = p_fot * Uo(:,iFte)
        ! al tiempo:
          p_fot = FFTW(nT,p_fot,+1,1/(nT*dt))
        ! remover efecto de la frecuencia imaginaria
          p_fot = p_fot * & 
          exp(-OMEI * Dt*((/(i,i=0,nT-1)/)))
        ! remover efecto de la velocidad imaginaria ?
          !p_fot = p_fot * & 
          !exp((1./2./Qq) * Dt*((/(i,i=0, nT-1)/))) 
          if (verbose .ge. 3)then
            write(titleN,"(i0,a,i0,a)") ix,"_",iz,".pdf"
            CALL SETFIL(trim(titleN))
            call qplot(real((/((i-1)*Dt,i=1,800)/),4),real(p_fot(1:800),4),800) 
          end if
      end do !ix
      end do !iz
      end do !imec
      close(6375)
      CALL chdir("..")
      end do !iFte
!     CALL chdir("..")
      end subroutine loadG_fotogramas
      

      subroutine plotSisGram(PSV,SH,guardarW)
      use resultVars , only : iPtini,iPtfin,allpoints
      use waveNumVars, only : NFREC,NPTSTIME
      use glovars, only : PlotFilledSabanas, comoFacDeAmpliDinamica
      use sourceVars, only : currentiFte
      use geometryvars, only : Xcoord_Voidonly, Xcoord_Incluonly
      integer :: iP,i,j
      character(LEN=100) :: yax,nam
      logical, intent(in) :: PSV,SH, guardarW
      complex*16, dimension(NPTSTIME) :: S,Sxx,Szx,Szz,Stt,Srt
      real*8 :: si,co
      if (guardarW) then
       ! archivo donde se guarda la solución en frecuencia para
       ! usarse en furturas corridas
       write(nam,'(a,I0,a)') "W",currentiFte,".bin"  
       OPEN(6373,FILE=trim(nam),STATUS='UNKNOWN', ACCESS='STREAM',ACTION='WRITE')
       write(nam,'(a,I0,a)') "Geom",currentiFte,".txt"  
       OPEN(9975,FILE=trim(nam),STATUS='UNKNOWN', ACTION='WRITE')
       if (allocated(Xcoord_Incluonly)) then
      do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont    !
      write(9975,'(EN22.4,2x,EN22.4,2x,EN22.4,2x,EN22.4)') & 
                       Xcoord_Incluonly(j,1,1),Xcoord_Incluonly(j,2,1) ,&
                       Xcoord_Incluonly(j,1,2),Xcoord_Incluonly(j,2,2) 
      end do
       end if!
       if (allocated(Xcoord_Voidonly)) then
      do j=1,size(Xcoord_Voidonly (:,1,1))!n_topo+1,n_topo+n_cont    !
      write(9975,'(EN22.4,2x,EN22.4,2x,EN22.4,2x,EN22.4)') & 
                       Xcoord_Voidonly(j,1,1), Xcoord_Voidonly(j,2,1) ,&
                       Xcoord_Voidonly(j,1,2), Xcoord_Voidonly(j,2,2) 
      end do
        end if
       close(9975)
       write(nam,'(a,I0,a)') "Secciones",currentiFte,".txt"  
       OPEN(9973,FILE=trim(nam),STATUS='UNKNOWN', ACTION='WRITE')
      end if !
      if (PSV) then
        do iP = iPtini,iPtfin
        !#< r ____________   W   _________________________________________ !#>
          write(yAx,'(a)') '$u_3$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%W
          end do; end if
!         print*,ip,"---------------------",allpoints(iP)%resp(:,currentiFte)%W; print*," "
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%W,'w--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,1,S,.false.)
          if (comoFacDeAmpliDinamica) &
          call W_to_t(allpoints(iP)%facAmpli(:,currentiFte)%W,'w--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,1,S,.true.)
        end do        
          call makeSabana('0_S-w__.pdf',.false.)
          if(PlotFilledSabanas) call makeSabana('0_S-w_f.pdf',.true.) ! filled traces
          
        !#< r ____________   U   _________________________________________ !#>
        do iP = iPtini,iPtfin
          write(yAx,'(a)') '$u_1$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%U
          end do; end if
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%U,'u--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,2,S,.false.)
          if (comoFacDeAmpliDinamica) &
          call W_to_t(allpoints(iP)%facAmpli(:,currentiFte)%U,'u--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,2,S,.true.)
        end do         
          call makeSabana('0_S-u__.pdf',.false.)
          if(PlotFilledSabanas) call makeSabana('0_S-u_f.pdf',.true.) ! filled traces
         
!       !#< r ____________   Tracciones   ________________________________ !#>
!       do iP = iPtini,iPtfin
!         write(yAx,'(a)') '$Tz_$ [m]'
!         call W_to_t(allpoints(iP)%resp(:,currentiFte)%Tz,'Tz-',yax,iP,& 
!                   allpoints(iP)%center%x,& 
!                   allpoints(iP)%center%z,0,S)
!         write(yAx,'(a)') '$Tx_$ [m]'
!         call W_to_t(allpoints(iP)%resp(:,currentiFte)%Tx,'Tx-',yax,iP,& 
!                   allpoints(iP)%center%x,& 
!                   allpoints(iP)%center%z,0,S)
!       end do
        
        !#< r ____________   sxx   _________________________________________ !#>
        do iP = iPtini,iPtfin
          write(yAx,'(a)') '$\sigma_{xx}$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%sxx
          end do; end if
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%sxx,'sxx',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,3,Sxx,.false.)
          if (comoFacDeAmpliDinamica) &
          call W_to_t(allpoints(iP)%facAmpli(:,currentiFte)%sxx,'sxx',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,3,Sxx,.true.)
!       end do
        !#< r ____________   szx   _________________________________________ !#>
!       do iP = iPtini,iPtfin
          write(yAx,'(a)') '$\sigma_{zx}$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%szx
          end do; end if
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%szx,'szx',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,4,Szx,.false.)
          if (comoFacDeAmpliDinamica) &
          call W_to_t(allpoints(iP)%facAmpli(:,currentiFte)%szx,'szx',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,4,Szx,.true.)
!       end do
        !#< r ____________   szz   _________________________________________ !#> 
!       do iP = iPtini,iPtfin
          write(yAx,'(a)') '$\sigma_{zz}$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%szz
          end do; end if
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%szz,'szz',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,5,Szz,.false.)
          if (comoFacDeAmpliDinamica) &
          call W_to_t(allpoints(iP)%facAmpli(:,currentiFte)%szz,'szz',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,5,Szz,.true.)
                    
          if (allpoints(iP)%isSeccion) then
             si = allpoints(iP)%sinT
             co = allpoints(iP)%cosT
             stt = si**2*Sxx + co**2*Szz - 2*si*co*Szx
             srt = si*co*(Szz-Sxx) + Szx * (co**2 - si**2)
              write(9973 ,'(EN22.4,2x,EN22.4,2x)', ADVANCE = "NO") & 
                 allpoints(iP)%center%x, allpoints(iP)%center%z  
            do i = 1, NPTSTIME ! stt !(tangencial)
              write(9973,'(EN22.4,2x)', ADVANCE = "NO") & 
              real(stt(i))
            end do!
            do i = 1, NPTSTIME ! srt !(cortante)
              write(9973,'(EN22.4,2x)', ADVANCE = "NO") & 
              real(srt(i))
            end do
              write(9973,'(a)') ""
          end if
        end do 
       end if !psv
        close(9973)
       if (SH) then
          do iP = iPtini,iPtfin
          write(yAx,'(a)') '$u_2$ [m]'
          if (guardarW) then ; do i = 1, NFREC+1
          write(6373) allpoints(iP)%resp(i,currentiFte)%V
          end do; end if
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%V,'v--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,0,S,.false.)    
          end do 
!         call chdir("..") 
          call makeSabana('0_S-v__.pdf',.false.) 
          call makeSabana('0_S-v_f.pdf',.true.)
!         call chdir("traces")
          do iP = iPtini,iPtfin
          write(yAx,'(a)') '$Ty_$ [m]'
          call W_to_t(allpoints(iP)%resp(:,currentiFte)%Ty,'Ty-',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z,0,S,.false.)     
          end do
!         call chdir("..")
       end if !sh  
       if (guardarW) close (6373)
      end subroutine plotSisGram
      
      subroutine loadW(PSV,SH)
      use resultVars , only : iPtini,iPtfin,allpoints
      use waveNumVars, only : NFREC
      use glovars, only : rutaOut
      use sourceVars, only : nFuentes
      integer :: iP,i,iFte
      logical, intent(in) :: PSV,SH
      character(LEN=100) :: nam
      call chdir(trim(adjustl(rutaOut)))
      do iFte = 1,nFuentes
      write(nam,'(a,I0)') "traces",iFte
      CALL chdir(trim(nam))
       write(nam,'(a,I0,a)') "W",iFte,".bin"
       OPEN(6372,FILE=trim(nam),STATUS='UNKNOWN', ACCESS='STREAM',ACTION='READ')
      if (PSV) then
        do iP = iPtini,iPtfin
           do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%W
           end do!
        end do!
        do iP = iPtini,iPtfin
           do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%U
           end do! 
        end do!
        do iP = iPtini,iPtfin
           do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%sxx
           end do! 
!       end do
!       do iP = iPtini,iPtfin
           do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%szx
           end do! 
!       end do
!       do iP = iPtini,iPtfin
           do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%szz
           end do! 
        end do
      end if!
      if (SH) then
          do iP = iPtini,iPtfin
          do i = 1, NFREC+1
             read(6372) allpoints(iP)%resp(i,iFte)%V
           end do!
          end do
      end if
      close (6372)
      CALL chdir("..")
      end do !ifte
      CALL chdir("..")
      end subroutine loadW
      

! Sismo/Foto- gramas 
      subroutine W_to_t(W,nombre,yAx,iP,x_i,z_i,icomp,Sout, soloFacAmpl)
      use waveNumVars, only : NFREC,DFREC, NPTSTIME, OMEI, t_vec
      use glovars
      use waveVars, only : dt,Uo,maxtime
      use ploteo10pesos
      use wavelets
      use resultvars, only : allpoints, Sabana, nSabanapts, nIpts, SabanaPlotIndividual
      use debugStuff
      use sourceVars, only : iFte=>currentiFte
      implicit none
      integer ,intent(in) :: iP,icomp
      complex*16, dimension(Nfrec+1), intent(in) :: W !espectro
      real*8, intent(in) :: x_i,z_i
      character(LEN=3)   :: nombre
      character(LEN=100) :: titleN,yAx, CTIT
      character(LEN=32)  :: name
      complex*16, dimension(NPTSTIME) :: S,Sout
      integer :: i,n_maxtime
      logical :: soloFacAmpl
      
      write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x_i,' , ',z_i,')'
 
      S = z0
      S(1:nfrec+1)= W(1:nfrec+1:+1)
      S(NPTSTIME-NFREC+2:NPTSTIME) = conjg(W(nfrec:2:-1))  
      
      ! (0) graficar en espectro sin corregir
        write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_a_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
        call plotXYcomp(S(1:int(nfrec)),real(DFREC,4), & 
         int(nfrec),titleN, 'frec[hz] ',yAx, CTIT ,1200,800,0.0)
        write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_a_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,I0)') 'Ampf_IP',iP,'_icomp',icomp
        OPEN(3202,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3202,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, & 
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(int(nfrec +1),1,S(1:int(nfrec +1)),name,3202)
        close (3202)
      
      if (soloFacAmpl) then  
      ! grafica simple:
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               '3_Ampfactor_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S(1:int(nfrec +1)),real(DFREC), & 
         int(nfrec +1),titleN, 'frec[hz] ',yAx, CTIT ,1200,800,0.0)
         
        write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_ampfac_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,I0)') 'Ampf_IP',iP,'_icomp',icomp
        OPEN(3212,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3212,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, & 
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(int(nfrec +1),1,S(1:int(nfrec +1)),name,3212)
        close (3212)
        return
      else
      S = S * t_vec ! t0 tiempo inicial
      end if
      S = S * Uo(:,iFte) ! conv con fucion de amplitud
      
      if (Verbose .ge. 4) call showMNmatrixZ(nptstime,1, S,"  S  ",6)
      
            
      !  (1) pasar al tiempo
         S = FFTW(NPTSTIME,S,+1,1/(NPTSTIME*dt)) !backward
         
      !  (2) remover efecto de la frecuencia imaginaria (DWN)
         S = S * exp(- OMEI * Dt*((/(i,i=0, NPTSTIME-1)/)))
         
!        !tiempo maximo para graficar
         n_maxtime = int(maxtime(iFte)/dt)
         if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
         if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
         S(n_maxtime+1: NPTSTIME) = z0;
         Sout = S
      ! guardar para hacer sabana o plotear
      if (allpoints(iP)%isSabana) then
         ! guardamos la sabana actual
!        print*,"saved sabana point, ",iP-(nIpts - nSabanapts)
         Sabana(iP-(nIpts - nSabanapts),1:NPTSTIME) = S
         if (SabanaPlotIndividual .eqv. .false.) return
      end if
      
      ! guardar componentes en el borde
      if (allpoints(iP)%atBou .and. icomp .ne. 0) then
!        !allocate(allpoints(iP)%S(NPTSTIME,5)) !W U sxx szx szz
         allpoints(iP)%S(1:NPTSTIME,icomp) = S
!        write(name,'(a,I0,a,I0,a)') 's_IP',iP,'_icomp',icomp,'.m'
!        OPEN(3211,FILE=trim(name),FORM="FORMATTED",ACTION='WRITE')
!        write(3211,'(a,a,EN26.9,a,EN26.9,a)') yax, & 
!        ' en (', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')'
!        call scripToMatlabMNmatrixZ(NPTSTIME,1,S,name,3211)
!        close (3211)
      end if
      
      !  (3) GRAFICAR en el tiempo
!#< g   imprimir sismograma  !#>
         write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_s_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S(1:n_maxtime) ,real(dt,4),n_maxtime,titleN, & 
         'time[sec]',yAx, CTIT ,1200,800, 0.)
        
        write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_s_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,I0)') 's_IP',iP,'_icomp',icomp
        OPEN(3214,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3214,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, & 
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(NPTSTIME,1,S(1:int(NPTSTIME)),name,3214)
        close (3214)
         
      !  (4) espectro correcto 
!       ! tapper ?
        
        ! fft 
        S = FFTW(NPTSTIME,S(1:int(NPTSTIME)),-1,dt) !forward
        
!#< g   imprimir espectro  !#>
      if (verbose .ge. 1) then
       if ((SabanaPlotIndividual .eqv. .false.) .and. & 
           (allpoints(iP)%isSabana .eqv. .true.)) then
       else
      ! grafica simple:
         write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_f_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S(1:int(NPTSTIME/2)),real((DFREC*NFREC)/(NPTSTIME/2),4), & 
         int(NPTSTIME/2),titleN, 'frec[hz] ',yAx, CTIT ,1200,800,0.0)
         
        write(titleN,'(I0,a,a,a,I0,a,I0,a,I0,a,I0,a)') & 
               iP,'_f_',nombre,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].m'
        write(name,'(a,I0,a,I0)') 'f_IP',iP,'_icomp',icomp
        OPEN(3212,FILE=trim(titleN),FORM="FORMATTED",ACTION='WRITE')
        write(3212,'(a,a,/,a,EN26.9,a,EN26.9,a,/,a)') '%',yAx, & 
         '%(', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')','%% '
        call scripToMatlabMNmatrixZ(int(NPTSTIME/2+1),1,S(1:int(NPTSTIME/2+1)),name,3212)
        close (3212)
        
!       if (allpoints(iP)%atBou .and. icomp .ne. 0) then
!        write(name,'(a,I0,a,I0,a)') 'f_IP',iP,'_icomp',icomp,'.m'
!        OPEN(3212,FILE=trim(name),FORM="FORMATTED",ACTION='WRITE')
!        write(3212,'(a,a,EN26.9,a,EN26.9,a)') yax, & 
!        ' en (', allpoints(iP)%center%x,',',allpoints(iP)%center%z,')'
!        call scripToMatlabMNmatrixZ(int(NPTSTIME/2),1,S(1:int(NPTSTIME/2)),name,3212)
!        close (3212)
!       end if  
       end if
!     ! grafica logaritmica
!           write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
!              'fL_',nombre,iP,'[', &
!              int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
!              int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
!           logflag = 'logx     '
!           write(xAx,'(a)') "frec[Hz]"
!           call plotSpectrum(S(1: NPTSTIME/2),real(DFREC,4), NPTSTIME/2, NPTSTIME/2, & 
!               titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
      end if
         
      end subroutine W_to_t
      
      subroutine makeSabana(nombre,filled)
      use waveNumVars, only : NPTSTIME,nfrec
      use glovars
      use waveVars, only : dt,maxtime
      use resultvars, only : Sabana, nSabanapts,sabZeroini,sabZerofin
      use sourceVars, only : iFte=>currentiFte
      use dislin   
      implicit none
      character(LEN=11)   :: nombre
      logical :: filled
      integer :: n_maxtime
      integer :: i,nPow10x,W,H
      real, dimension(NPTSTIME) :: x
      complex*16, dimension(:,:), pointer :: S
      real :: offset,escala
      character(LEN=100) :: dumb
      real maxY,minY,xstep,ystep
      character(LEN=1) :: imdone
      real, dimension(5,2) :: rec
      real, dimension(:), allocatable :: falda
      integer*4 :: lentitle
      character(LEN=60) :: CTIT
      
      if (nSabanapts .eq. 0) return
      !loop de escala y offset preguntando
      n_maxtime = int(maxtime(iFte)/dt)
      if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
      if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
      S => Sabana(1:nSabanapts,1:n_maxtime)
       
      if(sabZeroini .ne. 0) S(sabZeroini:sabZerofin,:) = z0
      
       ! primero usamos la escal para normalizar, luego para aumentar el constraste
       escala = 1.0/maxval(abs(real(S,4)))
       offset = 0.5 * maxval(abs(real(S,4)))*escala
       escala = escala * 2.0
       
       x = 0
      
      ! coordenada x (tiempo)
      ! para que se vea el texto en los ejes si es muy pequeño en formato F3.1
      nPow10x = 0
      if (Dt *(n_maxtime-1) < 0.6) then
        do i = 1,10
          if (Dt *(n_maxtime-1)*(10.0**i) > 1.0) then
            exit
          end if 
        end do 
        nPow10x = i
        
      elseif (Dt * (n_maxtime-1) > 6000.) then  
        do i = 1,10
          if (Dt *(n_maxtime-1)*(10.0**(-i)) < 1000.0) then
            exit
          end if
        end do
        nPow10x = -i
      end if
      
      DO i = 1,n_maxtime
        x(i) = real(Dt*(i-1)*(10.0**(nPow10x)),4)
      END DO
      
6720   W = 1200
       H = 800
       CALL METAFL('PDF')
       CALL SETFIL(trim(adjustl(nombre)))
       call filmod('DELETE')
       CALL SETPAG('DA4P')
       CALL PAGE(int(W+1200,4),int(H+350,4))
       CALL PAGMOD('NONE')
       CALL DISINI()
      CALL TEXMOD ('ON') ! latex!!
      call errmod ("all", "off") 
      CALL COMPLX 
      CALL HWFONT()
      CALL axspos (int(70,4) ,int(H+100,4))
      call axslen (int(W+1050,4) ,int(H,4))
      write(dumb,'(a,I0,a)') 't (sec) [ 10^',(nPow10x *(-1)),' ]'
       call name(trim(dumb),'X')
       call labdig(int(1,4),'X')
       call ticks (int(5,4) ,'X')
       call setgrf ('NAME','LINE','LINE','LINE')
       
       minY = -2.0*offset
       maxY = offset*nSabanapts+offset
       xstep = x(n_maxtime)/6.0
       ystep = (maxY-minY)/6.0
       
       call graf(real(x(1),4), & 
                real(x(n_maxtime),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
       
       do i = 1, nSabanapts
         call PENWID(real(1,4))
        if (filled) then ! sabana coloreada al estilo sismología
         if(.not. allocated(falda)) allocate(falda(n_maxtime))
         ! valores negativos van de gris
         call color('GRAY')
         call PENWID(real(0.0001,4))
         CALL SHDPAT (int(16,4)) ! filled shading
         falda = 0.0_4
         falda(1:n_maxtime) = real(S(i,1:n_maxtime),4)
         where(falda(1:n_maxtime) .gt. 0.0_4)
           falda(1:n_maxtime) = 0.0_4
         end where
         falda(1) = 0.0_4
         falda(n_maxtime) = 0.0_4
         call rlarea(real(x,4),real(falda(1:n_maxtime) * escala + (i-1)*offset,4),int(n_maxtime,4))
         
         ! valores positivos van de negro
         call SETRGB(0.05_4, 0.05_4, 0.05_4)
         call PENWID(real(0.0001,4))
         CALL SHDPAT (int(16,4)) ! filled shading
         falda = 0.0_4
         falda(1:n_maxtime) = real(S(i,1:n_maxtime),4)
         where(falda(1:n_maxtime) .lt. 0.0_4)
           falda(1:n_maxtime) = 0.0_4
         end where
         falda(1) = 0.0_4
         falda(n_maxtime) = 0.0_4
         call rlarea(real(x,4),real(falda(1:n_maxtime) * escala + (i-1)*offset,4),int(n_maxtime,4))
         
         call PENWID(real(0.1,4))
        end if 
         ! la curva
         call color('FORE')
         call curve(real(x,4) ,real(S(i,1:n_maxtime) * escala + (i-1)*offset,4) ,int(n_maxtime,4))
       end do
         call PENWID(real(1.0,4))
       ! cuadro blanco al final de la sabana
       call color('BACK') 
       CALL SHDPAT (int(16,4)) ! filled shading
      rec(1,1) = x(n_maxtime) - x(n_maxtime)*0.05
      rec(1,2) = minY
      rec(2,1) = x(n_maxtime) - x(n_maxtime)*0.05
      rec(2,2) = maxY
      rec(3,1) = x(n_maxtime)
      rec(3,2) = maxY
      rec(4,1) = x(n_maxtime)
      rec(4,2) = minY
      rec(5,1) = x(n_maxtime) - x(n_maxtime)*0.05
      rec(5,2) = minY
       CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
       
       ! traza de escala vertical
       call color('FORE')
       CALL HSYMBL(int(8,4)) !size of symbols
       CALL RLVEC (real(x(n_maxtime)*0.96,4), & 
                   real(minY + (maxY-minY)/2 + 2 ,4), & 
                   real(x(n_maxtime)*0.96,4), & 
                   real(minY + (maxY-minY)/2 - 2 ,4), int(1122,4)) 
     
      write(CTIT,'(EN12.4)') 4.0/escala
      lentitle = NLMESS(CTIT)
      CALL HEIGHT (int(19,4))
      CALL ANGLE(int(90,4))
      CALL MESSAG(CTIT,int((2250),4),int(580,4))
      
       
       call disfin() 
       return
       print *, char(7)
       write(6,'(a)', ADVANCE = "NO") & 
         'ok? [anykey] or replot with different scale [Y]: '
      
         read(5,*)imdone
         if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
            write(6,'(a,E10.1)') "Current scale,", escala
            write(6,'(a)', ADVANCE = "NO")'new scale = '
            read(5,*) escala
            go to 6720
         end if
         !
      if(associated(S)) nullify(S)
      end subroutine makeSabana
      

      subroutine F_K_exp(XF) !calc y graf de FK y CK a partir de sábana de sensores
      use waveNumVars, only : NFREC,dfrec,omei
      use soilvars, only : Qq
      use glovars
      use wavelets
      use resultvars, only : allpoints,nSabanapts,nIpts
      use dislin
      use debugStuff
      use sourceVars, only : currentiFte
      implicit none
      complex*16, dimension(nSabanapts,Nfrec+1,2) :: XF
      complex*16, dimension(nSabanapts,Nfrec+1,2) :: XFc
      complex*16, dimension(nfrec+1, nSabanapts/2,2) :: FC
      complex*16, dimension(nfrec+1, nSabanapts/2) :: uno
      integer :: ij, iP
      complex*16 :: come
      real*8 :: ome,dx
!     CHARACTER(len=32) :: arg
      logical :: alguno
      real, parameter :: p1 = 5. ! sharpness parameter
      real, parameter :: p2 = 1. ! sharpness parameter
      real, parameter :: p3 = 4 ! sharpness parameter
      ! guardar XF
!     print*, nIpts, nSabanapts
!     call winsiz(int(1200,4),int(1200,4))
      CALL SETPAG('DA4L')
      alguno = .false.
      do iP = 1,nIpts
       if (allpoints(iP)%isSabana) then
        alguno = .true.
        XF((iP-(nIpts-nSabanapts)) ,1:NFREC+1,1) = allpoints(ip)%resp(1:NFREC+1,currentiFte)%W
        XF((iP-(nIpts-nSabanapts)) ,1:NFREC+1,2) = allpoints(ip)%resp(1:NFREC+1,currentiFte)%U
!       print*,allpoints(ip)%center%x,"  ",(iP-(nIpts-nSabanapts))
       end if
      end do
      
      ! XF -> KF -------------------------------
      if (.not. alguno) return
      dx = (allpoints(nSabanapts)%center%x - allpoints(nIpts-nSabanapts+1)%center%x)
      dx = dx / (nSabanapts+1)
      do ij = 1,NFREC+1
      XF(1:nSabanapts,ij,1) = cshift(XF(1:nSabanapts,ij,1), nSabanapts/2)
      XF(1:nSabanapts,ij,2) = cshift(XF(1:nSabanapts,ij,2), nSabanapts/2)
      !escala
!     XF(1:nSabanapts,ij,1) = XF(1:nSabanapts,ij,1) * (sqrt(1.0*nSabanapts) * dx)
!     XF(1:nSabanapts,ij,2) = XF(1:nSabanapts,ij,2) * (sqrt(1.0*nSabanapts) * dx)
      
      call FORK(nSabanapts,XF(1:nSabanapts,ij,1),+1) 
      XF(1:nSabanapts,ij,1) = XF(1:nSabanapts,ij,1) * sqrt(1.0*nSabanapts) * dx ! factor de escala

      call FORK(nSabanapts,XF(1:nSabanapts,ij,2),+1) 
      XF(1:nSabanapts,ij,2) = XF(1:nSabanapts,ij,2) * sqrt(1.0*nSabanapts) * dx ! factor de escala

      
      XF(1:nSabanapts,ij,1) = cshift(XF(1:nSabanapts,ij,1), nSabanapts/2)
      XF(1:nSabanapts,ij,2) = cshift(XF(1:nSabanapts,ij,2), nSabanapts/2)
      
      end do
      
!      open(427,FILE= "outKF.m",action="write",status="replace")
!      write(arg,'(a)') "w_KF"
!      call scripToMatlabMNmatrixZ(nSabanapts,Nfrec+1,XF(:,:,1),arg,427)
!      write(arg,'(a)') "u_KF"
!      call scripToMatlabMNmatrixZ(nSabanapts,Nfrec+1,XF(:,:,2),arg,427)
!      close(427)
      
      ! comprimir
      XFc(1:nSabanapts,1:Nfrec+1,1) = real(log(1. + exp(p1)*abs(XF(1:nSabanapts,1:Nfrec+1,1))) / & 
           log(exp(p1)+1.),4)
      XFc(1:nSabanapts,1:Nfrec+1,1) = XFc(1:nSabanapts,1:Nfrec+1,1) / maxval(abs(XFc(1:nSabanapts,1:Nfrec+1,1)))
      CALL SETFIL("2_w_KF.pdf")
      CALL QPLCLR(real(abs(XFc(1:nSabanapts,1:NFREC+1,1)),4), nSabanapts,NFREC+1) 
      
      XFc(1:nSabanapts,1:Nfrec+1,2) = real(log(1. + exp(p1)*abs(XF(1:nSabanapts,1:Nfrec+1,2))) / & 
           log(exp(p1)+1.),4)
      XFc(1:nSabanapts,1:Nfrec+1,2) = XFc(1:nSabanapts,1:Nfrec+1,2) / maxval(abs(XFc(1:nSabanapts,1:Nfrec+1,2)))
      CALL SETFIL("2_u_KF.pdf")
      CALL QPLCLR(real(abs(XFc(1:nSabanapts,1:NFREC+1,2)),4), nSabanapts,NFREC+1) 
      
      uno = 1;
      
      ! velocidad  k positivo
      do ij = 1,NFREC+1
      ome = (ij*dfrec)*2*pi
      COME = CMPLX(OME, OMEI,8)!periodic sources damping
      COME = COME * cmplx(1.0, -1.0/2.0/Qq,8)
      FC(ij,nSabanapts/2:1:-1,1) = come/XF(nSabanapts/2+1:nSabanapts,ij,1)
      FC(ij,nSabanapts/2:1:-1,2) = come/XF(nSabanapts/2+1:nSabanapts,ij,2)
      end do
      
      ! comprimir
      FC(:,:,1) = real(log(1. + exp(p2)*abs(FC(:,:,1))) / & 
           log(exp(p2)+1.),4)
      FC(:,:,2) = real(log(1. + exp(p2)*abs(FC(:,:,2))) / & 
           log(exp(p2)+1.),4) 
      
      FC(:,:,1) = uno - FC(:,:,1)/maxval(abs(FC(:,:,1)))
      FC(:,:,2) = uno - FC(:,:,2)/maxval(abs(FC(:,:,2)))
      
      ! aumentar contraste un poquito
      FC(:,:,1) = real((abs(FC(:,:,1))**p3 * 0.4**p3)/(abs(FC(:,:,1))**p3 + 0.4**p3),4)
      FC(:,:,1) = FC(:,:,1) / maxval(abs(FC(:,:,1)))
      CALL SETFIL("3_w_CFp.pdf")
      CALL QPLCLR(real(abs(FC(:,:,1)),4), NFREC+1,nSabanapts/2) 
       
      FC(:,:,2) = real((abs(FC(:,:,2))**p3 * 0.4**p3)/(abs(FC(:,:,2))**p3 + 0.4**p3),4)
      FC(:,:,2) = FC(:,:,2) / maxval(abs(FC(:,:,2)))
      CALL SETFIL("3_u_CFp.pdf")
      CALL QPLCLR(real(abs(FC(:,:,2)),4), NFREC+1,nSabanapts/2) 
      
      ! velocidad  k negativo
      do ij = 1,NFREC+1
      ome = (ij*dfrec)*2*pi
      COME = CMPLX(OME, OMEI,8)!periodic sources damping
      COME = COME * cmplx(1.0, -1.0/2.0/Qq,8)
!     FC(ij,nSabanapts/2:1:-1,1) = come/XF(1:nSabanapts/2,ij,1)
!     FC(ij,nSabanapts/2:1:-1,2) = come/XF(1:nSabanapts/2,ij,2)
      FC(ij,1:nSabanapts/2,1) = come/XF(1:nSabanapts/2,ij,1)
      FC(ij,1:nSabanapts/2,2) = come/XF(1:nSabanapts/2,ij,2)
      end do
      
      ! comprimir
      FC(:,:,1) = real(log(1. + exp(p2)*abs(FC(:,:,1))) / & 
           log(exp(p2)+1.),4)
      FC(:,:,2) = real(log(1. + exp(p2)*abs(FC(:,:,2))) / & 
           log(exp(p2)+1.),4) 
      
      FC(:,:,1) = uno - FC(:,:,1)/maxval(abs(FC(:,:,1)))
      FC(:,:,2) = uno - FC(:,:,2)/maxval(abs(FC(:,:,2)))
      
      ! aumentar contraste un poquito
      FC(:,:,1) = real((abs(FC(:,:,1))**p3 * 0.4**p3)/(abs(FC(:,:,1))**p3 + 0.4**p3),4)
      FC(:,:,1) = FC(:,:,1) / maxval(abs(FC(:,:,1)))
      CALL SETFIL("3_w_CFn.pdf")
      CALL QPLCLR(real(abs(FC(:,:,1)),4), NFREC+1,nSabanapts/2) 
       
      FC(:,:,2) = real((abs(FC(:,:,2))**p3 * 0.4**p3)/(abs(FC(:,:,2))**p3 + 0.4**p3),4)
      FC(:,:,2) = FC(:,:,2) / maxval(abs(FC(:,:,2)))
      CALL SETFIL("3_u_CFn.pdf")
      CALL QPLCLR(real(abs(FC(:,:,2)),4), NFREC+1,nSabanapts/2) 
           
      end subroutine F_K_exp
      
     
      subroutine Churubusco(testPoints)
      use glovars, only : workBoundary,flip12,Printnum
      use DISLIN
      use peli, only : ypray => coords_Z, xpray => coords_X,& 
                     fotogramas,fotogramas_Region
      use meshVars, only : npixX,npixZ,MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ,nmarkX
      use waveVars, only : dt,maxtime
      use waveNumVars, only : NFREC, NPTSTIME
      use soilVars, only : Z,N,col=>layershadecolor, shadecolor_inc
      use geometryvars, only : nXI,Xcoord_ER, & 
                               Xcoord_Voidonly, Xcoord_Incluonly,Xcoord_flip_out,&
                               n_topo,n_cont,n_vall
      use resultvars, only : Punto,BouPoints,nbpts,allpoints,punEnlaFront,& 
                             nPtsolos,nIpts,& !nSabanapts,nSecciones,&
                             nBPt_topo,nBPt_cont,nBPt_vall
      use ploteo10pesos
      use sourceVars, only : iFte => currentiFte
      implicit none
      interface
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, & 
              plotReceptoresA, plotReceptoresB,plotFuente)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente
      end subroutine drawBoundary
      end interface
      real, dimension(:,:,:), allocatable :: xvmat,yvmat
      real :: maV1,maV2,minX,maxX,minY,maxY,xstep,zstep, & 
              encuadre, tlabel, madmax, escalaFlechas
      real*8 :: mxU,mxW
      real,dimension(nIpts*2) :: delX,delZ
      integer :: i,ii,j,n_maxtime,iT,k,fai,faf
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      integer*4 :: lentitle
      character(LEN=60) :: CTIT
      type (Punto), dimension(:), pointer :: BP
      real*8, dimension(:,:),allocatable :: rec
      logical :: testPoints
      integer :: nframes
      ! fotogramas tipo campo vectorial
      minx = MeshMinX
      maxx = MeshMaxX
      miny = MeshMinZ
      maxy = MeshMaxZ
      if (testPoints) then
      nframes = 1
      else
      
      !tiempo maximo para graficar
       n_maxtime = int(maxtime(iFte)/dt)
       if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
       if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
       write(Printnum,'(a,F5.2,a,F7.4,a,F6.2,a)')& 
       "maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," puntos"
       allocate(xvmat(npixX,npixZ,n_maxtime))
       allocate(yvmat(npixX,npixZ,n_maxtime))
       xstep = real(((maxX-minX) / nmarkX ))
      
!     CALL chdir("video")
      nframes = n_maxtime
      maV1 = maxVal(real(fotogramas(:,:,1:n_maxtime,1,iFte),4))
      maV2 = maxVal(real(fotogramas(:,:,1:n_maxtime,2,iFte),4))
      do i=1,npixZ
        do j=1,npixX
          xvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime,2,iFte),4) !U
          yvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime,1,iFte),4) !W
          
          ! para no imprimir lo que se ve muy pequeñito
          do iT = 1,n_maxtime
            if (abs(yvmat(j,i,iT)) .le. mav1*0.005) yvmat(j,i,iT) = 0.0
            if (abs(xvmat(j,i,iT)) .le. mav2*0.005) xvmat(j,i,iT) = 0.0
          end do
          
          !print*,j,i,yvmat(j,i,1)+xvmat(j,i,1)
          ! para no imprimir los picos raros
            if (abs(yvmat(j,i,1)) .gt. abs(mav1*0.1)) then
                 yvmat(j,i,1:n_maxtime) = 0.0
                 xvmat(j,i,1:n_maxtime) = 0.0
            end if!
            if (abs(xvmat(j,i,1)) .gt. abs(mav2*0.1)) then
                 yvmat(j,i,1:n_maxtime) = 0.0
                 xvmat(j,i,1:n_maxtime) = 0.0
            end if
        end do
      end do
      madmax = max(max(maxval(xvmat),maxval(yvmat)),max(maxval(abs(xvmat)),maxval(abs(yvmat))))
      escalaFlechas = real((xstep) / madmax)
      end if
      
      
      xstep = real(((maxX-minX) / 0 ))
      zstep = real(((maxY-minY) / 0 ))
      encuadre = (maxY-minY)/(maxX-minX)
!     encuadre = (ypray(npixZ)-ypray(1))/(xpray(npixX)-xpray(1))
      
      CALL METAFL('PNG')
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(3000,4),int(3000,4))
!     CALL PAGMOD('NONE')
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(1200,4)) !1200 ambos
      CALL SCRMOD('REVERS') !fondo blanco
      do i=1,nframes
       if (testPoints) then
      write(titleN,'(a)') '0___Sensors_movie.png'
       else
      write(titleN,'(a,I0,a)') 'foto_',i,'.png'
       end if
      CALL SETFIL(trim(titleN))
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
           !the position of an axis system.
      CALL axspos (int(300,4) ,int(2700,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600*encuadre,4)) !size of the axis system.
!     call name('X [m]','X')
!     call name('Z [m]','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(2,4),'Y')
      call ticks (int(1,4) ,'XY')
      call setgrf("NAME", "NAME", "LINE", "LINE") 
      call height(40) ! de los caracteres
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
     
      
      !estratigrafía --------------------------------------------
      call shdpat(int(16,4))                                    !
      ii = 1
      if (Z(0) .lt. 0.0) ii = 0
      do J=ii,N                                                  !
         call SETRGB(col(J), col(J), col(J))                    !
         call rlrec(real(minX,4),real(max(Z(J),minY),4),&       !
                    real(maxX-minX,4),real(Z(J+1)-max(Z(J),minY),4))!
         call color ('FORE')                                    !
         call rline(real(minX,4),real(max(Z(J),minY),4), &      !
                 real(maxX,4),real(max(Z(J),minY),4))           !
      end do                                                    !
      J = N+1                                                   !
      call SETRGB(col(J), col(J), col(J))                       !
      call rlrec(real(minX,4),real(Z(J),4),&                    !
                    real(maxX-minX,4),real(maxY-Z(J),4))        !
      call color ('FORE')                                       !
      call rline(real(minX,4),real(Z(J),4), &                   !
                 real(maxX,4),real(Z(J),4))                     !
      ! Borrar estratos en la cuenca                                  !
      call color ('BACK')                                             !
      call shdpat(int(16,4))                                          !
      if (abs(z(0)) .lt. 0.0001) then
      call rlrec(real(minX,4),real(minY,4),&                          !
                    real(maxX-minX,4),real(Z(1)-minY,4))              !
      end if!
      
      if (workboundary) then !#< r ----------------------------------!#>
      ! dibujar inclusión
      if (flip12) then ! Xcoord_flip_out
!     print*,"flip12"
       if (allocated(Xcoord_flip_out)) then
       if (size(Xcoord_flip_out(:,1,1)) .gt. 1) then
        if (allocated(rec)) deallocate(rec)
        allocate(rec(2* (size(Xcoord_flip_out(:,1,1))),2))              !
        ii=1                                                             !
        do j=1,size(Xcoord_flip_out (:,1,1))!n_topo+1,n_topo+n_cont    !
        rec(ii,1) = Xcoord_flip_out(j,1,1)                              !
        rec(ii,2) = Xcoord_flip_out(j,2,1)                              !
        rec(ii+1,1) = Xcoord_flip_out(j,1,2)                            !
        rec(ii+1,2) = Xcoord_flip_out(j,2,2)                            !
!     print*,j,rec(i,1),rec(i,2),rec(i+1,1),rec(i+1,2)
        ii=ii+2                                       !
        end do
        call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
        CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_flip_out(:,1,1)),4))  !
        deallocate(rec) 
       end if
       end if
      else ! flip12--------
       if (allocated(Xcoord_Incluonly)) then
       if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
        if (allocated(rec)) deallocate(rec)
        allocate(rec(2* (size(Xcoord_Incluonly(:,1,1))),2))              !
        ii=1                                                             !
        do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont     !
        rec(ii,1) = Xcoord_Incluonly(j,1,1)                     !
        rec(ii,2) = Xcoord_Incluonly(j,2,1)                     !
        rec(ii+1,1) = Xcoord_Incluonly(j,1,2)                 !
        rec(ii+1,2) = Xcoord_Incluonly(j,2,2)               !
        ii=ii+2                                                          !
        end do
      
        call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
        CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))
        deallocate(rec) 
       end if!
       if (size(Xcoord_Voidonly(:,1,1)) .gt. 1) then
        if (allocated(rec)) deallocate(rec)
        allocate(rec(2* (size(Xcoord_Voidonly(:,1,1))),2))              !
        ii = 1                                                          !
        do j=1,size(Xcoord_Voidonly(:,1,1))                             !
        rec(ii,1) = Xcoord_Voidonly(j,1,1)                              !
        rec(ii,2) = Xcoord_Voidonly(j,2,1)                              !
        rec(ii+1,1) = Xcoord_Voidonly(j,1,2)                            !
        rec(ii+1,2) = Xcoord_Voidonly(j,2,2)                            !
        ii=ii+2                                                         !
        end do                                                          !
        call color ('BACK')                                             !
        call shdpat(int(16,4))                                          !
        CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), &                !
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))              !
        deallocate(rec) 
       end if                                                !         !
      end if
      end if !flip12
      
      !#< r ## dibujar contorno naranja de geometría deformada     !#> 
      if (punEnlaFront .and. .not. testPoints) then 
!     print*,"nIpts = ",nIpts
!     print*,"nXi =   ",nXi
!     print*,"nPtsolos = ",nPtsolos
!     print*,"nSabana=",nSabanapts
!     print*,"nSeccio=",nSecciones
!     print*,"n_topo =",nBPt_topo
!     print*,"n_cont =",nBPt_cont
!     print*,"n_vall =",nBPt_vall
      
      call color ('ORANGE')                                           !
      call PENWID(real(4.0,4))
      ii = 1
      if (n_topo .gt. 0) then
      k = 100000
      fai = nPtsolos+1!nIpts-nXi-nSabanapts
      faf = nPtsolos+nBPt_topo!nIpts-nXi-nSabanapts+n_topo
!     print*,"hay n_topo =",fai,faf
!     fai = nIpts-nXi-nSabanapts-nSecciones
!     faf = nIpts-nXi-nSabanapts-nSecciones+n_topo
!     print*,"hay n_topo =",fai,faf
!     print*,j
      do j=fai,faf
        if (allpoints(j)%atBou) then
          mxU = maxval(abs(allpoints(j)%S(1:nframes,2)))
          mxW = maxval(abs(allpoints(j)%S(1:nframes,1)))
!         print*,"mxU=",mxU,"  mxW=",mxW
          k = min(ii,k)
          delX(ii) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2) / mxU,4)!U
          delZ(ii) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1) / mxW,4)!W
!       print*,"   ",ii,allpoints(j)%center%x,allpoints(j)%center%z," -> ",delX(ii),delZ(ii)
          ii = ii + 1
        end if
      end do
      if (k .ne. 100000) then
      delX(ii) = delX(1)
      delZ(ii) = delZ(1)
      do j = k,ii-1
        call rline(delX(j),delZ(j),delX(j+1),delZ(j+1))
      end do
      end if ! k
      ii = ii + 1
      end if ! hay n_topo
      if (n_cont .gt. 0) then
      k = 100000
      fai = nPtsolos+nBPt_topo+1!nIpts-nXi-nSabanapts+n_topo+1
      faf = nPtsolos+nBPt_topo+nBPt_cont!nIpts-nXi-nSabanapts+n_topo+n_cont
!     print*,"hay n_cont =",fai,faf
!     fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+1
!     faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont
!     print*,"hay n_cont =",fai,faf
      do j=fai,faf
!       print*,allpoints(j)%center
!       if (allpoints(j)%isOnInterface) cycle
        if (allpoints(j)%atBou) then
          mxU = maxval(abs(allpoints(j)%S(1:nframes,2)))
          mxW = maxval(abs(allpoints(j)%S(1:nframes,1)))
!         print*,"mxU=",mxU,"  mxW=",mxW
          k = min(ii,k)
          delX(ii) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2)/mxU,4)!x
          delZ(ii) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1)/mxW,4)!y
!       print*,"   ",ii,allpoints(j)%center%x,allpoints(j)%center%z," -> ",delX(ii),delZ(ii)
          ii = ii + 1
        end if
      end do
      if (k .ne. 100000) then
       delX(ii) = delX(k)
       delZ(ii) = delZ(k)
!     print*,"   ",ii,k," -> ",delX(ii),delZ(ii)
       do j = k,ii-1
        call rline(delX(j),delZ(j),delX(j+1),delZ(j+1))
       end do
      end if ! k
      ii = ii + 1
      end if! hay n_cont
      if (n_vall .gt. 0) then
      k = 100000
      fai = nPtsolos+nBPt_topo+nBPt_cont+1!nIpts-nXi-nSabanapts+n_topo+n_cont+1
      faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall!nIpts-nXi-nSabanapts+n_topo+n_cont+n_vall
!     print*,"hay n_vall =",fai,faf
!     fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+1
!     faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
!     print*,"hay n_vall =",fai,faf
      do j=fai,faf
!     print*,j
        if (allpoints(j)%atBou) then
          mxU = maxval(abs(allpoints(j)%S(1:nframes,2)))
          mxW = maxval(abs(allpoints(j)%S(1:nframes,1)))
!         print*,"mxU=",mxU,"  mxW=",mxW
          k = min(ii,k)
          delX(ii) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2)/mxU,4)!x
          delZ(ii) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1)/mxW,4)!y
!       print*,"   ",ii,allpoints(j)%center%x,allpoints(j)%center%z," -> ",delX(ii),delZ(ii)
          ii = ii + 1
        end if
      end do
      if (k .ne. 100000) then
      delX(ii) = delX(k)
      delZ(ii) = delZ(k)
!     print*,"   ",ii,k," -> ",delX(ii),delZ(ii)
      do j = k,ii-1
        call rline(delX(j),delZ(j),delX(j+1),delZ(j+1))
      end do
      end if !k
      ii = ii + 1
      end if!
      ! n_topo,n_cont,n_vall
!     do j = 1,ii-2
!       call rline(delX(j),delX(j+1),delZ(j),delZ(j+1))   !
!     end do
      else !(punEnlaFront .and. .not. testPoints)
        call color ('FORE')                                           !
        call PENWID(real(0.5,4)) 
        call marker(int(-1,4)) ! sin marcadores                       !
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do
      end if !(punEnlaFront .and. .not. testPoints)
      call color ('FORE')
      call PENWID(real(0.5,4))
      end if! WORKBOUNDARY
      if (testPoints) then
       ! indicar asignacion de regiones
      CALL HSYMBL(int(9,4)) !size of symbols
      do ii=1,npixZ
        do j=1,npixX
          if (fotogramas_Region(ii,j) .eq. 0) cycle
          if (fotogramas_Region(ii,j) .eq. 1) call color('RED')
          if (fotogramas_Region(ii,j) .eq. 2) call color('BLUE')
          CALL RLSYMB (16, xpray(j), ypray(ii))
        end do
      end do
      
      else !testPoints

      !campo de desplazamientos ———————————————————-------------------
      call color ('FORE')                                            !
      call vecclr(-1) ! (-2):color de las puntas de flecha activado  !
      CALL VECOPT(real(escalaFlechas,4),'SCALE')                     !
      CALL VECOPT(real(10.0,4),'ANGLE')                              !
      CALL VECOPT(real(0.9,4),'LENGTH')                              !
      call vecmat(xvmat(:,:,i), &                     !
                  yvmat(:,:,i), &                     !
                  npixX, npixZ,xpray,ypray,int(1201))                !
      
      call color ('FORE')
      tlabel = (i)*real(dt,4)
      write(CTIT,'(a,F9.5,a)') 't=',tlabel,' seg'
      lentitle = NLMESS(CTIT)
      CALL MESSAG(CTIT,int(300,4),int(2850,4))
      end if !testPoints
      call disfin
      
      end do ! i=1,n_maxtime
      
      if (testPoints) return
      
      write(titleN,'(a)') 'foto_0.png'
      write(extension,'(a)') 'PNG'
      BP => BouPoints
      call drawBoundary(BP,nbpts,titleN, extension,.false.,.false.,.false.,.true.)
      !  -framerate #   antes de -i para hacerlo más lento. Donde # es menor a 25 (default)
      write(titleN,'(a)')'ffmpeg -i foto_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p 0_video.mp4'
      write(Printnum,'(/,a,/)') trim(titleN)
      call system(trim(titleN))
!     call chdir("..")
!     call system('cp video/video.mp4 video.mp4')
!     call system('rm video/video.mp4')
      
      if (encuadre - 0.5 .le. 0.1) then
      write(titleN,'(a,I0,a,I0,a)') 'ffmpeg -i 0_video.mp4 -filter:v ''''crop=1200:',&
            int(encuadre*1200+150),':0:',1200-int(encuadre*1200+150),''''' 0_video_Crop.mp4'
      write(Printnum,'(/,a,/)') trim(titleN)
      call system(trim(titleN))
      end if
      end subroutine Churubusco    
      
      subroutine Hollywood(imec)
      use DISLIN
      use glovars, only : workboundary,verbose,Printnum
      use peli, only : ypray => coords_Z, xpray => coords_X,& 
                     fotogramas
      use meshVars, only : npixX,npixZ,nmarkZ,nmarkX
      use waveNumVars, only : NFREC, NPTSTIME
      use soilVars, only : Z,N!, shadecolor_inc
      use geometryvars, only : nXI,Xcoord_ER, Xcoord_Voidonly!, Xcoord_Incluonly
      use waveVars, only : dt,maxtime
      use resultvars, only : Punto,BouPoints,nbpts
      use sourceVars, only : iFte => currentiFte
      use ploteo10pesos !las rutinas externas no requieren interfaz 
      implicit none
      interface
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptoresA,plotReceptoresB,plotFuente)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente
      end subroutine drawBoundary
      end interface
      integer ,intent(in) :: imec
      real :: Sm(npixX,npixZ)
      character(LEN=3), dimension(3) :: nombre
      real     :: ColorRangeMaximumScale, tlabel
      character(LEN=1) :: imdone
      integer*4 :: lentitle
      character(LEN=3) :: extension
      character(LEN=60) :: CTIT
      real             :: maV,miV,p,Vstep,xstep,zstep,encuadre
      real, dimension(2) :: colorBounds
      real*8 :: maxY, minY, maxX, minX
      integer  :: i,j,iT,Iframe,Fframe, n_maxtime
      character(LEN=100) :: textoTimeStamp,path
      real, dimension(41)   :: ZLVRAY
      type (Punto), dimension(:), pointer :: BP
      real*8, dimension(:,:),allocatable :: rec
      
!     factor = sqrt(real(NPTSTIME))
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 'v--'
!     CALL chdir("video")
         
      !tiempo maximo para graficar
         n_maxtime = int(maxtime(iFte)/dt)
         if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
         if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
         write(Printnum,'(a,F5.2,a,F7.4,a,F6.2,a)')&
         "maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," puntos"
        
       ColorRangeMaximumScale = 0.1
  123  maV = maxVal(real(fotogramas(:,:,1:n_maxtime,imec,iFte),4))
       miV = minVal(real(fotogramas(:,:,1:n_maxtime,imec,iFte),4))
       maV = max(maV,abs(miV))
       miV = - maV
       
       IF (verbose .eq. 1)then
         write(Printnum,'(a)') " Used automatic scaling of amplitudes"
         imdone = "Y"
       ELSE
         print *, char(7)
         write(6,'(a,a,a,EN13.2,a,/,a,/,a,EN13.2,/,a)', ADVANCE = "NO") & 
         'Look at the seismograms for ', nombre(imec), & 
         '. Is the response too spiky (the max = ',maV,'? ', &
         'We can enhance detail by reducing the value for maximum color. ', &
         'We propose the maximum to be = ', & 
         ColorRangeMaximumScale * maV, &
         'Proceed [Y] or change it [else] ?'
         read(5,*)imdone
         
         if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
            write(6,'(a)') "keeping those plot limits"
         else
            write(6,'(a)', ADVANCE = "NO")'New maximum for plot = '
            read(5,*) p
            ColorRangeMaximumScale = real(p / maV,4)
         end if
       end if
         
       colorBounds(1) = maV * ColorRangeMaximumScale
       colorBounds(2) = miV * ColorRangeMaximumScale
        
          write(Printnum,'(a,a,a,a,E12.4,a,E12.4,a)') "colorbounds:", & 
          "(",nombre(imec),") [",colorBounds(2)," ; ",colorBounds(1),"]"
       
      minx = minval(xpray)
      maxx = maxval(xpray)
      miny = minval(ypray)
      maxy = maxval(ypray)
      xstep = real(abs(xpray(npixX)-xpray(1))/ nmarkX,4)
      zstep = real(abs(ypray(npixZ)-ypray(1))/ nmarkZ,4)
      encuadre = (ypray(npixZ)-ypray(1))/(xpray(npixX)-xpray(1))
      write(Printnum,'(a,F5.2)')"encuadre=",encuadre      
        
      Iframe = 1
      Fframe = n_maxtime 
      
      if(verbose .eq. 1) then
        imdone = "Y"
      else
        print *, char(7)
        write(6,'(a,I0,a)')'Look at the seismograms, I will plot ',Fframe,&
        'frames. Proceed [Y] or change frame range [else]'
        read(5,*)imdone
        
        if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        Write(6,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
        else
        write(6,'(a)', ADVANCE = "NO")'Start frame = '
        read(5,*)Iframe
        write(6,'(a)', ADVANCE = "NO")'Final frame = '
        read(5,*)Fframe
        Write(6,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
        end if
      end if!
      
      
      
      do iT=Iframe,Fframe !cada fotograma
      do i=1,npixX
      do j=1,npixZ
        Sm(i,j) = real(fotogramas(j,i,iT,imec,iFte),4)
      if (abs(Sm(i,j)) .le. mav*0.01) Sm(i,j) = 0.0
      end do
      end do
      write(textoTimeStamp,'(a,I0,a)') nombre(iMec), iT,'.png'
      
      maV = colorBounds(1)
      miV = colorBounds(2)
      Vstep = (maV-miV)/40.0
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      ! shaded contour plot
      CALL METAFL('PNG') !'PDF'
!     print*,'will print ', trim(textoTimeStamp)
      CALL SETFIL(trim(textoTimeStamp))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(3000,4),int(3000,4))
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(1200,4)) !1200,800
      CALL SCRMOD('REVERS') !fondo blanco
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL DISALF() !default font
           !the position of an axis system.
      CALL axspos (int(300,4) ,int(2700,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600*encuadre,4)) !size of the axis system.
      call name('X [m]','X')
      call name('Z [m]','Y')
      
      !number of decimal places for labels
      call labdig(int(2,4),'Z') 
      if (xstep .gt. 10) then; call labdig(int(-1,4),'X')
        else; call labdig(int(1,4),'X'); end if!
      if (zstep .gt. 10) then; call labdig(int(-1,4),'Y')
        else; call labdig(int(1,4),'Y'); end if!
      
      call setgrf("NAME", "NAME", "LINE", "LINE")      
      CALL LABELS ('EXP ','Z')
!     CALL SETVLT ('SPEC') !rainbow color table
      CALL SETVLT ('GREY') !Sh Grey color table
      call nobar
      call height(40) ! de los caracteres
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4),& 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
                 
      CALL CONSHD(real(xpray,4),int(npixX,4),real(ypray,4),int(npixZ,4), & 
                real(real(Sm(:,:)),4), real(ZLVRAY,4),int(41,4))
      
      !estratigrafía --------------------------------------------
      call shdpat(int(16,4))                                    !
      i = 1
      if (Z(0) .lt. 0.0) i = 0
      do J=i,N                                                  !
         call color ('FORE')                                    !
         call rline(real(minX,4),real(max(Z(J),minY),4), &      !
                 real(maxX,4),real(max(Z(J),minY),4))           !
      end do                                                    !
      J = N+1                                                   !
      call color ('FORE')                                       !
      call rline(real(minX,4),real(Z(J),4), &                   !
                 real(maxX,4),real(Z(J),4))                     !
                 
      ! Borrar estratos en la cuenca                                  !
      call color ('BACK')                                             !
      call shdpat(int(16,4))                                          !
      if (abs(z(0)) .lt. 0.0001) then
      call rlrec(real(minX,4),real(minY,4),&                          !
                    real(maxX-minX,4),real(Z(1)-minY,4))              !
      end if!
      if (workboundary) then !                                         !
!     if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
!     allocate(rec(2* (size(Xcoord_Incluonly(:,1,1))),2))             !
!     i=1                                                             !
!     do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont    !
!     rec(i,1) = Xcoord_Incluonly(j,1,1)                              !
!     rec(i,2) = Xcoord_Incluonly(j,2,1)                              !
!     rec(i+1,1) = Xcoord_Incluonly(j,1,2)                            !
!     rec(i+1,2) = Xcoord_Incluonly(j,2,2)                            !
!     i=i+2                                                           !
!     end do                                                          !
!     print*, "shadecolor_inc", shadecolor_inc
!     call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
!     CALL RLAREA(real(rec(:,1),4), & 
!                 real(rec(:,2),4), & 
!                 int(2*size(Xcoord_Incluonly(:,1,1)),4))  !
!     deallocate(rec) 
!     end if!
      if (size(Xcoord_Voidonly(:,1,1)) .gt. 1) then
      allocate(rec(2* (size(Xcoord_Voidonly(:,1,1))),2))              !
      i = 1                                                           !
      do j=1,size(Xcoord_Voidonly(:,1,1))                             !
      rec(i,1) = Xcoord_Voidonly(j,1,1)                               !
      rec(i,2) = Xcoord_Voidonly(j,2,1)                               !
      rec(i+1,1) = Xcoord_Voidonly(j,1,2)                             !
      rec(i+1,2) = Xcoord_Voidonly(j,2,2)                             !
      i=i+2                                                           !
      end do                                                          !
      call color ('BACK')                                             !
      call shdpat(int(16,4)) 
      CALL RLAREA(real(rec(:,1),4), & 
             real(rec(:,2),4), & 
              int(2*size(Xcoord_Voidonly(:,1,1)),4))                  !
      deallocate(rec)         
      end if                                        !
      ! dibujar topografia original                                   !
      call color ('FORE')                                             !
      call PENWID(real(0.5,4))                                        !
      call marker(int(-1,4)) ! sin marcadores                         !
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do                                                          !
      end if
      
      call color ('FORE')
      tlabel = (it)*real(dt,4)
      write(CTIT,'(a,F9.5,a)') 't=',tlabel,' seg'
      lentitle = NLMESS(CTIT)
      CALL MESSAG(CTIT,int(300,4),int(2850,4))
      
      CALL HEIGHT(int(150,4))
      CALL DISFIN 
!     stop "killed video 7609"
      end do !iT
      
      !now make the video with the frames
      write(path,'(a,a)') nombre(iMec),'0.png'
      write(extension,'(a)') 'PNG'
      BP => BouPoints
      call drawBoundary(BP,nbpts, path, extension,.false.,.false.,.false.,.true.)
      !  ffmpeg -framerate 15 -i foto_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4
      write(path,'(a,a,a)')'ffmpeg -i ',nombre(iMec), & 
                  '%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4'
      call system(trim(path))
!     call chdir("..")
      if (encuadre - 0.5 .le. 0.1) then
      ! ffmpeg -i video/video.mp4 -filter:v 'crop=1200:700:0:600' videoCr.mp4
!     write(path,'(a,a,a)') & 
!     'ffmpeg -i video/video.mp4 -filter:v ''''crop=1200:700:0:600'''' ',& 
!     nombre(iMec),'.mp4'
!     call system(trim(path))
      
      write(path,'(a,a,I0,a,I0,a,a,a)') 'ffmpeg -i video.mp4 -filter:v ','''crop=1200:',&
            int(encuadre*1200+150),':0:',1200-int(encuadre*1200+150),''' 0_',nombre(iMec),'.mp4'
      write(Printnum,'(a)') trim(path)
      call system(trim(path))
      else
      write(path,'(a,a,a)') 'cp video.mp4 0_',nombre(iMec),'.mp4'
!     call system(trim(path)) 
      end if
!     call system('rm video/video.mp4')
      
      if (verbose .ge. 2) then
      print *, char(7)
      write(6,'(a)') 'video is done. Do you want to change it?'
      write(6,'(a)', ADVANCE = "NO") 'replot [Y] , no I am done [else]: '  
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        go to 123
      end if
      print *, char(7)
      write(6,'(a)') 'Keep *.png files? [Y]'
      read(5,*)imdone
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        write(6,'(a)') 'kept the files'
      else 
!       call chdir("video")
!       call system('rm *.png')
!       call chdir("..")
      end if!
      
        write(6,'(a)') 'Keep *.txt files? [Y]'
        read(5,*)imdone
        if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
          write(6,'(a)') 'kept the files'
        else
!         call chdir("video")
!         call system('rm *.txt')
!         call chdir("..")
        end if
      end if !verbose 2
      end subroutine Hollywood           
      
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente)
      
      use DISLIN
      use soilVars, only : Z,N,col=>layershadecolor, shadecolor_inc
      use resultVars, only : Punto, IP => allpoints, nIpts ,overDeterminedSystem
      use sourceVars, only : Po,ifuente=>currentiFte!, tipoFuente, PW_pol
      use geometryvars, only : nXI,Xcoord_ER,normXI, & 
                               midPoint, Xcoord_Voidonly, Xcoord_Incluonly,Xcoord_flip_out
      use glovars, only : verbose, workBoundary,flip12
      use meshVars, only : MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ,nmarkZ,nmarkX
      
      implicit none
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptoresA, plotReceptoresB,plotFuente
      real :: maxY,minY,maxX,minX,xstep,zstep,encuadre
      integer :: i,J,sc
      real*8, dimension(:,:),allocatable :: rec
      real*8,dimension(2) :: theta,gamma
      if (verbose >= 4) Write(6,'(a)') "Will plot geometry..."
      sc=3
      ! Boundary boundaries
      maxX = MeshMaxX
      minX = MeshMinX
      maxY = MeshMaxZ
      minY = MeshMinZ
      
      xstep = real(((maxX-minX) / nmarkX ))
      zstep = real(((maxY-minY) / nmarkZ ))
      encuadre = (maxY-minY)/(maxX-minX)
      ! Dislin plotting routines 
      CALL METAFL(extension) ! define intended display  XWIN,PS,EPS,PDF
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(1200,4)) !1200 ambos
      CALL SCRMOD('REVERS') !fondo blanco
      
      CALL SETFIL(trim(titleN))
      call filmod('DELETE') ! para sobre escribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(3000,4),int(3000,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
      call errmod ("all", "off")
!     CALL PAGERA() ! plots a page border
!     CALL COMPLX ! sets a complex font
      call incmrk (int(1,4))
!     CALL HWFONT()
      CALL TRIPLX() !DISALF()
      CALL axspos (int(300,4) ,int(2700,4)) !the position of an axis system. Lower left corner
      call axslen (int(2600,4) ,int(2600*encuadre,4)) !size of the axis system. 
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(2,4),'Y')
      call ticks (int(1,4) ,'XY') 
      !            low X   left Y   upp X   right Y
      call setgrf("NAME", "NAME", "NONE", "LINE")
      CALL SETVLT ('SPEC')
      ! increment for labels 
      ! xini xfin yfirstvalue xstep ymin ymax yminvalueTicks yTicks
      call graf(real(minX,4),real(maxX,4),real(minX,4),& 
                real(xstep,4),real(maxY,4),real(minY,4),& 
                real(maxY,4),real(-zstep,4))
      
      !estratigrafía --------------------------------------------
      call shdpat(int(16,4)) !solid shade                       !
      i = 1
      if (Z(0) .lt. 0.0) i = 0
      do J=i,N                                                  !
         call SETRGB(col(J), col(J), col(J))                    !
         call rlrec(real(minX,4),real(max(Z(J),minY),4),&       !
                    real(maxX-minX,4),real(Z(J+1)-max(Z(J),minY),4))!
         call color ('FORE')                                    !
         call rline(real(minX,4),real(max(Z(J),minY),4), &      !
                 real(maxX,4),real(max(Z(J),minY),4))           !
      end do                                                    !
      J = N+1                                                   !
!     print*, "shade_layer",J," ",col(J)
      call SETRGB(col(J), col(J), col(J))                       !
      call rlrec(real(minX,4),real(Z(J),4),&                    !
                    real(maxX-minX,4),real(maxY-Z(J),4))        !
      call color ('FORE')                                       !
      call rline(real(minX,4),real(Z(J),4), &                   !
                 real(maxX,4),real(Z(J),4))                     !
      ! Borrar estratos en la cuenca                                  !
      call color ('BACK')                                             !
      call shdpat(int(16,4))                                          !
      if (abs(z(0)) .lt. 0.0001) then
      call rlrec(real(minX,4),real(minY,4),&                          !
                    real(maxX-minX,4),real(Z(1)-minY,4))              !
      end if!
      if (workboundary) then
      !topografía -----------------------------------------------------
      
      ! dibujar inclusión
      if (flip12) then ! Xcoord_flip_out
      if (allocated(Xcoord_flip_out)) then
      if (size(Xcoord_flip_out(:,1,1)) .gt. 1) then
!     do j = 1,8
!        print*,j,Xcoord_flip_out(j,1,1:2),Xcoord_flip_out(j,2,1:2)
!     end do
      allocate(rec(2* (size(Xcoord_flip_out(:,1,1))),2))              !
      i=1                                                             !
      do j=1,size(Xcoord_flip_out (:,1,1))!n_topo+1,n_topo+n_cont    !
      rec(i,1) = Xcoord_flip_out(j,1,1)                              !
      rec(i,2) = Xcoord_flip_out(j,2,1)                              !
      rec(i+1,1) = Xcoord_flip_out(j,1,2)                            !
      rec(i+1,2) = Xcoord_flip_out(j,2,2)                            !
!     print*,j,rec(i,1),rec(i,2),rec(i+1,1),rec(i+1,2)
      i=i+2                                       !
      end do
      call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
      CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_flip_out(:,1,1)),4))  !
      deallocate(rec) 
      end if
      end if
      else ! --------
      if (allocated(Xcoord_Incluonly)) then
      if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
      allocate(rec(2* (size(Xcoord_Incluonly(:,1,1))),2))             !
      i=1                                                             !
      do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont    !
      rec(i,1) = Xcoord_Incluonly(j,1,1)                              !
      rec(i,2) = Xcoord_Incluonly(j,2,1)                              !
      rec(i+1,1) = Xcoord_Incluonly(j,1,2)                            !
      rec(i+1,2) = Xcoord_Incluonly(j,2,2)                            !
      i=i+2                                                           !
      end do
      
!     print*, "shadecolor_inc", shadecolor_inc
      call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
      CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))  !
      deallocate(rec) !
      end if                                                    !
      end if
      end if!
      if (size(Xcoord_Voidonly(:,1,1)) .gt. 1) then
      allocate(rec(2* (size(Xcoord_Voidonly(:,1,1))),2))              !
      i = 1                                   
      do j=1,size(Xcoord_Voidonly(:,1,1))                             !
      rec(i,1) = Xcoord_Voidonly(j,1,1)                               !
      rec(i,2) = Xcoord_Voidonly(j,2,1)                               !
      rec(i+1,1) = Xcoord_Voidonly(j,1,2)                             !
      rec(i+1,2) = Xcoord_Voidonly(j,2,2)                             !
      i=i+2                                                           !
      end do                                                          !
      call color ('BACK')                                             !
      call shdpat(int(16,4)) 
      CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))  !
      deallocate(rec) 
      end if                                                          !
                                                      !                                                                !
      ! dibujar topografia original                                   !
      if (zoomGeom) then
      call color ('RED')                                             !
      call PENWID(real(3.0,4))                                        !
      else
      call color ('FORE')                                             !
      call PENWID(real(0.5,4))                                        !
      end if
      call marker(int(-1,4)) ! sin marcadores                         !
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do                                                          !
      
      if (zoomGeom) then
      !normales -------------------------------------------------------------
      call color('CYAN')                                                  !
      CALL HSYMBL(int(25,4)) !size of symbols                               ! 
      do j=1,nXI                                                            !
      CALL RLVEC (real(midPoint(j,1),4), real(midPoint(j,2),4), &           !
              real(midPoint(j,1)+normXI(j,1)* xstep*0.1,4), &                !
              real(midPoint(j,2)+normXI(j,2)* xstep*0.1,4), int(1001,4))     !
      end do                                                               !
      if (overDeterminedSystem) then                                    !
      do j=1,nipts                                                         !
        if (IP(j)%isOD) then                                               !
            call color('RED') 
            CALL RLSYMB (0, real(IP(j)%center%x,4), real(IP(j)%center%z,4))!
            CALL RLVEC (real(IP(j)%center%x,4), real(IP(j)%center%z,4), &  !
               real(IP(j)%center%x + IP(j)%normal%x * xstep*0.1,4), & !
               real(IP(j)%center%z + IP(j)%normal%z * xstep*0.1,4), &  !
              int(1001,4))                                                 !
          end if   
      end do                                                        !
      end if
      
      ! puntos centrales y gaussianos ------------------------------------
      call incmrk(-1) ! -1 : only symbols                                !
      CALL HSYMBL(int(10,4)) !size of symbols                            !
      call color('BACK')                                                 !
      do j=1,nBpts                                                       !
      !             cuadritos rellenos                                   !
      CALL RLSYMB (16, real(BP(j)%bord_A%x,4), real(BP(j)%bord_A%z,4))   !
      end do                                                             !
      do j=1,nBpts                                                       !
      CALL HSYMBL(int(5,4)) !size of symbols                             !
      call color('FORE')                                                 !
      !             circulitos negros                                    !
      CALL RLSYMB (17, real(BP(j)%center%x,4), real(BP(j)%center%z,4))   !
      call color('ORANGE')                                               !
      CALL HSYMBL(int(3,4)) !size of symbols                             !
      call marker(int(0,4)) !tachecitos                                  !
      call curve(real(BP(j)%Gq_xXx_coords(:,1),4), &                     !
                 real(BP(j)%Gq_xXx_coords(:,2),4), &                     !
                 int(size(BP(J)%Gq_xXx_coords(:,1)),4))                  !
      end do !j                                                          !
      end if ! zoomGeom
      end if !workboundary 
      
      if (plotFuente) then
      !fuente ------------------------------------------------------------
      call color('RED')                                                  !
!     do ifuente = 1,Nfuentes                          
      if (Po(ifuente)%tipoFuente .eq. 0) then !puntual
      CALL HSYMBL(int(40,4)) !size of symbols                            !
      CALL RLSYMB (8, real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4))!star     !
      CALL RLVEC (real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4), &            !
              real(Po(ifuente)%center%x + Po(ifuente)%normal%x * xstep*0.4,4), &           !
              real(Po(ifuente)%center%z + Po(ifuente)%normal%z * xstep*0.4,4), int(1101,4))!
      elseif (Po(ifuente)%tipoFuente .eq. 1) then !onda plana
      ! polarización
      if (Po(ifuente)%PW_pol .eq. 1) then !SV
        theta(1) = cos(Po(ifuente)%gamma)
        theta(2) = sin(Po(ifuente)%gamma)
      elseif (Po(ifuente)%PW_pol .eq. 2) then !P
        theta(1) = sin(Po(ifuente)%gamma)
        theta(2) = -cos(Po(ifuente)%gamma)
      end if
        gamma(1) = sin(Po(ifuente)%gamma)
        gamma(2) = -cos(Po(ifuente)%gamma)
      CALL HSYMBL(int(60,4)) !size of symbols
      CALL RLVEC (real(Po(ifuente)%center%x + theta(2) * xstep*0.05,4), & 
                  real(Po(ifuente)%center%z - theta(1) * xstep*0.05,4), &      !
              real(Po(ifuente)%center%x + theta(2) * xstep*0.05 & 
                                        + theta(1) * xstep*0.2,4), &           !
              real(Po(ifuente)%center%z - theta(1) * xstep*0.05 &
                                        + theta(2) * xstep*0.2,4), int(3101,4))!
      ! dirección de propagación
      CALL RLVEC (real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4), &!
              real(Po(ifuente)%center%x + gamma(1) * xstep*0.4,4), &           !
              real(Po(ifuente)%center%z + gamma(2) * xstep*0.4,4), int(1111,4))!
      ! frente de onda
      CALL RLVEC (real(Po(ifuente)%center%x - gamma(2) * xstep*0.25,4),&        !
                  real(Po(ifuente)%center%z + gamma(1) * xstep*0.25,4), &       !
              real(Po(ifuente)%center%x + gamma(2) * xstep*0.25,4), &           !
              real(Po(ifuente)%center%z - gamma(1) * xstep*0.25,4), int(1100,4))!
      elseif (Po(ifuente)%tipoFuente .eq. 2) then !segmento
        gamma(1) = sin(Po(ifuente)%gamma)
        gamma(2) = -cos(Po(ifuente)%gamma)
      CALL HSYMBL(int(60,4)) !size of symbols
      ! segmento
      CALL RLVEC (real(Po(ifuente)%center%x - gamma(2) * Po(ifuente)%length/2.,4),&        !
                  real(Po(ifuente)%center%z + gamma(1) * Po(ifuente)%length/2.,4), &       !
              real(Po(ifuente)%center%x + gamma(2) * Po(ifuente)%length/2.,4), &           !
              real(Po(ifuente)%center%z - gamma(1) * Po(ifuente)%length/2.,4), int(1100,4))!
      CALL RLVEC (real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4), &            !
              real(Po(ifuente)%center%x + Po(ifuente)%normal%x * xstep*0.3,4), &           !
              real(Po(ifuente)%center%z + Po(ifuente)%normal%z * xstep*0.3,4), int(1101,4))!
      end if !tipoFuente
!     end do
      end if!
      if (plotReceptoresA) then
      !receptores ----------------------------------------------------------
      CALL HSYMBL(int(25,4)) !size of symbols                              !
      if (overDeterminedSystem) then                                    !
      do j=1,nipts                                                         !
        if (IP(j)%isOD) then                                               !
            call color('RED') 
            CALL RLSYMB (0, real(IP(j)%center%x,4), real(IP(j)%center%z,4))!
            CALL RLVEC (real(IP(j)%center%x,4), real(IP(j)%center%z,4), &  !
               real(IP(j)%center%x + IP(j)%normal%x * xstep*0.2,4), & !
               real(IP(j)%center%z + IP(j)%normal%z * xstep*0.2,4), &  !
              int(1001,4))                                                 !
          end if   
      end do  !
      end if
      do j=1,nipts                                                         !
          if (IP(j)%isSeccion) cycle
          call color('BLUE') 
          CALL RLSYMB (2, real(IP(j)%center%x,4), real(IP(j)%center%z,4))  !
          if (IP(j)%atBou) then
            call color('CYAN') 
            CALL RLVEC (real(IP(j)%center%x,4), real(IP(j)%center%z,4), &  !
              real(IP(j)%center%x + IP(j)%normal%x * xstep*0.2,4), &       !
              real(IP(j)%center%z + IP(j)%normal%z * xstep*0.2,4), &       !
              int(1001,4))
          end if                                                           !
      end do                                                               !
      end if!
      if (plotReceptoresB) then
      !receptores ----------------------------------------------------------
      do j=1,nipts                
          if (IP(j)%isSeccion) then                                        !
          CALL HSYMBL(int(10,4)) !size of symbols                          !
          call color('BLUE') 
          CALL RLSYMB (2, real(IP(j)%center%x,4), real(IP(j)%center%z,4))  !
          CALL HSYMBL(int(25,4)) !size of symbols 
          call color('CYAN') 
          CALL RLVEC (real(IP(j)%center%x,4), real(IP(j)%center%z,4), &    !
              real(IP(j)%center%x + IP(j)%normal%x * xstep*0.2,4), &       !
              real(IP(j)%center%z + IP(j)%normal%z * xstep*0.2,4), &       !
              int(1001,4))                                                 !
          else                                                             !
          if (.not. IP(j)%atBou) then
          CALL HSYMBL(int(25,4)) !size of symbols                          !
          call color('BLUE')                                               !
          CALL RLSYMB (2, real(IP(j)%center%x,4), real(IP(j)%center%z,4))  !
          end if ; end if                                                  !
      end do                                                               !
      end if
      call disfin()
      end subroutine drawBoundary
      

      subroutine plot_at_eta (J, tt)
      use resultVars, only : Punto, allpoints, nIpts, nSabanapts
      use dislin
      use debugStuff
      use sourceVars, only : currentiFte
      use gloVars, only : Printnum
      implicit none
      integer :: iP,J
      character(LEN=100) :: tt,t2
      character(LEN=32) :: t3
      complex*16 :: BP
      real, dimension(nSabanapts) :: x
      complex*16, dimension(nSabanapts) :: y
      real*8 :: chNa
      logical :: warning,alguno
      warning = .false.;alguno = .false.
      
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4L')
      CALL SETFIL(trim(tt))
        do iP = 1,nIpts; if (allpoints(iP)%isSabana) then
           alguno = .true.
           if (tt(1:1) .eq. 'W') BP = allpoints(iP)%resp(J,currentiFte)%W
           if (tt(1:1) .eq. 'U') BP = allpoints(iP)%resp(J,currentiFte)%U
           if (tt(1:1) .eq. 'V') BP = allpoints(iP)%resp(J,currentiFte)%V
           chNa = abs(BP)
           if (chNa .eq. chNa+1) then; BP = 0; warning = .true.; end if
!          print*,iP,BP
           y(iP-(nIpts - nSabanapts)) = BP
           x(iP-(nIpts - nSabanapts)) = real(allpoints(ip)%center%x,4)
        end if;end do
        
        write(t2,'(a,a)') trim(tt(1:8)),".m"
        OPEN(739,FILE=trim(t2),FORM="FORMATTED",ACTION='WRITE')
        write(t3,'(a)') trim(tt(1:8))
        CALL scripToMatlabMNmatrixZ(nSabanapts,1,y(1:nSabanapts),t3,739)
!       call showMNmatrixZ(nSabanapts,1,y(1:nSabanapts),t2(1:5),739)
        close(739)
        
        write(t2,'(a,a)') trim(tt),"_.txt"
        if (alguno) then
        write(Printnum,'(a,a)') "plot_at_eta: plot ", trim(t2)
        if (warning) then 
        write(Printnum,'(a,a,/,a)')& 
        "plot_at_eta:",trim(tt), "   One or more values are either NaN or inf"
        else
!       call showMNmatrixZabs(nSabanapts,1, Sabana(1:nSabanapts,1),tt(1:5),6)
        write(t2,'(a,a)') "0__atEta_",trim(tt)
        CALL SETFIL(trim(t2))
        call disini()
        call errmod ("all", "off")
        CALL color('FORE')
        call solid()
        call QPLCRV(x(1:nSabanapts),real(abs(y(1:nSabanapts)),4), nSabanapts,'FIRST')
        
        CALL color('RED')
        call dash()
        call QPLCRV(x(1:nSabanapts),real(real(y(1:nSabanapts)),4), nSabanapts,'NEXT')
        
        CALL color('BLUE')
        call dash()
        call QPLCRV(x(1:nSabanapts),real(aimag(y(1:nSabanapts)),4), nSabanapts,'LAST')
        call disfin

!       write(t2,'(a,a)') "0__Abs_",trim(tt)
!       CALL SETFIL(trim(t2))
!       call disini()
!       call errmod ("all", "off")
!       CALL color('FORE')
!       call qplot(x(1:nSabanapts),real(abs(y(1:nSabanapts)),4), nSabanapts)
!       call disfin
!       write(t2,'(a,a)') "0__Rea_",trim(tt)
!       CALL SETFIL(trim(t2))
!       call disini()
!       call errmod ("all", "off")
!       CALL color('RED')
!       call qplot(x(1:nSabanapts),real(real(y(1:nSabanapts)),4), nSabanapts)
!       call disfin
!       write(t2,'(a,a)') "0__Ima_",trim(tt)
!       CALL SETFIL(trim(t2))
!       call disini()
!       call errmod ("all", "off")
!       CALL color('BLUE')
!       call qplot(x(1:nSabanapts),real(aimag(y(1:nSabanapts)),4), nSabanapts)
!       call disfin
        end if
        end if  
            
      end subroutine plot_at_eta
