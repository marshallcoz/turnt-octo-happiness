      PROGRAM DWNLIBEM 
      use gloVars
      use refSolMatrixVars, only : ak,B
      use waveNumVars
      use soilVars, only : alfa,beta, alfa0,beta0,minbeta,n,z,amu,lambda,rho,qq
      use debugStuff
      use resultVars
      use ploteo10pesos
      use sourceVars, only : Po,psv,sh,tipofuente
      use MeshVars, only : npixX
      use dislin
      use waveVars, only : t0
      use peli, only : fotogramas,fotogramas_region
      use GeometryVars, only : longitudcaracteristica_a
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
        
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptores
      end subroutine drawBoundary
        
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
      !#< blue
      call system('clear')
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .eq. 0) then 
        print*,"usage DWNIBEM.out [option]"
        print*,""
        print*,"     -r        : Correr programa normal."
        print*,"     -g        : Solo imprimir geometría."
        print*,"     -a        : Solo Imprimir funcion de amplitud."
        print*,"     -f [frec] : Solo Ejecutar la frecuencia [frec]."
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
      onlythisJ = .false. !#< b
      IF (LEN_TRIM(arg) .ne. 0) then
        if (trim(arg) .eq. '-f') then
          CALL get_command_argument(2, arg)
          IF (LEN_TRIM(arg) .ne. 0) then
            read(arg,*) frecIni
            if (frecIni .lt. 1) stop "-f lower than 1"
            if (frecIni .gt. NFREC+1) stop "-f higher than NFREC+1 "
            frecEnd = frecIni
            onlythisJ = .true.
          else
            print*," "; print*," "
            print*,"no single frequency specified on argument -f"
            stop
          end if
        end if
      end if !#>
      call checarWisdom(2*nfrec,2*nmax,NPTSTIME) ! FFTw
        nIpts=0; nMpts=0; nBpts = 0; iPtfin = 0; mPtfin = 0
      call getsource(skipdir,PSV,SH)
      if (PSV) write(PrintNum,*) "  this a P-SV case"
      if (SH)  write(PrintNum,*) "  this a SH case"
      if (workBoundary .eqv. .true.) call getTopography
      call getInquirePoints
      if(.not. onlythisJ) call sourceAmplitudeFunction
      call getVideoPoints
      
        Npts = nIpts + nMpts
        write(PrintNum,'(a,I0)') '   Number of fixed receptors: ',Npts
        allocate (allpoints(Npts))
        allpoints(iPtini:iPtfin)= inqPoints(iPtini:iPtfin); deallocate(inqPoints)
      
      call setInqPointsRegions
      call setVideoPointsRegions
 !#< blue
      call chdir(trim(adjustl(rutaOut)))
         write(titleN,'(a)') '0___OriginalGeometry.pdf'
         write(extension,'(a)') 'PDF'
         BP => BouPoints
         call drawBoundary(BP,nbpts,titleN,extension,.false.,.false.)
         
         write(titleN,'(a)') '0___Sensors.pdf'
         write(extension,'(a)') 'PDF'
         BP => BouPoints
         call drawBoundary(BP,nbpts,titleN,extension,.false.,.true.)
         
         if (workBoundary .eqv. .true.) then 
         write(titleN,'(a)') '0___Inclusion.pdf'
         write(extension,'(a)') 'PDF'
         BP => BouPoints
         call drawBoundary(BP,nbpts,titleN,extension,.true.,.false.)
         end if
      call chdir("..")
              
      call get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then 
        if (trim(arg) .eq. '-g') stop "argumento -g : Geometría"
      end if
  !#>
! predimensionamientos
      write(PrintNum,'(a)') &
       "---------------------------------------------------------------------------------"
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
          call loadG_fotogramas
          if (PSV) call Churubusco
          if (SH) call Hollywood(3)
          stop "done"
        end if ! load 
      end if! argument
      
      if (SH) then 
        i = 2*N+1
        if (Z(0) .lt. 0.0) i = i + 1 !HS arriba
      end if !
        if (PSV) then 
        i = 4*N+2
        if (Z(0) .lt. 0.0) i = i + 2 !HS arriba
      end if
      ik = 1 ;  if (tipofuente .eq. 1) ik = 0 !onda plana 
        allocate (Ak(i,i,ik:2*nmax)); allocate (B(i,ik:2*nmax)) 
        allocate(ipivA(i)); allocate(workA((i)*(i)))
   !   to store the strata diffracted displacement: W,U,V,...
      do iP=iPtini,iPtfin
        allocate(allpoints(iP)%W(NFREC+1,3)); allpoints(iP)%W = Z0 !W,U,V
        allocate(allpoints(ip)%resp(NFREC+1)) ! para prescindir de W
        allpoints(ip)%resp(1:NFREC+1)%U = z0
        allpoints(ip)%resp(1:NFREC+1)%V = z0
        allpoints(ip)%resp(1:NFREC+1)%W = z0
        allpoints(ip)%resp(1:NFREC+1)%Tx = z0
        allpoints(ip)%resp(1:NFREC+1)%Ty = z0
        allpoints(ip)%resp(1:NFREC+1)%Tz = z0
         ! Plot integrand FK of Green(l) function
        if(plotFKS) then; if(allpoints(iP)%guardarFK)then
         allocate(allpoints(iP)%FK(NFREC+1,NMAX,3)) !W,U,V
         allpoints(iP)%FK = Z0
        end if;end if
      end do
       if (makeVideo) then
         allpoints(mPtini:mPtfin) = moviePoints; deallocate(moviePoints)
        do iP=mPtini,mPtfin
         allocate(allpoints(iP)%Wmov(NFREC+1,3,npixX)) !W,U,V
         allpoints(iP)%Wmov = z0
        end do
       end if
      ! un indice de punto para todos
      do i=1,Npts
        allpoints(i)%pointIndex = i
      end do ! i
      !#< b
      if (verbose .ge. 2) then;    print*,"inquire points, center,region,layer,movie"
        do iP=1,nPts
          print*,iP,allpoints(iP)%center,allpoints(ip)%region,& 
          allpoints(ip)%layer,allpoints(ip)%guardarMovieSiblings
        end do
      end if !#>
      
      allocate(k_vec(2*nmax))
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
      !#< g
      if (developerfeature .ne. 0) then
      allocate(developerAUXvec(NFREC,2))
      end if !#>
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
      !#< g
      if (developerfeature .ne. 0) then 
      developerAUXvec(J,1) = amu(N+1) * (come/beta(N+1))**2. 
!     if (J .lt. 5) developerAUXvec(J,1) = developerAUXvec(J,1) + 2*UI*OMEI 
      !<- trampita para que no sea un numero muy chico
      developerAUXvec(J,2) = longitudcaracteristica_a * come / alfa(N+1) !ok
!     print*,J,come,developerAUXvec(J,1)
!    cycle
      print*,"J=",J
      print*,"   lamb_alf_0=",alfa(0)/frec
      print*,"   lamb_bet_0=",beta(0)/frec
      print*,"   lamb_alf_1=",alfa(1)/frec
      print*,"   lamb_bet_1=",beta(1)/frec
      print*,"   lamb_alf_2=",alfa(2)/frec
      print*,"   lamb_bet_2=",beta(2)/frec
      print*,"   lamb_alf_3=",alfa(3)/frec
      print*,"   lamb_bet_3=",beta(3)/frec
      print*,"-------------------------"
      end if!#> !developerfeature .ne. 0
      !#< b
      if(verbose .le. 1)then
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),60)
        write(6,'(A)', ADVANCE = "NO") repeat(char(8),17) !eta
        if (workBoundary) write(6,'(A)', ADVANCE = "NO") repeat(char(8),10) !nbp
        write(6,'(A)', ADVANCE = "NO") "["
        write(6,'(A,A)', ADVANCE = "NO") & 
        repeat("X",int((58.0/NFREC)*(NFREC+1-J))),&
        repeat("_",58-int((58.0/NFREC)*(NFREC+1-J)))
        write(6,'(A)', ADVANCE = "NO") "]" 
      else
        call ETIME(tarray, result)
        write(PrintNum,'(A,I0,A,EN18.2,A,EN18.2,A)', ADVANCE = "YES") &
        'w(',J,') | ',FREC," | ",result,"sec"
      end if
      if (onlythisJ) write(PrintNum,*) "Frec = ",COME, " eta = ",&
      (ome*longitudcaracteristica_a)/(pi*beta0(N+1))!#> 
!        print*,"come=",come
!        write(6,'(A,EN10.2,2x)', ADVANCE = "NO") & 
!        "eta(ome)= ", (ome*longitudcaracteristica_a)/(pi*beta0(N+1))
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
         end if
       end if!
       if (i .eq. 1) then
        if ((overDeterminedSystem) .and. & 
            (OD_Jini .le. J) .and. &
            (J .le. OD_Jend)) then
         allocate(ibemMat(ik * (l+n_OD),ik * l));allocate(trac0vec  (ik * (l+n_OD)))
         allocate(IPIVbem(1));                   allocate(auxGvector(ik * (l+n_OD)))
        else
         allocate(ibemMat(ik * l,ik * l));allocate(trac0vec  (ik * l))
         allocate(IPIVbem(ik * l));       allocate(auxGvector(ik * l))
        end if
       end if
       ibemMat = z0 ; trac0vec = z0 ; auxGvector = z0
      end if!workbou
     
!     call ETIME(tarray, result)
!     write(6,'(a,f10.3,a)') "Start= ",result,"seconds" 
!     lastresult = result
!     print*,"antes de hacer e invertir matriz A (todas k)"!,result
         pt_ipivA => ipivA
         pt_workA => workA
      if (PSV) then
         tam = size(Ak,1)
       do ik = 1,vecNK(J) ! k positivo (Aorig -> nmax+1)
         pointAp => Ak(1:tam,1:tam,ik)
         pt_k => k_vec(ik)
         call gloMat_PSV(pointAp,pt_k,pt_cOME_i)
         call inverseA(pointAp,pt_ipivA,pt_workA,tam)
       end do ! ik
      
!      call parImpar_gloMat
!      open(421,FILE= "outA.m",action="write",status="replace")
!      do ik=1,2*NMAX
!      write(arg,'(a,I0,a)') "Aorig{",ik,"}"
!      call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),arg,421)
!      end do
!      close(421)
!      stop     
       
       k0 = vecNK(J)
       call intrplr_gloMat(k0,15,pt_cOME_i,pt_ipivA,pt_workA)         
       call parImpar_gloMat ! k negativo
       
!     call ETIME(tarray, result)
!     write(6,'(a,f10.3,a)') "Elapsed ",result,"seconds"       
!     write(6,'(a,f10.3,a)') "Delta time =",result-lastresult,"seconds"
      
      
!      open(421,FILE= "outA.m",action="write",status="replace")
!      do ik=1,2*NMAX
!      write(arg,'(a,I0,a)') "A{",ik,"}"
!      call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),arg,421)
!      end do
!      close(421)
!      print*,"printed outA.m"
!      stop

      end if!psv
      if (SH) then
         Ak = Z0
      Do ik=1,nmax+1
      ! k positivo
         pointAp => Ak(1:size(Ak,1),1:size(Ak,2),ik)
         pt_k => k_vec(ik)
         call globalmatrix_SH(pointAp,pt_k,pt_cOME_i)
         call inverseA(pointAp,pt_ipivA,pt_workA,size(Ak,1))
      end do ! ik
      ! k negativo (es simétrico)
         Ak(1:size(Ak,1),1:size(Ak,2),Nmax+2:2*nmax) = &
         Ak(1:size(Ak,1),1:size(Ak,2),Nmax:2:-1)
      end if!sh
     
      do dir= 1,3 !x,y,z direction of force application
        if(dir .eq. 2) then
         if(skipdir(dir)) cycle
        else ! 1 o 3
         if(skipdir(1) .and. skipdir(3)) cycle
        end if! dir
       if(.not. skipdir(dir)) then 
         call diffField_at_iz(0,dir,J,cOME)
       end if
       ! fill IBEM matrix
       if (workboundary) then
         if(.not. skipdir(dir)) then 
          if (Po(1)%region .eq. 2) then
           call termIndepR(dir,cOME)
          end if
         end if
      !********(campo difractda en medio estratificado)*********************
         do iz = 1,nZs !por cada Z de fuentes en tabla de segmentos tipo 0 y 1
           if (thereisavirtualsourceat(iz)) then ! (ipxi = 1,n_top_sub)
             call diffField_at_iz(iz,dir,J,cOME)
           end if 
         end do !iz
      !********(campo refractado en inclusión y FF de inclusion (columnas xi=d2))***  
         if (n_con_sub .gt. 0) then 
           do iPxi = n_top_sub +1, n_top_sub + n_con_sub
             call reffField_by_(iPxi,dir,cOME)
           end do ! iPxi
         end if !n_cont
      !********(frontera libre en inclusión (columnas de xi=d1))*****************
         if (n_val_sub .gt. 0) then
           do iPxi = n_top_sub + n_con_sub +1, n_top_sub + n_con_sub + n_val_sub
             call reffField_by_(iPxi,dir,cOME)
           end do ! iPxi
         end if !n_vall
      !*****(funcions de Green de desplazamiento en la región R)*****************
         call GreenReg_R(J,dir,cOME)
       end if !workboundary
      end do !dir
     
      if (workboundary) then  !; print*,"solve ibem"
!     call showMNmatrixZabs(2*nBpts,1, trac0vec,"  t0 ",PrintNum)
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
!     if (verbose .ge. 1) call showMNmatrixZ(Mi,Ni, ibemMat ," mat ",6)
!     if (verbose .ge. 1) call showMNmatrixZ(Mi,1 , trac0vec,"  b  ",6) ;stop
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
      call zgesv(Mi,1,ibemMat,Mi,IPIVbem,trac0vec,Mi,info)
      
      if(info .ne. 0) stop "problem with ibem system"
      if (any(isnan(real(trac0vec)))) stop "487 valió madres el ibem"
      end if !#< r -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- !#>
      !#< b
!     call showMNmatrixZ(Mi,1, trac0vec,"phi  ",6)
      if (verbose .ge. 2) then
         call chdir(trim(adjustl(rutaOut)))
         CALL chdir("phi")
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
          if (allpoints(iP_x)%region .eq. 1) then !'estr'
            iPhi_I = 1
!           iPhi_F = (n_top_sub + n_con_sub) * 2
            ipxi_I = 1
            ipxi_F =  n_top_sub + n_con_sub
          else if (allpoints(iP_x)%region .eq. 2) then !'incl'
            iPhi_I = (n_top_sub + n_con_sub) * 2 + 1
!           iPhi_F = (n_top_sub + 2* n_con_sub + n_val_sub) * 2
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
          if(i .eq. 1) allpoints(iP_x)%resp(J)%W = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J)%W
          if(i .eq. 2) allpoints(iP_x)%resp(J)%U = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J)%U
          
          l = i + 3 !componente de la tracción {Tz,Tx}
          auxGvector = z0
          iPhi = iPhi_I 
          do iPxi = ipxi_I,ipxi_F ! recopilamos  G_ij
            auxGvector(iPhi)   = boupoints(iPxi)%G(iP_X,l,1) ! por fzas horizontales:
            auxGvector(iPhi+1) = boupoints(iPxi)%G(iP_X,l,3) ! por fzas verticales:
            iPhi = iPhi + 2 
          end do
          if(l .eq. 4) allpoints(iP_x)%resp(J)%Tz = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J)%Tz
          if(l .eq. 5) allpoints(iP_x)%resp(J)%Tx = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                       allpoints(iP_x)%resp(J)%Tx
      end do !i
      end do !iP_x
      
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
          fotogramas(iP_x-nIpts,m,J,i) = sum(trac0vec(1:Mi) * auxGvector(1:Mi)) + &
                                         allpoints(iP_x)%Wmov(J,i,m)
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
      !#< b
      if (verbose .ge. 3) call showMNmatrixZ(Mi,Ni, ibemMat," mat ",6)
      if (verbose .ge. 3) call showMNmatrixZ(Mi,1, trac0vec,"  b  ",6)
      !#>
      call zgesv(Mi,1,ibemMat,Mi,IPIVbem,trac0vec,Mi,info)
      if(info .ne. 0) stop "problem with ibem system"
      
      if (any(isnan(real(trac0vec)))) then ! NAN is not equal even to itself
      stop "2120 valio madres el ibem"; end if!
      !#< b
      if (verbose .ge. 2) then  
         call chdir(trim(adjustl(rutaOut))) 
         CALL chdir("phi")
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
      do i= 3,6,3!dirección de desplazamient V,Ty !SH
         auxGvector = z0 !(1:ik)
!        print*, iPhi_I,iPhi_F, " -- ", ipxi_I,ipxi_F
         iPxi = ipxi_I
         do iPhi = iPhi_I,iPhi_F
           auxGvector(iPhi) = boupoints(iPxi)%G(iP_X,i,2)
           iPxi = iPxi + 1
         end do !iPhi         
          if (verbose .ge. 4) call showMNmatrixZ(Mi,1,auxGvector," auxG",6)
          
          if (i .eq. 3) allpoints(iP_x)%W(J,i) = sum(trac0vec * auxGvector) + &
          allpoints(iP_x)%W(J,i) 
          
          if (i .eq. 3) allpoints(iP_x)%resp(J)%V = allpoints(iP_x)%resp(J)%V + &
                        sum(trac0vec * auxGvector)
          if (i .eq. 6) allpoints(iP_x)%resp(J)%Ty = allpoints(iP_x)%resp(J)%Ty + &
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
          fotogramas(iP_x-nIpts,m,J,i) = sum(trac0vec * auxGvector) + &
          allpoints(iP_x)%Wmov(J,i,m) 
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
          fotogramas(iP_x-nIpts,m,J,i) = allpoints(iP_x)%Wmov(J,i,m) 
      end do !i
      end do !iP_x
      end do !m
      end if !PSV
      if (SH) then
      do m = 1,npixX
      do iP_x = mPtini,mPtfin
          fotogramas(iP_x-nIpts,m,J,3) = allpoints(iP_x)%Wmov(J,3,m) 
      end do !iP_x
      end do !m
      end if !SH
      end if !makevideo
      end if !workboundary
      if (onlythisJ) then
          print*,"";print*,"Rultados en puntos receptores: ------"
          do iP_x = 1,nIpts
           if (.not. allpoints(iP_x)%isSabana) then
            print*,"ip",ip_x,"[",allpoints(iP_x)%center%x,",",allpoints(iP_x)%center%z,"]", &
            "n = [",allpoints(iP_x)%normal%x,",",allpoints(iP_x)%normal%z,"]"
            print*,"   W=",allpoints(iP_x)%resp(J)%W
            print*,"   U=",allpoints(iP_x)%resp(J)%U
            print*,"   Tz=",allpoints(iP_x)%resp(J)%Tz
            print*,"   Tx=",allpoints(iP_x)%resp(J)%Tx 
           end if
          end do !iP_x
        end if!
!     call ETIME(tarray, result)
!     write(PrintNum,'(a,f10.3,a)') "Elapsed ",result,"seconds"
!     print*,"deltaT=",result-lastresult
!     stop 2218
      END DO ! J: frequency loop
      
!     deallocate(B);deallocate(IPIV)
      if(verbose >= 1) write(PrintNum,'(a)')" done"
      ! w->time & plot 
      call ETIME(tarray, result)
      if (result .ge. 60) then
      write(PrintNum,'(a,f10.3,a)') "Elapsed ",result/60,"minutes"
      write(6,'(a,f10.3,a)') "Elapsed ",result/60,"minutes"
      else
      write(PrintNum,'(a,f10.3,a)') "Elapsed ",result,"seconds"
      write(6,'(a,f10.3,a)') "Elapsed ",result,"seconds"
      end if
      
      write(6,'(a)')"Printing seismograms, etc..."
      if (plotFKS .and. (tipoFuente .ne. 1)) then  
           write(xAx,'(a)')"frec [Hz]"
           write(yAx,'(a)')" K [1/m] "
           call chdir(trim(adjustl(rutaOut)))
         do iP = iPtini,iPtfin
           if (allpoints(iP)% guardarFK) then
             if (SH) then
             write(tt,'(a,I0)')"FK_",iP
             call plotFK(allpoints(iP)%FK(1:NFREC+1,1:NMAX,1:3), &
                         real(allpoints(iP)%center%x,4), & 
                         real(allpoints(iP)%center%z,4), & 
                         tt,xAx,yAx,PrintNum,3,3, onlythisJ, frecEnd)
             end if !
             if (PSV) then
             write(tt,'(a,I0)')"FK_",iP
             call plotFK(allpoints(iP)%FK(1:NFREC+1,1:NMAX,1:3), &
                         real(allpoints(iP)%center%x,4), & 
                         real(allpoints(iP)%center%z,4), & 
                         tt,xAx,yAx,PrintNum,1,2, onlythisJ, frecEnd)
             end if
           end if
         end do
           CALL chdir("..")
      end if
      ! sabana en frecuencia eta
      if (verbose .eq. 2) then 
        call chdir(trim(adjustl(rutaOut)))
        CALL chdir("phi")
        write(tt,'(a,I0,a)') "W__eta",J,".pdf"
        call plot_at_eta(frecIni,tt)
        write(tt,'(a,I0,a)') "U__eta",J,".pdf"
        call plot_at_eta(frecIni,tt)
        CALL chdir("..");CALL chdir("..")
        if ( onlythisJ )then 
          print*,""
          print*,"at normalized frequency eta=",abs(come*1.0/(pi*beta(N+1)))
          Stop "-f argum finish"
        end if
      end if
      
      ! mostrar sismogramas en los puntos de interés
           call chdir(trim(adjustl(rutaOut)))
           call chdir("traces")
       if (PSV) then !#< g  
        if (developerfeature .eq. 1) then
        write(yAx,'(a)') 'denominador'
          call W_to_t(developerAUXvec(1:nfrec,1),'r__',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)
        do iP = 2,4,2
        ! tx = sxx = \sigma_{\theta \theta}
          write(yAx,'(a)') '$|\sigma_{\theta \theta}^{*}|$'
          call W_to_t(allpoints(iP)%resp(:)%Tx,'Tx-',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)
        end do!
        do iP = 1,3,2
        ! tz = szz = \sigma_{\theta \theta}
          write(yAx,'(a)') '$|\sigma_{\theta \theta}^{*}|$'
          call W_to_t(allpoints(iP)%resp(:)%Tz,'Tz-',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)
        end do
        end if !#>
        
                !w
        do iP = iPtini,iPtfin
          write(yAx,'(a)') '$u_3$ [m]'
          call W_to_t(allpoints(iP)%resp(:)%W,'w--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)
          call chdir("..")       
          call makeSabana('1_S-w__.pdf',.false.)
          call makeSabana('1_S-w_f.pdf',.true.) ! filled traces
          call chdir("traces")
        end do 
        !u
        do iP = iPtini,iPtfin
          write(yAx,'(a)') '$u_1$ [m]'
          call W_to_t(allpoints(iP)%resp(:)%U,'u--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)
          call chdir("..")         
          call makeSabana('1_S-u__.pdf',.false.)
          call makeSabana('1_S-u_f.pdf',.true.) ! filled traces
          call chdir("traces")
        end do 
        !Tractions
!       do iP = iPtini,iPtfin
!         write(yAx,'(a)') '$Tz_$ [m]'
!         call W_to_t(allpoints(iP)%resp(:)%Tz,'Tz-',yax,iP,& 
!                   allpoints(iP)%center%x,& 
!                   allpoints(iP)%center%z,PrintNum)
!         write(yAx,'(a)') '$Tx_$ [m]'
!         call W_to_t(allpoints(iP)%resp(:)%Tx,'Tx-',yax,iP,& 
!                   allpoints(iP)%center%x,& 
!                   allpoints(iP)%center%z,PrintNum)
!       end do
       call chdir("..") 
       end if !psv
        
       if (SH) then
          do iP = iPtini,iPtfin
          write(yAx,'(a)') '$u_2$ [m]'
          call W_to_t(allpoints(iP)%resp(:)%V,'v--',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)    
          end do 
          call chdir("..") 
          call makeSabana('1_S-v__.pdf',.false.) 
          call makeSabana('1_S-v_f.pdf',.true.)
          
          call chdir("traces")
          do iP = iPtini,iPtfin
          write(yAx,'(a)') '$Ty_$ [m]'
          call W_to_t(allpoints(iP)%resp(:)%Ty,'Ty-',yax,iP,& 
                    allpoints(iP)%center%x,& 
                    allpoints(iP)%center%z)     
          end do
          call chdir("..")
       end if !sh         
            
      if (makeVideo) then 
        call crepa_four_fotogramas
        if (PSV) call Churubusco
        if (SH) call Hollywood(3)
      end if
      call vaciarWisdom
      Write(PrintNum,'(a)') ' done '
      Write(6,'(a)') ' done '
      close(PrintNum)
      
!     if (borrarDirectorio .eqv. .false.) then
!     call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
!     call chdir("..")
!     write(path,'(a,a)') "mv outs outs_",time
!     call system(trim(adjustl(path)))
!     end if
      END program
      
      
      ! SET UP & READ 
      subroutine setupdirectories
      use glovars, only : borrarDirectorio, rutaOut,PrintNum
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
      write(path,'(a,a)') "mkdir outs_",time
      call system(trim(adjustl(path)))
      write(rutaOut,'(a,a)') "outs_",time
      call chdir(trim(adjustl(rutaOut)))
      else
      call system('mkdir outs')
      CALL chdir("outs",status)
      write(rutaOut,'(a)') "outs"
      end if!
      if (status .eq. 0) call chdir("..") !workdir
      write(path,'(a,a,a)') 'cp -r ins ',trim(adjustl(rutaOut)),'/insbackup'
      call system(trim(adjustl(path)))
      call chdir(trim(adjustl(rutaOut)),status) !outs
      if (status .eq. 0) call system("rm *.*")
      if (PrintNum /= 6) open(PrintNum,FILE= "logfile.txt")
      call date_and_time(TIME=time); write(PrintNum,'(a,a)') "hhmmss.sss = ",time
      call system('mkdir phi')
      CALL chdir("phi",status)
      if (status .eq. 0) call system("rm *.*")
      if (status .eq. 0) call chdir("..")
      call system('mkdir traces')
      CALL chdir("traces",status)
      if (status .eq. 0) call system("rm *.txt")
      if (status .eq. 0) call system("rm *.pdf")
      if (status .eq. 0) call chdir("..")
      call system('mkdir subdivs')
      CALL chdir("subdivs",status)
      if (status .eq. 0) call system("rm *.*")
      if (status .eq. 0) call chdir("..")
      call system('mkdir video')
      call chdir("video",status)
      if (status .eq. 0) call system("rm *.png")
      if (status .eq. 0) call system("rm *.txt")
      if (status .eq. 0) call chdir("..")
      call system('cp ../../DWNLIBEM.f90 DWNLIBEM.f90')
      call chdir("..") 
      write(PrintNum,'(///)') 
      end subroutine setupdirectories
      
      subroutine getMainInput
      use glovars
      use GeometryVars , only : staywiththefinersubdivision,finersubdivisionJ,&
       longitudcaracteristica_a, fraccionDeSmallestWL_segm_de_esquina
      use resultvars, only : overDeterminedSystem,OD_Jini,OD_Jend
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
      READ(35,'(L1)') workBoundary!; print*,"boundary? ",workBoundary
      READ(35,'(L1)') flip12; if(.not. workBoundary) flip12 = .false.
      READ(35,'(L1)') plotFKS!; print*,"plotFK?",plotFKS
      READ(35,'(L1)') saveG!; print*,"Save Green funcs?", saveG
      READ(35,*) multSubdiv!; print*,"division multiple = ", multSubdiv
      READ(35,*) staywiththefinersubdivision, finersubdivisionJ
      READ(35,*) overDeterminedSystem,OD_Jini,OD_Jend
      READ(35,*) cKbeta!; print*,"multBminIntegrando = ", cKbeta
      read(35,*) periodicdamper!; print*,"Periodic sources damping factor = ", periodicdamper
      READ(35,*) useAzimi
      read(35,*) developerfeature
      read(35,*) borrarDirectorio
      read(35,*) PrintNum
      read(35,*) longitudcaracteristica_a
      read(35,*) fraccionDeSmallestWL_segm_de_esquina
      close(35)    
      CALL chdir("..")
      end subroutine getMainInput
      
      subroutine getSoilProps
      use soilVars; use waveNumVars; use waveVars, only : dt, maxtime; use fitting
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
       ' -inf.  - ',Z(1),ALFA0(0),BETA0(0),real(AMU(0)),  RHO(0), real(LAMBDA(0)),real(anu(J))!#>
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
          Z(J),' - ',Z(J+1),ALFA0(J),BETA0(J),real(AMU(J)),& 
          RHO(J), real(LAMBDA(J)),real(anu(J)) !#>
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
         ' -inf.  - ',Z(1),ALFA0(0),BETA0(0),real(AMU(0)),&
         RHO(0), real(LAMBDA(0)),real(anu(0)) !#>
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
         Z(1),' -  inf. ',ALFA0(N+1),BETA0(N+1),real(AMU(N+1)),&
         RHO(N+1), real(LAMBDA(N+1)),real(anu(N+1)) !#>
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
      READ(7,*)DFREC,NFREC,NPTSTIME,nK,Qq,TW,maxtime
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
      ! el resto hasta nmax+1 se interpola
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
       write(outpf,'(a,F9.7,/,a,F15.5)') "   dt = ",Dt,"   Tmax=",Dt* NPTSTIME
       write(outpf,'(a,F8.1)') '   Atenuation Q = ',Qq
       write(outpf,'(a,I0,a,F8.4,a,F12.4,a,/)') & 
           '   N. frequencies: ',NFREC,'  @',DFREC,'Hertz :. Fmax = ', & 
           NFREC*DFREC,'Hertz'
      
       write(outpf,'(a)') &
       "--- DWN -------------------------------------------------------------------------"
       write(outpf,'(a,F9.7)') "   DK = ",DK
       write(outpf,'(a,I0,a,I0,a,I0,a,I0,a,I0,a,I0,a)') "   nk a 3pt spline with (",& 
          1,",",int(0.55*nK),"), (",& 
          int(nfrec * frac)+1,",",int(0.6*nk),"), (",& 
          NFREC+1,",",nk,")"
       write(outpf,'(a,I0)') '   nk average = ',int(sum(vecNK)/(NFREC+1))
       write(outpf,'(a,I0)') '   nmax (each sign): ',NMAX
       write(outpf,'(a,F12.7)') "   delta X = ", real(pi / (nMax * DK),4)
       write(outpf,'(a,EN14.2E2,a)') "   L = ",2*pi/DK, "m"
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
      use resultVars, only : inqPoints, nIpts, iPtini,iPtfin, & 
                       nSabanapts, Sabana,sabZeroini,sabZerofin,&
                       SabanaPlotIndividual, sabanaBajarAFrontera, & 
                       n_OD,overDeterminedSystem
      use waveNumVars, only : NPTSTIME
      implicit none
      integer :: i,j,k,thelayeris,iIndex, nnsabanas,thisnsab
      logical :: lexist, tellisoninterface, ths_isoninterface,wtfk
      integer :: auxGuardarFK, ths_layer, sabanabajarmax
      real :: xini,deltax,zini,deltaz,dx,dz, SbanadeltaZ
      real*8 :: escalax,escalay,offsetx,offsety
      logical :: cn,adentroOafuera
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
      READ(7,*) nnsabanas, nSabanapts, & 
                SabanaPlotIndividual, SbanadeltaZ, sabanabajarmax
      READ(7,*) !__________________________________________
      READ(7,*) escalax,escalay
      READ(7,*) offsetx,offsety
      iPtini = 1
      if (nnsabanas .eq. 0) nSabanapts=0
      iPtfin = nIpts + nSabanapts
      allocate(inqPoints(nIpts + nSabanapts))
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
      READ(7,*) !  X        Z          nx       nz     guardarFK
      do i=1, nIpts 
         READ(7,*) inqPoints(i)%center%x, inqPoints(i)%center%z, &
             inqPoints(i)%normal%x, inqPoints(i)%normal%z, auxGuardarFK, & 
             inqPoints(i)%isOD,inqPoints(i)%tipoFrontera
           inqPoints(i)%center%x = inqPoints(i)%center%x * escalax + offsetx
           inqPoints(i)%center%z = inqPoints(i)%center%z * escalay + offsety
           IF (inqPoints(i)%isOD) n_OD = n_OD +1
      !encontrar el layer en el que estan o 0 si está sobre la interfaz
           inqPoints(i)%layer = thelayeris(real(inqPoints(i)%center%z,8))
           if (auxGuardarFK .ge. 1 ) inqPoints(i)%guardarFK = .true.
           inqPoints(i)%isOnInterface = & 
                         tellisoninterface(real(inqPoints(i)%center%z,8))                         
      end do
      if (.not. overDeterminedSystem) n_OD = 0        ! en caso de que no hay puntos
      if (n_OD .eq. 0) overDeterminedSystem = .false. ! de coloc. adicionales  
      iIndex = nIpts
      if (nSabanapts .gt. 0) then
      allocate(Sabana(nSabanapts, NPTSTIME))
      read(7,*) !Sabanapoints -------------------------
      read(7,*) escalax,escalay
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
      close(7)
      nIpts = nIpts + nSabanapts
      end if !sabanas
      CALL chdir("..")
      end subroutine getInquirePoints   
      
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
      real*8 ::  errT = 0.005_8
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
      
      subroutine getsource(skipdir,PSV,SH)
      use wavevars, only: Escala,Ts,Tp, ampfunction, sigGaus, t0
      use sourceVars, only: Po, tipoFuente, PW_pol, nFuentes
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
      
      logical :: skipdir(3),PSV,SH
      integer :: thelayeris,efsource,i,regi
      logical :: lexist, tellisoninterface, intfsource
      real    :: xfsource,zfsource,l
      real*8  :: nxfsource,nyfsource,nzfsource, PW_theta
      type (Punto), pointer :: BPi
      CALL chdir("ins")
      inquire(file="source.txt",exist=lexist)
      if (lexist) then
        OPEN(77,FILE="source.txt",FORM="FORMATTED")
      else
        write(6,'(a)') 'There is a missing input file. '
        stop 'Check "source.txt" on Working directory' 
      end if
      READ(77,*);READ(77,*);READ(77,*);READ(77,*)
      READ(77,*) tipoFuente;READ(77,*);READ(77,*)
      READ(77,*) nFuentes;READ(77,*) 
      allocate(Po(nFuentes))
      do i=1,nFuentes
       READ(77,*) xfsource, zfsource, nxfsource,&
                 nyfsource, nzfsource, PW_theta, l, regi
       Po(i)%center%x = xfsource
       Po(i)%center%z = zfsource
       Po(i)%normal%x = nxfsource 
       Po(i)%normal%y = nyfsource 
       Po(i)%normal%z = nzfsource 
       Po(i)%cosT = cos(PW_theta*pi/180.0)
       Po(i)%sinT = sin(PW_theta*pi/180.0)
       Po(i)%gamma = PW_theta*pi/180.0
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
      end do
      READ(77,*); READ(77,*) Escala; READ(77,*) ampfunction; READ(77,*) t0
      READ(77,*) Ts; READ(77,*) Tp; READ(77,*) sigGaus
      READ(77,*) PW_pol; close(77); CALL chdir("..")
      do i = 1,nFuentes
      write(Printnum,'(a)') &
       "---------------------------------------------------------------------------------"
      write(PrintNum,'(/,a,F8.2,a,F8.2,a,2x,a,F9.2,a,F9.2,a,F9.2,a,I0,a,I0)') & 
      "   Source: (",Po(i)%center%x,",",Po(i)%center%z,")", &
      "n=[",Po(i)%normal%x,",",Po(i)%normal%y,",",Po(i)%normal%z,& 
      "] r= ",Po(i)%region," e=",Po(i)%layer
      end do
      
      skipdir = .true.
      SH = .false. ; PSV = .false.
      do i=1,Nfuentes
       if (abs(Po(i)%normal%x) .gt. 0.001_8) skipdir(1) = .false.
       if (abs(Po(i)%normal%y) .gt. 0.001_8) skipdir(2) = .false.
       if (abs(Po(i)%normal%z) .gt. 0.001_8) skipdir(3) = .false.
      end do
      if (skipdir(2) .eqv. .false.) SH = .true.
      if (skipdir(1) .eqv. .false.) PSV = .true.
      if (skipdir(3) .eqv. .false.) PSV = .true.
      end subroutine getsource
      
      subroutine getVideoPoints
      use resultVars, only : moviePoints, nMpts, & 
                             iPtfin,mPtini,mPtfin
      use peli, only : coords_Z,coords_X,fotogramas,fotogramas_Region
      use meshVars
      use gloVars,only: makevideo,z0
      use waveNumVars, only : NPTSTIME
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
      READ(7,*) npixZ,npixX
      READ(7,*) !      Vertical (start,end)
      READ(7,*) MeshMinZ,MeshMaxZ,firstZ
      READ(7,*) !      Horizontal (start,end)
      READ(7,*) MeshMinX,MeshMaxX
      close(7)
      CALL chdir("..")
      if (makeVideo .eqv. .false.) return 
        allocate(coords_X(npixX))
        allocate(coords_Z(npixZ))
        allocate(fotogramas(npixZ,npixX,NPTSTIME,3))
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
          moviePoints(iz)%center%x = 0.0
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
      boxVoid_minY, midPoint,normXI, origGeom, surf_poly
      use gloVars
      use fitting
      use soilVars, only : Z,N,RHO,BETA0,ALFA0,shadecolor_inc,LAMBDA0,ANU,AMU0
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
      huboCambios = .false.
      allocate(auxA(1)); allocate(auxB(1))
      
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
      shadecolor_inc = 0.8_4
      else
      shadecolor_inc = real(0.6-(maxBeta-beta0(N+2))*((0.6-0.89)/(maxBeta-minBeta)),4)
      end if
      
!     go to 384
      nXIoriginal = nXI
      ! cortar en intersección con estratos y determinar estrato de cada segemento
      if (verbose >= 4) write(PrintNum,*)"begin slice with layers"
      iXI = 1
      DO !para cada segmento
      ! we check if there is a change of medium along segment iXI
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
      
      if (iXI .eq. nXI) exit ! salidr de Do de iXI
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
      use glovars, only : makeVideo, flip12
      use soilVars, only : Z
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
      if (makeVideo .eqv. .false.) return
       do i=1,npixZ
        do j=1,npixX
          fotogramas_Region(i,j) = reg(1)!'estr'
          cn = adentroOafuera(coords_X(j), coords_Z(i),'incl')
          if (cn .eqv. .true.) then
            fotogramas_Region(i,j) = reg(2)!'incl' 
          end if
          cn = adentroOafuera(coords_X(j), coords_Z(i),'void')
          if (cn .eqv. .true.) then
            fotogramas_Region(i,j) = reg(0)!'void' 
          end if
          if (abs(z(0)) .gt. 0.0001) then ! si no es cero
            if (coords_Z(i) .lt. 0.0) then
              fotogramas_Region(i,j) = reg(0)!'void'
            end if
          end if
        end do
       end do
      end subroutine setVideoPointsRegions    
      
      subroutine setInqPointsRegions
      use resultVars, only : allpoints,nIpts,nPts
      use soilVars, only : N
      use glovars, only : verbose, Printnum, flip12
      implicit none
      integer :: i
      integer, dimension(0:2) :: reg
      logical :: cn,adentroOafuera
      !!1 estr, 2 incl, 0 void
      reg(0) = 0; reg(1)= 1; reg(2) = 2
      if (flip12) then
        reg(0) = 0; reg(1)= 2; reg(2) = 1
      end if
      allpoints(1:nPts)%region = reg(1)!'estr'
      if (verbose .ge. 2) &
        write(PrintNum,*) "center,region,layer"
      do i = 1, nIpts
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
        if (allpoints(i)%region .eq. 2) allpoints(i)%layer = N+2
        if (verbose .ge. 2) then
        write(PrintNum,*)i,"[",allpoints(i)%center%x, & 
        ",",allpoints(i)%center%z,"] is ",allpoints(i)%region,allpoints(i)%layer
        end if
      end do!;stop "setInqPointsRegions"
      end subroutine setInqPointsRegions
      
      function adentroOafuera(x,y,region)
       use geometryvars, only : Xcoord_Voidonly,Xcoord_Incluonly, & 
            boxIncl_maxX,boxIncl_maxY,boxIncl_minX,boxIncl_minY, &
            boxVoid_maxX,boxVoid_maxY,boxVoid_minX,boxVoid_minY
                         
       implicit none
       logical :: adentroOafuera
       character(LEN=4) :: region
       real*4,intent(in) :: x,y
       real*16 :: vt
       integer :: cn,k,nXI
       real*8 :: XcooBox_maxX,XcooBox_minX,XcooBox_maxY,XcooBox_minY
       real*8, dimension(:,:,:), pointer :: Xcoord_ER
       nullify(Xcoord_ER)
       
       adentroOafuera = .false.
       if (region .eq. 'void') then
       XcooBox_maxX = boxVoid_maxX
       XcooBox_maxY = boxVoid_maxY
       XcooBox_minX = boxVoid_minX
       XcooBox_minY = boxVoid_minY
       if (allocated( Xcoord_Voidonly)) then
       Xcoord_ER => Xcoord_Voidonly
       else
       return
       end if
       else if (region .eq. 'incl') then
       XcooBox_maxX = boxIncl_maxX
       XcooBox_maxY = boxIncl_maxY
       XcooBox_minX = boxIncl_minX
       XcooBox_minY = boxIncl_minY
       if (allocated( Xcoord_Incluonly)) then
       Xcoord_ER => Xcoord_Incluonly
       else
       return
       end if
       else 
       return
       end if
!      print*,XcooBox_maxX,XcooBox_maxY,XcooBox_minX,XcooBox_minY
       nXI = size(Xcoord_ER(:,1,1))
          if (XcooBox_minX .le. x) then
          if (x .le. XcooBox_maxX) then
          if (XcooBox_minY .le. y) then
          if (y .le. XcooBox_maxY) then
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
            adentroOafuera = .true. !el punto pertenece a la ragion de adentro
          end if
          
          end if
          end if
          end if
          end if
       
      end function adentroOafuera
      
! pointer table 
      subroutine preparePointerTable(pota,firstTime,smallestWL)
      use resultVars, only : nPts,allpoints,nBpts,BouPoints, & 
                             fixedpota,nZs,Punto,n_top_sub,n_con_sub
      use debugstuff
      use Gquadrature, only : Gquad_n
      use glovars, only : verbose,workBoundary, rutaOut
!     use interpolPlaneWavesOfStratification
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
         if (allpoints(ip)%region .ne. 1) cycle !nada mas agregar receptores en medio estratificado
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
       if(verbose .ge. 4) then
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
         
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptores
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
!     real*8     :: errT = 0.01_8
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
!         if (resi .gt. errT) then
!          nsubdivs = 2*ndivsiguales + 1 + 1
!          allocate(subdiv(ixi)%x(nsubdivs));allocate(subdiv(ixi)%z(nsubdivs))
!          subdiv(ixi)%x(ndivsiguales+1) = origGeom(iXI)%bord_A%x + deltaX*ndivsiguales
!          subdiv(ixi)%z(ndivsiguales+1) = origGeom(iXI)%bord_A%z + deltaZ*ndivsiguales
!          subdiv(ixi)%x(ndivsiguales+2) = origGeom(iXI)%bord_B%x - deltaX*ndivsiguales
!          subdiv(ixi)%z(ndivsiguales+2) = origGeom(iXI)%bord_B%z - deltaZ*ndivsiguales 
!         else !por fin si solo por la mitad
           nsubdivs = 2*ndivsiguales + 1
           allocate(subdiv(ixi)%x(nsubdivs));allocate(subdiv(ixi)%z(nsubdivs))
           subdiv(ixi)%x(ndivsiguales+1) = (origGeom(iXI)%bord_A%x + deltaX*ndivsiguales + & 
                                            origGeom(iXI)%bord_B%x - deltaX*ndivsiguales)/2.
           subdiv(ixi)%z(ndivsiguales+1) = (origGeom(iXI)%bord_A%z + deltaZ*ndivsiguales + &
                                            origGeom(iXI)%bord_B%z - deltaZ*ndivsiguales)/2.
!         end if
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
      BouPoints(1+bou_conta_ini:bou_conta_fin)%normal%y = 0.0_8
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
        allocate(boupoints(ixi)%G(nIpts,6,3)); boupoints(ixi)%G = z0
        !                               `--- 1 a 6 elemento mecánico
        !                                          1 : W      4 : Tz
        !                                          2 : U      5 : Tx
        !                                          3 : V      6 : Ty
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
       call drawBoundary(BP,nBpts,txt, extension,.true.,.false.)
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
                             Gqu_t => Gqu_t_8, & 
                             Gqu_A => Gqu_A_8
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
      use wavelets !las funciones: ricker, fork
      use waveVars, only : Dt,Uo, ampfunction,Escala,dt,maxtime
      use waveNumVars, only : DFREC,nfrec, NPTSTIME
      use gloVars, only : ve => verbose,Ur,rutaOut,z0,Printnum
      use ploteo10pesos
      implicit none
      character(LEN=9)          :: logflag
      character(LEN=100)        :: titleN,xAx,yAx,CTIT
      CHARACTER(len=32) :: arg
      !complex*16 :: FFTW!(NPTSTIME)
      integer :: n_maxtime
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then 
        if (trim(arg) .eq. '-a') ve = 6
      end if
      
       n_maxtime = int(maxtime/dt)
       if(maxtime .lt. dt) n_maxtime = 2*nfrec
       if(maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
      
!     real :: factor
      ! Amplitude function of incident wave
      ! prepare the signal we will use as incident wave amplitude.
      ! el tamaño del ricker es 2*NFREC porque incluirá el conjugado
      allocate(Uo(NPTSTIME))
      Uo=z0
      if(ampfunction .eq. 1) then
        call ricker ! Ricker wavelet saved on Uo
          if (ve .ge. 1) then
             write(Printnum,'(a)')'   Incident wave amplitude function: Ricker'
            call chdir(trim(adjustl(rutaOut)))
            write(titleN,'(a)') 'WaveAmplitude-ricker_time.pdf'
            write(CTIT,'(a)') 'WaveAmplitude of Ricker wavelet'
            xAx = 'time[sec]'
            write(yAx,'(a)') 'amplitude'
            call plotXYcomp(Uo(1:n_maxtime),real(Dt,4), & 
                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
            CALL chdir("..")
          end if
        ! forward
        Uo = FFTW(NPTSTIME,Uo,-1,Dt)        
          if (ve .ge. 1) then
            call chdir(trim(adjustl(rutaOut)))
            write(titleN,'(a)') 'WaveAmplitude-ricker_frec.pdf'
            xAx = 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'      
            logflag = 'logx     '      
!           logflag = 'none     '
            call plotSpectrum(Uo,real(DFREC,4), size(Uo),int(size(Uo)/2), & 
            titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
            CALL chdir("..")
          end if
          
     !*************************** test *********************************************
!         Uo = Uo * (sqrt(1.0*NPTSTIME) * dfrec) !backward
!         call fork(size(Uo),Uo,+1,verbose,outpf)
!           CALL chdir("outs")
!           write(titleN,'(a)') 'WaveAmplitude-ricker_time_BACK.pdf'
!           xAx = 'time[sec]'
!           write(yAx,'(a)') 'amplitude'      
!           call plotXYcomp(Uo,real(Dt,4),size(Uo),titleN,xAx,yAx,CTIT,1200,800,0.0)
!           CALL chdir("..")
!           stop 3903
     !********************************************************************************
      elseif(ampfunction .eq. 2) then ! Gaussian
        call gaussian
          if (ve .ge. 1) then
           write(Printnum,'(a)')'   Incident wave amplitude function: Gaussian'
            call chdir(trim(adjustl(rutaOut)))
            write(titleN,'(a)') 'WaveAmplitude-Gaussian_frec.pdf'
            write(CTIT,'(a)') 'WaveAmplitude of Gaussian wavelet'
            xAx = 'frec[Hz] '
            write(yAx,'(a)') 'amplitude'
            call plotXYcomp(Uo(1:n_maxtime),real(DFREC,4), & 
                 n_maxtime,titleN,xAx,yAx,CTIT,1200,800,0.0)
            CALL chdir("..")
          end if
      else  ! DIRAC -----------------------------------------
       if (ve .ge. 1) write(Printnum,'(a)')'   Incident wave amplitude function: Dirac delta'
        Uo=UR
      end if !tipoPulso
      Uo = Uo*Escala !escala de la señal
      write(Printnum,'(a)') &
       "---------------------------------------------------------------------------------"
      if (ve .eq. 6) stop "argumento -a"
      end subroutine sourceAmplitudeFunction
      subroutine diffField_at_iz(i_zF,dir_j,J,cOME_in)
!#define ver 1
      ! esta función es llamada con cada profundidad donde hay
      ! por lo menos una fuente.
      use gloVars, only: z0, plotFKS,UI,UR,outpf => PrintNum
      use resultVars, only : pota,Punto,nZs,MecaElem,FFres
      use refSolMatrixVars, only : B,Ak
      use waveNumVars, only : NMAX,k_vec,dk,vecNK!,DFREC
      use wavelets !fork 
      use dislin
      use sourceVars, only: tipofuente, nFuentes, PW_pol
      use soilvars, only:N,Z,alfa0,beta0!,alfa,beta
      use, intrinsic :: iso_c_binding!, only : C_INT
      implicit none
      interface
        include 'interfaz.f'
        function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k)
          use soilvars, only:N
          complex*16, dimension(1:5) :: PSVdiffByStrata
          real*8, intent(in)           :: z_i,k
          complex*16, intent(in)       :: cOME_i  
          integer, intent(in)          :: e
          complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
        end function PSVdiffByStrata
      end interface
      integer, intent(in) :: i_zF,dir_j,J
      complex*16, intent(in),target  :: cOME_in
      logical,pointer :: intf
      integer, pointer :: ef
      real*8,target :: k
      real*8, pointer :: zf,xf,pt_k
      complex*16, dimension(:,:), allocatable, target :: auxK,savedAuxK
      complex*16, target  :: cOME
      complex*16, pointer :: pt_cOME_i
      integer, dimension(:),allocatable,target :: ipivA
      integer, dimension(:),pointer :: pt_ipivA
      complex*16, dimension(:),allocatable,target :: workA
      complex*16, dimension(:),pointer :: pt_workA
      complex*16,dimension(:,:),pointer :: pointA
      type(Punto), pointer :: pXi,p_X
      logical :: auxLogic ,porLoMenosUnoEsEstr
      integer :: ik,tam,itabla_z,itabla_x,iMec,mecS,mecE,&
                 nXis,n_Xs,iXi,dj,i_Fuente,i_FuenteFinal,po,ne
      type(MecaElem)  :: Meca_diff, SHdiffByStrata
      allocate(auxK(2*nmax,5)); allocate(savedAuxK(2*nmax,5))
#ifdef ver
      real :: result,lastresult
      real, dimension(2) :: tarray
#endif
      cOME = cOME_in 
      if (i_zF .eq. 0) then
         itabla_x = 3 ! En la tabla (pota) de índices: la fuente real -> (0,3)
         i_FuenteFinal = nFuentes ! Cantidad de fuentes reales
         if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then ! onda plana incidente ······················p
            ! con incidencia de onda plana no usamos atenuación                  ·l
            cOME = real(cOME_in) * UR!                                           ·a
         end if! ·································································n
      else; itabla_x = 2 + pota(i_zF,1) + 1 !    la primer fuente virtual
            i_FuenteFinal = 1; end if       !        a esa profundidad
      ! ------- para cada una de las fuentes en la tabla -------------------------
      do i_Fuente = 1,i_FuenteFinal
       call asociar(pXi,i_Fuente,i_zF,itabla_x) ! asociar apuntador a fuente [pXi]
       xf=>pXi%center%x;zf=>pXi%center%z;ef=>pXi%layer;intf=>pXi%isOnInterface
#ifdef ver
       call ETIME(tarray, result)
       print*,"pXI%center%x",pXi%center%x,"pXI%center%z",pXi%center%z,& 
       "isboundary",pXi%isboundary, result
       lastresult = result
#endif
      ! Si es la fuente real y es una onda plana no se usa el DWN. Se calcula para
      ! el número de onda horizontal asociado al ángulo de incidencia ············
      if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then ! onda plana incidente   ·
          if (dir_j .eq. 2) then! SH                                             ·
            if(PW_pol .eq. 3) k = real(cOME/beta0(N+1)*sin(pXi%gamma))!          ·o
            tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1!                       ·n
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·d
            pt_k => k; pt_come_i => cOME!                                        ·a
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call globalmatrix_SH(pointA,pt_k,pt_come_i)!                         ·p
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·l
            B(:,0) = matmul(Ak(:,:,0),B(:,0))!                                   ·a
          else!  P-SV                                                            ·n
            if(PW_pol .eq. 1) k = real(cOME/beta0(N+1)*sin(pXi%gamma))!          ·a
            if(PW_pol .eq. 2) k = real(cOME/alfa0(N+1)*sin(pXi%gamma))!          ·
            tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2!                       ·
            pointA => Ak(1:tam,1:tam,0) !indice 0 reservado para onda plana      ·
            pt_k => k; pt_come_i => cOME!                                        ·
            allocate(ipivA(tam)); allocate(workA((tam)*(tam)))!                  ·
            pt_ipivA => ipivA; pt_workA => workA!                                ·
            call gloMat_PSV(pointA,pt_k,pt_come_i)!                              ·
            call inverseA(pointA,pt_ipivA,pt_workA,tam)!                         ·
            call PSVvectorB_ondaplana(B(:,0),cOME,pxi%gamma)!                    ·
            B(:,0) = matmul(Ak(1:tam,1:tam,0),B(:,0))!                           ·
          end if!                                                                ·
      else ! · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · · 
      ! La fuente es cilíndrica, puede tratarse de la fuente real o una virtual ·· 
        po = min(int(vecNK(J)*1.25),nmax); ne = 2*nmax-(po-2)                    !
        B(:,po:ne) = 0                                                           !
        if (dir_j .eq. 2) then                                                   !
         tam = 2*N+1; if (Z(0) .lt. 0.0) tam = tam + 1                           !o
         do ik = 1,po                                                            !n
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik))             !d
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !a
         end do                                                                  !                                     
         do ik = ne,2*NMAX                                                       !c
           call SHvectorB_force(i_zF,B(:,ik),tam,pXi,cOME,k_vec(ik))             !i
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !l
         end do                                                                  !i
        else!  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .     !n        
         tam = 4*N+2; if (Z(0) .lt. 0.0) tam = tam + 2                           !d
         if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) return         
         do ik = 1,po                                                            !r
         call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik))        !i
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !c
         end do                                                                  !a
         do ik = ne,2*NMAX                                                       !
           call PSVvectorB_force(i_zF,B(:,ik),tam,pXi,dir_j,cOME,k_vec(ik))      !
           B(:,ik) = matmul(Ak(:,:,ik),B(:,ik))                                  !
         end do                                                                  !
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
#ifdef ver 
       print*,"  itabla_z",itabla_z,"de",nZs,"******************"
       print*,"   pota(itabla_z,:)=",pota(itabla_z,:)
#endif
      !(en la tabla todos los puntos son receptores)
       auxLogic = porLoMenosUnoEsEstr(itabla_z)
       if (auxLogic .eqv. .false.) cycle
        ! (el primer receptor de la tabla --------| a esa profundidad) 
            call asociar(p_X, i_Fuente, itabla_z, 3)
            if (p_X%layer .eq. N+2) cycle
            dj = dir_j; if(dj .eq. 3) dj = 2 
        ! ¿calcular sólo G, sólo T o ambos?
            if (dir_j .eq. 2) then; mecS = 1; mecE = 3 !V,s32,s12
            else;                   mecS = 1; mecE = 5 !W,U,s33,s31,s11
            end if
#ifdef ver 
       call ETIME(tarray, result)
       print*,"  p_x%center%z",p_x%center%z,"[m]  ",result-lastresult
       lastresult = result
#endif
            
      ! ... elementos mecánicos a la profundidad p_X%center%z ........
      ! .... usando los coeficientes de las ondas en el estrato ......
      savedauxk = z0
      if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then ! onda plana·············
          if (dir_j .eq. 2) then!                                               ·
            if(PW_pol .eq. 3) k = real(cOME/beta0(N+1)*sin(pXi%gamma))!         ·
                 Meca_diff = SHdiffByStrata(B(:,0), &!                          ·o
                             p_X%center%z, p_X%layer, & !                       ·n
                             cOME,k,mecS,mecE,outpf) !                          ·d
                 savedauxK(1,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE) !          ·a
          else !                                                                ·
            if(PW_pol .eq. 1) k = real(cOME/beta0(N+1)*sin(pXi%gamma))!         ·p
            if(PW_pol .eq. 2) k = real(cOME/alfa0(N+1)*sin(pXi%gamma))!         ·l
                 savedauxk(1,1:5) = PSVdiffByStrata(B(:,0), &!                  ·a
                              p_X%center%z, p_X%layer,cOME,k)!                  ·a
          end if!                                                               ·
      else ! onda plana incidente / onda cilíndrica circular ····················
        do ik = 1,po                                                            !
          if (dir_j .eq. 2) then !SH                                            !
              Meca_diff = SHdiffByStrata(B(:,ik), &                             !
                             p_X%center%z, p_X%layer, &                         !o
                             cOME,k_vec(ik),mecS,mecE,outpf)                    !n
              savedauxK(ik,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE)              !d
          else !PSV                                                             !a
              savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                    !c
                              p_X%center%z, p_X%layer,cOME,k_vec(ik))           !i
          end if                                                                !l
        end do ! ik                                                             !i
        do ik = ne,2*Nmax                                                       !n
          if (dir_j .eq. 2) then !SH                                            !d
              Meca_diff = SHdiffByStrata(B(:,ik), &                             !r
                             p_X%center%z, p_X%layer, &                         !i
                             cOME,k_vec(ik),mecS,mecE,outpf)                    !c
             savedauxK(ik,mecS:mecE) = Meca_diff%Rw_SH(mecS:mecE)               !a
          else !PSV                                                             !
             savedauxk(ik,1:5) = PSVdiffByStrata(B(:,ik), &                     !
                                 p_X%center%z, p_X%layer,cOME,k_vec(ik))        !
          end if                                                                !
        end do ! ik                                                             !
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
            call asociar(pXi, i_Fuente, i_zF, itabla_x)
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
          call asociar(p_x, i_Fuente, itabla_z, 2+ itabla_x)
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
          auxK(1:2*nmax,mecS:mecE) = savedAuxK(1:2*nmax,mecS:mecE)
            ! agregar información fase horizontal de fuente y receptor 
          do imec = mecS,mecE !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
            if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then !             ·
              if (dir_j .eq. 2) then!                                       ·
!               if(PW_pol .eq. 3) kx = cOME/beta(N+1)*sin(pXi%gamma)!       ·p
                if(PW_pol .eq. 3) k = real(cOME/beta0(N+1)*sin(pXi%gamma))! ·p
              else!                                                         ·l
!               if(PW_pol .eq. 1) kx = cOME/beta(N+1)*sin(pXi%gamma)!       ·a
!               if(PW_pol .eq. 2) kx = cOME/alfa(N+1)*sin(pXi%gamma)!       ·n
                if(PW_pol .eq. 1) k = real(cOME/beta0(N+1)*sin(pXi%gamma))! ·a
                if(PW_pol .eq. 2) k = real(cOME/alfa0(N+1)*sin(pXi%gamma))! ·n
              end if  !                                                     ·a
                auxk(1,imec) = auxk(1,imec) * &           ! onda plana      ·
                exp(-UI*k*(p_x%center%x - xf))            !                 ·
                CYCLE ! imec                              !                 ·
            end if ! ························································
            do ik = 1,po                                                    !
                auxk(ik,imec) = auxk(ik,imec) * &                           !
                exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_x%center%x - xf), 8))!c
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
                    if( imec .le. 2) then                                   !
                      p_x%FK(J,1:nmax,imec) = & 
                      p_x%FK(J,1:nmax,imec) + auxK(1:nmax,imec) 
                    end if
                end if;end if;end if;end if!-------------------------       !
#ifdef ver
       call ETIME(tarray, result)
       print*,"        guardado el integrando",result-lastresult            !
       lastresult = result
#endif
      ! K -> X  .........................................................   !
         if (p_x%guardarMovieSiblings .eqv. .false.) then
             auxK(:,iMec) = FFTW(2*nmax,auxK(:,iMec),+1,dk) !backward       !
#ifdef ver
       call ETIME(tarray, result)                                           !
       print*,"        fork non movie", result-lastresult                   !
       lastresult = result                                                  !
#endif
         end if                                                             !
      end do !imec !.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
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
      end do !i_Fuente
#ifdef ver
       stop "2921 end of diffField_at_iz"
#endif
      end subroutine diffField_at_iz
      
! G_stra - matrix pointAp,pt_k,pt_cOME_i
      subroutine gloMat_PSV(this_A,k,cOME_i)
      ! Calcular para +k. Una vez invertida la matrix, hay paridades para -k.
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
      use debugStuff
      implicit none
      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      complex*16, intent(in),pointer     :: cOME_i 
      
      real*8     :: k2
      complex*16, dimension(2,4) :: subMatD,subMatS, &
                                   subMatD0,subMatS0
      complex*16, dimension(4,4) :: diagMat
      complex*16 :: gamma,nu,xi,ega,enu,ck
      integer    :: i,iR,iC,e,bord
      this_A = 0 
      ck = -k*UR
      iR= 0;iC= 0;i=1
      if (Z(0) .lt. 0.0) then !Half-Space por arriba de z=0 
        i = 0;iR= -2; iC= -2
      end if
      
      DO e = i,N+1
          gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
          nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z crece.
          if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
          xi = k**2.0 - nu**2.0 
          ! en fortran los elementos se indican por columnas:
          subMatD0 = RESHAPE((/ -gamma, ck, ck,nu,& 
                                 gamma, ck, ck,-nu /), &
                           (/ 2,4 /))
          subMatD0 = subMatD0 * UI
          k2 = 2.0*k
          subMatS0 = RESHAPE((/ xi,-k2*gamma,-k2*nu,-xi,& 
                               xi,k2*gamma,k2*nu,-xi /),&
                           (/2,4/)) 
          subMatS0 = amu(e) * subMatS0 
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
          diagMat = RESHAPE((/ UR, Z0, Z0, Z0, & 
                               Z0, UR, Z0, Z0, & 
                               Z0, Z0, ega, Z0, & 
                               Z0, Z0, Z0, enu /), &
                           (/ 4,4 /))
        else !bord .eq. 1 -----------------------------------
          if (e /= 0) then !(radiation condition upper HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e))) 
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))       
          else
            ega = Z0
            enu = Z0
          end if !<----<----<----<----<----<----<----<----<--
          diagMat = RESHAPE((/ ega, Z0, Z0, Z0, & 
                               Z0, enu, Z0, Z0, & 
                               Z0, Z0, UR, Z0, & 
                               Z0, Z0, Z0, UR /), &
                           (/ 4,4 /))
        end if
          ! desplazamientos
          subMatD = matmul(subMatD0,diagMat)
          ! esfuerzos
          subMatS = matmul(subMatS0,diagMat)
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
      
!     ! Calcular para +k y -k al mismo tiempo. Se ahorran algunas operaciones
!     subroutine globalmatrix_PSV(this_A,this_An,k,cOME_i)
!     use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
!     use gloVars, only : UI,UR,Z0
!     use debugStuff
!     implicit none
!     
!     complex*16,intent(inout),dimension(:,:),pointer :: this_A,this_An
!     real*8,     intent(in),pointer     :: k
!     complex*16, intent(in),pointer     :: cOME_i  
!     
!     real*8     :: z_i
!     complex*16, dimension(2,4) :: subMatDp,  subMatDn, & 
!                                   subMatSp,  subMatSn, &
!                                  subMatD0p, subMatD0n, & 
!                                  subMatS0p, subMatS0n
!     complex*16, dimension(4,4) :: diagMat
!     complex*16 :: gamma,nu,xi
!     complex*16 :: egammaN,enuN,egammaP,enuP
!     integer    :: i,iR,iC,e,bord
!     complex*16 :: mink_zp,mink_zn
!         iR= 0
!         iC= 0 
!         mink_zp = -k*UR ! porque armamos matrices para +k y para -k al mismo tiempo.
!         mink_zn =  k*UR
!     !la matriz global se llena por columnas con submatrices de cada estrato
!     i=1
!     if (Z(0) .lt. 0.0) then 
!       i = 0 !Half-Space por arriba de z=0 
!       iR= -2
!       iC= -2
!     end if
!     DO e = i,N+1
!         ! algunas valores constantes para todo el estrato
!         gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
!         nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
!         ! Se debe cumplir que la parte imaginaria del número de onda 
!         ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
!         ! nencial a medida que z es mayor que cero.
!         if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
!         if(aimag(nu).gt.0.0)nu= conjg(nu)
!         
!         xi = k**2.0 - nu**2.0 
!         ! en fortran los elementos se indican por columnas:
!         subMatD0p = RESHAPE((/ -gamma, mink_zp, mink_zp,nu,& 
!                                gamma, mink_zp, mink_zp,-nu /), &
!                          (/ 2,4 /))
!         subMatD0p = subMatD0p * UI
!         
!         subMatD0n = RESHAPE((/ -gamma, mink_zn, mink_zn,nu,& 
!                                gamma, mink_zn, mink_zn,-nu /), &
!                          (/ 2,4 /))
!         subMatD0n = subMatD0n * UI
!         
!         subMatS0p = RESHAPE((/ xi,-2.0*k*gamma,-2.0*k*nu,-xi,& 
!                              xi,2.0*k*gamma,2.0*k*nu,-xi /),&
!                          (/2,4/)) 
!         subMatS0p = amu(e) * subMatS0p 
!         
!         subMatS0n = RESHAPE((/ xi, 2.0*k*gamma, 2.0*k*nu,-xi,& 
!                              xi,-2.0*k*gamma,-2.0*k*nu,-xi /),&
!                          (/2,4/)) 
!         subMatS0n = amu(e) * subMatS0n 
!         
!         ! la profundidad z de la frontera superior del estrato
!!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!!         z_f = Z(e+1) ! e=1  ->  z = z1 
!       do bord = 0,1
!         if ((e .eq. 0) .and. (bord .eq. 0)) cycle
!         if (e+bord > N+1) then ! si 1+0;1+1;2+0;[2+1] > 2
!           exit
!         end if                           
!         ! la profundidad z de la frontera superior del estrato
!               z_i = Z(e+bord)   ! e=1 , bord=0  ->  z = z0 = 0
!                                 ! e=1 , bord=1  ->  z = Z1 = h1
!         !downward waves
!         !  bord=0  ... 0 | bord=1  ... h(e)=Z(e+1)-Z(e)
!         !        -> 1    |       -> exp(-UI * nu * h(e))
!         if (e /= 0) then !(radiation condition upper HS)
!           egammaN = exp(-UI * gamma * (z_i-Z(e))) 
!           enuN = exp(-UI * nu * (z_i-Z(e)))       
!         else
!           egammaN = Z0
!           enuN = Z0
!         end if
!         !upward waves    
!         if (e /= N+1) then !(radiation condition)
!           egammaP = exp(UI * gamma * (z_i-Z(e+1))) !
!           enuP = exp(UI * nu * (z_i-Z(e+1)))
!         else
!           egammaP = Z0
!           enuP = Z0
!         end if
!           !la matrix diagonal
!        diagMat = RESHAPE((/ egammaN, Z0, Z0, Z0, & 
!                             Z0,    enuN, Z0, Z0, & 
!                             Z0, Z0, egammaP, Z0, & 
!                             Z0, Z0, Z0, enuP /), &
!                          (/ 4,4 /))
!         ! desplazamientos
!         subMatDp = matmul(subMatD0p,diagMat)
!         subMatDn = matmul(subMatD0n,diagMat)
!         ! esfuerzos
!         subMatSp = matmul(subMatS0p,diagMat)
!         subMatSn = matmul(subMatS0n,diagMat)
!       !ensamble de la macro columna i
!         !evaluadas en el borde SUPERIOR del layer i
!         if (bord == 0 .AND. e /= 0 ) then 
!          if (e /= N+1) then !(radiation condition)
!           if (e/=i) then !(only stress bound.cond. in the surface
!            this_A( iR-1 : iR   , iC+1 : iC+4 ) = -subMatDp
!           end if
!            this_A( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatSp
!          else
!          if (e/=i) then !(only stress bound.cond. in the surface
!            this_A( iR-1 : iR   , iC+1 : iC+2 ) = -subMatDp(:,1:2)
!          end if
!            this_A( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatSp(:,1:2)
!           ! exit
!          end if
!         end if
!         
!         !evaluadas en el borde INFERIOR del layer i
!         if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
!          if (e /= 0) then !(radiation condition upward)
!           ! ambas ondas
!           this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatDp
!           this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatSp
!          else
!           ! solo onda hacia arriba
!           this_A( iR+3 : iR+4 , iC+3 : iC+4 ) = subMatDp(:,3:4)
!           this_A( iR+5 : iR+6 , iC+3 : iC+4 ) = subMatSp(:,3:4)
!          end if
!         end if
!        !
!     if (associated(this_An)) then
!         if (bord == 0 .AND. e /= 0 ) then 
!          if (e /= N+1) then !(radiation condition)
!           if (e/=i) then !(only stress bound.cond. in the surface
!           this_An( iR-1 : iR   , iC+1 : iC+4 ) = -subMatDn
!           end if
!           this_An( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatSn
!          else
!          if (e/=i) then !(only stress bound.cond. in the surface
!           this_An( iR-1 : iR   , iC+1 : iC+2 ) = -subMatDn(:,1:2)
!          end if
!           this_An( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatSn(:,1:2)
!          end if
!         end if
!         !evaluadas en el borde INFERIOR del layer i
!         if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
!          if (e /= 0) then !(radiation condition upward)
!           ! ambas ondas
!           this_An( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatDn
!           this_An( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatSn
!          else
!           ! solo onda hacia arriba
!           this_An( iR+3 : iR+4 , iC+3 : iC+4 ) = subMatDn(:,3:4)
!           this_An( iR+5 : iR+6 , iC+3 : iC+4 ) = subMatSn(:,3:4)
!          end if
!         end if
!     end if
!       end do !bord loop del borde i superior o nferior
!         iR= iR+4 
!         iC= iC+4
!     END DO !{e} loop de las macro columnas para cada estrato
!     end subroutine globalmatrix_PSV
      
      subroutine globalmatrix_SH(this_A,k,cOME_i)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : UI,UR,Z0
      use debugStuff  
      implicit none
      
      complex*16,    intent(inout), dimension(:,:),pointer :: this_A
      real*8,     intent(in),pointer     :: k
      complex*16, intent(in),pointer     :: cOME_i  
      
      real*8     :: z_i
      complex*16, dimension(1,2) :: subMatD, subMatS, subMatD0, subMatS0
      complex*16, dimension(2,2) :: diagMat
      complex*16 :: nu,enuN,enuP
      integer    :: i,iR,iC,e,bord
          iR= 0
          iC= 0 
      !la matriz global se llena por columnas con submatrices de cada estrato
      i = 1
      if (Z(0) .lt. 0.0) then 
        i = 0 !HS arriba (no free surface)
        iR= -1
        iC= -1
      end if
      DO e = i,N+1!cada estrato
          ! algunas valores constantes para todo el estrato
          nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
!         ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z es mayor que cero.
          if(aimag(nu).gt.0.0)nu= conjg(nu)

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
!         call showMNmatrixZ(2*N+2,2*N+2,this_A,"Ash  ",6) 
!         stop 5248
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
      
      ! los elementos de la inversa de la matrix global
      ! tienen tendencia asintótica lineal después para 
      ! k muy grande. Se interpola linealmente
      subroutine inter_gloMat(k0,k1,k2)
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
!     use waveNumVars, only : NMAX
      use fitting
      implicit none
      integer      ::  k0, k1, k2
      complex*16   :: vk0,vk1,vk2
      integer :: r,c,ik,tam
      real*8,dimension(k2-k0+1) :: vecX
      tam = size(Ak,1)
      print*,k0,k1,k2!;stop
      vecX = (/(ik,ik=k0,k2)/)
      do r=1,tam
      do c=1,tam
      vk0 = Ak(r,c,k0)
      vk1 = Ak(r,c,k1)
      vk2 = Ak(r,c,k2)
!     print*,""
!     print*,1.0_8*k0, vk0
!     print*,1.0_8*k1, vk1
!     print*,1.0_8*k2, vk2
      if (abs(vk0) .lt. 1E-020) then
      Ak(r,c,k0+1:k2) = 0
      cycle
      end if
      
!     Ak(r,c,k0+1:k2) =  Zpolyfit(vx,vy,d) ! complex*16, dimension(d+1)
!     
      call spline3p(vecX, & ! vector x [real*8]
                    Ak(r,c,k0:k2), & ! vector interpolado [complex*16]
                    k1-k0, & ! n puntos [x0 ... x1-1]
                    k2-k1+1, & !        [x1 ...  x2 ]
                    1.0_8*k0, vk0, &  !x0,y0
                    1.0_8*k1, vk1, &  !x1,y1
                    1.0_8*k2, vk2 )   !x2,y2
      
!     do ik=1,k2-k0+1
!     print*,ik,vecX(ik),Ak(r,c,int(vecX(ik)))
!     end do
      end do !c
      end do !r
!     stop "inter_gloMat"
      end subroutine inter_gloMat
      
      
      subroutine intrplr_gloMat(k0,n,pt_cOME_i,pt_ipivA,pt_workA)
      use refSolMatrixVars, only : Ak !Ak(#,#,ik:2*nmax)
      use waveNumVars, only : k_vec, NMAX
      use fitting
      implicit none
      interface
      include 'interfaz.f'
      end interface
      integer :: i,k0,k1,k2,n,tam,r,c,ik
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
         pointAp => Ak(1:tam,1:tam,idat(i))
         pt_k => k_vec(idat(i))
         call gloMat_PSV(pointAp,pt_k,pt_cOME_i)
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
      
      
      subroutine intrplr_zpoly_gloMat(k0,n,pt_cOME_i,pt_ipivA,pt_workA)
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
      complex*16, pointer :: pt_cOME_i
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
         call gloMat_PSV(pointAp,pt_k,pt_cOME_i)
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
!       print*,ik,ikp;cycle
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),"antes")
        do r=1,tam
        do c=1,tam
        Ak(r,c,ik) = Ak(r,c,ikp) * Comal(r,c)
        end do
        end do
!     call scripToMatlabMNmatrixZ(tam,tam,Ak(1:tam,1:tam,ik),"despu")
!     stop
      end do
!     stop "parImpar_gloMat"
      end subroutine parImpar_gloMat
      
! G_stra - term indep
      subroutine PSVvectorB_force(i_zF,this_B,tam,pXi,direction,cOME,k)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO,NPAR
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only: tipofuente
      implicit none
      
      integer, intent(in) :: i_zF,tam
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
      complex*16 :: omeAlf,omeBet
      
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: egamz,enuz,gamz,nuz
      complex*16, dimension(2) :: sincGamma, sincNu
      
                                  !  1 para Greeni1 (horizontal),
                                  !  3 para Greeni3 (vertical)
      complex*16, dimension(2) :: G11,G31,G33
      complex*16, dimension(2) :: s111,s331,s131,s113,s333,s313
      
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
      end if   
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1)  z_loc(2) = Z(e_f+1) - z_f !upward (+)
      
      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if (i_zF .eq. 0) el_tipo_de_fuente = tipofuente 
      
      G11=Z0;G31=Z0;G33=Z0
      s111=Z0;s331=Z0;s131=Z0;s113=Z0;s333=Z0;s313=Z0
      
      DEN = 4.0*PI*RHO(e_f)*cOME**2.0
      omeAlf = cOME**2.0/ALFA(e_f)**2.0
      omeBet = cOME**2.0/BETA(e_f)**2.0
      L2M = LAMBDA(e_f) + 2.0*AMU(e_f)
      
          ! algunas valores constantes para todo el estrato          
          gamma = sqrt(omeAlf - k**2.0)
          nu = sqrt(omeBet - k**2.0)
          ! Se debe cumplir que la parte imaginaria del número de onda 
          ! vertical debe ser menor que cero. La parte imaginaria contri-
          ! buye a tener ondas planas inhomogéneas con decaimiento expo-
          ! nencial a medida que z crece.
          if(aimag(gamma).gt.0.0_8)gamma = conjg(gamma)
          if(aimag(nu).gt.0.0_8)nu=conjg(nu)
          
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
!       !sincGamma(2) = egamz(2)*(sin(k*A)/(k*A)*(2*A)) !no se inclina
        
        argum = nu*A*seno + UR * (k*A*cose)
        sincNu(2) = sincmod(argum,nuz(2))*(2*A)
!       !sincNu(2) = enuz(2) *(sin(k*A)/(k*A)*(2*A)) !no se inclina
        
        ! caso (receptor arriba de la fuente) : interfaz de arriba (1)
        argum = -gamma*A*seno + UR * (k*A*cose)
        sincGamma(1) = sincmod(argum,gamz(1))*(2*A)
!       !sincGamma(1) = egamz(1) *(sin(k*A)/(k*A)*(2*A)) !no se inclina
        
        argum = -nu*A*seno + UR * (k*A*cose)
        sincNu(1) = sincmod(argum,nuz(1))*(2*A)
!       !sincNu(1) = enuz(1) *(sin(k*A)/(k*A)*(2*A)) !no se inclina
        
      end if !fuente tipo 0 o 1
      
      do iIf = 1,2    
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if

      G31(iIf) = -UI/DEN * SGNz*k*(sincGamma(iIf) & 
                           -sincNu(iIf))
                           
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
      
      if (fisInterf) then
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
      
      subroutine PSVvectorB_ondaplana(this_B,come,gamma)
      use soilvars, only : n,lambda0,amu0,alfa0,beta0,Z!,beta,alfa
      use glovars, only:UI,z0
      use sourceVars, only: PW_pol
      use debugStuff
      implicit none
      complex*16, intent(inout), dimension(1:4*N+2) :: this_B
      complex*16, intent(in)    :: come ! no trae amortiguamiento
      real*8, intent(in) :: gamma
      integer :: i,e
      real*8,dimension(1:2) :: theta
      complex*16 :: kx,kz,U,W,c
      real*8 :: z_loc
      !     Colocamos la onda incindente en la interfaz
      !     con el semiespacio de abajo.
      e = N+1 
      z_loc = 0! (Z(N+1)- Z(N+1)) 
      if (PW_pol .eq. 1) then
        c = beta0(N+1) !SV
        theta(1) = cos(gamma)
        theta(2) = sin(gamma)
      elseif (PW_pol .eq. 2) then 
        c = alfa0(N+1) !P
        theta(1) = sin(gamma)
        theta(2) = -cos(gamma)
      end if
      
      kx = come/c * sin(gamma)
      kz = come/c * cos(gamma)
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
                            ( W * kz * (LAMBDA0(e) + 2.0*AMU0(e))) & !
                          - ( U * kx * LAMBDA0(e))) ! szz           !
      this_B(1+4*(e-1)+1 + i) = UI * AMU0(e) &                      !
                          * ( kz * U - kx * W ) ! szx              !
      !                   sxx = UI * ( &                           !
      !                   - ( U * kx * (LAMBDA(e) + 2.0*AMU(e))) & !
      !                   + ( W * kz * LAMBDA(e)))                 !
      !............................................................!
!     call showMNmatrixZ(4*N+2,1, this_B ,"  B  ",6)
      end subroutine PSVvectorB_ondaplana
      
      subroutine SHvectorB_force(i_zF,this_B,tam,pXi,cOME,k)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA,RHO,NPAR
      use gloVars, only : UR,UI,PI,Z0
!     use debugStuff
      use resultvars, only : Punto
      use sourceVars, only : tipofuente!,PW_theta
      implicit none
      
      integer, intent(in) :: i_zF,tam
      complex*16,    intent(inout), dimension(tam) :: this_B
      real*8,     intent(in)    :: k
      complex*16, intent(in)    :: cOME
      type(Punto),intent(in),target    :: pXi
      
      integer, pointer :: e_f
      real*8, pointer :: z_f,cose,seno
      logical, pointer :: fisInterf
      integer :: iIf,nInterf, el_tipo_de_fuente
      real    :: SGNz
      complex*16 :: nu,DEN,omeBet, argum,sincmod
      real*8     :: errT = 0.0001_8
      ! una para cada interfaz (la de arriba [1] y la de abajo [2])
      real*8,     dimension(2) :: z_loc
      complex*16, dimension(2) :: G22,S22,enuz,nuz,sincNu
      real*8 :: a
      
      this_B = Z0
      sincNu = Z0
      e_f => pXi%layer
      z_f => pXi%center%z
      fisInterf => pXi%isOnInterface
      nInterf = 2
      if (fisInterf) then
         nInterf = 1
      end if 
      z_loc(1:2) = 0.0_8
      if (e_f .ne. 0) z_loc(1) = Z(e_f) - z_f !downward (-)
      if (e_f .ne. N+1) z_loc(2) = Z(e_f+1) - z_f !upward (+)
      
      el_tipo_de_fuente = 2 ! fuente segmento (para ibem)
      if (i_zF .eq. 0) el_tipo_de_fuente = tipofuente !(puntual u onda plana)
      !  0  puntual 
      !  2  segmento
      
      G22=Z0;S22=z0 
      DEN = (4.0*AMU(e_f)*UI*PI)
      omeBet = cOME**2.0/BETA(e_f)**2.0
      nu = sqrt(omeBet - k**2.0)
      if(aimag(nu).gt.0.0_8)nu= conjg(nu)
      
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
      
        a = pxi%length * 0.5_8 ! 2a=lenght
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
      do iIf = 1, nInterf
          if (abs(z_loc(iIf)) > errT ) then
            SGNz = real(z_loc(iIf) / ABS(z_loc(iIf)),4)
          else
            SGNz = 0.0
          end if
          
          G22(iIf) = UR/DEN * sincNu(iIf) / nu
          S22(iIf) = -UR / (4.0_8*pi) * sincNu(iIf) * SGNz
          
          
          if (fisInterf) then
            S22(1) = + cmplx(1.0 / (2.0 * PI ),0.0,8)
            G22(1) = 0
          end if
      end do !iIf
      
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
      end subroutine SHvectorB_force
      
      FUNCTION SINCMOD(ARG,ARGZ)  
      use glovars, only : UI    
      COMPLEX*16 :: SINCMOD,ARG, ARGZ
      SINCMOD=EXP(-UI*ARGZ)
      IF(ABS(ARG).LE.0.0001_8)RETURN
      SINCMOD=(EXP(UI*(ARG-ARGZ))-EXP(-UI*(ARG+ARGZ)) )/2.0_8/UI/ARG
      END function SINCMOD
! G_stra - coefficients
      ! coeficientes de las ondas planas emitidias en cada interfaz 
      ! para representar el campo difractado por estratigrafía
      function PSVdiffByStrata(coefOndas_PSV,z_i,e,cOME_i,k)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : PrintNum,verbose,UI,UR,Z0
      use resultVars, only : MecaElem
      implicit none
      
      complex*16, dimension(1:5) :: PSVdiffByStrata
!     type (MecaElem)              :: PSVdiffByStrata
      real*8, intent(in)           ::  z_i,k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e!,outpf,mecStart,mecEnd
      complex*16, dimension(1:4*N+2),intent(in) :: coefOndas_PSV
      complex*16 :: gamma,nu,xi,eta
      complex*16 :: egammaN,enuN,egammaP,enuP
      complex*16, dimension(2,4) :: subMatD
      complex*16, dimension(3,4) :: subMatS
      complex*16, dimension(4,4) :: diagMat
      complex*16, dimension(1:4) :: coeffsPSV
      complex*16, dimension(1:2) :: resD
      complex*16, dimension(1:3) :: resS
      integer :: i !#< b
      if (verbose > 4) then
       write(PrintNum,'(a,F7.3,a,F12.7,a,F10.2,a,I0)') & 
                    "PSVdiffByStrata at w:", & 
                    real(cOME_i),"k=",k," z_i{",z_i,"} e=",e
      end if !#> 
      PSVdiffByStrata = 0
!     PSVdiffByStrata%Rw(1:5) = z0
      gamma = sqrt(cOME_i**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
      nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
      if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
      if(aimag(nu).gt.0.0)nu= conjg(nu)
      xi = k**2.0 - nu**2.0
      eta = 2.0*gamma**2.0 - cOME_i**2.0 / BETA(e)**2.0
          
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
         diagMat = RESHAPE((/ egammaN, Z0, Z0, Z0, & 
                              Z0,    enuN, Z0, Z0, & 
                              Z0, Z0, egammaP, Z0, & 
                              Z0, Z0, Z0, enuP /), &
                           (/ 4,4 /))
      !coeficientes de las ondas en el estrato
      i = 0
      if (z(0) .lt. 0.0) i = 2
      
        if (e .eq. 0) then! semiespacio de arriba
          coeffsPSV(1:2) = (/z0,z0/)
          coeffsPSV(3:4) = coefOndas_PSV(1 : 2)
        elseif (e .eq. N+1) then ! semiespacio de abajo
          coeffsPSV(1:2) = coefOndas_PSV(4*(e-1)+1+i: 4*(e-1)+2+i)
          coeffsPSV(3:4) = (/z0,z0/)
        else! estrato
          coeffsPSV(1:4) = coefOndas_PSV(4*(e-1)+1+i : 4*(e-1)+4+i)
        end if      
 
      ! desplazamientos
        ! {W}
        ! {U}
        subMatD = RESHAPE((/ -gamma,-k*UR,-k*UR,nu, & 
                              gamma,-k*UR,-k*UR,-nu /), &
                           (/ 2,4 /))

        subMatD = UI * subMatD
        subMatD = matmul(subMatD,diagMat)
        ! Pdown Sdown Pup Sup
        resD = matmul(subMatD, coeffsPSV)
        
        PSVdiffByStrata(1) = resD(1) !W
        PSVdiffByStrata(2) = resD(2) !U
       
      ! esfuerzos
        subMatS = RESHAPE((/ xi,      -2.0*k*gamma,     eta,     &
                          -2.0*k*nu,     -xi,        2.0*k*nu,   &
                           xi,       2.0*k*gamma,     eta,       &
                           2.0*k*nu,     -xi,       -2.0*k*nu /),&
                           (/3,4/))      
     
        subMatS = amu(e) * subMatS
        subMatS = matmul(subMatS,diagMat)
        resS = matmul(subMatS, coeffsPSV)
        
        PSVdiffByStrata(3) = resS(1) !s33
        PSVdiffByStrata(4) = resS(2) !s31
        PSVdiffByStrata(5) = resS(3) !s11 
      end function PSVdiffByStrata
      
      function SHdiffByStrata(coefOndas_SH,  & 
                               z_i,e,cOME_i,k,mecStart,mecEnd,outpf)
      use soilVars !N,Z,AMU,BETA,ALFA,LAMBDA
      use gloVars, only : verbose,UI,UR,Z0
      use resultVars, only : MecaElem
      implicit none
      type (MecaElem)              :: SHdiffByStrata
      real*8, intent(in)           :: z_i
      real*8, intent(in)           :: k
      complex*16, intent(in)       :: cOME_i  
      integer, intent(in)          :: e,outpf,mecStart,mecEnd
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
      if (verbose > 4) then
       write(outpf,'(a,F7.3,a,F12.7,a,F10.2,a,I0)') & 
                    "SHdiffByStrata at w:", & 
                    real(cOME_i),"k=",k," z_i{",z_i,"} e=",e
      end if
      SHdiffByStrata%Rw_sh(1:3) = z0
      ! algunas valores constantes para todo el estrato
      nu = sqrt(cOME_i**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
      if(aimag(nu).gt.0.0)nu= conjg(nu)
          !downward waves
          if (e /= 0) then !(radiation condition upper HS)
            enuN = exp(-UI * nu * abs(z_i-Z(e)))
          else
            enuN = Z0
          end if
          !upward waves 
          if (e /= N+1) then !(radiation condition)
            enuP = exp(-UI * nu * abs(z_i-Z(e+1)))
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
      if (mecStart .eq. 1)then
        subMatDsh = RESHAPE((/ UR,UR /), (/ 1,2 /))
        subMatDsh = matmul(subMatDsh,diagMatSH)
        resDsh = matmul(subMatDsh, coeffsSH)
        SHdiffByStrata%Rw_sh(1) = resDsh(1,1) !V
      end if ! desplazamientos
      
      ! esfuerzos
      if (mecEnd .eq. 3) then
        subMatSsh = RESHAPE((/ -UI*nu,-UI*k,UI*nu,-UI*k/),(/2,2/)) 
        subMatSsh = amu(e) * subMatSsh
        subMatSsh = matmul(subMatSsh,diagMatSH)
        resSsh = matmul(subMatSsh, coeffsSH)
        SHdiffByStrata%Rw_sh(2) = resSsh(1,1) !s32
        SHdiffByStrata%Rw_sh(3) = resSsh(2,1) !s12
      end if ! esfuerzos
          
      end function SHdiffByStrata
      
! G_full space 
      Subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
      ! Sanchez-Sesma y Campillo, 1991.  Mismo resultado que:
      ! Kaussel, Fundamental solutions in elastodynamics... pag 38
      use soilVars ,only : alfa0,beta0,Lambda0,AMU0,alfa,beta,amu,lambda,rho,N,Z
      use gloVars, only:UI,UR,one,z0
      use hank !     use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: tipofuente, PW_pol
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
      complex*16 :: ome
      integer :: i,j
      integer, pointer :: e
      real*8 :: nX(2)
      integer :: iGq,nGq
      real*8, pointer :: xf,zf,GqC
      real*8 :: deltaij
      real*8,dimension(2) :: theta
      integer :: el_tipo_de_fuente
      integer, target :: estrato
      logical :: shouldI,XinoEstaEnInterfaz,usarGreenex!,estratosIguales
      
      FF%W=z0;FF%U=z0;FF%Tz=z0;FF%Tx=z0
!     estratosIguales = .false.
      XinoEstaEnInterfaz = .false.
      usarGreenex = .false.
      shouldI = .false.
      j = dir_j ! ***********  dir_j = 3  (vertical)
      if(j .eq. 3) j = 2 ! para coincidir con los indicies
      nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
      xf => one; zf => one
      
      if (i_zF .eq. 0) then ! es la fuente real
       el_tipo_de_fuente = tipofuente !(puntual:0 u onda plana:1)
!      if (p_X%region .eq. pXI%region) then 
         shouldI = .true.
         if (pXi%region .eq. 2) then
          XinoEstaEnInterfaz = .true.
          estrato = N+2
          e => estrato
         else 
          if (p_x%layer .eq. pXi%layer) shouldI = .true.
          if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
          estrato = p_x%layer
          e => p_x%layer
         end if
!      end if
      else ! una fuente virtual
        el_tipo_de_fuente = 2 !(fuente segmento)
        XinoEstaEnInterfaz = .true.
        if (i_zF .eq. -1) then !(campo refractado en inclusion)
          shouldI = .true.
          estrato = N+2
          e => estrato
        else !(campo refractado en medio estratificado)
          if (p_x%layer .eq. pXi%layer) then 
            shouldI = .true.
            estrato = p_x%layer
            e => p_x%layer
          end if
        end if
      end if
      !
      if (shouldI) then
      
!     if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
!     if (pXi%isOnInterface .eqv. .false.)XinoEstaEnInterfaz = .true.
!     j = dir_j ! ***********  dir_j = 3  (vertical)
!     if(j .eq. 3) j = 2 ! para coincidir con los indicies
!     
!     el_tipo_de_fuente = 2 !(fuente segmento)
!     if (i_zF .eq. 0) el_tipo_de_fuente = tipofuente !(puntual:0 u onda plana:1) 
!     if (i_zF .eq. -1) then !(campo refractado en inclusion)
!       estratosIguales = .true.
!       XinoEstaEnInterfaz = .true.
!     end if!  
!     
!     if ((i_zF .eq. 0) .and. (pXi%region .eq. 2)) then
!       if (p_X%region .eq. pXI%region) estratosIguales = .true.
!       XinoEstaEnInterfaz = .true.
!       estrato = N+2
!       e => estrato
!       nx(1) = p_X%normal%x; nx(2) = p_X%normal%z
!     end if
!     
!     if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then 
!       if (p_x%layer   .eq. pXi%layer) estratosIguales = .true.
!     end if!
!     
!     if (estratosIguales .eqv. .true.) then !should I
!     estrato = p_x%layer
!     e => p_x%layer
!     nx(1) = p_X%normal%x; nx(2) = p_X%normal%z 
!     xf => one; zf => one !para que no chiste el compilador
!     if (i_zF .eq. -1) then !(en la inclusión)
!       estrato = N+2
!       e => estrato
!       nx(1) = p_X%normal%x; nx(2) = p_X%normal%z !las normales sin voltear
!     end if!
        
      if ((i_zF .eq. 0) .and. (el_tipo_de_fuente .eq. 1)) then ! onda plana !
!      if (p_x%isOnInterface .eqv. .true.) return  ! creo                   !
       ome = real(cOME) * UR                                                !
       if (PW_pol .eq. 1) then                                              !
        c = beta0(N+1) !SV                                                  !
        theta(1) = cos(pxi%gamma)                                           !
        theta(2) = sin(pxi%gamma)                                           !
       elseif (PW_pol .eq. 2) then                                          !
        c = alfa0(N+1) !P                                                   !
        theta(1) = sin(pxi%gamma)                                           !
        theta(2) = -cos(pxi%gamma)                                          !
       end if                                                               !
        kx = ome/c * sin(pxi%gamma)                                         !
        kz = ome/c * cos(pxi%gamma)                                         !
        FF%U = (theta(1))* exp(UI * kz * (p_x%center%z - Z(N+1))) &         !
              * exp(-UI * kx * (p_x%center%x - pXi%center%x))               !
        FF%W = (theta(2))* exp(UI * kz * (p_x%center%z - Z(N+1))) &         !
              * exp(-UI * kx * (p_x%center%x - pXi%center%x))               !
        szz = UI * ( &                                                      !
                    ( FF%W * kz * (LAMBDA0(e) + 2.0* AMU0(e))) &            !
                  - ( FF%U * kx * LAMBDA0(e)))                              !
        szx = UI * AMU0(e) * ( kz * FF%U - kx * FF%W )                      !
        sxx = UI * ( &                                                      !
                  - ( FF%U * kx * (LAMBDA0(e) + 2.0*AMU0(e))) &             !
                  + ( FF%W * kz * LAMBDA0(e)))                              !
        FF%Tx = sxx * nx(1) + szx * nx(2)                                   !
        FF%Tz = szx * nx(1) + szz * nx(2)                                   !
        return                                                              !
      end if! fin onda plana -´-´-´-´-´-´-´-´-´-´-´-´-´-´--´-´ fin onda plana          
      if (XinoEstaEnInterfaz .eqv. .true.) then 
        xf => pXi%center%x
        zf => pXi%center%z
        r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.) 
        
        ! para funciones de Greeen en frontera con cont. de desplazamiento
        if ((pXi%isboundary .eqv. .true.) .and. & 
            (p_x%boundaryIndex .eq. pXi%boundaryIndex)) usarGreenex = .true.

        ! para campo cercano por fuente virtual
        if ((pXi%isboundary .eqv. .true.) .and. & 
            (abs(r) .lt. pXi%length/2)) usarGreenex = .true.
        ! si fuente real y cilindrica entonces no
        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false. 
        
        if (p_X%isOD) usarGreenex = .false.
      if (usarGreenex .eqv. .false.) then
      if (el_tipo_de_fuente .eq. 0) then !; print*,"fuente real y cilindrica"
        nGq = 1
        GqC => ONE
      else !.eq. 2 ; print*,"fuente segmento (real o virtual para IBEM)"
        if (pXi%isBoundary) nGq = Gquad_n ! IBEM
        if (pXi%isSourceSegmentForce) nGq = Gquad_n ! fuente segmento
        if (nGq .ne. Gquad_n) stop "chin 6600"
      end if
      do iGq = 1,nGq
        if (nGq .gt. 1) then
          xf => pXi%Gq_xXx_coords(iGq,1)
          zf => pXi%Gq_xXx_coords(iGq,2)
          GqC => pXi%Gq_xXx_C(iGq) 
        end if
      r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
      !print*,iGq,xf,zf,r
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
      end if !mece5
      end do ! iGq
      
      else !...............................                          
        call greenexPSV(FF,dir_j,p_X,pXi,estrato,cOME) ! W,U
      end if ! greenex o gaussiana
      end if !on the interface
      end if !should I?
      end subroutine FFpsv

      subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
      use soilVars ,only : amu,beta,N,z
      use gloVars, only:UR,UI,z0,ONE
      use hank !use specfun
      use resultvars, only : Punto,FFres
      use sourceVars, only: tipofuente
      use Gquadrature, only : Gquad_n
      implicit none
      type (FFres), intent(out) :: FF !V,Ty
      integer,    intent(in)     :: mecS,mecE, i_zF
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16, intent(in)     :: cOME
      
      real*8 :: r,gamma(2)!,longonda
      complex*16 :: omeS
      complex*16 :: H0s,H1s,kx,kz,szy,sxy
      integer, pointer :: e
      integer :: iGq,nGq
      real*8 ::  nX(2)
      real*8, pointer :: xf,zf,C
      integer :: el_tipo_de_fuente
      integer, target :: NmasDos
      logical :: estratosIguales,noEstaEnInterfaz,usarGreenex
      complex*16 :: greenexG22
!     integer :: NM
!     complex*16, dimension(0:1) :: CBJ,CDJ,CBY,CDY
      
      FF%V=z0;FF%Ty=z0
      estratosIguales = .false.
      noEstaEnInterfaz = .false.
      usarGreenex = .false.
      if (p_x%layer .eq. pXi%layer) estratosIguales = .true.
      if (pXi%isOnInterface .eqv. .false.)noEstaEnInterfaz = .true.
      el_tipo_de_fuente = 2 ! fuente segmento para IBEM
      if (i_zF .eq. 0) el_tipo_de_fuente = tipofuente !(puntual u onda plana) 
      if (i_zF .eq. -2) then !(campo refractado en frontera d2R y d1R) 
        estratosIguales = .true.
        noEstaEnInterfaz = .true.
      end if! 
      if (i_zF .eq. -1) then !(campo refractado en frontera d1R si voltear normal) 
        estratosIguales = .true.
        noEstaEnInterfaz = .true.
      end if! 
      if (estratosIguales .eqv. .true.) then !should I
      e => p_x%layer
      nx(1) = p_X%normal%x;nx(2) = p_X%normal%z
      xf => one;zf => one
      if (i_zF .eq. -2) then
        NmasDos = N+2
        e => NmasDos !(en la inclusión)
        nx(1) = -p_X%normal%x; nx(2) = -p_X%normal%z !las normales volteadas
      end if!
      if (i_zF .eq. -1) then
        NmasDos = N+2
        e => NmasDos !(en la inclusión)
        nx(1) = p_X%normal%x; nx(2) = p_X%normal%z !las normales sin voltear
      end if!
      
      if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then ! es una onda plana incidente
        if (p_x%isOnInterface .eqv. .true.) return  ! creo
        kx = cOME/beta(N+1)*sin(pxi%gamma)
        kz = cOME/beta(N+1)*cos(pxi%gamma)
        ! SV:
        FF%V = (UR)* exp(UI * kz * (p_x%center%z - Z(N+1))) & 
              * exp(-UI * kx * (p_x%center%x - pXi%center%x))
!       szy, sxy
        sxy = amu(N+1) * FF%V * (-UI * kx)
        szy = amu(N+1) * FF%V * (UI * kz)
        FF%Ty = sxy * nx(1) + szy * nx(2)
        return
      end if! onda plana
      if (noEstaEnInterfaz .eqv. .true.) then !..............................
        xf => pXi%center%x
        zf => pXi%center%z
        r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
!       if (abs(r) .lt. pXi%length/2) usarGreenex = .true.
        if (p_x% boundaryIndex .eq. pXi%boundaryIndex) usarGreenex = .true.
!       if (abs(r) .gt. 0.001_8) usarGreenex = .false.
        if (el_tipo_de_fuente .eq. 0) usarGreenex = .false.
!       if (p_X%region .ne. 2) usarGreenex = .false. ! 'incl'
        if (p_X% tipoFrontera .ne. 1) usarGreenex = .false. ! si no es frontera de continuidad
!       usarGreenex = .false. !*+++++++++++++++++++++++++++++++++++++++++

      if (el_tipo_de_fuente .eq. 0) then                                !
        C => ONE                                                        ! 
        nGq = 1                                                         !
      else !; print*,"fuente segmento (real o virtual para IBEM)"       !
        if (pXi%isBoundary) nGq = Gquad_n ! IBEM                        !
        if (pXi%isSourceSegmentForce) nGq = Gquad_n ! fuente segmento   !
      end if                                                            !
      do iGq = 1,nGq !.............................................     !
      if (nGq .gt. 1) then                                        !     !
        xf => pXi%Gq_xXx_coords(iGq,1)                            !     !
        zf => pXi%Gq_xXx_coords(iGq,2)                            !     !
        C => pXi%Gq_xXx_C(iGq)                                    !     !
      end if                                                      !     !
      r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)     !     !
      gamma(1) = (p_X%center%x - xf) / r ! gamma x                !     !
      gamma(2) = (p_X%center%z - zf) / r ! gamma z                !     !
      omeS = cOME * r / beta(e)                                   !     !
      call hankels(omeS,H0s,H1s) ! Hankel de segunda especie      !     !  
!     call CJYNB(1,omeS,NM,CBJ,CDJ,CBY,CDY)
!     H0s = CBJ(0)-UI*CBY(0)
!     H1s = CBJ(1)-UI*CBY(1) 
      if(mecS .eq. 1) FF%V = FF%V + (-UI/(4.*amu(e))*H0s) * C     !     !
      if(mecE .eq. 3) then                                        !     !
!     FF%s32 = FF%s32 + (UI*cOME*(p_x%center%z-zf)) &             !     !
!                        /(4.0_8*beta(e)*r)*H1s * C               !     !
!     FF%s12 = FF%s12 + (UI*cOME*(p_x%center%x-xf))/ &            !     !
!                         (4.0_8*beta(e)*r)*H1s * C               !     !
      FF%Ty = FF%Ty + (UI*cOME/beta(e)/4.0*H1s* &                 !     !
               (gamma(1)*nx(1) + gamma(2)*nx(2))) * C             !     !
      end if !mec3                                                !     !
      end do ! iGq ................................................     !
        
      if (usarGreenex .eqv. .true.) then !...............................
        xf => pXi%center%x
        zf => pXi%center%z
        r = sqrt((p_x%center%x-xf)**2. + (p_x%center%z-zf)**2.)
           FF%V = greenexG22(cOME/beta(e),r,pXi%length,e) 
      end if ! greenex o gaussiana .....................................!
      end if ! should I? ....................................................
      end if
      end subroutine FFsh
      
      subroutine greenexPSV(FF,dir_j,p_X,pXi, e,cOME)
      use glovars, only : UR,UI,pi,z0
      use resultvars, only : Punto,FFres
      use soilVars ,only : alfa,beta,rho
      implicit none
      type (FFres), intent(out) :: FF !U y W dados dir_j 1 o 3
      integer,intent(in) :: dir_j,e
      type(Punto),intent(in), pointer :: p_X,pXi
      complex*16,intent(in) :: cOME !frecuencia compleja (rad/s)
      complex*16 :: AQA,AKA,BEALF,BEALF2,ALFBE
      complex*16 :: argk,ark2,h0kr,argq,arq2,h0qr,AA,BB,h2kr,h2qr
      complex*16 :: GN11,GN31,GN13,GN33
      real*8 :: G1,G3,depi,RNM,DSM,c1,c2,c13,c23,geu
!     stop "Greeneer"
      FF%W=z0;FF%U=z0     
      GEU=0.5772156649_8
      
      !come = ur * real(comein)
      
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
         
         if (dir_j .eq. 1) then
           FF%U = GN11
           FF%W = GN13
         else if (dir_j .eq. 3) then
           FF%U = GN31
           FF%W = GN33
         end if
!     print*,dir_j,RNM,p_x%center,pXI%center,FF%U,FF%W
      end subroutine greenexPSV
      
      function greenexG22(k,rij,dr,e)
             
      !funcion de Green G22 analitica en la fuente
      ! se usan las series ascendentes de las funciones de Bessel
      use glovars, only : UR,UI,pi
      use soilVars ,only : amu
      implicit none
      complex*16 :: k ! numero de onda
      real*8 :: rij ! distancia (muy pequeña o 0)
      real*8 :: dr ! longitud del segmento de integración
      integer :: e
      
      real*8 :: c1,c2,c13,c23,geu,depi
      complex*16 :: argk,ark2,h0kr,greenexG22
      
      DEPI=2.0/PI
      c1=1.0+rij/dr*2.0
      c2=1.0-rij/dr*2.0
      c13=c1**3.0
      c23=c2**3.0
      geu=0.5772156649_8
      argk= k * dr
      ark2=argk**2.0
      
      h0kr = ur*(1.0-ark2*(c13+c23)/96.0) &
          - ui*depi*( geu-1.0 + 0.5*c1*log(argk*c1/4.0) &
                              + 0.5*c2*log(argk*c2/4.0) &
          +(4.0/3.0-geu)*ark2*(c13+c23)/96.0 &
            - ark2*c13*log(argk*c1/4.0)/96.0 &
            - ark2*c23*log(argk*c2/4.0)/96.0)
            
      greenexG22 = - ui/(4.0*amu(e)) * h0kr * rij
      
      end function greenexG22
                  
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
      use sourceVars, only: tipofuente
      implicit none
      interface
         subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
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
      nf(1) = pXi%normal%x
      nf(2) = pXi%normal%y
      nf(3) = pXi%normal%z
      ! si es onda plana no importa nf conque sea unitario.
      if (tipofuente .eq. 1) nf(1:3) = 1 ! fuente real y onda plana
      
      if (pXi%region .eq. 2) then 
        auxK(1,mecS:mecE) = 0
        if (p_x%boundaryIndex .le. n_top_sub) return
        nf(1:3) = -nf(1:3)
      end if !
      
      if (dir_j .eq. 2) then ! SH
        call FFsh(0,FF,p_X,pXi,cOME,1,3)
      !  | Ty |
        trac0vec(p_x%boundaryIndex) = & 
        trac0vec(p_x%boundaryIndex) - (& 
          (TractionSH(auxk(1,2),auxk(1,3),p_x%normal) + FF%Ty) * nf(dir_j))
           
        if (p_X%tipoFrontera .eq.1) then !los desplazamientos
      !  |   V   |
        trac0vec(p_x%boundaryIndex+ n_con_sub) = & 
        trac0vec(p_x%boundaryIndex+ n_con_sub) - & 
          (auxk(1,1) + FF%V)* nf(dir_j)
        end if
      else !PSV 
        call FFpsv(0,FF,dir_j,p_X,pXi,cOME,1,5)
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
      use sourceVars, only : Po,nFuentes,tipofuente
      use resultvars, only : trac0vec,Punto,FFres,& 
            boupoints,n_con_sub,n_top_sub,n_val_sub
      implicit none
      interface
        subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface
      integer, intent(in) :: dir_j
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      type(Punto), pointer :: p_X
      real*8 :: nf(3)
      integer :: iPo,ip_X
      
      do iPo =1,nFuentes
!     print*,",ipo",ipo
      nf(1) = Po(iPo)%normal%x
      nf(2) = Po(iPo)%normal%y
      nf(3) = Po(iPo)%normal%z
!     print*,nf
      if (tipofuente .eq. 1) nf(1:3) = 1 ! fuente real y onda plana
      
      if (dir_j .eq. 2) then ! SH
        stop "termIndepR _ falta"
      else !PSV
        do ip_X = n_top_sub + n_con_sub + 1, n_top_sub + n_con_sub + n_val_sub
          nullify(p_X)
          p_X => boupoints(ip_X)
!         print*,ip_x,p_x%center,p_x%tipoFrontera
          call FFpsv(0,FF,dir_j,p_X,Po(iPo),cOME,3,5)
      !  | Tx |
      !  | Tz |
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 0))+ 2* n_con_sub) - (FF%Tx * nf(dir_j)) 
          
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) = &
        trac0vec((p_x%boundaryIndex *2 - (1 - 1))+ 2* n_con_sub) - (FF%Tz * nf(dir_j))
        end do
      end if !dir_j
      end do !iPo
      end subroutine termIndepR
!#< r IBEM - reg E        !#>
      subroutine fill_ibemMat(i_zF,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres,ibemMat,n_con_sub
      use glovars, only : z0,UR
      implicit none
      interface
        subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface
      
      type(Punto), pointer :: pXi,p_X
      integer :: mecS,mecE,i_zF,dj
      integer, intent(in) :: dir_j
      complex*16, dimension(1,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      complex*16 :: TractionPSV, TractionSH
      if (p_X%tipoFrontera .gt. 1) return !o sea =2: (tracciones libres en inclusión)
!     print*,"ibemmat"
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
      
      end subroutine fill_ibemMat
            
      subroutine fill_diffbyStrata(i_zf,J,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto,FFres, nIpts
      use waveNumVars, only : NMAX,k_vec,dk,vecNK
      use meshvars, only: npixX,MeshMaxX,MeshMinX 
      use wavelets !fork
      use soilvars, only:alfa0,beta0,N
      use sourceVars, only: tipofuente, PW_pol
      use glovars, only : UI
      use peli, only : fotogramas_Region
      implicit none
      interface   
       subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface 
      type(Punto), pointer :: p_X,pXi,p_Xmov
      type(Punto),target :: p_xaux
      integer :: i,ik,i_zf,iMec,mecS,mecE,J,mecaElemEnd,po,ne
      integer, intent(in) :: dir_j
      complex*16, dimension(2*nmax,mecS:mecE), target :: auxK
      complex*16, dimension(2*nmax,mecS:mecE) :: auxKmo
      real*8 :: mov_x
      complex*16 :: kx
      complex*16, intent(in)  :: cOME
      type(FFres),target :: FF
      real*8 :: nf(3)
      nf(1) = pXi%normal%x
      nf(2) = pXi%normal%y
      nf(3) = pXi%normal%z
      mecaElemEnd = 2 !PSV
      if (dir_j .eq. 2) mecaElemEnd = 1 !SH
      if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) nf(1:3) = 1.0 !fuente real y onda plana
      if ((i_zf .eq.0) .and. (abs(nf(dir_j)) .lt. 0.0001)) return !fuente real fuerza; componente nulo
      ! diffracted field due to the stratification U,V,W or G
      if (p_x%guardarMovieSiblings) then
               p_xaux%center%x = p_x%center%x
               p_xaux%center%z = p_x%center%z
               p_xaux%normal = p_x%normal
               p_xaux%isOnInterface = p_x%isOnInterface
               p_xaux%layer = p_x%layer
               p_Xmov => p_xaux
        po = min(int(vecNK(J)*1.25),nmax); ne = 2*nmax-(po-2)
        do i = 1,npixX ! para cada hermanito
          mov_x = MeshMinX + (MeshMaxX - MeshMinX)/(npixX-1) * (i-1)! la coordenada x
          p_Xmov%center%x = mov_x
          
          ! si no es de la region 0 no vale la pena calcularlo
          if(fotogramas_Region(p_x%pointIndex-nIpts,i) .ne. 1) then !'estr'
!           print*,"cycled [",p_Xmov%center%x,",",p_xaux%center%z,"]"
            cycle
          end if
          
         do imec = 1,mecaElemEnd
           ! si es incidencia de una onda plana para una k y cycle imec
           if ((i_zF .eq. 0) .and. (tipofuente .eq. 1)) then ! onda plana
             if (dir_j .eq. 2) then
               if(PW_pol .eq. 3) kx = real(cOME)/beta0(N+1)*sin(pXi%gamma)
             else
               if(PW_pol .eq. 1) kx = real(cOME)/beta0(N+1)*sin(pXi%gamma)
               if(PW_pol .eq. 2) kx = real(cOME)/alfa0(N+1)*sin(pXi%gamma)
             end if
             auxkmo(1,imec) = auxk(1,imec) * & 
             exp(-UI*kx*(p_Xmov%center%x))
             cycle !imec 
           end if ! onda plana-------------------------------------------
           ! incidencia de una onda cilindrica:
           auxKmo(po:ne,imec) = 0
           do ik = 1,po!2*Nmax
            auxKmo(ik,imec) = auxk(ik,imec) * &
            exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_Xmov%center%x), 8))
           end do!  ik
           do ik = ne,2*Nmax
            auxKmo(ik,imec) = auxk(ik,imec) * &
            exp(cmplx(0.0_8, (-1.0_8)*k_vec(ik)*(p_Xmov%center%x), 8))
           end do!  ik
           !K->X
            !auxKmo(:,imec) = auxKmo(:,imec) * (sqrt(2.0*nmax) * dk) !backward
            !call fork(2*nmax, auxKmo(:,iMec),+1,0,6)
            auxKmo(:,imec) = FFTW(2*nmax,auxKmo(:,iMec),+1,dk)
         end do! imec
          
          ! almacenar
          if (dir_j .eq. 2) then ! SH
          call FFsh(i_zf,FF,p_Xmov,pXi,cOME,1,1)
            if(i_zf .eq.0) then
              p_x%Wmov(J,3,i) = & 
              p_x%Wmov(J,3,i) + (auxKmo(1,1) + FF%V) * nf(dir_j) ! V
            else
              pXi%Gmov(p_x%pointIndex-nIpts,3,dir_j,i) = & 
              pXi%Gmov(p_x%pointIndex-nIpts,3,dir_j,i) + auxKmo(1,1) + FF%V ! V
            end if
          else !PSV
          call FFpsv(i_zF,FF,dir_j,p_Xmov,pXi,cOME,1,1)
            if(i_zf .eq.0) then
              p_x%Wmov(J,1,i) = p_x%Wmov(J,1,i) + (auxKmo(1,1) + FF%W) * nf(dir_j) !W
              p_x%Wmov(J,2,i) = p_x%Wmov(J,2,i) + (auxKmo(1,2) + FF%U) * nf(dir_j) !U
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
       call FFsh(i_zf,FF,p_X,pXi,cOME,1,3) !incidencia directa
       if(i_zf .eq.0) then
!         p_x%W(J,3) = p_x%W(J,3) + (auxk(1,1) + FF%V) * nf(dir_j) ! V
          p_x%resp(J)%V = p_x%resp(J)%V + (auxk(1,1) + FF%V) * nf(dir_j) ! V
          !                     s12                       s32
          p_x%resp(J)%Ty = ((auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + &
                             FF%Ty) * nf(dir_j) + p_x%resp(J)%Ty
       else
          pXi%G(p_x%pointIndex,3,dir_j) = & 
          pXi%G(p_x%pointIndex,3,dir_j) + auxK(1,1) + FF%V ! V
          pXi%G(p_x%pointIndex,6,dir_j) = & 
          pXi%G(p_x%pointIndex,6,dir_j) + & 
          (auxk(1,3)* p_x%normal%x + auxk(1,2)* p_x%normal%z) + FF%Ty ! Ty
       end if
      else !PSV
       call FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,1,5)  !incidencia directa
       
       if(i_zf .eq.0) then
        if (pXi%region .ne. p_x%region) return
         p_x%resp(J)%W = p_x%resp(J)%W + (auxK(1,1) + FF%W) * nf(dir_j) !W
         p_x%resp(J)%U = p_x%resp(J)%U + (auxK(1,2) + FF%U) * nf(dir_j) !U
         p_x%resp(J)%Tz = p_x%resp(J)%Tz + &
         ((auxk(1,4)* p_x%normal%x + auxk(1,3)* p_x%normal%z) + FF%Tz) * nf(dir_j)
         !      s31                       s33
         p_x%resp(J)%Tx = p_x%resp(J)%Tx + &
         ((auxk(1,5)* p_x%normal%x + auxk(1,4)* p_x%normal%z) + FF%Tx) * nf(dir_j)
         !      s11                       s31
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
       end if
      end if !dir_j
      end if ! guardarMovieSiblings
      end subroutine fill_diffbyStrata
                

      function TractionPSV(RW,normal,l)
      use resultVars, only: Punto3d
      implicit none
      complex*16 :: TractionPSV
      complex*16, intent(in) :: RW(3:5)
      type(Punto3d), intent(in) :: normal
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
      use resultVars, only: Punto3d
      implicit none
      complex*16 :: TractionSH
      complex*16, intent(in) :: s32,s12 
      type(Punto3d), intent(in) :: normal
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
         subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface
      integer, intent(in) :: ipXi,dir_j ! ip_X : n_topo+1,n_topo+n_cont
      complex*16, intent(in),target  :: cOME
      type(Punto), pointer :: pXi,p_X
      integer :: ip_X,dj
      type(FFres),target :: FF
!     real*8 :: r
!     complex*16 :: greenexG22
      
      nullify(pXi)
      pXi => boupoints(ipXi)
      
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
      use sourceVars, only : Po,nFuentes
      use soilvars, only : N
      implicit none
      interface
      subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface
      integer, intent(in) :: Jfrec,dir_j
      complex*16, intent(in)  :: cOME
      integer :: iP_x,iPXi
      type(Punto), pointer :: p_X,pXi,p_Xmov
      integer :: i,j,iPo
      type(Punto),target :: p_xaux
      type(FFres),target :: FF
      real*8 :: nf(3)
      ! ciclar en todos los puntos y averiguar se le corresponde
      do iP_x = 1, nIpts ! (todos los puntos receptores menos las fronteras)
        p_X => allpoints(iP_x)
        if (p_x%isboundary) stop "GreenReg_R: p_x is boundary "
        if (p_x%region .ne. 2) cycle !'incl'
        ! para todas las fuentes en la región R
      do iPXi = n_top_sub+1,nBpts
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
        end if !dir_j
      end do !iPXi
      
      ! y la fuente real
      do iPo = 1,nfuentes
      if (Po(iPo)%region .eq. 2) then
      pXi => Po(iPo)
      nf(1) = pXi%normal%x
      nf(2) = pXi%normal%y
      nf(3) = pXi%normal%z
      if (dir_j .eq. 2) then ! SH
       call FFsh(0,FF,p_X,pXi,cOME,1,3)
          p_x%resp(Jfrec)%V = p_x%resp(Jfrec)%V + (FF%V) * nf(dir_j) ! V
          p_x%resp(Jfrec)%Ty = p_x%resp(Jfrec)%Ty + (FF%Ty) * nf(dir_j)
      else !PSV
       call FFpsv(0,FF,dir_j,p_X,pXi,cOME,1,5) 
!        if (pXi%region .ne. p_x%region) cycle
         p_x%resp(Jfrec)%W = p_x%resp(Jfrec)%W + (FF%W) * nf(dir_j) !W
         p_x%resp(Jfrec)%U = p_x%resp(Jfrec)%U + (FF%U) * nf(dir_j) !U
         p_x%resp(Jfrec)%Tz = p_x%resp(Jfrec)%Tz + (FF%Tz) * nf(dir_j)
         p_x%resp(Jfrec)%Tx = p_x%resp(Jfrec)%Tx + (FF%Tx) * nf(dir_j)
      end if !dir_j
      end if !Po 2
      end do !iPo
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
            p_xaux%normal%x = 1.0_8 !
            p_xaux%normal%y = 1.0_8 !  no se usan 
            p_xaux%normal%z = 1.0_8 !
            p_xaux%region = 2; p_xaux%layer = N+2 !'incl'
            p_xaux%boundaryindex = -10000; p_xaux%isOD = .false.
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
                nf(2) = pXi%normal%y
                nf(3) = pXi%normal%z
                if (abs(nf(dir_j)) .lt. 0.0001) cycle
                if (dir_j .eq. 2) then ! SH
                  call FFsh(0,FF,p_Xmov,pXi,cOME,1,1)
         p_Xmov%Wmov(Jfrec,3,j) = p_Xmov%Wmov(Jfrec,3,j) + FF%V * nf(dir_j) !V
                else !PSV
                  call FFpsv(0,FF,dir_j,p_Xmov,pXi,cOME,1,1) 
         p_x%Wmov(Jfrec,1,j) = p_x%Wmov(Jfrec,1,j) + FF%W * nf(dir_j) !W
         p_x%Wmov(Jfrec,2,j) = p_x%Wmov(Jfrec,2,j) + FF%U * nf(dir_j) !U
                end if !dir_j
              end if ! Po 2
            end do ! iPo
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
      implicit none
      interface
      subroutine FFpsv(i_zF,FF,dir_j,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF
         integer,    intent(in)     :: mecS,mecE, i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME
         integer, intent(in) :: dir_j
         end subroutine FFpsv
         
         subroutine FFsh(i_zF,FF,p_X,pXi,cOME,mecS,mecE)
         use resultvars, only : Punto,FFres
         type (FFres), intent(out) :: FF !V,Ty
         integer,    intent(in)     :: mecS,mecE,i_zF
         type(Punto),intent(in), pointer :: p_X,pXi
         complex*16, intent(in)     :: cOME 
         end subroutine FFsh
      end interface
      integer, intent(in) :: J
      integer :: ren,renStep,col
      logical, intent(in) :: PSV
      type(Punto), pointer :: pXi,p_X
      integer :: iP_x,iPxi,ipxi_I,ipxi_F,dj,dir_j!,iPhi_I,iPhi_F
      type(FFres),target :: FF
!     print*,"ren=",ren !primer renglon del primer elemento de sobredeterminado
      do iP_x = 1,nIpts  !cada receptor X
      p_X => allpoints(iP_x) 
!     col = 1
      if (.not. p_X%isOD) cycle!; print*,"iP_X = ",iP_x,"  iP_XCen=",p_X%center
      if (PSV) then
       if (p_X% tipoFrontera .eq.0) then!; print*,"frontera tipo trontera libre en medio estrat"
          ipxi_I = 1
          ipxi_F = n_top_sub
          do iPxi = ipxi_I,ipxi_F!; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi)!; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2!; print*,"dir_j ",dir_j," dj ",dj
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
            col = pXi%boundaryIndex *2 -(2 - dj)
              ibemMat(ren,   col) = pXi%G(p_X%pointIndex,5,dir_j) !Tx
              ibemMat(ren+1, col) = pXi%G(p_X%pointIndex,4,dir_j) !Tz
            end do !dj
          end do !iPxi
            !  | Tx |
            !  | Tz |
              trac0vec(ren  ) = p_x%resp(J)%Tx
              trac0vec(ren+1) = p_x%resp(J)%Tz
          renStep = 2
       else if (p_X%tipoFrontera .eq.1) then !; print*,"frontera tipo continuidad"
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
            col = pXi%boundaryIndex *2 -(2 - dj)+ 2* n_con_sub
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
              trac0vec(ren  ) = p_x%resp(J)%Tx
              trac0vec(ren+1) = p_x%resp(J)%Tz
              trac0vec(ren+2) = p_x%resp(J)%W
              trac0vec(ren+3) = p_x%resp(J)%U
          renStep = 4
       else if (p_X%tipoFrontera .eq.2) then!; print*,"frontera tipo frontera libre en inclusión elástica"
          ipxi_I = n_top_sub + 1
          ipxi_F = n_top_sub + n_con_sub + n_val_sub
          do iPxi = ipxi_I,ipxi_F!; print*,"iPxi = ",iPxi
            pXI => boupoints(iPxi)!; print*,"pxiCen=",pXI%center
            do dir_j=1,3,2 !por fuerza hozintal 1 y vertical 3->2
            dj = dir_j; if (dj .eq. 3) dj = 2!; print*,"dir_j ",dir_j," dj ",dj
            call FFpsv(-1,FF,dir_j,p_X,pXi,cOME,3,5)
            !  |  Txx Txz  |
            !  |  Tzx Tzz  |
            col = pXi%boundaryIndex *2 -(2 - dj)+ 2* n_con_sub
              ibemMat(ren,   col) =  FF%Tx !Tx
              ibemMat(ren+1, col) =  FF%Tz !Tz
            end do !dj
          end do !iPxi
            !  | Tx |
            !  | Tz |
              trac0vec(ren  ) = p_x%resp(J)%Tx
              trac0vec(ren+1) = p_x%resp(J)%Tz
          renStep = 2
       end if !%tipoFrontera
      else !PSV / SH
         stop "overDetermineSystem: Falta SH sobreDet"
      end if !PSV/SH
      ren = ren + renStep
      end do !iP_x
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

! Sismo/Foto- gramas 
      subroutine W_to_t(W,nombre,yAx,iP,x_i,z_i)
      use waveNumVars, only : NFREC,DFREC, NPTSTIME, OMEI, t_vec
      use glovars
      use waveVars, only : dt,Uo,maxtime
      use ploteo10pesos
      use wavelets
      use resultvars, only : allpoints, Sabana, nSabanapts, nIpts, SabanaPlotIndividual
      use debugStuff
      implicit none
      integer ,intent(in) :: iP
      complex*16, dimension(Nfrec+1), intent(in) :: W !espectro
      real*8, intent(in) :: x_i,z_i
      character(LEN=3)   :: nombre
      character(LEN=100) :: titleN,xAx,yAx, CTIT
      
      complex*16, dimension(NPTSTIME) :: S
      
      character(LEN=9)   :: logflag
      integer :: i,n_maxtime
      real*8 :: this_maxtime
      
      if (Verbose .ge. 4) print*,"at W_to_t"
      if (developerfeature .eq. 1) then
       if (nombre(1:1) .eq. "T") then
        write(xAx,'(a)') "$a \omega c_p^{-1}$"
        S = z0
        S(1:nfrec)= W(1:nfrec:+1)
        write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),']_sin.pdf'
        logflag = 'none     '
        call plotSpectrum(S(1:100),real(developerAUXvec(2,2),4),& 
             100,100,titleN,xAx,yAx,logflag,1200,800,real(developerAUXvec(2,2)*100,4))
        !guardar en texto
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),']_sin.txt'
         
         OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         write(3211,'(I0,A)') 1," espectro antes de convolución con funcion de amplitud"
         write(3211,'(a)') "dfrec="
         write(3211,'(F15.8)') DFREC
         do i = 1, nfrec
          write(3211,'(I0,2x,ES14.5E2,2x,ES14.5E2,2x,ES14.5E2,2x,ES14.5E2)') & 
          i,real(developerAUXvec(i,2),4),real(aimag(developerAUXvec(i,2)),4),real(S(i)),aimag(S(i))
         end do
         close (3211)
       
        return
       
!       S(1:nfrec) = S(1:nfrec) / developerAUXvec(1:nfrec,1)
!       write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
!              'f_',nombre,iP,'[', &
!              int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
!              int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
!        logflag = 'none     '
!       call plotSpectrum(S(1:100),real(developerAUXvec(2,2),4),& 
!            100,100,titleN,xAx,yAx,logflag,1200,800,real(developerAUXvec(2,2)*100,4))
!       !guardar en texto
!        write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
!              'f_',nombre,iP,'[', &
!              int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
!              int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
!        
!        OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
!        write(3211,'(I0,A)') 1," espectro antes de convolución con funcion de amplitud"
!        write(3211,'(a)') "dfrec="
!        write(3211,'(F15.8)') DFREC
!        do i = 1, nfrec
!         write(3211,'(ES14.5E2,2x,ES14.5E2,2x,ES14.5E2,2x,ES14.5E2)') & 
!         real(developerAUXvec(i,2),4),real(aimag(developerAUXvec(i,2)),4),real(S(i)),aimag(S(i))
!        end do
!        close (3211) 
!         return
       end if!
       if (nombre(1:1) .eq. "r") then
          write(xAx,'(a)') "$a \omega c_p^{-1}$"
          S = z0
          S(1:nfrec)= developerAUXvec(1:nfrec,1)
          write(titleN,'(a)') 'denominador.pdf'
          logflag = 'none     '
          call plotSpectrum(S(1:100),real(developerAUXvec(2,2),4),& 
             100,100,titleN,xAx,yAx,logflag,1200,800,real(developerAUXvec(2,2)*100,4))
          !guardar en texto
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
         
         OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         write(3211,'(I0,A)') 1," denominador"
         write(3211,'(a)') "dfrec="
         write(3211,'(F15.8)') DFREC
         do i = 1, nfrec
          write(3211,'(ES14.5E2,2x,ES14.5E2,2x,ES14.5E2,2x,ES14.5E2)') & 
          real(developerAUXvec(i,2),4),real(aimag(developerAUXvec(i,2)),4),real(S(i)),aimag(S(i))
         end do
         close (3211) 
          return
        end if
      end if
      S = z0
      S(1:nfrec+1)= W(1:nfrec+1:+1)
      S(NPTSTIME-NFREC+2:NPTSTIME) = conjg(W(nfrec:2:-1)) 
      S = S * t_vec 
      
      ! conv con fucion de amplitud
      S = S * Uo
!     CALL SETFIL("SUr.pdf")
!     call qplot(real((/((i-1)*dfrec,i=1,nptstime)/),4),real(S,4),nptstime)
!     CALL SETFIL("SUi.pdf")
!     call qplot(real((/((i-1)*dfrec,i=1,nptstime)/),4),real(aimag(S),4),nptstime)
!     stop
      if (Verbose .ge. 4) call showMNmatrixZ(nptstime,1, S,"  S  ",6)
!     print*,"6751"
      write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x_i,' , ',z_i,')'
      if (verbose .ge. 1) then
       if ((SabanaPlotIndividual .eqv. .false.) .and. & 
           (allpoints(iP)%isSabana .eqv. .true.)) go to 527
      ! grafica simple:
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'f_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
        ! plotXYcomp(y_in,Dt,n,titleN,xAx,yAx,CTIT,W,H,ma)
         call plotXYcomp(S(1: NFREC),real(DFREC,4), NFREC,titleN, & 
         'frec[hz] ',yAx, CTIT ,1200,800,0.0)
      end if!
      if (Verbose .ge. 2) then
      ! grafica logaritmica
            write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'fL_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
            logflag = 'logx     '
            write(xAx,'(a)') "frec[Hz]"
            call plotSpectrum(S(1: NPTSTIME/2),real(DFREC,4), NPTSTIME/2, NPTSTIME/2, & 
                titleN,xAx,yAx,logflag,1200,800,real(DFREC*(NFREC+1),4))
      end if!
      
      !  (1) pasar al tiempo
!527     S = S * (sqrt(1.0*NPTSTIME) * dfrec) !backward
!        call fork(NPTSTIME,S,+1,verbose,outpf)
 527     S = FFTW(NPTSTIME,S,+1,1/(NPTSTIME*dt)) !backward
         
      !  (2.1) remover efecto de la frecuencia imaginaria
         S = S * exp(- OMEI * Dt*((/(i,i=0, NPTSTIME-1)/)))
      !  (2.2) remover efecto de la velocidad imaginaria  ?
!        S = S * exp((2./2./Qq) * Dt*((/(i,i=0, NPTSTIME-1)/)))
         
         !tiempo maximo para graficar
         n_maxtime = int(maxtime/dt)
         if(maxtime .lt. dt) n_maxtime = 2*nfrec
         if(maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
         
      if (verbose .ge. 4) print*,"maxtime = ", this_maxtime," segs :: @",dt, & 
                                  " : ",n_maxtime," puntos"
      !  (3) plot the damn thing
      if (verbose .ge. 3) then
      !guardar en texto
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               'S_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].txt'
!        print*,titleN
         OPEN(3211,FILE=titleN,FORM="FORMATTED",ACTION='WRITE')
         write(3211,'(F15.8)') DT
         do i = 1,size(S)
          write(3211,'(ES14.5E2,2x,ES14.5E2)') real(S(i)),aimag(S(i))
         end do
         close (3211) 
      end if
      
      ! guardar para hacer sabana o plotear
      if (allpoints(iP)%isSabana) then
         ! guardamos la sabana actual
!        print*,"saved sabana point, ",iP-(nIpts - nSabanapts)
         Sabana(iP-(nIpts - nSabanapts),1:NPTSTIME) = S
         if (SabanaPlotIndividual .eqv. .false.) return
      end if
         
!        if (plotmaxy .lt. maxval(real(S))) plotmaxy = maxval(real(S,4))
         write(titleN,'(a,a,I0,a,I0,a,I0,a,I0,a,I0,a)') & 
               '2_S_',nombre,iP,'[', &
               int(x_i),'.',abs(int((x_i-int(x_i))*10)),';', & 
               int(z_i),'.',abs(int((z_i-int(z_i))*10)),'].pdf'
         call plotXYcomp(S(1:n_maxtime) ,real(dt,4),n_maxtime,titleN, & 
         'time[sec]',yAx, CTIT ,1200,800, 0.)
      
!     end do !imec
      end subroutine W_to_t
      
      subroutine makeSabana(nombre,filled)
      use waveNumVars, only : NPTSTIME,nfrec
      use glovars
      use waveVars, only : dt,maxtime
      use resultvars, only : Sabana, nSabanapts,sabZeroini,sabZerofin
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
      n_maxtime = int(maxtime/dt)
      if(maxtime .lt. dt) n_maxtime = 2*nfrec
      if(maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
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
      
      subroutine crepa_four_fotogramas
      use glovars, only : saveG,verbose
      use peli, only : fotogramas 
      use sourceVars, only : SH,PSV
      use waveNumVars, only : t_vec,omei, nF => NFREC, nT => NPTSTIME!,dfrec
      use meshVars, only : npixX,npixZ
      use waveVars, only : Uo,Dt
!     use soilVars, only : Qq
      use wavelets ! FORK
      use dislin
      implicit none
      integer :: i,ix,iz,imec,mecS,mecE
      complex*16, dimension(:), pointer :: p_fot
      character(LEN=100) :: titleN
!     real*8 :: factor
      if (verbose .ge. 1) print*,"frec -> time"
!     factor = sqrt(1.0*nT)
      mecS = 3; mecE = 2; if (PSV) mecS =1; if (SH) mecE =3
      if (saveG .eqv. .true.) then
       CALL chdir("video")
       OPEN(6374,FILE="G.bin",STATUS='UNKNOWN', ACCESS='STREAM',ACTION='WRITE')
       write(6374) mecS
       write(6374) mecE
       write(6374) npixZ
       write(6374) npixX
       write(6374) nT
      end if
      do imec = mecS,mecE
      do iz=1,npixZ
        do ix=1,npixX
         !fotogramas(iz,ix,1:Nfrec+1,imec)  ok
          fotogramas(iz,ix,nT-nF+2:nT,imec) = conjg(fotogramas(iz,ix,nF:2:-1,imec))
          p_fot => fotogramas(iz,ix,1:nT,imec)
          p_fot = p_fot * t_vec
          if (saveG .eqv. .true.) then 
           do i=1,nT
            write(6374) p_fot(i)
           end do
          end if
          p_fot = p_fot * Uo
        ! al tiempo:
          !p_fot = p_fot * (sqrt(1.0*nT) * dfrec) !backward
          !call fork(nT,p_fot,+1,1,6)
          p_fot = FFTW(nT,p_fot,+1,1/(nT*dt))
        ! remover efecto de la frecuencia imaginaria
          p_fot = p_fot * & 
          exp(-OMEI * Dt*((/(i,i=0,nT-1)/)))
        ! remover efecto de la velocidad imaginaria
          !p_fot = p_fot * & 
          !exp((1./2./Qq) * Dt*((/(i,i=0, nT-1)/))) 
          if (verbose .ge. 2) then
            write(titleN,"(i0,a,i0,a)") ix,"_",iz,".pdf"
            CALL SETFIL(trim(titleN))
            call qplot(real((/((i-1)*Dt,i=1,800)/),4),real(p_fot(1:800),4),800)         
          end if
        end do !ix
      end do !iz
      end do !imec
      if (saveG .eqv. .true.) THEN 
      close(6374)
      print*,"saved G function"
      CALL chdir("..")
      end if
      end subroutine crepa_four_fotogramas
      
      subroutine loadG_fotogramas
      use glovars, only:verbose,rutaOut
      use peli, only : fotogramas 
      use waveNumVars, only : omei, nT => NPTSTIME!,dfrec
      use meshVars, only : npixX,npixZ
      use waveVars, only : Uo,Dt
!     use soilVars, only : Qq
      use wavelets ! FORK
      use dislin
      implicit none
      integer :: i,ix,iz,imec,mecS,mecE
      complex*16, dimension(:), pointer :: p_fot
!     real*8 :: factor
      character(LEN=100) :: titleN
      call chdir(trim(adjustl(rutaOut)))
      CALL chdir("video")
       OPEN(6375,FILE="G.bin",STATUS='UNKNOWN', ACCESS='STREAM',ACTION='READ')
       read(6375) mecS
       read(6375) mecE
       read(6375) npixZ
       read(6375) npixX
       read(6375) nT
!     factor = sqrt(1.0*nT)
      do imec = mecS,mecE
      do iz=1,npixZ
      do ix=1,npixX
      do i=1,nT
         read(6375) fotogramas(iz,ix,i,imec)
      end do
          p_fot => fotogramas(iz,ix,1:nT,imec)
          p_fot = p_fot * Uo
        ! al tiempo:
          !p_fot = p_fot * (sqrt(1.0*nT) * dfrec) !backward
          !call fork(nT,p_fot,+1,1,6)
          p_fot = FFTW(nT,p_fot,+1,1/(nT*dt))
        ! remover efecto de la frecuencia imaginaria
          p_fot = p_fot * & 
          exp(-OMEI * Dt*((/(i,i=0,nT-1)/)))
        ! remover efecto de la velocidad imaginaria
          !p_fot = p_fot * & 
          !exp((1./2./Qq) * Dt*((/(i,i=0, nT-1)/))) 
          if (verbose .ge. 2)then
            write(titleN,"(i0,a,i0,a)") ix,"_",iz,".pdf"
            CALL SETFIL(trim(titleN))
            call qplot(real((/((i-1)*Dt,i=1,800)/),4),real(p_fot(1:800),4),800) 
          end if
      end do !ix
      end do !iz
      end do !imec
      close(6375)
      CALL chdir("..")
      end subroutine loadG_fotogramas
      

      subroutine Churubusco
      use glovars, only : verbose, workBoundary,flip12
      use DISLIN
      use peli, only : ypray => coords_Z, xpray => coords_X,& 
                     fotogramas
      use meshVars, only : npixX,npixZ,MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ
      use waveVars, only : dt,maxtime
      use waveNumVars, only : NFREC, NPTSTIME
      use soilVars, only : Z,N,col=>layershadecolor, shadecolor_inc
      use geometryvars, only : nXI,Xcoord_ER, & 
                               n_cont, Xcoord_Voidonly, Xcoord_Incluonly,Xcoord_flip_out
      use resultvars, only : Punto,BouPoints,nbpts
      use ploteo10pesos
      implicit none
      interface
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptores
      end subroutine drawBoundary
      end interface
      real, dimension(:,:,:), allocatable :: xvmat,yvmat
      real :: maV1,maV2,minX,maxX,minY,maxY,xstep,zstep,encuadre, tlabel, madmax, escalaFlechas
      integer :: i,ii,j,n_maxtime,iT
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      integer*4 :: lentitle
      character(LEN=60) :: CTIT
      type (Punto), dimension(:), pointer :: BP
      real*8, dimension(:,:),allocatable :: rec
      ! fotogramas tipo campo vectorial
      if (verbose >= 1) print*, "Will make a movie..."
      CALL chdir("video")
      
      !tiempo maximo para graficar
       n_maxtime = int(maxtime/dt)
       if(maxtime .lt. dt) n_maxtime = 2*nfrec
       if(maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
       print*,"maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," puntos"
       
      allocate(xvmat(npixX,npixZ,n_maxtime))
      allocate(yvmat(npixX,npixZ,n_maxtime))
      
      maV1 = maxVal(real(fotogramas(:,:,1:n_maxtime,1),4))
      maV2 = maxVal(real(fotogramas(:,:,1:n_maxtime,2),4))
      do i=1,npixZ
        do j=1,npixX
          xvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime,2),4) !U
          yvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime,1),4) !W
          
          do iT = 1,n_maxtime
            if (abs(yvmat(j,i,iT)) .le. mav1*0.0025) yvmat(j,i,iT) = 0.0
            if (abs(xvmat(j,i,iT)) .le. mav2*0.0025) xvmat(j,i,iT) = 0.0
          end do
        end do
      end do
      
      minx = MeshMinX
      maxx = MeshMaxX
      miny = MeshMinZ
      maxy = MeshMaxZ
      xstep = real(((maxX-minX) / 5.0 ))
      zstep = real(((maxY-minY) / 6. ))
      encuadre = (maxY-minY)/(maxX-minX)
!     xstep = real(abs(xpray(npixX)-xpray(1))/4.0,4)
!     zstep = real(abs(ypray(npixZ)-ypray(1))/10.0,4)
      madmax = max(max(maxval(xvmat),maxval(yvmat)),max(maxval(abs(xvmat)),maxval(abs(yvmat))))
      escalaFlechas = real((xstep * 0.6) / madmax)
!     encuadre = (ypray(npixZ)-ypray(1))/(xpray(npixX)-xpray(1))
      print*,"encuadre=",encuadre
      CALL METAFL('PNG')
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(3000,4),int(3000,4))
!     CALL PAGMOD('NONE')
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(1200,4)) !1200 ambos
      CALL SCRMOD('REVERS') !fondo blanco
      do i=1,n_maxtime
      write(titleN,'(a,I0,a)') 'foto_',i,'.png'
      CALL SETFIL(trim(titleN))
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL DISALF() !default font
           !the position of an axis system.
      CALL axspos (int(300,4) ,int(2700,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600*encuadre,4)) !size of the axis system.
      call name('X [m]','X')
      call name('Z [m]','Y')
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
      if (workboundary) then !
      ! dibujar inclusión
      if (flip12) then ! Xcoord_flip_out
      if (allocated(Xcoord_flip_out)) then
      if (size(Xcoord_flip_out(:,1,1)) .gt. 1) then
!     do j = 1,8
!        print*,j,Xcoord_flip_out(j,1,1:2),Xcoord_flip_out(j,2,1:2)
!     end do
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
      else ! --------
      if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
      if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
      allocate(rec(2* (size(Xcoord_Incluonly(:,1,1))),2))              !
      ii=1                                                             !
      do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont     !
      rec(ii,1) = Xcoord_Incluonly(j,1,1)                              !
      rec(ii,2) = Xcoord_Incluonly(j,2,1)                              !
      rec(ii+1,1) = Xcoord_Incluonly(j,1,2)                            !
      rec(ii+1,2) = Xcoord_Incluonly(j,2,2)                            !
      ii=ii+2                                                          !
      end do
      
      call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
      CALL RLAREA(real(rec(:,1),4),real(rec(:,2),4),int(2*n_cont,4))  !
      deallocate(rec) 
      end if!
      if (size(Xcoord_Voidonly(:,1,1)) .gt. 1) then
      ! Borrar lo que está en el aire !
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
      CALL RLAREA(real(rec(:,1),4),real(rec(:,2),4), &                !
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))              !
      deallocate(rec) 
      end if                                                !         !
      end if
      end if
      ! dibujar contorno de topografia original                       !
      call color ('FORE')                                             !
      call PENWID(real(0.5,4))                                        !
      call marker(int(-1,4)) ! sin marcadores                         !
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do                                                          !
      end if
      
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
      
      call disfin
      
!     stop "6708 killed video"
      
      end do ! i=1,n_maxtime
      write(titleN,'(a)') 'foto_0.png'
      write(extension,'(a)') 'PNG'
      BP => BouPoints
      call drawBoundary(BP,nbpts,titleN, extension,.false.,.false.)
      
      write(titleN,'(a)')'ffmpeg -i foto_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4'
      call system(trim(titleN))
      call chdir("..")
      call system('cp video/video.mp4 video.mp4')
      call system('rm video/video.mp4')
      
      if (encuadre - 0.5 .le. 0.1) then
      write(titleN,'(a,I0,a,I0,a)') 'ffmpeg -i video.mp4 -filter:v ''''crop=1200:',&
                            int(encuadre*1200+150),':0:',1200-int(encuadre*1200+150),''''' video_Crop.mp4'
      call system(trim(titleN))
      end if
      end subroutine Churubusco    
      
      subroutine Hollywood(imec)
      use DISLIN
      use glovars, only : workboundary,verbose
      use peli, only : ypray => coords_Z, xpray => coords_X,& 
                     fotogramas
      use meshVars, only : npixX,npixZ
      use waveNumVars, only : NFREC, NPTSTIME
      use soilVars, only : Z,N, shadecolor_inc
      use geometryvars, only : nXI,Xcoord_ER, Xcoord_Voidonly, Xcoord_Incluonly
      use waveVars, only : dt,maxtime
      use resultvars, only : Punto,BouPoints,nbpts
      use ploteo10pesos !las rutinas externas no requieren interfaz 
      implicit none
      interface
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
      use resultVars, only : Punto
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptores
      end subroutine drawBoundary
      end interface
      integer ,intent(in) :: imec
      real :: Sm(npixX,npixZ)
      character(LEN=3), dimension(3) :: nombre
      real     :: factor, ColorRangeMaximumScale, tlabel
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
      
      factor = sqrt(real(NPTSTIME))
      nombre(1)= 'w--'
      nombre(2)= 'u--'
      nombre(3)= 'v--'
      CALL chdir("video")
         
      !tiempo maximo para graficar
         n_maxtime = int(maxtime/dt)
         if(maxtime .lt. dt) n_maxtime = 2*nfrec
         if(maxtime .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
         print*,"maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," puntos"
        
       ColorRangeMaximumScale = 0.1
  123  maV = maxVal(real(fotogramas(:,:,1:n_maxtime,imec),4))
       miV = minVal(real(fotogramas(:,:,1:n_maxtime,imec),4))
       maV = max(maV,abs(miV))
       miV = - maV
       
         print *, char(7)
         write(6,'(a,a,a,EN13.2,a,/,a,/,a,EN13.2,/,a)', ADVANCE = "NO") & 
         'Look at the seismograms for ', nombre(imec), & 
         '. Is the response too spiky (the max = ',maV,'? ', &
         'We can enhance detail by reducing the value for maximum color. ', &
         'We propose the maximum to be = ', & 
         ColorRangeMaximumScale * maV, &
         'Proceed [Y] or change it [else] ?'
       
       IF (verbose .eq. 1)then
         print*," YES (automatic)"
         imdone = "Y"
       ELSE
         read(5,*)imdone
       end if
       
         if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
            write(6,'(a)') "keeping those plot limits"
         else
            write(6,'(a)', ADVANCE = "NO")'New maximum for plot = '
            read(5,*) p
            ColorRangeMaximumScale = real(p / maV,4)
         end if
       colorBounds(1) = maV * ColorRangeMaximumScale
       colorBounds(2) = miV * ColorRangeMaximumScale
       
          write(6,'(a,a,a,a,E12.4,a,E12.4,a)') "colorbounds:", & 
          "(",nombre(imec),") [",colorBounds(2)," ; ",colorBounds(1),"]"
       
      minx = minval(xpray)
      maxx = maxval(xpray)
      miny = minval(ypray)
      maxy = maxval(ypray)
      xstep = real(abs(xpray(npixX)-xpray(1))/3.0,4)
      zstep = real(abs(ypray(npixZ)-ypray(1))/10.0,4)
      encuadre = (ypray(npixZ)-ypray(1))/(xpray(npixX)-xpray(1))
      print*,"encuadre=",encuadre      
        
      Iframe = 1
      Fframe = n_maxtime 
      
      if(verbose .eq. 1) then
        imdone = "Y"
      else
        print *, char(7)
        write(6,'(a,I0,a)')'Look at the seismograms, I will plot ',Fframe,&
        'frames. Proceed [Y] or change frame range [else]'
        read(5,*)imdone
      end if!
      
      if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
        Write(6,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
      else
        write(6,'(a)', ADVANCE = "NO")'Start frame = '
        read(5,*)Iframe
        write(6,'(a)', ADVANCE = "NO")'Final frame = '
        read(5,*)Fframe
        Write(6,'(a,I5,a)') ' plotting',Fframe-Iframe+1,' photograms'
      end if
      
      do iT=Iframe,Fframe !cada fotograma
      do i=1,npixX
      do j=1,npixZ
        Sm(i,j) = real(fotogramas(j,i,iT,imec),4)
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
      if (workboundary) then ! si es una topografía                                         !
      if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
      allocate(rec(2* (size(Xcoord_Incluonly(:,1,1))),2))             !
      i=1                                                             !
      do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont    !
      rec(i,1) = Xcoord_Incluonly(j,1,1)                              !
      rec(i,2) = Xcoord_Incluonly(j,2,1)                              !
      rec(i+1,1) = Xcoord_Incluonly(j,1,2)                            !
      rec(i+1,2) = Xcoord_Incluonly(j,2,2)                            !
      i=i+2                                                           !
      end do                                                          !
!     print*, "shadecolor_inc", shadecolor_inc
      call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
      CALL RLAREA(real(rec(:,1),4), & 
                  real(rec(:,2),4), & 
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))  !
      deallocate(rec) 
      end if!
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
      call drawBoundary(BP,nbpts, path, extension,.false.,.false.)
      
      write(path,'(a,a,a)')'ffmpeg -i ',nombre(iMec), & 
                  '%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4'
      call system(trim(path))
      call chdir("..")
      if (encuadre - 0.5 .le. 0.1) then
      write(path,'(a,a,a)') & 
      'ffmpeg -i video/video.mp4 -filter:v ''''crop=1200:700:0:600'''' ',& 
      nombre(iMec),'.mp4'
      call system(trim(path))
      else
      write(path,'(a,a,a)') 'cp video/video.mp4 ',nombre(iMec),'.mp4'
      call system(trim(path))
      end if
      call system('rm video/video.mp4')
      
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
        call chdir("video")
        call system('rm *.png')
        call chdir("..")
      end if!
      
        write(6,'(a)') 'Keep *.txt files? [Y]'
        read(5,*)imdone
        if(imdone .eq. 'Y' .or. imdone .eq. 'y') then
          write(6,'(a)') 'kept the files'
        else
          call chdir("video")
          call system('rm *.txt')
          call chdir("..")
        end if
      end if !verbose 2
      end subroutine Hollywood           
      
      subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
      
      ! my interface
!     interface
!     subroutine drawBoundary(BP, nbpts, titleN, extension, zoomGeom, plotReceptores)
!     use resultVars, only : Punto
!     type (Punto), dimension(:), pointer :: BP
!     integer, intent(in) :: nbpts
!     character(LEN=100) :: titleN
!     character(LEN=3) :: extension
!     logical, intent(in) :: zoomGeom, plotReceptores
!     end subroutine drawBoundary
!     end interface
      
      
      use DISLIN
      use soilVars, only : Z,N,col=>layershadecolor, shadecolor_inc
      use resultVars, only : Punto, IP => allpoints, nIpts ,overDeterminedSystem
      use sourceVars, only : Po,NFuentes, tipoFuente, PW_pol
      use geometryvars, only : nXI,Xcoord_ER,normXI, & 
                               midPoint, Xcoord_Voidonly, Xcoord_Incluonly,Xcoord_flip_out
      use glovars, only : verbose, workBoundary,flip12
      use meshVars, only : MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ
      
      implicit none
      type (Punto), dimension(:), pointer :: BP
      integer, intent(in) :: nbpts
      character(LEN=100) :: titleN
      character(LEN=3) :: extension
      logical, intent(in) :: zoomGeom, plotReceptores
      real :: maxY,minY,maxX,minX,xstep,zstep,encuadre
      integer :: i,J,ifuente,sc
      real*8, dimension(:,:),allocatable :: rec
      real*8,dimension(2) :: theta,gamma
      if (verbose >= 4) Write(6,'(a)') "Will plot geometry..."
      sc=3
      ! Boundary boundaries
      maxX = MeshMaxX
      minX = MeshMinX
      maxY = MeshMaxZ
      minY = MeshMinZ
      
      
      xstep = real(((maxX-minX) / 5.0 ))
      zstep = real(((maxY-minY) / 10. ))
      encuadre = (maxY-minY)/(maxX-minX)

      if (zoomGeom) then
      maxX = max(MeshMaxX,maxval(real(Xcoord_ER(:,1,:),4)))
      minX = min(MeshMinX,minval(real(Xcoord_ER(:,1,:),4)))
      maxY = max(MeshMaxZ,maxval(real(Xcoord_ER(:,2,:),4)))
      minY = min(MeshMinZ,minval(real(Xcoord_ER(:,2,:),4)))
      end if
      
      ! Dislin plotting routines 
      CALL METAFL(extension) ! define intended display  XWIN,PS,EPS,PDF
      call imgfmt('RGB')
      call winsiz(int(1200,4),int(1200,4)) !1200 ambos
      CALL SCRMOD('REVERS') !fondo blanco
      
      CALL SETFIL(trim(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
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
      call labdig(int(1,4),'Y')
      call ticks (int(10,4) ,'XY') 
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
      call PENWID(real(1.0,4))                                        !
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
      call color('ORANGE')                                                  !
      CALL HSYMBL(int(25,4)) !size of symbols                               ! 
      do j=1,nXI                                                            !
      CALL RLVEC (real(midPoint(j,1),4), real(midPoint(j,2),4), &           !
              real(midPoint(j,1)+normXI(j,1)*xstep*0.1,4), &                !
              real(midPoint(j,2)+normXI(j,2)*xstep*0.1,4), int(1001,4))     !
      end do                                                                !
      
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
      
      if (.not. zoomGeom) then
      !fuente ------------------------------------------------------------
      call color('RED')                                                  !
      do ifuente = 1,Nfuentes                          
      if (tipoFuente .eq. 0) then !puntual
      CALL HSYMBL(int(40,4)) !size of symbols                            !
      CALL RLSYMB (8, real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4))!star     !
      CALL RLVEC (real(Po(ifuente)%center%x,4), real(Po(ifuente)%center%z,4), &            !
              real(Po(ifuente)%center%x + Po(ifuente)%normal%x * xstep*0.3,4), &           !
              real(Po(ifuente)%center%z + Po(ifuente)%normal%z * xstep*0.3,4), int(1101,4))!
      elseif (tipoFuente .eq. 1) then !onda plana
      ! polarización
      if (PW_pol .eq. 1) then !SV
        theta(1) = cos(Po(ifuente)%gamma)
        theta(2) = sin(Po(ifuente)%gamma)
      elseif (PW_pol .eq. 2) then !P
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
      elseif (tipoFuente .eq. 2) then !segmento
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
      end do
      end if!
      if (plotReceptores) then
      !receptores ----------------------------------------------------------
      CALL HSYMBL(int(25,4)) !size of symbols                              !
      call color('BLUE')                                                   !
      do j=1,nipts                                                         !
        if (IP(j)%isOD) then                                               !
          if (overDeterminedSystem) then                                   !
            CALL RLSYMB (0, real(IP(j)%center%x,4), real(IP(j)%center%z,4))!
            CALL RLVEC (real(IP(j)%center%x,4), real(IP(j)%center%z,4), &  !
              real(IP(j)%center%x + IP(j)%normal%x * xstep*0.1,4), &       !
              real(IP(j)%center%z + IP(j)%normal%z * xstep*0.1,4), &       !
              int(1001,4))                                                 !
          end if                                                           !
        else                                                               !
          CALL RLSYMB (2, real(IP(j)%center%x,4), real(IP(j)%center%z,4))  !
        end if                                                             !
      end do                                                               !
      end if
      call disfin()
      end subroutine drawBoundary
      

      subroutine plot_at_eta (J, tt)
      use resultVars, only : Punto, allpoints, nIpts, nSabanapts
      use dislin
      use debugStuff
      implicit none
      integer :: iP,J
      character(LEN=100) :: tt,t2
      complex*16 :: BP
      real, dimension(nSabanapts) :: x
      complex*16, dimension(nSabanapts) :: y
      real*8 :: chNa
      logical :: warning
      warning = .false.
      print*,"plot ", trim(tt)
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4L')
      CALL SETFIL(trim(tt))
        do iP = 1,nIpts; if (allpoints(iP)%isSabana) then
           if (tt(1:1) .eq. 'W') BP = allpoints(iP)%resp(J)%W
           if (tt(1:1) .eq. 'U') BP = allpoints(iP)%resp(J)%U
           chNa = abs(BP)
           if (chNa .eq. chNa+1) then; BP = 0; warning = .true.; end if
!          print*,iP,BP
           y(iP-(nIpts - nSabanapts)) = BP
           x(iP-(nIpts - nSabanapts)) = real(allpoints(ip)%center%x,4)
        end if;end do
        write(t2,'(a,a)') trim(tt),"_.txt"
        OPEN(739,FILE=trim(t2),FORM="FORMATTED",ACTION='WRITE')
        call showMNmatrixZ(nSabanapts,1,y(1:nSabanapts),t2(1:5),739)
        close(739)
        
        if (warning) then 
        print*,trim(tt), ": One or more values are either NaN or inf"
        else
        
!       call showMNmatrixZabs(nSabanapts,1, Sabana(1:nSabanapts,1),tt(1:5),6)
        write(t2,'(a,a)') "0__Abs_",trim(tt)
        CALL SETFIL(trim(t2))
        CALL color('BLACK')
        call qplot(x(1:nSabanapts),real(abs(y(1:nSabanapts)),4), nSabanapts)
        
        write(t2,'(a,a)') "0__Rea_",trim(tt)
        CALL SETFIL(trim(t2))
        CALL color('RED')
        call qplot(x(1:nSabanapts),real(real(y(1:nSabanapts)),4), nSabanapts)
        
        write(t2,'(a,a)') "0__Ima_",trim(tt)
        CALL SETFIL(trim(t2))
        CALL color('BLUE')
        call qplot(x(1:nSabanapts),real(aimag(y(1:nSabanapts)),4), nSabanapts)
        end if
                
!       CALL color('BLACK')
!       call qplcrv(x(1:nSabanapts),real(abs(Sabana(1:nSabanapts,1)),4), nSabanapts,'FIRST')
!       call color('RED')
!       call qplcrv(x(1:nSabanapts),real(real(Sabana(1:nSabanapts,1)),4), nSabanapts,'NEXT')
!       call color('BLUE')
!       call qplcrv(x(1:nSabanapts),real(aimag(Sabana(1:nSabanapts,1)),4), nSabanapts,'LAST')
        
      
      end subroutine plot_at_eta
