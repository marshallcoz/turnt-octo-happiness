      module cine
      contains
      subroutine CINETECA
      use DISLIN
      use glovars, only : makevideo
      use peli, only : ypray => coords_Z, xpray => coords_X,fotogramas
      use meshVars, only : npixX,npixZ, MeshVecLen,MeshVecLen2, &
                           MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ
      use waveVars, only : dt,maxtime
      use waveNumVars, only : NFREC, NPTSTIME
      use soilVars, only : Z,N,col=>layershadecolor, shadecolor_inc
      use geometryvars, only : nXI,Xcoord_ER, &
                               Xcoord_Voidonly, Xcoord_Incluonly,&
                               n_topo,n_cont,n_vall
      use resultvars, only : Punto,allpoints,nIpts,nPtsolos,nSabanapts,&
                              nSecciones,nBPt_topo,nBPt_cont,nBPt_vall
      use ploteo10pesos
      use sourceVars, only : iFte => currentiFte
      implicit none
      real, dimension(:,:,:), allocatable :: xvmat,yvmat
      real :: maV1,maV2,minX,maxX,minY,maxY,xstep,zstep, &
      tlabel, escalaFlechas,escalaFlechas2,centrox,centroz,radio, madmax
      real,dimension(nIpts*2,NPTSTIME,2) :: delX,delZ
      integer :: i,ii,iii,j,j2,n_maxtime,iT,k,fai,faf
      character(LEN=100) :: titleN
      integer*4 :: lentitle
      character(LEN=60) :: CTIT
      real*8, dimension(:,:),allocatable :: recIncluonly,recVoidonly
      real, dimension(5,2) :: rec
      real*8 :: rat
      real*8 :: getstt,maxstt,getsrt,maxsrt
      real*8, dimension(:,:), allocatable :: stt,srt
      integer :: nframes
      real*8 :: si,co,szz,szx,sxx

      !tiempo maximo para graficar
       n_maxtime = int(maxtime(iFte)/dt)
       if(maxtime(iFte) .lt. dt) n_maxtime = 2*nfrec
       if(maxtime(iFte) .gt. NPTSTIME * real(dt,4)) n_maxtime = NPTSTIME
       print*,"maxtime = ",maxtime," segs :: @",dt," : ",n_maxtime," fotogramas"
       allocate(xvmat(npixX,npixZ,n_maxtime))
       allocate(yvmat(npixX,npixZ,n_maxtime))

      nframes = n_maxtime
      if (makevideo) then
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
      madmax = max(mav1,mav2)
      escalaFlechas = real(MeshVecLen / madmax)
      print*,"Escalaflechas = ",escalaFlechas

      end if
      fai = nIpts-nXi-nSabanapts-nSecciones
      faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
      print*,"nIpts = ",nIpts
      print*,"nXi =   ",nXi
      print*,"nPtsolos = ",nPtsolos
      print*,"nSabana=",nSabanapts
      print*,"nSeccio=",nSecciones
      print*,"n_topo =",n_topo
      print*,"n_cont =",n_cont
      print*,"n_vall =",n_vall
      mav1 = -100000.0
      mav2 = -100000.0
      do j=fai,faf
        if (allpoints(j)%atBou) then
          mav1 = max(mav1,maxval(real(allpoints(j)%S(1:n_maxtime,1),4)))!y
          mav2 = max(mav2,maxval(real(allpoints(j)%S(1:n_maxtime,2),4)))!x
        end if
      end do
      madmax = max(mav1,mav2)
      !escalaFlechas = real(MeshVecLen / madmax)
      escalaFlechas2 = real(MeshVecLen2 / madmax)
      print*,"Escalaflechas2 = ",escalaFlechas2
      print*,"madmax=",madmax

      maxX = maxval(real(Xcoord_ER(1:nXI,1,:),4))
      minX = minval(real(Xcoord_ER(1:nXI,1,:),4))
      maxY = maxval(real(Xcoord_ER(1:nXI,2,:),4))
      minY = minval(real(Xcoord_ER(1:nXI,2,:),4))
      centroX = (maxX+minX)/2
      centroZ = (maxY+minY)/2
      radio = max((maxY-minY)/2,(maxX-minX)/2)*1.5
      maxx = centroX + radio
      minx = centroX - radio
      maxy = centroz + radio
      miny = centroz - radio
      print*,minX,maxX,minY,maxY
      xstep = real(((maxX-minX) / 0 ))
      zstep = real(((maxY-minY) / 0 ))

      if (allocated(Xcoord_Incluonly)) then
       if (size(Xcoord_Incluonly(:,1,1)) .gt. 1) then
        if (allocated(recIncluonly)) deallocate(recIncluonly)
        allocate(recIncluonly(2* (size(Xcoord_Incluonly(:,1,1))),2))
        ii=1                                                             !
        do j=1,size(Xcoord_Incluonly (:,1,1))!n_topo+1,n_topo+n_cont     !
        recIncluonly(ii,1) = Xcoord_Incluonly(j,1,1)                     !
        recIncluonly(ii,2) = Xcoord_Incluonly(j,2,1)                     !
        recIncluonly(ii+1,1) = Xcoord_Incluonly(j,1,2)                 !
        recIncluonly(ii+1,2) = Xcoord_Incluonly(j,2,2)               !
        ii=ii+2                                                          !
        end do
!       if (allocated(recIncluonly)) then
!       call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
!       CALL RLAREA(real(recIncluonly(:,1),4), &
!                 real(recIncluonly(:,2),4), &
!                 int(2*size(Xcoord_Incluonly(:,1,1)),4))
!       end if
       end if!
       if (size(Xcoord_Voidonly(:,1,1)) .gt. 1) then
        if (allocated(recVoidonly)) deallocate(recVoidonly)
        allocate(recVoidonly(2* (size(Xcoord_Voidonly(:,1,1))),2))              !
        ii = 1                                                          !
        do j=1,size(Xcoord_Voidonly(:,1,1))                             !
        recVoidonly(ii,1) = Xcoord_Voidonly(j,1,1)                              !
        recVoidonly(ii,2) = Xcoord_Voidonly(j,2,1)                              !
        recVoidonly(ii+1,1) = Xcoord_Voidonly(j,1,2)                            !
        recVoidonly(ii+1,2) = Xcoord_Voidonly(j,2,2)                            !
        ii=ii+2                                                         !
        end do                                                          !
!       if (allocated(recVoidonly)) then
!       call color ('BACK')                                             !
!       call shdpat(int(16,4))                                          !
!       CALL RLAREA(real(recVoidonly(:,1),4), &
!                 real(recVoidonly(:,2),4), &                !
!                 int(2*size(Xcoord_Voidonly(:,1,1)),4))              !
!       end if
       end if
      end if

      allocate(stt(nframes,nIpts))
      allocate(srt(nframes,nIpts))
      maxstt = -1000000;
      maxsrt = -1000000;
      do i=1,nframes
      do j=nIpts-nXi-nSabanapts-nSecciones,nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
      if (allpoints(j)%atBou) then

      si = allpoints(j)%sinT
      co = allpoints(j)%cosT
      !W U sxx szx szz
      sxx = real(allpoints(j)%S(i,3),8)
      szx = real(allpoints(j)%S(i,4),8)
      szz = real(allpoints(j)%S(i,5),8)
      getstt = si**2*sxx + co**2*szz - 2*si*co*szx
      stt(i,j) = getstt
      maxstt = max(maxstt,stt(i,j))
      getsrt = si*co*(szz - sxx) + szx * (co**2 - si**2)
      srt(i,j) = getsrt
      maxsrt = max(maxsrt,srt(i,j))
      end if
      end do! j
      end do! i

      ! ******************** envolvente


      ! ******************** pelicula
!     maxstt = maxval(stt(1:nframes,nIpts-nXi-nSabanapts:nIpts-nXi-nSabanapts+n_topo+n_cont+n_vall))
!     print*,"maxstt=",maxstt
      stt = stt/maxstt * MeshVecLen2 !;print*,"new max=",maxval(stt(:,:))
      srt = srt/maxsrt * MeshVecLen2
      CALL METAFL('PNG')
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(7800,4),int(3000,4))
!     CALL PAGMOD('NONE')
      call imgfmt('RGB')
      call winsiz(int(3120,4),int(1200,4))
      CALL SCRMOD('REVERS') !fondo blanco
      do i=1,nframes
      !#< r ############################################ DESPLAZAMIENTOS !#>
      write(titleN,'(a,I0,a)') 'mecElem_',i,'.png'
      CALL SETFIL(trim(titleN))
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
      CALL TEXMOD ('ON') ! latex!!
           !the position of an axis system.
      CALL axspos (int(0,4) ,int(2600,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600,4)) !size of the axis system.
!     call name('X [m]','X')
!     call name('Z [m]','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(2,4),'Y')
      call ticks (int(1,4) ,'XY')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call graf(real(MeshMinX,4),real(MeshMaxX,4),real(MeshMinX,4),real(xstep,4), &
                 real(MeshMaxZ,4),real(MeshMinZ,4),real(MeshMaxZ,4),real(-zstep,4))

      !estratigrafía --------------------------------------------
      call shdpat(int(16,4))                                    !
      ii = 1
      if (Z(0) .lt. 0.0) ii = 0
      do J=ii,N                                                 !
         call SETRGB(col(J), col(J), col(J))                    !
         call rlrec(real(MeshMinX,4),real(max(Z(J),MeshMinZ),4),&       !
                    real(MeshMaxX-MeshMinX,4),real(Z(J+1)-max(Z(J),MeshMinZ),4))!
         call color ('FORE')                                    !
         call rline(real(MeshMinX,4),real(max(Z(J),MeshMinZ),4), &      !
                 real(MeshMaxX,4),real(max(Z(J),MeshMinZ),4))           !
      end do                                                    !
      J = N+1                                                   !
      call SETRGB(col(J), col(J), col(J))                       !
      call rlrec(real(MeshMinX,4),real(Z(J),4),&                    !
                    real(MeshMaxX-MeshMinX,4),real(MeshMaxZ-Z(J),4))        !
      call color ('FORE')                                       !
      call rline(real(MeshMinX,4),real(Z(J),4), &                   !
                 real(MeshMaxX,4),real(Z(J),4))                     !
      ! Borrar estratos en la cuenca                                  !
      call color ('BACK')                                             !
      call shdpat(int(16,4))                                          !
      if (abs(z(0)) .lt. 0.0001) then
      call rlrec(real(MeshMinX,4),real(MeshMinZ,4),&                          !
                    real(MeshMaxX-MeshMinX,4),real(Z(1)-MeshMinZ,4))              !
      end if!
      !#< r ##################                    topografia original    !#>
      call marker(int(-1,4))
        if (allocated(recIncluonly)) then
        call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
        call shdpat(int(16,4))
        CALL RLAREA(real(recIncluonly(:,1),4), &
                  real(recIncluonly(:,2),4), &
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))
        end if!
        if (allocated(recVoidonly)) then
        call color ('BACK')                                             !
        call shdpat(int(16,4))                                          !
        CALL RLAREA(real(recVoidonly(:,1),4), &
                  real(recVoidonly(:,2),4), &                !
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))              !
        end if
      call color ('FORE')                                           !
      call PENWID(real(5.0,4))
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do

      ii = 1
      if (n_topo .gt. 0) then !#< r ######  topografía medio estratificado !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo
      fai = nPtsolos+1!nIpts-nXi-nSabanapts
      faf = nPtsolos+nBPt_topo!nIpts-nXi-nSabanapts+n_topo
      !print*,"hay n_topo =",fai,faf
      do j=fai,faf
        if (allpoints(j)%atBou) then
          delX(ii,i,1) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2),4)!U
          delZ(ii,i,1) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1),4)!W
          delX(ii,i,2) = real(allpoints(j)%center%x + escalaFlechas2 * allpoints(j)%S(i,2),4)!U
          delZ(ii,i,2) = real(allpoints(j)%center%z + escalaFlechas2 * allpoints(j)%S(i,1),4)!W
          ii = ii + 1
        end if
      end do
      delX(ii,i,1) = delX(1,i,1)
      delZ(ii,i,1) = delZ(1,i,1)
      delX(ii,i,2) = delX(1,i,2)
      delZ(ii,i,2) = delZ(1,i,2)
!     call color ('BACK')
!     call shdpat(int(16,4))
!     CALL RLAREA(real(delX(1:ii),4), &
!                 real(delZ(1:ii),4), &
!                 int(ii,4))

      do j = 1,ii-1 ! para cada elementocall !odograma
      call PENWID(real(2.0,4))
      call color ('RED')
      do iii = 1,i-1 ! de un tiempo a otro
        call rline(delX(j,iii,2),delZ(j,iii,2),delX(j,iii+1,2),delZ(j,iii+1,2))
      end do!
      end do
      call PENWID(real(20.0,4))
      call color ('ORANGE')
      do j = 1,ii-1
        call rline(delX(j,i,1),delZ(j,i,1),delX(j+1,i,1),delZ(j+1,i,1))
      end do
      ii = ii + 1
      end if!
      if (n_cont .gt. 0) then !#< r ##############  frontera de continuidad !#>
      k = 100000
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont
      fai = nPtsolos+nBPt_topo+1!nIpts-nXi-nSabanapts+n_topo+1
      faf = nPtsolos+nBPt_topo+nBPt_cont!nIpts-nXi-nSabanapts+n_topo+n_cont
      !print*,"hay n_cont =",fai,faf
      do j=fai,faf
!     print*,j
        if (allpoints(j)%atBou) then
          k = min(ii,k)
          delX(ii,i,1) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2),4)!x
          delZ(ii,i,1) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1),4)!y
          delX(ii,i,2) = real(allpoints(j)%center%x + escalaFlechas2 * allpoints(j)%S(i,2),4)!x
          delZ(ii,i,2) = real(allpoints(j)%center%z + escalaFlechas2 * allpoints(j)%S(i,1),4)!y
          ii = ii + 1
        end if
      end do
      delX(ii,i,1) = delX(k,i,1)
      delZ(ii,i,1) = delZ(k,i,1)
      delX(ii,i,2) = delX(k,i,2)
      delZ(ii,i,2) = delZ(k,i,2)
!     print*,"   ",ii,k," -> ",delX(ii),delZ(ii)
!     call SETRGB(shadecolor_inc, shadecolor_inc, shadecolor_inc)     !
!     CALL RLAREA(real(delX(k:ii),4), &
!                 real(delZ(k:ii),4), &
!                 int(ii-k+1,4))
      do j = 1,ii-1 ! para cada elementocall !odograma
      call PENWID(real(2.0,4))
      call color ('RED')
      do iii = 1,i-1 ! de un tiempo a otro
        call rline(delX(j,iii,2),delZ(j,iii,2),delX(j,iii+1,2),delZ(j,iii+1,2))
      end do!
      end do
      call PENWID(real(20.0,4))
      call color ('ORANGE')
      do j = k,ii-1
        call rline(delX(j,i,1),delZ(j,i,1),delX(j+1,i,1),delZ(j+1,i,1))
      end do
      ii = ii + 1
      end if!
      if (n_vall .gt. 0) then !#< r ########  frontera libre en incusión !#>
      k = 100000
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
      fai = nPtsolos+nBPt_topo+nBPt_cont+1!nIpts-nXi-nSabanapts+n_topo+n_cont+1
      faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall!nIpts-nXi-nSabanapts+n_topo+n_cont+n_vall
      !print*,"hay n_vall =",fai,faf
      do j=fai,faf
!     print*,j
        if (allpoints(j)%atBou) then
          k = min(ii,k)
          delX(ii,i,1) = real(allpoints(j)%center%x + escalaFlechas * allpoints(j)%S(i,2),4)!x
          delZ(ii,i,1) = real(allpoints(j)%center%z + escalaFlechas * allpoints(j)%S(i,1),4)!y
          delX(ii,i,2) = real(allpoints(j)%center%x + escalaFlechas2 * allpoints(j)%S(i,2),4)!x
          delZ(ii,i,2) = real(allpoints(j)%center%z + escalaFlechas2 * allpoints(j)%S(i,1),4)!y
          ii = ii + 1
        end if
      end do
      delX(ii,i,1) = delX(k,i,1)
      delZ(ii,i,1) = delZ(k,i,1)
      delX(ii,i,2) = delX(k,i,2)
      delZ(ii,i,2) = delZ(k,i,2)
!     call color ('BACK')
!     call shdpat(int(16,4))
!     CALL RLAREA(real(delX(k:ii),4), &
!                 real(delZ(k:ii),4), &
!                 int(ii-k+1,4))
      do j = 1,ii-1 ! para cada elementocall !odograma
      call PENWID(real(2.0,4))
      call color ('RED')
      do iii = 1,i-1 ! de un tiempo a otro
        call rline(delX(j,iii,2),delZ(j,iii,2),delX(j,iii+1,2),delZ(j,iii+1,2))
      end do!
      end do
      call PENWID(real(20.0,4))
      call color ('ORANGE')
      do j = k,ii-1
        call rline(delX(j,i,1),delZ(j,i,1),delX(j+1,i,1),delZ(j+1,i,1))
      end do
      ii = ii + 1
      end if!
      call PENWID(real(1.0,4))
      if (makevideo) then
      !#< r campo de desplazamientos ———————————————————-------------!#>
      call color ('FORE')                                            !
      call vecclr(-1) ! (-2):color de las puntas de flecha activado  !
      CALL VECOPT(real(escalaFlechas,4),'SCALE')                     !
      CALL VECOPT(real(10.0,4),'ANGLE')                              !
      CALL VECOPT(real(0.9,4),'LENGTH')                              !
      call vecmat(xvmat(:,:,i), &
                  yvmat(:,:,i), & ! porque los desplazamientos negativos van para arriba
                  npixX, npixZ,xpray,ypray,int(1201))
      end if ! workboundary

      !#< r cota  MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ !#>  MeshVecLen,madmax
      ! 0,2700 ! Lower left corner
      call PENWID(real(1.0,4))
      call color ('BACK')                                             !
      call shdpat(int(16,4))
      call rlrec(real(MeshMinX+(MeshMaxX-MeshMinX)*0.05,4),&
                 real(MeshMaxZ-(MeshMaxZ-MeshMinZ)*0.05,4),&
                 real(MeshVecLen,4),&
                 real((MeshMaxZ-MeshMinZ)*0.03,4))
      call color ('FORE')
      call rlrec(real(MeshMinX+(MeshMaxX-MeshMinX)*0.05 + 0.5*MeshVecLen,4),&
                 real(MeshMaxZ-(MeshMaxZ-MeshMinZ)*0.035 ,4),&
                 real(0.5*MeshVecLen,4),&
                 real((MeshMaxZ-MeshMinZ)*0.015,4))
      call rlrec(real(MeshMinX+(MeshMaxX-MeshMinX)*0.05 + 0.25*MeshVecLen,4),&
                 real(MeshMaxZ-(MeshMaxZ-MeshMinZ)*0.05 ,4),&
                 real(0.25*MeshVecLen,4),&
                 real((MeshMaxZ-MeshMinZ)*0.015,4))
      call rlrec(real(MeshMinX+(MeshMaxX-MeshMinX)*0.05,4),&
                 real(MeshMaxZ-(MeshMaxZ-MeshMinZ)*0.035 ,4),&
                 real(0.25*MeshVecLen,4),&
                 real((MeshMaxZ-MeshMinZ)*0.015,4))
      call PENWID(real(5.0,4))
      write(CTIT,'(a,ES8.2E2,a)') 'max. |$u_{i}$| = ',madmax,' '
      CALL HEIGHT (int(80,4))
      CALL MESSAG(CTIT,int((30),4),int(2660,4))

      call color ('FORE')
      call height(80) ! de los caracteres
      tlabel = (i)*real(dt,4)
      write(CTIT,'(a,F9.5,a)') '$t= ',tlabel,' seg$'
      lentitle = NLMESS(CTIT)
      CALL MESSAG(CTIT,int(7800-lentitle-100,4),int(2850,4))

      call height(80) ! de los caracteres
      write(CTIT,'(a)') '$u_{i}$'
      lentitle = NLMESS(CTIT)
      call color ('FORE')
      CALL MESSAG(CTIT,int(15,4),int(15,4))
      CALL ENDGRF
      !#< r ################################################### ESFUERZOS TANGENCIALES !#>
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
      CALL TEXMOD ('ON') ! latex!!
           !the position of an axis system.
      CALL axspos (int(2600,4) ,int(2600,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600,4)) !size of the axis system.
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(2,4),'Y')
      call ticks (int(1,4) ,'XY')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call PENWID(real(1.0,4))
      call color ('FORE')
      call height(40) ! de los caracteres
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), &
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
      call marker(int(-1,4))
      !#< r ##################                    topografia original    !#>
      ii = 1
      if (Z(0) .lt. 0.0) ii = 0
      call color ('FORE')                                       !
      do J=ii,N                                                 !
         call rline(real(minX,4),real(max(Z(J),minY),4), &      !
                 real(maxX,4),real(max(Z(J),minY),4))           !
      end do                                                    !
      J = N+1                                                   !
      call rline(real(minX,4),real(Z(J),4), &                   !
                 real(maxX,4),real(Z(J),4))                     !

        if (allocated(recIncluonly)) then
        call SETRGB(0.8, 0.8, 0.8)
        call shdpat(int(16,4))
        CALL RLAREA(real(recIncluonly(:,1),4), &
                  real(recIncluonly(:,2),4), &
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))
        end if!
        if (allocated(recVoidonly)) then
        call color ('BACK')
        call shdpat(int(16,4))
        CALL RLAREA(real(recVoidonly(:,1),4), &
                  real(recVoidonly(:,2),4), &
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))
        end if
        ! normales
        !fai = nIpts-nXi-nSabanapts-nSecciones
        !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
        fai = nPtsolos+1
        faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall
        call color ('CYAN')
        call PENWID(real(2,4))
        CALL HSYMBL(int(60,4))
        do j=fai,faf
        if (allpoints(j)%atBou) then
            call RLVEC(real(allpoints(j)%center%x,4),&
                       real(allpoints(j)%center%z,4),&
                       real(allpoints(j)%center%x+allpoints(j)%normal%x*MeshVecLen2/6,4),&
                       real(allpoints(j)%center%z+allpoints(j)%normal%z*MeshVecLen2/6,4),&
                       int(1111,4))
        end if
        end do!j

      !#< r ##################   TANGENCIALES     !#>
      call mypat(45,5,3,1)
      call PENWID(real(2.5,4))
      if (n_topo .gt. 0) then !#< r ######  topografía medio estratificado !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo
      fai = nPtsolos+1
      faf = nPtsolos+nBPt_topo!nIpts-nXi-nSabanapts+n_topo
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + stt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + stt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
      IF ((stt(i,j) .gt. 0) .and. (stt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .gt. 0) .and. (stt(i,j2) .lt. 0)) then
        rat = stt(i,j)/(stt(i,j)-stt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .gt. 0)) then
        rat = stt(i,j2)/(stt(i,j2)-stt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      if (n_cont .gt. 0) then !#< r ##############  frontera de continuidad !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont
      fai = nPtsolos+nBPt_topo+1
      faf = nPtsolos+nBPt_topo+nBPt_cont
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + stt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + stt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
!     IF (stt(i,j) .gt. 0) then
!     call color ('RED')
!     else; call color ('BLUE'); end if
!     CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
!       end if
      IF ((stt(i,j) .gt. 0) .and. (stt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .gt. 0) .and. (stt(i,j2) .lt. 0)) then
        rat = stt(i,j)/(stt(i,j)-stt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .gt. 0)) then
        rat = stt(i,j2)/(stt(i,j2)-stt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      if (n_vall .gt. 0) then !#< r ########  frontera libre en incusión !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
      fai = nPtsolos+nBPt_topo+nBPt_cont+1
      faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + stt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + stt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
!     IF (stt(i,j) .gt. 0) then
!     call color ('RED')
!     else; call color ('BLUE'); end if
!     CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
!       end if
      IF ((stt(i,j) .gt. 0) .and. (stt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((stt(i,j) .gt. 0) .and. (stt(i,j2) .lt. 0)) then
        rat = stt(i,j)/(stt(i,j)-stt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((stt(i,j) .lt. 0) .and. (stt(i,j2) .gt. 0)) then
        rat = stt(i,j2)/(stt(i,j2)-stt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + stt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + stt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      call color ('FORE')
      call PENWID(real(5.0,4))
      do j=1,nXI
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), &
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))
      end do
      ! stt = stt/maxstt * MeshVecLen
      call shdpat(int(16,4))
      call PENWID(real(1.0,4))
      call color ('GRAY')
      call rlrec(real(minX+(maxX-minX)*0.045,4),&
                 real(maxY-(maxY-minY)*0.055,4),&
                 real(MeshVecLen2*1.1,4),&
                 real((maxY-minY)*0.03*1.3,4))

      call color ('BACK')
      call rlrec(real(minX+(maxX-minX)*0.05,4),&
                 real(maxY-(maxY-minY)*0.05,4),&
                 real(MeshVecLen2,4),&
                 real((maxY-minY)*0.03,4))
      call color ('FORE')
      call rlrec(real(minX+(maxX-minX)*0.05 + 0.5*MeshVecLen2,4),&
                 real(maxY-(maxY-minY)*0.035 ,4),&
                 real(0.5*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))
      call rlrec(real(minX+(maxX-minX)*0.05 + 0.25*MeshVecLen2,4),&
                 real(maxY-(maxY-minY)*0.05 ,4),&
                 real(0.25*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))
      call rlrec(real(minX+(maxX-minX)*0.05,4),&
                 real(maxY-(maxY-minY)*0.035 ,4),&
                 real(0.25*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))

      call PENWID(real(5.0,4))
      write(CTIT,'(a,ES8.2E2)') 'max. |$\sigma_{\theta \theta}$| = ',maxstt
      lentitle = NLMESS(CTIT)
      CALL HEIGHT (int(80,4))
      CALL MESSAG(CTIT,int((2630),4),int(2660,4))

      call height(80) ! de los caracteres
      write(CTIT,'(a)') '$\sigma_{\theta \theta}$'
      lentitle = NLMESS(CTIT)
      call color ('FORE')
      CALL MESSAG(CTIT,int(2615,4),int(15,4))

      CALL ENDGRF
      !#< r ########################################################### CORTANTES !#>
      call color ('FORE')
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
      CALL TEXMOD ('ON') ! latex!!
           !the position of an axis system.
      CALL axspos (int(5200,4) ,int(2600,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600,4)) !size of the axis system.
!     call name('X [m]','X')
!     call name('Z [m]','Y')
!     call labdig(int(1,4),'X') !number of decimal places for labels
!     call labdig(int(2,4),'Y')
!     call ticks (int(1,4) ,'XY')
      call PENWID(real(1.0,4))
      call color ('FORE')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), &
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
      call marker(int(-1,4))
      CALL HSYMBL(int(25,4)) !size of symbols
      !#< r ##################                    topografia original    !#>
      ii = 1
      if (Z(0) .lt. 0.0) ii = 0
      call color ('FORE')                                       !
      do J=ii,N                                                 !
         call rline(real(minX,4),real(max(Z(J),minY),4), &      !
                 real(maxX,4),real(max(Z(J),minY),4))           !
      end do                                                    !
      J = N+1                                                   !
      call rline(real(minX,4),real(Z(J),4), &                   !
                 real(maxX,4),real(Z(J),4))                     !

        if (allocated(recIncluonly)) then
        call SETRGB(0.8, 0.8, 0.8)     !
        call shdpat(int(16,4))
        CALL RLAREA(real(recIncluonly(:,1),4), &
                  real(recIncluonly(:,2),4), &
                  int(2*size(Xcoord_Incluonly(:,1,1)),4))
        end if!
        if (allocated(recVoidonly)) then
        call color ('BACK')                                             !
        call shdpat(int(16,4))                                          !
        CALL RLAREA(real(recVoidonly(:,1),4), &
                  real(recVoidonly(:,2),4), &                !
                  int(2*size(Xcoord_Voidonly(:,1,1)),4))              !
        end if
      call color ('FORE')                                           !
      call PENWID(real(5.0,4))
      do j=1,nXI                                                      !
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), & !
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))   !
      end do
      ! normales
        !fai = nIpts-nXi-nSabanapts-nSecciones
        !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
        fai = nPtsolos+1
        faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall
        call color ('CYAN')
        call PENWID(real(2,4))
        CALL HSYMBL(int(60,4))
        do j=fai,faf
        if (allpoints(j)%atBou) then
            call RLVEC(real(allpoints(j)%center%x,4),&
                       real(allpoints(j)%center%z,4),&
                       real(allpoints(j)%center%x+allpoints(j)%normal%x*MeshVecLen2/6,4),&
                       real(allpoints(j)%center%z+allpoints(j)%normal%z*MeshVecLen2/6,4),&
                       int(1111,4))
        end if
        end do!j
      !#< r ##################   CORTANTES     !#>
      call mypat(45,5,3,1)
      call PENWID(real(2.5,4))
      if (n_topo .gt. 0) then !#< r ######  topografía medio estratificado !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo
      fai = nPtsolos+1
      faf = nPtsolos+nBPt_topo
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + srt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + srt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
      IF ((srt(i,j) .gt. 0) .and. (srt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .gt. 0) .and. (srt(i,j2) .lt. 0)) then
        rat = srt(i,j)/(srt(i,j)-srt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .gt. 0)) then
        rat = srt(i,j2)/(srt(i,j2)-srt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      if (n_cont .gt. 0) then !#< r ##############  frontera de continuidad !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont
      fai = nPtsolos+nBPt_topo+1
      faf = nPtsolos+nBPt_topo+nBPt_cont
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + srt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + srt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
!     IF (srt(i,j) .gt. 0) then
!     call color ('RED')
!     else; call color ('BLUE'); end if
!     CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
!       end if
      IF ((srt(i,j) .gt. 0) .and. (srt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .gt. 0) .and. (srt(i,j2) .lt. 0)) then
        rat = srt(i,j)/(srt(i,j)-srt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .gt. 0)) then
        rat = srt(i,j2)/(srt(i,j2)-srt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      if (n_vall .gt. 0) then !#< r ########  frontera libre en incusión !#>
      !fai = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+1
      !faf = nIpts-nXi-nSabanapts-nSecciones+n_topo+n_cont+n_vall
      fai = nPtsolos+nBPt_topo+nBPt_cont+1
      faf = nPtsolos+nBPt_topo+nBPt_cont+nBPt_vall
      do j=fai,faf
        if (allpoints(j)%atBou) then
      rec(1,1) = real(allpoints(j)%center%x ,4)
      rec(1,2) = real(allpoints(j)%center%z ,4)

      rec(2,1) = real(allpoints(j)%center%x + srt(i,j)*allpoints(j)%normal%x ,4)
      rec(2,2) = real(allpoints(j)%center%z + srt(i,j)*allpoints(j)%normal%z ,4)
      j2 = j+1
      if (j .eq. faf) j2 = fai

      rec(3,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
      rec(3,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)

      rec(4,1) = real(allpoints(j2)%center%x ,4)
      rec(4,2) = real(allpoints(j2)%center%z ,4)

      rec(5,1) = real(allpoints(j)%center%x ,4)
      rec(5,2) = real(allpoints(j)%center%z ,4)
!     IF (srt(i,j) .gt. 0) then
!     call color ('RED')
!     else; call color ('BLUE'); end if
!     CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
!       end if
      IF ((srt(i,j) .gt. 0) .and. (srt(i,j2) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((srt(i,j) .gt. 0) .and. (srt(i,j2) .lt. 0)) then
        rat = srt(i,j)/(srt(i,j)-srt(i,j2))
        rec(3,1) = real((allpoints(j)%center%x*(1-rat) + allpoints(j2)%center%x*rat),4)
        rec(3,2) = real((allpoints(j)%center%z*(1-rat) + allpoints(j2)%center%z*rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((srt(i,j) .lt. 0) .and. (srt(i,j2) .gt. 0)) then
        rat = srt(i,j2)/(srt(i,j2)-srt(i,j))
        rec(3,1) = real((allpoints(j)%center%x*rat + allpoints(j2)%center%x*(1-rat)),4)
        rec(3,2) = real((allpoints(j)%center%z*rat + allpoints(j2)%center%z*(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(allpoints(j2)%center%x ,4)
        rec(1,2) = real(allpoints(j2)%center%z ,4)
        rec(2,1) = real(allpoints(j2)%center%x + srt(i,j2)*allpoints(j2)%normal%x ,4)
        rec(2,2) = real(allpoints(j2)%center%z + srt(i,j2)*allpoints(j2)%normal%z ,4)
        rec(4,1) = real(allpoints(j2)%center%x ,4)
        rec(4,2) = real(allpoints(j2)%center%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      end if
      end do
      end if!
      call color ('FORE')
      call PENWID(real(5.0,4))
      do j=1,nXI
      call rline(real(Xcoord_ER(j,1,1),4),real(Xcoord_ER(j,2,1),4), &
                 real(Xcoord_ER(j,1,2),4),real(Xcoord_ER(j,2,2),4))
      end do
      ! srt = srt/maxstt * MeshVecLen
      call shdpat(int(16,4))
      call PENWID(real(1.0,4))
      call color ('GRAY')
      call rlrec(real(minX+(maxX-minX)*0.045,4),&
                 real(maxY-(maxY-minY)*0.055,4),&
                 real(MeshVecLen2*1.1,4),&
                 real((maxY-minY)*0.03*1.3,4))

      call color ('BACK')
      call rlrec(real(minX+(maxX-minX)*0.05,4),&
                 real(maxY-(maxY-minY)*0.05,4),&
                 real(MeshVecLen2,4),&
                 real((maxY-minY)*0.03,4))
      call color ('FORE')
      call rlrec(real(minX+(maxX-minX)*0.05 + 0.5*MeshVecLen2,4),&
                 real(maxY-(maxY-minY)*0.035 ,4),&
                 real(0.5*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))
      call rlrec(real(minX+(maxX-minX)*0.05 + 0.25*MeshVecLen2,4),&
                 real(maxY-(maxY-minY)*0.05 ,4),&
                 real(0.25*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))
      call rlrec(real(minX+(maxX-minX)*0.05,4),&
                 real(maxY-(maxY-minY)*0.035 ,4),&
                 real(0.25*MeshVecLen2,4),&
                 real((maxY-minY)*0.015,4))

      ! -.-.-.-.-.-.-.-..
      call PENWID(real(5.0,4))
      write(CTIT,'(a,ES8.2E2)') 'max. |$\sigma_{r \theta}$| = ',maxsrt
      lentitle = NLMESS(CTIT)
      CALL HEIGHT (int(80,4))
      CALL MESSAG(CTIT,int((5230),4),int(2660,4))

      call height(80) ! de los caracteres
      write(CTIT,'(a)') '$\sigma_{r \theta}$'
      lentitle = NLMESS(CTIT)
      call color ('FORE')
      CALL MESSAG(CTIT,int(5215,4),int(15,4))
      CALL ENDGRF
      call disfin

      end do ! i=1,n_maxtime

      !  -framerate #   antes de -i para hacerlo más lento. Donde # es menor a 25 (default)
      write(titleN,'(a)')'ffmpeg -i mecElem_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p 0_MecElemvideo.mp4'
      print*,trim(titleN)
      call system(trim(titleN))

!     if (encuadre - 0.5 .le. 0.1) then
!     write(titleN,'(a,I0,a,I0,a)') 'ffmpeg -i 0_MecElemvideo.mp4 -filter:v ''''crop=1200:',&
!           int(encuadre*1200+150),':0:',1200-int(encuadre*1200+150),''''' 0_MecElemvideoCrop.mp4'
!     print*,trim(titleN)
!     call system(trim(titleN))
!     end if
      end subroutine CINETECA

!      function getstt(i,j)
!      ! esfuerzo tangencial sigma theta theata
!      use resultvars, only : allpoints
!      implicit none
!      integer :: i,j
!      real*8 :: si,co,getstt,szz,szx,sxx
!      si = allpoints(j)%sinT
!      co = allpoints(j)%cosT
!      !W U sxx szx szz
!      sxx = real(allpoints(j)%S(i,3),8)
!      szx = real(allpoints(j)%S(i,4),8)
!      szz = real(allpoints(j)%S(i,5),8)
!
!      getstt = si**2*sxx + co**2*szz - 2*si*co*szx
!!     print*,sxx,szx,szz,getstt
!      end function getstt

      end module cine
