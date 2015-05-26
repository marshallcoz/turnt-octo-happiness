      module datos
      save
!     real*8, dimension(2), parameter :: rho = (/2.0, 2.0/)
!     real*8, dimension(2), parameter :: nu = (/0.333, 0.333/)
!     real*8, dimension(2), parameter :: bet = (/500.0, 750.0/) !1 medio, 2 liner
      real*8, dimension(2), parameter :: rho = (/2700.0, 7850.0/)
      real*8, dimension(2), parameter :: nu = (/0.25, 0.3/)
      real*8, dimension(2), parameter :: bet = (/3333.333, 3397.23/) !1 medio, 2 liner
      real*8, dimension(2) :: alf
      integer,parameter :: nRes = 120
      real*8, dimension(2) :: radios ! a, b 
      real, parameter :: Qq = 10000.0, Ts = 0.19, Tp = 0.06
      real*8, parameter :: DFREC = 0.66666!0.33
      integer, parameter :: NFREC = 150, nplot = 150, NPTSTIME = 2048
      integer, parameter :: nmax = 500, nfracs = 1
      integer, parameter :: frameInicial = 1, nframes = 168
      real*8, parameter :: Dns = 1.
      complex*16, parameter :: UI = cmplx(0.0d0,1.0d0,8), &
                               UR = cmplx(1.0d0,0.0d0,8), &
                               Z0 = cmplx(0.0d0,0.0d0,8)
      real*8, parameter :: PI = real(4.0d0*ATAN(1.0d0),8)
      logical, parameter :: imprimirEspectros = .true.
      logical, parameter :: hacerSnapshots = .false.
      real :: ventana
      real, parameter :: MeshVecLen = 1.0, giro = 0!-PI/2
      contains
      subroutine set_radios
!     radios(1) = 5.00_8 ! a   : in the liner
!     radios(2) = 5.60_8 ! b   : in the medium
      radios(1) = 5.45_8 ! a   : in the liner
      radios(2) = 5.50_8 ! b   : in the medium
      ventana = 6.5
      end subroutine set_radios
      end module datos
      module vars_func_of_w
      save
      complex*16 :: cOME
      real*8 :: FREC,OME,OMEI
      complex*16, dimension(2) :: lambda,amu
      
      complex*16, dimension(2,2),target :: VEL !( velocidad, region)
      !    cp(1)  cp(2) 
      !    cs(2)  cs(2)
      complex*16, dimension(:), pointer :: cp,cs 
      
      ! compressional and shear wave numbers: 
      complex*16, dimension(2,2),target :: w_c  ! w/c ( velocidad, region) !(alfa:beta,reg1:reg2)
      !    w/cp(1)  w/cp(2) 
      !    w/cs(2)  w/cs(2)
      complex*16, dimension(:), pointer :: alfa,beta 
      
      !ratios of the medium and the liner
      real*8 :: eta
      complex*16 :: muR,gammaP,gammaS
      
      !wavenumbers (squared)
      complex*16, dimension(2) :: s2,p2
      
      ! delta t
      real*8 :: Dt
      end module vars_func_of_w
      module RES 
      use datos, only : NFREC,nRes,NPTSTIME
      type receptor
      real*8 :: x,z,r,th,n(2) ! coordenadas
      integer :: ir
      integer :: reg ! 1: region de afuera
                     ! 2: inclusion
      complex*16,dimension(NFREC) :: & 
        s_rr,s_tt,s_rt,u_r,u_t
        
      complex*16,dimension(nptstime,5) :: S !s_tt,u_r,u_t,u_x,u_z
      end type receptor
      
      type(receptor), dimension(nRes),target :: Rw
      complex*16,dimension(NFREC) :: abscisa
      complex*16, dimension(NPTSTIME) :: Uo
      contains
      
      subroutine setup_resu
      use datos, only : pi,radios,z0
      integer :: i,iside,iini(2),ifin(2),thisradio(2),ir
      real*8 :: delth
      iini(1) = 1             ; ifin(1) = int(nRes/2) ; thisradio(1) = 2!b  !reg 1 in the medium
      iini(2) = 1+int(nRes/2) ; ifin(2) = nRes        ; thisradio(2) = 1!a  !reg 2 in the lining
      print*,"putnos receptroes ------------------"
      print*,"      r      th      x      z      region"
      delth = (2.0_8*pi) / int(nRes/2)
      do iside = 1,2 ! medio,cilindro
      ir = 0
      do i = iini(iside),ifin(iside)
        Rw(i)%r = radios(thisradio(iside))
        RW(i)%ir = thisradio(iside)
        Rw(i)%th = delth * ir ; ir = ir + 1
        Rw(i)%reg = iside
        write(6,'(I0,F8.2,F8.2,3x,I0)') & 
        i,RW(i)%r,RW(i)%th,RW(i)%reg
        call initRW(i)
      end do
      end do
!       Rw(1)%r = radios(2) !r=b in the medium at the boundary of the liner
!       Rw(1)%th = pi!/2.
!       Rw(1)%reg = 1 ! in the medium
!      
!       Rw(2)%r = radios(1) !r=a at the inner side of the cylinder
!       Rw(2)%th = pi!/2.
!       Rw(2)%reg = 2 ! in the liner
!     stop "setup_resu"
      end subroutine setup_resu
      
      subroutine initRW(i)
      use datos, only : z0
      implicit none 
      integer :: i
          Rw(i)%s_rr = z0
          Rw(i)%s_tt = z0
          Rw(i)%s_rt = z0
          Rw(i)%u_r = z0
          Rw(i)%u_t = z0
      end subroutine initRW
      end module RES
       
      module debug
      contains
      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(A,E9.2,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E9.2,A)',advance='no')      AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      end module debug
      module plotter 
      contains
      subroutine plot(i_Px,i_nom,reg, abscisa ,pt_RES,nplot)
      use dislin
      implicit none
      integer :: i_Px,i_nom,reg,nplot
      character(len=4),dimension(5) :: nom
      complex*16,dimension(nplot) :: abscisa,pt_RES
      character(len=100) :: tx
      nom(1) = "s_rr"
      nom(2) = "s_tt"
      nom(3) = "s_rt"
      nom(4) = "u_r_"
      nom(5) = "u_t_"
      ! graficar
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4P')
      write(tx,'(a,i0,a,i0,a)') "RES_",i_Px,nom(i_nom),reg,".pdf"
      CALL SETFIL(trim(adjustl(tx)))
      
      call qplcrv(real(abscisa,4), & 
                 real(real(pt_RES),4), nplot,'FIRST')
      call color('RED')
      call qplcrv(real(abscisa,4), & 
                 real(real(pt_RES),4), nplot,'NEXT')
      call color('BLUE')
      call qplcrv(real(abscisa,4), & 
                 real(aimag(pt_RES),4), nplot,'LAST')
      
      write(tx,'(a,i0,a,i0,a)') "RES_",i_Px,nom(i_nom),reg,"abs.pdf"
      CALL SETFIL(trim(adjustl(tx)))
      CALL color('BLACK')
      call qplot(real(abscisa,4), & 
                 real(abs(pt_RES),4), nplot)
      end subroutine plot
      
      subroutine plotSpectrum(y_in,Df,full_n,n,titleN,xAx,yAx,logflag,CTIT,W,H,maxfrec)
      ! (Uo,DFREC,size(Uo),size(Uo)/2.0,titleN,xAx,yAx,logflag,1200,800,maxfrec)
      use DISLIN
      implicit none
      real, intent(in)                              :: Df,maxfrec
      integer, intent(in)                           :: full_n,n,H,W
      character(LEN=100), intent(in)                 :: xAx
      character(LEN=100), intent(in)                 :: yAx
      character(LEN=9)                             :: logflag
      character(LEN=100)                            :: titleN,CTIT
      COMPLEX*16, DIMENSION(full_n), intent(in) :: y_in
      complex,    dimension(:), allocatable :: y
      real,       dimension(:), allocatable :: x
      real maxY,minY,maxYc,minYc,xstep,ystep
      integer :: i
      character(LEN=100) :: dumb
      CHARACTER(LEN=6)  :: CBUF
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      allocate(x(n))
      allocate(y(n))
      DO i = 1,n
        x(i) = Df*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
!       !write(6,*) x(i), y(i)
      END DO
      
!     
      minY=MINVAL(real(y(:)),1)
!     !write(6,*)"MinReal Val= ",minY
      maxY=MAXVAL(real(y(:)),1)
!     !write(6,*)"MaxReal Val= ",maxY
      minYc=MINVAL(aimag(y(:)),1)
!     !write(6,*)"MinComplex Val= ",minYc
      maxYc=MAXVAL(aimag(y(:)),1)
!     !write(6,*)"MaxComplex Val= ",maxYc
      minY =MIN(minY,minYc)
      maxY =MAX(maxY,maxYc)
!     !write(6,*)"MinY Val= ",minY
!     !write(6,*)"MaxY Val= ",maxY
      
      logflag = trim(adjustl(logflag))
      if (trim(logflag) == 'logx') then
       logflag = 'X'
       ! los ceros no son muy populares:
       x(1)=x(2)/2.
      elseif (trim(logflag) == 'logy') then
       logflag = 'Y'
!     ! minY =MIN(minY,minYc,-0.1)
!     ! maxY =MAX(maxY,maxYc, 0.1)
      elseif (trim(logflag) == 'loglog') then
       logflag = 'XY'
       ! los ceros no son muy populares:
       x(1)=x(2)/2.
!     ! minY =MIN(minY,minYc,-0.1)
!     ! maxY =MAX(maxY,maxYc, 0.1)
      else
       logflag = 'none'
      end if
      
      
      
      ! Dislin plotting routines 
      CALL METAFL('PDF') !define intended display  XWIN,PS,EPS,PDF
!     !write(titleN,'(a,a)') trim(titleN),'.eps' 
!     !titleN = trim(titleN)
!     !titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
      call errmod ("all", "off")
!     !CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL TEXMOD ('ON') ! latex!!
      CALL HWFONT()
      CALL axspos (int(370,4) ,int(H+100,4)) !Axis system pos. Lower left corner
      call axslen (int(W,4) ,int(H,4)) !size of the axis system.
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(5,4) ,'XY') 
!     !call titlin ( titleN , 1 )
      xstep = x(n) /6. ! increment for labels 
      ystep = (maxY-minY)/6.0
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      
      !
      if (logflag == 'Y') then
        CALL AXSSCL('LOG',trim(logflag))
!       !print*,x(n)
!       !print*,log10(x(n))
        call graf(real(x(1),4),real(x(n),4) &
             ,real(x(1),4) ,real(xstep,4) &
             ,real(-7.0,4),real(log10(maxY),4) &
             ,real(-7.0,4),real(1.0,4))  
             
      elseif (logflag == 'X' .or. logflag == 'XY') then
        CALL AXSSCL('LOG',trim(logflag))
        call graf(real(-2.0,4),real(log10(x(n)),4) &
             ,real(-2.0,4) ,real(1.0,4) &
             ,real(minY,4),real(maxY,4) &
             ,real(minY,4),real(ystep,4))
      else
        call graf(real(0.0,4),real(x(n),4), & 
                  real(0.0,4),real(max(1.0,xstep),4), &
                  real(minY,4),real(maxY,4), & 
                  real(minY,4),real(ystep,4))
      end if
      
      
      call color ('RED')
      call curve(real(abs(x),4) ,real(y,4) ,int(n,4))
      call color('BLUE')
      call curve(real(abs(x),4), real(aimag(y),4), int(n,4))
      call color ('FORE') 
      call curve(real(abs(x),4), & 
                    real(sqrt(real(y)**2.0 +aimag(y)**2.0),4), int(n,4))
!     
!     maxfrec
      call color ('FORE') 
      call rline(real(maxfrec,4),real(minY,4), & 
                 real(maxfrec,4),real(maxY,4))
!     
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
!     
!     
      call legini(CBUF,int(3,4),int(20,4))
!     !nx = nxposn(x(n) + x(n) / 20.)
!     !ny = nyposn(minY + (maxY-minY)*0.7)
!     !print*,nx
!     !print*,ny
      call legpos(int(1600,4),int(320,4))
      write(dumb,'(a)') 'Re(z) '
!     !print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      write(dumb,'(a)') 'Im(z) '
!     !print*,dumb
      call leglin(CBUF,dumb,int(2,4))
      write(dumb,'(a)') 'Abs(z)'
!     !print*,dumb
      call leglin(CBUF,dumb,int(3,4))
      
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(3,4))
      CALL MESSAG (CTIT, int(1600,4), int(100,4))
      call disfin()
      
!     print*,'plotted ',trim(titleN)
!     !print*,''
      deallocate(x,y)
      end subroutine plotSpectrum
      
                            
      subroutine plotXYcomp(y_in,Dt,n,titleN,xAx,yAx,CTIT,W,H,ma)
      ! (Uo,Dt,size(Uo),'FIGURE_NAME.pdf','time[sec]','amplitude',1200,800) 
      USE DISLIN
!     use glovars, only : pi
      implicit none
      real, intent(in)                              :: Dt,ma
      integer, intent(in)                           :: n,H,W
      character(LEN=9), intent(in)                  :: xAx
      character(LEN=100), intent(in)                :: yAx
      character(LEN=100)                            :: titleN
      COMPLEX*16, DIMENSION(n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      real             :: yr,yi
      character(LEN=50) :: tx
      
      real, dimension(n) :: x
      real maxY,minY,xstep,ystep,maxYc,minYc,val
      integer :: i,nPow10x,nPow10y,signo
!     integer :: Modo,nx,ny
      character(LEN=100) :: dumb
      CHARACTER(LEN=30) :: CBUF
      character(LEN=100) :: CTIT
!     integer*4 :: lentitle
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      
!     allocate(x(n))
!     print*,size(y)
      
      ! para que se vea el texto en los ejes si es muy pequeño en formato F3.1
      nPow10x = 0
      if (Dt *(n-1) < 0.6) then
        do i = 1,10
          if (Dt *(n-1)*(10.0**i) > 1.0) then
            exit
          end if 
        end do 
        nPow10x = i
        
      elseif (Dt * (n-1) > 6000.) then  
        do i = 1,10
          if (Dt *(n-1)*(10.0**(-i)) < 1000.0) then
            exit
          end if
        end do
        nPow10x = -i
      end if
      DO i = 1,n
        x(i) = Dt*(i-1)*(10.0**(nPow10x))
        write(tx,'(EN18.2,2x,EN18.2)') real(y_in(i)),aimag(y_in(i))
        read(tx,*) yr,yi
        y(i) = cmplx(yr,yi)
      END DO
      if (abs(ma) .lt. 0.001_4) then
      minY=minval(real(y(:)),1)
!     write(6,*)"MinReal Val= ",minY
      maxY=MAXVAL(real(y(:)),1)
!     write(6,*)"MaxReal Val= ",maxY
      minYc=MINVAL(aimag(y(:)),1)
!     write(6,*)"MinComplex Val= ",minYc
      maxYc=MAXVAL(aimag(y(:)),1)
!     write(6,*)"MaxComplex Val= ",maxYc
      minY =MIN(minY,minYc)
      maxY =MAX(maxY,maxYc)
      else
      maxy = ma
      miny = -ma
      end if
      val = max(abs(miny),abs(maxy))
      signo = 1
      if (miny*maxy < 0) then
        if (val - maxy > 0.0001) signo = -1
      else
        if (maxy < 0.0)  signo = -1
      end if
      nPow10y = 0
      do i = 1,10
!       if (
      end do
!     print*,"plotting"
! Dislin plotting routines ...
      CALL METAFL('PDF') !define intended display XWIN,PS,EPS,PDF
!     write(titleN,'(a,a)') trim(titleN),'.eps' 
!     titleN = trim(titleN)
!     titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
!     print*,"file: ",trim(adjustl(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize
      call errmod ("all", "off") 
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL TEXMOD ('ON') ! latex!!
      CALL HWFONT()
      CALL axspos (int(380,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W+650,4) ,int(H,4)) !size of the axis system.
      if (nPow10x .ne. 0) then
       write(dumb,'(a,a,I0)') trim(xAx),'x10^',(nPow10x *(-1))
       call name(trim(dumb),'X')
      else
       call name(trim(xAx),'X') 
      end if
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(2,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(10,4) ,'X') 
      call ticks (int(5,4) ,'Y') 
!     call titlin ( titleN , 1 )
      xstep = x(n)/6.0 ! incremen for labels
      ystep = (maxY-minY)/6.0
      
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0 !solo 3 etiquetas
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      
      
      
      
!     print*,"maxmin= ",maxy,miny
      call graf(real(x(1),4), & 
                real(x(n)+x(2),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
      
!     call title 
      call color ('GREEN') 
      call curve(real(x,4) ,real(y,4) ,int(n,4))
      call color('BLUE')
      call curve(real(x,4), real(aimag(y),4), int(n,4))
      call color ('FORE') 
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
      call legini(CBUF,int(2,4),int(20,4))
 !     nx = nxposn(x(n)*n + x(n)*n / 20.)
 !     ny = nyposn(minY + (maxY-minY)*0.7)
 !     print*,nx
 !     print*,ny
      if (maxval(abs(aimag(y))) .gt. 0.05*maxval(abs(real(y)))) then
      call legpos(int(1840,4),int(720,4))
      write(dumb,'(a)') 'Re(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      write(dumb,'(a)') 'Im(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(2,4))
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(2,4))
      end if
      
!     write(CTIT,'(a,ES11.4E2,a)') 'dt=',Dt,' seg'
!     lentitle = NLMESS(CTIT)
!     CALL MESSAG(CTIT,int((1700),4),int(200,4))
      CALL MESSAG (CTIT, int(900,4), int(110,4))
!     call errmod ("protocol", "off") !suppress dislin info
      call disfin()      
      
!     print*,'plotted ',trim(titleN)
      end subroutine plotXYcomp
      end module plotter
      module Hank
      use datos, only : nmax, nfracs
      use RES, only : nRes
      type argvelocidad
        complex*16, dimension (2) :: c ! 1 para alfa; 2 para beta
      end type argvelocidad
      
      type argradiovelocidad
        type (argvelocidad), dimension(2) :: r ! 1 para a;  2 para b; ///nomas 3... para receptores
      end type argradiovelocidad
      
      type Bessel
        type(argradiovelocidad),dimension(1:4,-1:nmax*nfracs+nfracs) :: JYH1H2 ! -1,0,1,2,...,nmax
      ! i = 1 :: Jn
      ! i = 2 :: Yn
      ! i = 3 :: H(1)n
      ! i = 4 :: H(2)n
      end type Bessel
      
      type(Bessel), dimension(2),target :: Bess ! 1 para afuera; 2 para recubrimiento
      
      contains
      subroutine makeBessels(minNS,imprimo)
      USE SPECFUN
      use vars_func_of_w, only : w_c !( velocidad, region)
      use datos, only : nmax,nfracs,Dns,radios,ui,z0
      use RES, only : nRes,Rw
      use dislin
      use debug
      use plotter
      implicit none
      logical :: imprimo
      integer :: reg,c,r,imec,i,ns,minNS
      complex*16 :: z!,H20,H21,H22
      real*8 :: v,vm
      complex*16, dimension(0:nmax) :: CBJ,CDJ,CBY,CDY
      character(Len=3), dimension(4) :: nom
      character(len=100) :: tx
      real*4,dimension(nmax) :: x, pt_RES 
      nom(1) = "__J"
      nom(2) = "__Y"
      nom(3) = "_H1"
      nom(4) = "_H2"
      minNS = 100000
      do reg=1,2  ! afuera 1 (medio) y adentro 2 (lining)
      do r=1,2 !a,b
!     if (r .ge. 3) then
!       if (Rw(r-2)%reg .ne. reg) cycle
!     end if
      do c=1,2 !alfa y beta 
      Z = w_c(c,reg) * radios(r)  !(alfa:beta,reg1:reg2)
      do i = 1,nfracs
      V = Dns * real(i-1,8) + real(nmax,8)
!     print*,V,Z
      Bess(reg)%JYH1H2(1,i-1:(nmax*nfracs)+(i-1):nfracs)%r(r)%c(c) = z0
      call CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
      ns = int(vm)
      minNS = min(minNS,ns)
      ! J
      Bess(reg)%JYH1H2(1,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns)
      ! Y
      Bess(reg)%JYH1H2(2,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBY(0:ns)
      ! H1
      Bess(reg)%JYH1H2(3,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns) + UI*CBY(0:ns)
      ! H2
      Bess(reg)%JYH1H2(4,i-1:(ns*nfracs)+(i-1):nfracs)%r(r)%c(c) = CBJ(0:ns) - UI*CBY(0:ns)
      end do ! i
      Bess(reg)%JYH1H2(1:4,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(1:4,1)%r(r)%c(c)
      
!     print*,"z=",Z
!     do i = 0,ns
!     print*,i,",J=",Bess(reg)%JYH1H2(1,i)%r(r)%c(c), & 
!              ",Y=",Bess(reg)%JYH1H2(2,i)%r(r)%c(c)
!     end do ! i
!     stop 588
      
      
      
!     ! Hankel of the second kind 
!     ! n = 0,1,2
!     z = w_c(c,reg) * radios(r) ! w_c(c,reg) * radios(r)
!     call hankels(z,H20,H21)
!     H22 = -H20 + 2./z * H21
!     Bess(reg)%JYH1H2(4,-1)%r(r)%c(c) = -H21 !-1
!     Bess(reg)%JYH1H2(4,0)%r(r)%c(c) = H20   !-0
!     Bess(reg)%JYH1H2(4,1)%r(r)%c(c) = H21   !+1
!     Bess(reg)%JYH1H2(4,2)%r(r)%c(c) = H22   !+2
!     ! n = 3,4,...nmax
!     do n=2,nmax
!     Bess(reg)%JYH1H2(4,n+1)%r(r)%c(c) = & 
!     2* n / z * Bess(reg)%JYH1H2(4,n)%r(r)%c(c) - Bess(reg)%JYH1H2(4,n-1)%r(r)%c(c)
!     end do
!     ! Hankel of the first kind
!     ! n = 0,1,2
!     ! tomando de la continuaci´on anal´itica modificamos el argumento
!     z = w_c(c,reg) * radios(r) * exp(-pi*UI) 
!     call hankels(z,H20,H21)
!     H22 = -H20 + 2./z * H21
!     Bess(reg)%JYH1H2(3,-1)%r(r)%c(c) = -H21
!     Bess(reg)%JYH1H2(3,0)%r(r)%c(c) = H20
!     Bess(reg)%JYH1H2(3,1)%r(r)%c(c) = H21
!     Bess(reg)%JYH1H2(3,2)%r(r)%c(c) = H22
!     do n=2,nmax
!     Bess(reg)%JYH1H2(3,n+1)%r(r)%c(c) = & 
!     2 * n / z * Bess(reg)%JYH1H2(3,n)%r(r)%c(c) - Bess(reg)%JYH1H2(3,n-1)%r(r)%c(c)
!     end do
!     ! continuación analítica
!     do n=-1,nmax+1
!     Bess(reg)%JYH1H2(3,n)%r(r)%c(c) = exp(-real(n)*pi*UI) * Bess(reg)%JYH1H2(3,n)%r(r)%c(c) ! (-) ?
!     end do
!     
!     z = w_c(c,reg) * radios(r) ! alfa(c) * radios(r)
!     call BESJC(z,nmax,Bess(reg)%JYH1H2(1,0:nmax+1)%r(r)%c(c))
!     Bess(reg)%JYH1H2(1,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(1,1)%r(r)%c(c) 
!     z = w_c(c,reg) * radios(r) ! alfa(c) * radios(r)
!     call BESYC(z,nmax,Bess(reg)%JYH1H2(2,0:nmax+1)%r(r)%c(c))
!     Bess(reg)%JYH1H2(2,-1)%r(r)%c(c) = - Bess(reg)%JYH1H2(2,1)%r(r)%c(c) 
!     
      if (imprimo .eqv. .true.) then
      do imec = 1,4
        write(tx,'(a,a,i0,a,i0,a,i0,a)') nom(imec),"_",reg,"_",r,"_",c,".pdf"
        print*,trim(tx)," | z=",radios(r) * w_c(c,reg)
        call showMNmatrixZ(NMAX,1,Bess(reg)%JYH1H2(imec,0:NMAX)%r(r)%c(c),"-----",6)
          
          x = (/(real(i*dns,4),i=0,nmax-1)/)
      CALL METAFL('PDF')
      call filmod('DELETE')
      CALL SETPAG('DA4P')
      CALL SETFIL(trim(adjustl(tx)))
      
          pt_RES = real(Bess(reg)%JYH1H2(imec,0:NMAX-1)%r(r)%c(c),4)
      call qplcrv(real(x,4),pt_RES, nmax,'FIRST')
      call color('RED')
      call qplcrv(real(x,4),pt_RES, nmax,'NEXT')
      call color('BLUE')
      pt_RES = real(aimag (Bess(reg)%JYH1H2(imec,0:NMAX-1)%r(r)%c(c)),4)
      call qplcrv(real(x,4),pt_RES, nmax,'LAST')
      
          
      end do
      end if
      end do !c
      end do !r
      end do !reg
      end subroutine makeBessels
      end module hank

      module fft
      contains
      SUBROUTINE FORK(LX,CX,SIGNI)
      implicit none
      integer, intent(in) :: LX,SIGNI
      COMPLEX*16 :: CARG,CW,CTEMP 
      complex*16,intent(inout) :: CX(LX)
      real*8, parameter :: pi = 4.*ATAN(1.)
      real*8 :: SC
      integer :: i,j,m,istep,l
      J=1
      SC=DSQRT(real(1.0,8)/real(LX,8))
      DO 30 I=1,LX
      IF(I > J)GO TO 10
      CTEMP=CX(J)*cmplx(SC,0.0,8)
      CX(J)=CX(I)*cmplx(SC,0.0,8)
      CX(I)=CTEMP
   10 M=LX/2
   20 IF(J <= M)GO TO 30
      J=J-M
      M=M/2
      IF(M >= 1)GO TO 20
   30 J=J+M
      L=1
   40 ISTEP=2*L
      DO 50 M=1,L
      CARG=cmplx(0.0,(pi*real(SIGNI*(M-1)))/real(L),8)  
      CW=EXP(CARG)
      DO 50 I=M,LX,ISTEP
      CTEMP=CW*CX(I+L)
      CX(I+L)=CX(I)-CTEMP
   50 CX(I)=CX(I)+CTEMP
      L=ISTEP
      IF(L < LX)GO TO 40
      RETURN
      END subroutine fork
      end module fft
        module specfun
        contains
        SUBROUTINE CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
!
!       ===========================================================
!       Purpose: Compute Bessel functions Jv(z), Yv(z) and their 
!                derivatives for a complex argument
!       Input :  z --- Complex argument
!                v --- Order of Jv(z) and Yv(z)
!                      ( v = n+v0, n = 0,1,2,..., 0 Û v0 < 1 )
!       Output:  CBJ(n) --- Jn+v0(z)
!                CDJ(n) --- Jn+v0'(z)
!                CBY(n) --- Yn+v0(z)
!                CDY(n) --- Yn+v0'(z)
!                VM --- Highest order computed
!       Routines called:
!            (1) GAMMA for computing the gamma function
!            (2) MSTA1 and MSTA2 for computing the starting
!                point for backward recurrence
!       ===========================================================
!
        IMPLICIT DOUBLE PRECISION (A,B,G,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:*),CDJ(0:*),CBY(0:*),CDY(0:*)
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PV0=PI*V0
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           IF (V0.EQ.0.0) THEN
              CBJ(0)=(1.0D0,0.0D0)
              CDJ(1)=(0.5D0,0.0D0)
           ELSE
              CDJ(0)=(1.0D+300,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0D0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CJV0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+V0))
              CJV0=CJV0+CR
              IF (CDABS(CR).LT.CDABS(CJV0)*1.0D-15) GO TO 20
15         CONTINUE
20         VG=1.0D0+V0
           CALL GAMMA(VG,GA)
           CA=(0.5D0*Z1)**V0/GA
           CJV0=CJV0*CA
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           VV=4.0D0*V0*V0
           CPZ=(1.0D0,0.0D0)
           CRP=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-&
                          (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
25            CPZ=CPZ+CRP
           CQZ=(1.0D0,0.0D0)
           CRQ=(1.0D0,0.0D0)
           DO 30 K=1,K0
              CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-&
                          (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
30            CQZ=CQZ+CRQ
           CQZ=0.125D0*(VV-1.0)*CQZ/Z1
           ZK=Z1-(0.5D0*V0+0.25D0)*PI
           CA0=CDSQRT(RP2/Z1)
           CCK=CDCOS(ZK)
           CSK=CDSIN(ZK)
           CJV0=CA0*(CPZ*CCK-CQZ*CSK)
           CYV0=CA0*(CPZ*CSK+CQZ*CCK)
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              CJVN=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 35 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K-V0))
                 CJVN=CJVN+CR
                 IF (CDABS(CR).LT.CDABS(CJVN)*1.0D-15) GO TO 40
35            CONTINUE
40            VG=1.0D0-V0
              CALL GAMMA(VG,GB)
              CB=(2.0D0/Z1)**V0/GB
              CJU0=CJVN*CB
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 45 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
45               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
           ENDIF
        ENDIF
        IF (N.EQ.0) N=1
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 50 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1-CF2
           IF (K.LE.N) CBJ(K)=CF
           CF2=CF1
50         CF1=CF
        CS=CJV0/CF
        DO 55 K=0,N
55         CBJ(K)=CS*CBJ(K)
        IF (REAL(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
           ENDIF
           DO 60 K=0,N
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBJ(K)=CDEXP(-PI*(K+V0)*CI)*CBJ(K)
              ELSE IF (DIMAG(Z).GT.0.0D0) THEN
                 CBJ(K)=CDEXP(PI*(K+V0)*CI)*CBJ(K)
              ENDIF
60         CONTINUE
           Z1=Z1
        ENDIF
        CBY(0)=CYV0
        DO 65 K=1,N
           CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           CBY(K)=CYY
65      CONTINUE
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 70 K=1,N
70         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 75 K=1,N
75         CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
        VM=N+V0
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,&
              -0.6558780715202538D0, -0.420026350340952D-1,&
               0.1665386113822915D0,-.421977345555443D-1,&
               -.96219715278770D-2, .72189432466630D-2,&
               -.11651675918591D-2, -.2152416741149D-3,&
               .1280502823882D-3, -.201348547807D-4,&
               -.12504934821D-5, .11330272320D-5,&
               -.2056338417D-6, .61160950D-8,&          
                .50020075D-8, -.11812746D-8,&          
                .1043427D-9, .77823D-11,&          
               -.36968D-11, .51D-12,&          
               -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward  
!                recurrence such that the magnitude of    
!                Jn(x) at that point is about 10^(-MP)
!       Input :  x     --- Argument of Jn(x)
!                MP    --- Value of magnitude
!       Output:  MSTA1 --- Starting point   
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
!
!       ===================================================
!       Purpose: Determine the starting point for backward
!                recurrence such that all Jn(x) has MP
!                significant digits
!       Input :  x  --- Argument of Jn(x)
!                n  --- Order of Jn(x)
!                MP --- Significant digit
!       Output:  MSTA2 --- Starting point
!       ===================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
        end module specfun
      PROGRAM tunelRecub
! Solución analítica de la difración de una onda P incidente 
! en una cavidad cilíndrica circular con recubrimiento
      use datos; use vars_func_of_w; use Hank
      use debug; use RES; use plotter; use fft
      implicit none
      integer :: J,n,et,info,i,ii,ir,minNS
      real*8, pointer :: r,th
      integer, pointer :: reg
      complex*16, dimension(6,6) :: M
      complex*16, dimension(6,2) :: B !P , SV
      integer, dimension(6) :: ipiv 
      complex*16,dimension(nplot) :: pt_RES
      ! terms defined in the appendix (functions)
      complex*16 :: e11,e12,e21,e22,e41,e42,e71,e72,e81,e82,sig0
      real*8,dimension(2) :: BEALF
      character(LEN=4)          :: nom(5)
      character(LEN=9)          :: logflag
      character(LEN=100)        :: titleN,xAx,yAx,CTIT
      
!     complex*16,dimension(:,:),allocatable :: AF
!     complex*16,dimension(:,:),allocatable :: Xbem
!     complex*16,dimension(:),allocatable :: WORKbem
!     real*8,dimension(:),allocatable :: Rbem,Cbem,RWORK
!     real*8 :: RCOND,FERR,BERR
!     character*1 :: EQUED
      
      call set_radios !a,b
      call setup_resu !puntos receptores
      Dt = (1.0) / (real(NPTSTIME) * DFREC)
      eta = radios(2)/radios(1) !b/a
      BEALF(1:2)=SQRT((0.5-NU(1:2))/(1.0-NU(1:2))) !IF POISSON RATIO IS GIVEN
      alf(1:2) = bet(1:2)/BEALF(1:2)
      print*,"alf", alf(1),alf(2)
      print*,"bet", bet(1),bet(2)
      print*,"mu1=",(RHO(1) * bet(1)**2.)
      print*,"mu2=",(RHO(2) * bet(2)**2.)
      print*,"muR=",(RHO(1) * bet(1)**2.)/(RHO(2) * bet(2)**2.)
      print*,"gamm=alf1/alf2=",alf(1) / alf(2)
      print*,"nu1=",NU(1)
      print*,"nu2=",NU(2)
      print*,"eta= b/a =",eta 
      print*,"DFREC=",DFREC 
      print*,"Dt=",Dt
      print*,"Fmax= [Hz]", DFREC* NFREC
      print*,"Tmax= [s]",Dt*NPTSTIME
      write(6,'(A)', ADVANCE = "NO") " "
      write(6,'(A)', ADVANCE = "NO") repeat("_",58)
      write(6,'(A)', ADVANCE = "NO") " "
      print*,""
      write(6,'(A)', ADVANCE = "NO") repeat(" ",60)
      
!     allocate(AF(6,6))
!     allocate(Xbem(6,1))
!     allocate(Rbem(6))
!     allocate(Cbem(6))
!     allocate(WORKbem(2*6))
!     allocate(RWORK(2*6))
      do J=1,NFREC !*********************************************
      FREC=DFREC*real(J); if (J .eq. 1)  FREC = 0.5_8 * DFREC ! Hz
      OME=2.0*PI*FREC !rad/s
!     VEL(1,1:2) = cmplx(alf(1:2)*& 
!                  (1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
!     VEL(2,1:2) = cmplx(bet(1:2)*& 
!                  (1.+1./pi/Qq*log(ome/2./pi)),-1./2./Qq,8)
      COME = OME * CMPLX(1., -1./(2.*Qq),8)
      VEL(1,1:2) = cmplx(alf(1:2),(/0.,0./),8)
      VEL(2,1:2) = cmplx(bet(1:2),(/0.,0./),8)
      cp(1:2) => VEL(1,1:2) ! dilatación
      cs(1:2) => VEL(2,1:2) ! corte
      w_c(1:2,1:2) = cOME / VEL(1:2,1:2) !(alfa:beta,reg1:reg2)
      beta(1:2) => w_c(2,1:2) !shear wave number
      alfa(1:2) => w_c(1,1:2) !compressional wave number
      p2(1:2) = (alfa(1:2))**2.0 !compressional wave number (square)
      s2(1:2) = (beta(1:2))**2.0 !shear wave number (square)
      aMU(1:2) = RHO(1:2) * cs(1:2)**2.
      Lambda(1:2) = RHO(1:2)* cp(1:2)**2. - real(2.)* aMU(1:2)
      muR = amu(1)/amu(2) !shear moduli ratio
      gammaP = z0! spacing variable  2.5D
      gammaS = z0! spacing variable
      call makeBessels(minNS,.false.) ! n = -1,0,1,2,...,vm (imprimir)
      do n=0,minNS!nmax*nfracs ! ensamblar matriz 4.26 y terminos independientes
      M(1:6,1:6) = z0; B(1:6,1) = z0; iPIV = 6
      ! sigma_{rr}1 = sigma_{rr}2   @ r = b
      M(1,1) = - muR * e11(3,1,1,2,n)
      M(1,2) = - muR * e12(3,2,1,2,n)
      M(1,3) = e11(3,1,2,2,n)
      M(1,4) = e11(4,1,2,2,n)
      M(1,5) = e12(3,2,2,2,n)
      M(1,6) = e12(4,2,2,2,n) 
      B(1,1) = (1.0) * et(n) * UI**n * muR * e11(1,1,1,2,n)
      ! sigma_{r th}1 = sigma_{r th}2   @ r = b
      M(2,1) = - muR * e41(3,1,1,2,n)
      M(2,2) = - muR * e42(3,2,1,2,n)
      M(2,3) = e41(3,1,2,2,n)
      M(2,4) = e41(4,1,2,2,n)
      M(2,5) = e42(3,2,2,2,n)
      M(2,6) = e42(4,2,2,2,n)
      B(2,1) = (1.0) * et(n) * UI**n * muR * e41(1,1,1,2,n)

      ! u_{r}1 = u_{r}2   @ r = b
      M(3,1) = - e71(3,1,1,2,n) 
      M(3,2) = - e72(3,2,1,2,n) 
      M(3,3) = e71(3,1,2,2,n)
      M(3,4) = e71(4,1,2,2,n) 
      M(3,5) = e72(3,2,2,2,n)
      M(3,6) = e72(4,2,2,2,n)
      B(3,1) = (1.0) * et(n) * UI**n * e71(1,1,1,2,n)
      ! u_{th}1 = u_{th}2   @ r = b
      M(4,1) = - e81(3,1,1,2,n) 
      M(4,2) = - e82(3,2,1,2,n) 
      M(4,3) = e81(3,1,2,2,n) 
      M(4,4) = e81(4,1,2,2,n) 
      M(4,5) = e82(3,2,2,2,n)
      M(4,6) = e82(4,2,2,2,n)
      B(4,1) = (1.0) * et(n) * UI**n * e81(1,1,1,2,n)
      ! sigma_{rr}2 = 0   @ r = a
      M(5,1) = z0
      M(5,2) = z0
      M(5,3) = e11(3,1,2,1,n)
      M(5,4) = e11(4,1,2,1,n)
      M(5,5) = e12(3,2,2,1,n)
      M(5,6) = e12(4,2,2,1,n)
      B(5,1) = z0
      ! sigma_{r th}2 = 0   @ r = a
      M(6,1) = z0
      M(6,2) = z0
      M(6,3) = e41(3,1,2,1,n)
      M(6,4) = e41(4,1,2,1,n)
      M(6,5) = e42(3,2,2,1,n)
      M(6,6) = e42(4,2,2,1,n)
      B(6,1) = z0

!#< r  solve system       . . !#>
      ! driver simple
      call zgesv(6,1,M(1:6,1:6),6,IPIV,B(:,1),6,info)
      
      ! driver experto:
!     call ZGESVX('E','N',6,1,M,6,AF,6,IPIV,EQUED,Rbem, &
!     Cbem,B,6,Xbem,6,RCOND,FERR,BERR,WORKbem,RWORK,INFO)
      
      if(info .ne. 0) then
        write(6,'(A,I0,a,I0)', ADVANCE = "NO") &
        "/",n, " info =",info
        if (info .gt. 6) write(6,'(A)', ADVANCE = "NO") " working precision"
        if (info .le. 6) write(6,'(A)', ADVANCE = "NO") "factor U is exactly singular"
!       stop 1120
        exit
!     else if (abs(B(1,1)) .lt. 0.00000001) then  !NaN
!      write(6,'(A,I0,a)', ADVANCE = "NO") &
!     "trim at",n, " por chiquito abs(B(1)) < 10^{-8}"
!      exit
      else
        write(6,'(A)', ADVANCE = "NO") "."
      end if
!     call showMNmatrixZ(6,1,B,"  A  ",6)     
      !#< r elementos mecanicos !#>
      do i = 1,nRes
!      r => Rw(i)%r
       ir = Rw(i)%ir
       th => Rw(i)%th
       reg => Rw(i)%reg
        if (reg .eq. 1) then  !reg 1 in the medium
          Rw(i)%s_rr(J) = Rw(i)%s_rr(J) + &
          (1. * et(n) * UI**n * e11(1,1,1,ir,n) + &
                       B(1,1) * e11(3,1,1,ir,n) + &
                       B(2,1) * e12(3,2,1,ir,n)) * (cos(n * th))
          Rw(i)%s_tt(J) = Rw(i)%s_tt(J) + &
          (1. * et(n) * UI**n * e21(1,1,1,ir,n) + & 
                       B(1,1) * e21(3,1,1,ir,n) + & 
                       B(2,1) * e22(3,2,1,ir,n)) * (cos(n * th))
          Rw(i)%s_rt(J) = Rw(i)%s_rt(J) + &
          (1. * et(n) * UI**n * e41(1,1,1,ir,n) + & 
                       B(1,1) * e41(3,1,1,ir,n) + & 
                       B(2,1) * e42(3,2,1,ir,n)) * (sin(n * th))
          Rw(i)%u_r(J) = Rw(i)%u_r(J) + &
          (1. * et(n) * UI**n * e71(1,1,1,ir,n) + & 
                       B(1,1) * e71(3,1,1,ir,n) + & 
                       B(2,1) * e72(3,2,1,ir,n)) * (cos(n * th))
          Rw(i)%u_t(J) = Rw(i)%u_t(J) + &
          (1. * et(n) * UI**n * e81(1,1,1,ir,n) + & 
                       B(1,1) * e81(3,1,1,ir,n) + & 
                       B(2,1) * e82(3,2,1,ir,n)) * (sin(n * th))
        else if (reg .eq. 2) then !reg 2 in the lining
          Rw(i)%s_rr(J) = Rw(i)%s_rr(J) + &
                      (B(3,1) * e11(3,1,2,ir,n) + &
                       B(4,1) * e11(4,1,2,ir,n) + &
                       B(5,1) * e12(3,2,2,ir,n) + &
                       B(6,1) * e12(4,2,2,ir,n)) * (cos(n * th))
          Rw(i)%s_tt(J) = Rw(i)%s_tt(J) + &
                      (B(3,1) * e21(3,1,2,ir,n) + &
                       B(4,1) * e21(4,1,2,ir,n) + &
                       B(5,1) * e22(3,2,2,ir,n) + &
                       B(6,1) * e22(4,2,2,ir,n)) * (cos(n * th))
          Rw(i)%s_rt(J) = Rw(i)%s_rt(J) + &
                      (B(3,1) * e41(3,1,2,ir,n) + &
                       B(4,1) * e41(4,1,2,ir,n) + &
                       B(5,1) * e42(3,2,2,ir,n) + &
                       B(6,1) * e42(4,2,2,ir,n)) * (sin(n * th))
          Rw(i)%u_r(J) = Rw(i)%u_r(J) + &
                      (B(3,1) * e71(3,1,2,ir,n) + &
                       B(4,1) * e71(4,1,2,ir,n) + &
                       B(5,1) * e72(3,2,2,ir,n) + &
                       B(6,1) * e72(4,2,2,ir,n)) * (cos(n * th))
          Rw(i)%u_t(J) = Rw(i)%u_t(J) + &
                      (B(3,1) * e81(3,1,2,ir,n) + &
                       B(4,1) * e81(4,1,2,ir,n) + &
                       B(5,1) * e82(3,2,2,ir,n) + &
                       B(6,1) * e82(4,2,2,ir,n)) * (sin(n * th))
        end if
      end do! i:nRes
      end do !n
      
      ! términos fuera de la suma:
      do i = 1,nRes
       r => Rw(i)%r
!      th => Rw(i)%th
       reg => Rw(i)%reg
       sig0 = amu(1) * beta(1)**2. !eq 3.15   para hacerlo factor de amplificación
          Rw(i)%s_rr(J) = Rw(i)%s_rr(J) * 2. *amu(reg) / r**2. 
          Rw(i)%s_tt(J) = Rw(i)%s_tt(J) * 2. *amu(reg) / r**2. !/ sig0
          Rw(i)%s_rt(J) = Rw(i)%s_rt(J) * 2. *amu(reg) / r**2. 
          Rw(i)%u_r(J) = Rw(i)%u_r(J) / r 
          Rw(i)%u_t(J) = Rw(i)%u_t(J) / r 
      end do! i:nRes

!       write(6,'(A)', ADVANCE = "NO") repeat(char(8),60)
        print*,""
        write(6,'(A)', ADVANCE = "NO") "["
        write(6,'(A,A)', ADVANCE = "NO") & 
        repeat("X",58-int((58.0/NFREC)*(NFREC+1-J))), &
        repeat("_",int((58.0/NFREC)*(NFREC+1-J)))
        write(6,'(A,I0)', ADVANCE = "NO") "]",minNS
      abscisa(J) = dfrec * 2 * pi * J * radios(1) / cp(1)
      end do !J   
      ! plot curves
      ! end program tunelRecub
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      !#< r print results !#>
      if (imprimirEspectros) then
      call system('mkdir outEspectrosOriginales')
      CALL chdir("outEspectrosOriginales")
      do i = 1,nRes
        reg => Rw(i)%reg 
!         pt_RES = Rw(i)%s_rr(1:nplot)
!         call plot(i,1,reg, abscisa ,pt_RES,nplot)
          pt_RES = Rw(i)%s_tt(1:nplot)
          call plot(i,2,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)%s_rt(1:nplot)
!         call plot(i,3,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_r(1:nplot)
!         call plot(i,4,reg, abscisa ,pt_RES,nplot)
!         pt_RES = Rw(i)% u_t(1:nplot)
!         call plot(i,5,reg, abscisa ,pt_RES,nplot)

      end do! i:nRes  
      CALL chdir("..")  
      end if
      

      !#< b func. amplitud !#>
      call ricker(Uo)
      write(titleN,'(a)') 'WaveAmplitude-ricker_time.pdf'
      write(CTIT,'(a)') 'WaveAmplitude of Ricker wavelet'
      xAx = 'time[sec]'
      write(yAx,'(a)') 'amplitude'
      call plotXYcomp(Uo(1:nframes),real(Dt,4), & 
                 nframes,titleN,xAx,yAx,CTIT,1200,800,0.0)
      
      call fork(size(Uo),Uo,-1) !forward
      Uo = Uo * sqrt(real(size(Uo))) * dt ! factor de escala
      
!     call fork(size(Uo),Uo,+1) !backward
!     Uo = Uo / sqrt(real(size(Uo))) / dt ! factor de escala
      
      
      write(titleN,'(a)') 'WaveAmplitude-ricker_frec.pdf'
      xAx = 'frec[Hz] '
      write(yAx,'(a)') 'amplitude'      
      logflag = 'logx     '      
!     logflag = 'none     '
      call plotSpectrum(Uo(:),real(DFREC,4), size(Uo(:)),int(size(Uo(:))/2), & 
            titleN,xAx,yAx,logflag,CTIT,1200,800,real(DFREC*(NFREC+1),4))

      !snapshots 
      if (imprimirEspectros) &
      call system('mkdir outEspectros')
      call system('mkdir outTraces')
        
      nom(1) = 's_tt'
      nom(2) = 'u__r'
      nom(3) = 'u__t'
      nom(4) = 'u__x'
      nom(5) = 'u__z'
      write(yAx,'(a)') 'amplitude'
      do i=1,nRes
        ! Coordenadas cartesianas
        Rw(i)%th = Rw(i)%th + giro
        Rw(i)%x = Rw(i)%r * cos(Rw(i)%th)
        Rw(i)%z = Rw(i)%r * sin(Rw(i)%th)
        Rw(i)%n(1) = cos(Rw(i)%th)
        Rw(i)%n(2) = sin(Rw(i)%th)
        write(CTIT,'(a,F7.2,a,F7.2,a)')'(', Rw(i)%x,' , ',Rw(i)%z,')'
        Rw(i)%S = 0
        
        !s_tt
        Rw(i)%S(1:nfrec,1) = Rw(i)%s_tt(1:nfrec)
        Rw(i)%S(NPTSTIME-NFREC+2:NPTSTIME,1) = conjg(Rw(i)%s_tt(nfrec:2:-1)) 
        
        !u_r
        Rw(i)%S(1:nfrec,2) = Rw(i)%u_r(1:nfrec)
        Rw(i)%S(NPTSTIME-NFREC+2:NPTSTIME,2) = conjg(Rw(i)%u_r(nfrec:2:-1)) 
        
        !u_t
        Rw(i)%S(1:nfrec,3) = Rw(i)%u_t(1:nfrec)
        Rw(i)%S(NPTSTIME-NFREC+2:NPTSTIME,3) = conjg(Rw(i)%u_t(nfrec:2:-1)) 
        
        !u_x
        Rw(i)%S(:,4) = cos(Rw(i)%th) * Rw(i)%S(:,2) - &
                       sin(Rw(i)%th) * Rw(i)%S(:,3)
        
        !u_z
        Rw(i)%S(:,5) = sin(Rw(i)%th) * Rw(i)%S(:,2) + &
                       cos(Rw(i)%th) * Rw(i)%S(:,3)
      
      
      print*,Rw(i)%r,Rw(i)%th,abs(Rw(i)%S(NFREC,1))
       
        
      do ii = 1,5  !#< g     al tiempo       !#>
      Rw(i)%S(:,ii) = Rw(i)%S(:,ii) * Uo
      if (imprimirEspectros) then
      write(titleN,'(a,a,I0,a)') 'f_',nom(ii),i,'.pdf'
      
      write(xAx,'(a)') 'Hz[sec]'
      CALL chdir("outEspectros")
!     call plotXYcomp(Rw(i)%S(1:nfrec,ii),real(DFREC,4), & 
!                nfrec,titleN,xAx,yAx,CTIT,1200,800,0.0)
      call plotSpectrum(Rw(i)%S(:,ii),real(DFREC,4), NPTSTIME, nfrec, & 
            titleN,xAx,yAx,logflag,CTIT,1200,800,real(DFREC*(NFREC+1),4))
      CALL chdir("..")
      end if
      
      call fork(NPTSTIME,Rw(i)%S(:,ii),+1) !backward
      Rw(i)%S(:,ii) = Rw(i)%S(:,ii) / sqrt(1.0*NPTSTIME) / dt ! factor de escala
     
      write(titleN,'(a,I0,a)') nom(ii),i,'.pdf'
      write(xAx,'(a)') 'time[sec]'
      CALL chdir("outTraces")
      call plotXYcomp(Rw(i)%S(1:nframes,ii),real(Dt,4), & 
                 nframes,titleN,xAx,yAx,CTIT,1200,800,0.0)
      CALL chdir("..")
        end do
      end do ! i:nRes
      if (hacerSnapshots) then
      call system('mkdir outSnapshots')
      CALL chdir("outSnapshots")
      call cineteca
      end if
      end program tunelRecub
      function et(n)
      integer,intent(in) :: n
      integer :: et
      et = 2
      if (n .eq. 0) et = 1
      end function et
      
      function e11(i,c,reg,r,n)
      use vars_func_of_w, only : alfa,s2,gammaP
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e11
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 1 !alfa
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(1)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(1)
      e11 = (n**2 + n - s2(reg) * radios(r)**2 / 2. + & 
      gammaP**2 * radios(r)**2) * Bess_n - alfa(reg)* radios(r) * Bess_n_1
      end function e11
      function e12(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e12
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 2 !beta
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(2)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(2)
      e12 = n * (-(n+1)* Bess_n + beta(reg) * radios(r) * Bess_n_1)
      end function e12
      function e21(i,c,reg,r,n)
      use vars_func_of_w, only : alfa,s2,p2
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e21
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 1 !alfa
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(1)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(1)
      e21 = - (n**2 + n + s2(reg)*radios(r)**2 / 2. - p2(reg)*radios(r)**2) * &
             Bess_n + alfa(reg) * radios(r) * Bess_n_1
      end function e21
      function e22(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e22
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 2 !beta
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(2)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(2)
      e22 = n * ((n+1)*Bess_n - beta(reg)*radios(r)*Bess_n_1)
      end function e22
      function e41(i,c,reg,r,n)
      use vars_func_of_w, only : alfa
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e41
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 1 !alfa
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(1)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(1)
      e41 = - n * (-(n+1) * Bess_n + alfa(reg)*radios(r)*Bess_n_1)
      end function e41
      function e42(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e42
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 2 !beta
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(2)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(2)
      e42 = - (n**2 + n - beta(reg)**2 * radios(r)**2 / 2.) * Bess_n + &
            beta(reg) * radios(r) * Bess_n_1
      end function e42
      function e71(i,c,reg,r,n)
      use vars_func_of_w, only : alfa
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e71
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 1 !alfa
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(1)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(1)
      e71 = alfa(reg) * radios(r) * Bess_n_1 - n * Bess_n
      end function e71
      function e72(i,c,reg,r,n)
      use Hank, only : Bess
      implicit none
      complex*16 :: e72
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n
!     c = 2 !beta
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(2)
      e72 = n * Bess_n
      end function e72
      function e81(i,c,reg,r,n)
      use Hank, only : Bess
      implicit none
      complex*16 :: e81
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n
!     c = 1 !alfa
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(1)
      e81 = - n * Bess_n
      end function e81
      function e82(i,c,reg,r,n)
      use vars_func_of_w, only : beta
      use datos, only : radios
      use Hank, only : Bess
      implicit none
      complex*16 :: e82
      integer :: i,c,reg,r,n
      complex*16, pointer :: Bess_n,Bess_n_1
!     c = 2 !beta
      Bess_n => Bess(reg)%JYH1H2(i,n)%r(r)%c(2)
      Bess_n_1 => Bess(reg)%JYH1H2(i,n-1)%r(r)%c(2)
      e82 = - (beta(reg) * radios(r) * Bess_n_1 & 
              - n * Bess_n)
      end function e82
      subroutine ricker(Uo)
      use datos, only : PI,UR,NPTSTIME,Ts,Tp
      use vars_func_of_w, only : Dt
      implicit none
      integer :: i
      complex*16, dimension(NPTSTIME) :: Uo
      real*8 :: A,B
      Uo = 0;
      A = nearest(pi*(-Ts) / Tp,1.0)
      A = nearest(real((A*A-0.5)* exp(- A * A), 8),1.0)
      if (abs(A) .lt. 0.0001) A = 0
      Uo(1) = A * UR
      
      do i = 2, NPTSTIME/2+1
      ! NEAREST(X, S) returns the processor-representable number 
      ! nearest to X in the direction indicated by the sign of S
        A = nearest(pi*(Dt*(i*1.0_8-1.0_8)-Ts) / Tp,1.0)
        A = nearest(A * A,1.0)
        B = nearest(exp(-A),1.0)
        A = nearest((A - 0.5_8) * B ,1.0)
        if (abs(A) .lt. 0.001) A = 0
        Uo(i) = A * UR
      end do
      end subroutine ricker
      
      subroutine CINETECA 
      use DISLIN
      use datos
      use RES
      use vars_func_of_w, only : dt
      implicit none
      real :: maV1,maV2,minX,maxX,minY,maxY,xstep,zstep,EsC,madmax
      character(LEN=100) :: titleN
      integer*4 :: lentitle
      character(LEN=60) :: CTIT
      real, dimension(5,2) :: rec
      real*8 :: rat
      integer :: i,iii,j,j2,iside,iini(2),ifin(2)
      
      ! area de la gráfica --------------------
      mav1 = -10000.
      do j=1,nRes
          mav1 = max(mav1,maxval(real(RW(j)%S(frameInicial:nframes,1),4)))
      end do
      EsC = real(MeshVecLen / mav1)
      do j=1,nRes
      RW(j)%S(frameInicial:nframes,1) = RW(j)%S(frameInicial:nframes,1) * EsC
      end do
      print*,"";print*,"mav1=",mav1,"  EsC=",EsC
      
      mav1 = -10000. ; mav2 = -10000.
      do j=1,nRes
          mav1 = max(mav1,maxval(real(RW(j)%S(frameInicial:nframes,4),4)))!x
          mav2 = max(mav2,maxval(real(RW(j)%S(frameInicial:nframes,5),4)))!z
      end do
      madmax = max(mav1,mav2)
      EsC = real(MeshVecLen / madmax)
      do j=1,nRes
      RW(j)%S(frameInicial:nframes,4) = RW(j)%S(frameInicial:nframes,4) * EsC
      RW(j)%S(frameInicial:nframes,5) = RW(j)%S(frameInicial:nframes,5) * EsC
      end do
      maxx = ventana
      minx = -maxx
      maxy = maxx
      miny = minx
      print*,"BOX=",minX,maxX,minY,maxY
      xstep = real(((maxX-minX) / 0 ))
      zstep = real(((maxY-minY) / 0 ))
      ! indicies de las fronteras
      iini(1) = 1             ; ifin(1) = int(nRes/2) !b  !reg 1 in the medium
      iini(2) = 1+int(nRes/2) ; ifin(2) = nRes        !a  !reg 2 in the lining
      
      
      CALL METAFL('PNG')
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(7800,4),int(3000,4))
      call imgfmt('RGB')
      call winsiz(int(3120,4),int(1200,4))
      CALL SCRMOD('REVERS') !fondo blanco
      do i=frameInicial,nframes
      write(6,'(a)',ADVANCE = "NO") "X"
!     !#< r ############################################ DESPLAZAMIENTOS !#>
      write(titleN,'(a,I0,a)') 'mecElem_',i,'.png'
      CALL SETFIL(trim(titleN))
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
      CALL TEXMOD ('ON') ! latex!!
!          !the position of an axis system.
      CALL axspos (int(0,4) ,int(2600,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600,4)) !size of the axis system.
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(2,4),'Y')
      call ticks (int(1,4) ,'XY')
      call PENWID(real(4.0,4))
      call color ('FORE')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call height(40) ! de los caracteres
      call graf(real(minX,4),real(MaxX,4),real(MinX,4),real(xstep,4), &
                 real(MaxY,4),real(MinY,4),real(MaxY,4),real(-zstep,4))
      
      do iside = 1,2
      do j = iini(iside),ifin(iside)
      j2 = j+1
      if (j .eq. ifin(iside)) j2 = iini(iside)
      
      call PENWID(real(20.0,4))
      call color ('FORE')
      call rline(real(Rw(j)%x,4), real(Rw(j)%z,4), &
                 real(Rw(j2)%x,4),real(Rw(j2)%z,4))
      call PENWID(real(16.0,4))
      call color ('ORANGE')
      call rline(real(Rw(j)%x+Rw(j)%S(i,4),4), real(Rw(j)%z+Rw(j)%S(i,5),4), &
                 real(Rw(j2)%x+Rw(j2)%S(i,4),4),real(Rw(j2)%z+Rw(j2)%S(i,5),4))
      end do !j
      end do !iside
                 
      call color ('FORE')
      call height(80) ! de los caracteres
      write(CTIT,'(a,F9.5,a)') '$t= ',(i)*real(dt,4),' seg$'
      lentitle = NLMESS(CTIT)
      call PENWID(real(4.0,4))
      CALL MESSAG(CTIT,int(7800-lentitle-100,4),int(2850,4))
      CALL ENDGRF
!     !#< r ################################################### ESFUERZOS TANGENCIALES !#>
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
      call PENWID(real(4.0,4))
      call color ('FORE')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call height(40) ! de los caracteres
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), &
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
      call marker(int(-1,4))
      call mypat(45,5,3,1)
      call PENWID(real(2.5,4))
      
      do iside = 1,2
      do j = iini(iside),ifin(iside)
      
      rec(1,1) = real(Rw(j)%x ,4)
      rec(1,2) = real(Rw(j)%z ,4)
      rec(2,1) = real(Rw(j)%x + Rw(j)%S(i,1) * Rw(j)%n(1) ,4)
      rec(2,2) = real(Rw(j)%z + Rw(j)%S(i,1) * Rw(j)%n(2) ,4)
      j2 = j+1
      if (j .eq. ifin(iside)) j2 = iini(iside)
      
      call color ('FORE')
      call rline(real(Rw(j)%x,4), real(Rw(j)%z,4), &
                 real(Rw(j2)%x,4),real(Rw(j2)%z,4))
                 
      rec(3,1) = real(Rw(j2)%x + Rw(j2)%S(i,1)* Rw(j2)%n(1) ,4)
      rec(3,2) = real(Rw(j2)%z + Rw(j2)%S(i,1)* Rw(j2)%n(2) ,4)
      rec(4,1) = real(Rw(j2)%x ,4)
      rec(4,2) = real(Rw(j2)%z ,4)
      rec(5,1) = real(Rw(j)%x ,4)
      rec(5,2) = real(Rw(j)%z ,4)
      IF ((real(Rw(j)%S(i,1),4) .gt. 0) .and. (real(Rw(j2)%S(i,1),4) .gt. 0)) then
        call color ('RED')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((real(Rw(j)%S(i,1),4) .lt. 0) .and. (real(Rw(j2)%S(i,1),4) .lt. 0)) then
        call color ('BLUE')
        CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      else if ((real(Rw(j)%S(i,1),4) .gt. 0) .and. (real(Rw(j2)%S(i,1),4) .lt. 0)) then
        rat = real(Rw(j)%S(i,1)) / real(Rw(j)%S(i,1) - Rw(j2)%S(i,1))
        rec(3,1) = real((Rw(j)%x *(1-rat) + Rw(j2)%x * rat),4)
        rec(3,2) = real((Rw(j)%z *(1-rat) + Rw(j2)%z * rat),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(Rw(j2)%x ,4)
        rec(1,2) = real(Rw(j2)%z ,4)
        rec(2,1) = real(Rw(j2)%x + Rw(j2)%S(i,1) * Rw(j2)%n(1) ,4)
        rec(2,2) = real(Rw(j2)%z + Rw(j2)%S(i,1) * Rw(j2)%n(2) ,4)
        rec(4,1) = real(Rw(j2)%x ,4)
        rec(4,2) = real(Rw(j2)%z ,4)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      else if ((real(Rw(j)%S(i,1),4) .lt. 0) .and. (real(Rw(j2)%S(i,1),4) .gt. 0)) then
        rat = real(Rw(j2)%S(i,1))/ real(Rw(j2)%S(i,1)-Rw(j)%S(i,1))
        rec(3,1) = real((Rw(j)%x * rat + Rw(j2)%x *(1-rat)),4)
        rec(3,2) = real((Rw(j)%z * rat + Rw(j2)%z *(1-rat)),4)
        rec(4,1) = rec(1,1)
        rec(4,2) = rec(1,2)
        call color ('BLUE')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
        rec(1,1) = real(Rw(j2)%x ,4)
        rec(1,2) = real(Rw(j2)%z ,4)
        rec(2,1) = real(Rw(j2)%x + Rw(j2)%S(i,1)* Rw(j2)%n(1) ,4)
        rec(2,2) = real(Rw(j2)%z + Rw(j2)%S(i,1)* Rw(j2)%n(2) ,4)
        rec(4,1) = real(Rw(j2)%x ,4)
        rec(4,2) = real(Rw(j2)%z ,4)
        call color ('RED')
        CALL RLAREA(rec(1:4,1),rec(1:4,2),int(4,4))
      end if
      
      end do !j
      end do !iside
      
      CALL ENDGRF
!     !#< r ########################################################### ODOGRAMA !#>
      call color ('FORE')
      CALL DISINI()
      call errmod ("all", "off")
      call incmrk (int(1,4))
      CALL TRIPLX()! CALL DISALF() !default font
      CALL TEXMOD ('ON') ! latex!!
           !the position of an axis system.
      CALL axspos (int(5200,4) ,int(2600,4)) ! Lower left corner
      call axslen (int(2600,4), int(2600,4)) !size of the axis system.
      call PENWID(real(4.0,4))
      call color ('FORE')
      call setgrf("NAME", "NAME", "LINE", "LINE")
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), &
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
      
      do iside = 1,2
      do j = iini(iside),ifin(iside)
      j2 = j+1
      if (j .eq. ifin(iside)) j2 = iini(iside)
      
      call PENWID(real(2.5,4))
      call color ('FORE')
      call rline(real(Rw(j)%x,4), real(Rw(j)%z,4), &
                 real(Rw(j2)%x,4),real(Rw(j2)%z,4))
      call PENWID(real(3.0,4))
      call color ('FORE')
      do iii = 1,i-1 ! de un tiempo a otro
      call rline(real(Rw(j)%x+Rw(j)%S(iii,4),4),   real(Rw(j)%z+Rw(j)%S(iii,5),4), &
                 real(Rw(j)%x+Rw(j)%S(iii+1,4),4), real(Rw(j)%z+Rw(j)%S(iii+1,5),4))
      end do
      call color ('ORANGE')
      CALL RLSYMB (5, real(Rw(j)%x+Rw(j)%S(iii,4),4),   real(Rw(j)%z+Rw(j)%S(iii,5),4))   !
      
      end do !j
      end do !iside
      CALL ENDGRF
      call disfin
      end do ! i=1,n_maxtime
 !     !  -framerate #   antes de -i para hacerlo más lento. Donde # es menor a 25 (default)
      write(titleN,'(a)')'ffmpeg -y -i mecElem_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p 0_MecElemvideo.mp4'
      print*,trim(titleN)
      call system(trim(titleN))
      end subroutine CINETECA

