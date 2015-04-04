      module gloVars 
      save
      ! verbose = 0   ! no output 
      !         = 1   ! main calls, properties and images
      !         = 2   ! 1 + counters in loops and subfunctions
      !         = 3   ! 2 + matrix values
      integer    :: verbose 
      logical    :: makeVideo
      logical    :: workBoundary
      logical    :: flip12
      logical    :: plotFKS
      logical    :: saveG
      integer    :: multSubdiv ! De 6 en adelante 
      real       :: cKbeta ! 1.1 - 3 (default 1.5)
      real       :: periodicdamper
      integer    :: developerfeature
      logical    :: borrarDirectorio
      logical    :: useAzimi
      integer    :: PrintNum 
      complex*16, dimension(:,:), allocatable :: developerAUXvec !sigma0
      complex*16, parameter :: UI = cmplx(0.0d0,1.0d0,8), &
                               UR = cmplx(1.0d0,0.0d0,8), &
                               Z0 = cmplx(0.0d0,0.0d0,8)
      real*8, target :: ONE = 1.0_8
      real*8, parameter :: PI = real(4.0d0*ATAN(1.0d0),8)
      CHARACTER(len=400) :: rutaOut
      end module gloVars 
      

      module soilVars
      save
      integer ::  N !number of layers from z>= 0. HALF-SPACE at N+1
      real*8 :: minBeta
      real*8,     dimension(:), allocatable :: ALFA0,BETA0,AMU0,LAMBDA0,ANU,Z,RHO
      complex*16 ,dimension(:), allocatable :: ALFA ,BETA ,AMU, LAMBDA
      real*8 :: Qq
      real*4, dimension(:),  allocatable :: layershadecolor
      complex*16 :: ALFA_inc,BETA_inc,AMU_inc,LAMBDA_inc
      real*4 :: shadecolor_inc
      end module soilVars
      

      module waveVars
      real, save :: Escala
      real*8, save :: Dt,maxtime  !segundos
      real, save :: t0
      integer, save :: ampfunction ! 0 dirac; 1 ricker
      complex*16, dimension(:), allocatable :: Uo
      real, save :: Ts 
      real, save :: Tp 
      real, save :: sigGaus
      end module waveVars
      

      module waveNumVars
      use, intrinsic :: iso_c_binding
      !include 'fftw3.f03'
      !frequency loop vars:
      integer,save      :: NFREC,NPTSTIME
      integer,dimension (:),allocatable :: vecNK
      complex*16,target :: cOME  
      real*8   ,save    :: FREC,DFREC,OME,OMEI,TW, smallestWL
      !Discrete Wave-number:
      real*8   ,save    :: DK    ! delta k
      real*8, dimension(:), allocatable,target :: k_vec
      integer,save      :: NK,NMAX
      complex*16, dimension(:), allocatable :: t_vec
      logical,save      :: trimKplease
      type(C_PTR),save  :: planNmaxF,planNmaxB,& 
                           planNfrecF,planNfrecB,& 
                           planNtimeF,planNtimeB
      end module waveNumVars
      

      module refSolMatrixVars
        complex*16, save, allocatable :: B(:,:)
        complex*16, dimension(:,:,:), allocatable,target :: Ak
        integer, dimension(:), allocatable :: IPIV
        
        complex*16, dimension(:,:,:,:), allocatable,target :: BparaGa,BparaNu
!       integer :: info    
      end module refSolMatrixVars
            

      module resultVars
       type FFres
         complex*16 :: U,V,W,Tx,Ty,Tz!,s33,s31,s11!,s32,s12
       end type FFres
      
!     use gloVars, only: dp
       type Punto2d
        real*8 :: x,z
       end type Punto2d
       
       type Punto3d
        real*8 :: x,y,z
       end type Punto3d
       
       type MecaElem
        type (Punto2d) :: center
        complex*16, dimension(5) :: Rw !u1,u3,s33,s31,s11
        complex*16, dimension(3) :: Rw_SH !u2,s32,s12
       end type MecaElem
       
!      type tipo_qlmk
!        logical :: shouldUseQuadrature
!      end type tipo_qlmk
       
       type Punto              
        type(Punto2d) :: center
        type(Punto3d) :: normal
        type(Punto2d) :: bord_A,bord_B !1x, 2y (a line on the xz plane)
        real*8  :: length
        logical :: segmentoDeEsquina
        real*8  :: gamma
        real*8  :: cosT,sinT
        integer :: layer
        logical :: isBoundary
        logical :: isSourceSegmentForce 
        logical :: isOnInterface
        logical :: guardarFK
        logical :: isSabana
        logical :: isOD
        logical :: guardarMovieSiblings
        integer :: boundaryIndex
        integer :: pointIndex
        integer :: region !1 estr, 2 incl, 0 void
        integer :: tipoFrontera !para diferenciar los puntos de colocación
        ! 0 TE^0 + TE^d = 0 
        ! 1 TE^0 + TE^d = TR^r; uE^0 + uE^d = uR^r
        ! 2 TR^r = 0

      !                       ,--- f: 1...nfrec+1
      !                       | ,--- k: 1...NMAX+1 / 2*NMAX
      !                       | | ,--- iMec: 1:3 (desps) W,U,V
      !                       | | | 
        complex*16, dimension(:,:,:), allocatable :: FK
      !                     
        complex*16, dimension(:,:,:)   , allocatable :: G
        complex*16, dimension(:,:,:,:)   , allocatable :: Gmov
      
        ! espectro campo total inquirePoints : 
      !                        ,--- f: 1...nfrec+1
      !                        | ,--- iMec: 1:2 y 3
        complex*16, dimension (:,:), allocatable :: W
        type(FFres),dimension (:), allocatable :: resp
        complex*16, dimension (:,:,:), allocatable :: Wmov
      
      !                 ,--- xXx (indice punto integracion Gaussiana)
      !                 | ,--- (1,2) -> (x,z)
      !                 | |
        real*8, dimension(:,:), allocatable :: Gq_xXx_coords
        real*8, dimension(:), allocatable :: Gq_xXx_C
      
       end type Punto   
      ! bondary elements:     ,--- POINT index / x (receptor)
      type (Punto), dimension(:), allocatable, save, target :: allpoints
      type (Punto), dimension(:), allocatable, save :: inqPoints
      type (Punto), dimension(:), allocatable, save :: moviePoints
      type (Punto), dimension(:), allocatable, save, target :: BouPoints !xi
      
        logical :: overDeterminedSystem
        integer :: OD_Jini,OD_Jend
        logical :: SabanaPlotIndividual,sabanaBajarAFrontera
        integer, save :: nIpts, nSabanapts, nMpts, nBpts, nPts,&
                       iPtini,iPtfin,mPtini,mPtfin, & 
                       n_top_sub,n_con_sub,n_val_sub,n_OD
        complex*16, dimension(:,:), allocatable :: ibemMat,copyibemmat
        complex*16, dimension(:), allocatable :: trac0vec 
        integer, dimension(:), allocatable :: IPIVbem
          complex*16, allocatable, save, target :: Sabana(:,:) !(punto,traza)
          integer :: sabZeroini,sabZerofin
        integer, allocatable, dimension(:,:), save :: fixedPoTa,pota
        integer :: nZs ! depths at pota
        complex*16, allocatable, dimension(:,:,:) :: XF
      end module resultVars
      

      module peli
        real, dimension(:),allocatable :: coords_Z ,coords_X
        complex*16, dimension(:,:,:,:), allocatable,target :: fotogramas
        integer, dimension(:,:), allocatable :: fotogramas_Region
      end module peli
      

      module sourceVars
        use resultVars, only : Punto
        type (Punto),dimension(:),allocatable, save, target :: Po
        logical, save :: SH,PSV
        integer :: tipofuente
        integer :: nFuentes
        integer :: PW_pol
      end module sourceVars
      

      module meshVars
        integer, save :: npixX,npixZ,nmarkZ,nmarkX
        real, save :: MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ,MeshVecLen
      end module meshVars
                   

      module Gquadrature
        real :: WLmulti ! una o dos longitudes de onda
        integer, parameter :: Gquad_n = 8
      ! ##############################
      ! # Para cambiar Num de puntos #
      ! # actualizar encabezado de   #
      ! # rutina punGa               #
      ! ##############################
      
      ! courtesy of http://pomax.github.io/bezierinfo/legendre-gauss.html
      
      real, parameter, dimension(8) :: Gqu_t_8 = & 
      (/ -0.9602898564975363, &
         -0.7966664774136267, &
         -0.5255324099163290, &
         -0.1834346424956498, &
          0.1834346424956498, &
          0.5255324099163290, &
          0.7966664774136267, &
          0.9602898564975363 /)
      
      real, parameter, dimension(8) :: Gqu_A_8 = &
      (/ 0.1012285362903763, &
         0.2223810344533745, &
         0.3137066458778873, &
         0.3626837833783620, &
         0.3626837833783620, &
         0.3137066458778873, &
         0.2223810344533745, &
         0.1012285362903763 /)
      
      real, parameter, dimension(14) :: Gqu_t_14 = & 
      (/ -0.9862838086968123, &
         -0.9284348836635735, &
         -0.8272013150697650, &
         -0.6872929048116855, &
         -0.5152486363581541, &
         -0.3191123689278897, &
         -0.1080549487073437, &
          0.1080549487073437, &
          0.3191123689278897, &
          0.5152486363581541, &
          0.6872929048116855, &
          0.8272013150697650, &
          0.9284348836635735, &
          0.9862838086968123 /)
      
      real, parameter, dimension(14) :: Gqu_A_14 = & 
      (/  0.0351194603317519, &
          0.0801580871597602, &
          0.1215185706879032, &
          0.1572031671581935, &
          0.1855383974779378, &
          0.2051984637212956, &
          0.2152638534631578, &
          0.2152638534631578, &
          0.2051984637212956, &
          0.1855383974779378, &
          0.1572031671581935, &
          0.1215185706879032, &
          0.0801580871597602, &
          0.0351194603317519 /)
      
      real, parameter, dimension(30) :: Gqu_t_30 = & 
      (/  -0.9968934840746495, &
       -0.9836681232797472, &
       -0.9600218649683075, &
       -0.9262000474292743, &
       -0.8825605357920527, &
       -0.8295657623827684, &
       -0.7677774321048262, &
       -0.6978504947933158, &
       -0.6205261829892429, &
       -0.5366241481420199, &
       -0.4470337695380892, &
       -0.3527047255308781, &
       -0.2546369261678899, &
       -0.1538699136085835, &
       -0.0514718425553177, &
       0.0514718425553177, &
       0.1538699136085835, &
       0.2546369261678899, &
       0.3527047255308781, &
       0.4470337695380892, &
       0.5366241481420199, &
       0.6205261829892429, &
       0.6978504947933158, &
       0.7677774321048262, &
       0.8295657623827684, &
       0.8825605357920527, &
       0.9262000474292743, &
       0.9600218649683075, &
       0.9836681232797472, &
       0.9968934840746495/)
      
      real, parameter, dimension(30) :: Gqu_A_30 = & 
      (/ 0.0079681924961666, &
       0.0184664683110910, &
       0.0287847078833234, &
       0.0387991925696271, &
       0.0484026728305941, &
       0.0574931562176191, &
       0.0659742298821805, &
       0.0737559747377052, &
       0.0807558952294202, &
       0.0868997872010830, &
       0.0921225222377861, &
       0.0963687371746443, &
       0.0995934205867953, &
       0.1017623897484055, &
       0.1028526528935588, &
       0.1028526528935588, &
       0.1017623897484055, &
       0.0995934205867953, &
       0.0963687371746443, &
       0.0921225222377861, &
       0.0868997872010830, &
       0.0807558952294202, &
       0.0737559747377052, &
       0.0659742298821805, &
       0.0574931562176191, &
       0.0484026728305941, &
       0.0387991925696271, &
       0.0287847078833234, &
       0.0184664683110910, &
       0.0079681924961666/)
      end module Gquadrature
      

      module GeometryVars
      use resultVars, only : Punto
      integer, save :: n_topo,n_cont,n_vall
      integer, save :: nXI !total number of original segment nodes 
      real*8, dimension(:,:,:), allocatable, target, save :: Xcoord_ER
      real*8, dimension(:,:,:), allocatable, target, save :: Xcoord_Voidonly
      real*8, dimension(:,:,:), allocatable, target, save :: Xcoord_Incluonly
      real*8, dimension(:,:,:), allocatable, target, save :: Xcoord_flip_out
      real*8, save :: boxIncl_maxX,boxIncl_maxY,boxIncl_minX,boxIncl_minY
      real*8, save :: boxVoid_maxX,boxVoid_maxY,boxVoid_minX,boxVoid_minY
      real,  dimension(1+1) :: surf_poly 
      real*8, allocatable, save :: midPoint(:,:)
      real*8, allocatable, save :: normXI(:,:)
      type (Punto), dimension(:), allocatable, save, target :: origGeom,origGeom_R
      Type segemntedcoords
        real*8, dimension(:), allocatable :: x,z !cantidad de nodos de la subdivisión
      end Type segemntedcoords
      logical :: staywiththefinersubdivision
      integer :: finersubdivisionJ, fraccionDeSmallestWL_segm_de_esquina
      real :: longitudcaracteristica_a
      integer :: N_de_regdionesR, N_de_segmentosR
      integer,dimension(:),allocatable :: Xcoord_Incluonly_e
      end module GeometryVars
      module wavelets
      contains 
      
      !  The Ricker wavelet on the time domain saved on   Uo
      subroutine ricker
      use gloVars, only : PI,UR
      use waveVars, only : Uo,Ts,Tp,Dt
      use waveNumVars, only : NPTSTIME
      implicit none
      integer :: i
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
        if (abs(A) .lt. 0.0001) A = 0
        Uo(i) = A * UR
      end do
      end subroutine ricker
      
      subroutine gaussian
      use waveVars, only : Uo,sigGaus
      use waveNumVars, only : NPTSTIME!,nfrec
      implicit none
      integer :: i
      real*8 :: f,s
      
      f = 0.0
      !positivos
      s = real(sigGaus/100.0 * NPTSTIME/2,8)
      do i=1, NPTSTIME/2+1
        f = real(i-1,8) ! Hz
        Uo(i) = cmplx(exp(-0.5*(f/s)**2.),0.,8)
      end do
      !negativos
      Uo(NPTSTIME/2+2:NPTSTIME) = conjg(Uo(NPTSTIME/2:2:-1))
      
      end subroutine gaussian
      
      function FFTW(n,Uin,direccion,escala)
      use waveNumVars,only:planNmaxF,planNmaxB,& 
                           planNfrecF,planNfrecB,& 
                           planNtimeF,planNtimeB,&
                           nfrec,nmax,NPTSTIME
      use, intrinsic :: iso_c_binding
      include 'fftw3.f03'
      integer, intent(in) :: n,direccion
      complex*16 :: Uin(n)
      complex*16 :: FFTW(n)
      real*8,intent(in) :: escala
      type(C_PTR) :: plan
      ! El plano creado en checarWisdom a partir del wisdom se reusa
      if (direccion .lt. 1) then ! forward -1
      if (n .eq. 2*NFREC) plan = planNfrecF
      if (n .eq. 2*NMAX) plan = planNmaxF
      if (n .eq. NPTSTIME) plan = planNtimeF
      else !backward +1
      if (n .eq. 2*NFREC) plan = planNfrecB
      if (n .eq. 2*NMAX) plan = planNmaxB
      if (n .eq. NPTSTIME) plan = planNtimeB
      end if
      ! ejecutar con plan de reuso
      if (C_ASSOCIATED (plan) .eqv. .false.) stop "fail FFTW_WISDOM_ONLY"
      call fftw_execute_dft(plan,Uin,FFTW)
      ! hacia adelante
      !Uof = Uof * Dt
      ! hacia atrás
      !Uot = Uot / (NPTSTIME*dt) ! *DFREC
      FFTW = FFTW * escala
      end function FFTW
      
      SUBROUTINE FORK(LX,CX,SIGNI,verbose,outpf)
      implicit none
      integer, intent(in) :: outpf
      integer, intent(in) :: LX,SIGNI,verbose
      COMPLEX*16 :: CARG,CW,CTEMP 
      complex*16,intent(inout) :: CX(LX)
      real*8, parameter :: pi = 4.*ATAN(1.)
      real*8 :: SC
      integer :: i,j,m,istep,l
      if (verbose >= 4) then
        write(outpf,'(a,I4,a)')'FFT on ',LX,' length vector'
      end if
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
      end module
      !

      module specfun
!     integer :: NM
!     complex*16 :: z
!     real*8 :: v,vm
!     complex*16, dimension(0:nmax) :: CBJ,CDJ,CBY,CDY
        
        private
        
        public :: CJYNA
        public :: CJYNB
        public :: CJYVB
        
        contains
        
        SUBROUTINE CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)
!
!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Rouitines called:
!            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
!            (2) MSTA1 and MSTA2 to calculate the starting 
!                point for backward recurrence
!       =======================================================
!
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 5 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
5             CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
        CBJ(0)=CBJ0
        CBJ(1)=CBJ1
        CBY(0)=CBY0
        CBY(1)=CBY1
        CDJ(0)=CDJ0
        CDJ(1)=CDJ1
        CDY(0)=CDY0
        CDY(1)=CDY1
        IF (N.LE.1) RETURN
        IF (N.LT.INT(0.25*A0)) THEN
           CJ0=CBJ0
           CJ1=CBJ1
           DO 70 K=2,N
              CJK=2.0D0*(K-1.0D0)/Z*CJ1-CJ0
              CBJ(K)=CJK
              CJ0=CJ1
70            CJ1=CJK
        ELSE
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CBJ0).GT.CDABS(CBJ1)) THEN
              CS=CBJ0/CF
           ELSE
              CS=CBJ1/CF2
           ENDIF
           DO 80 K=0,NM
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        DO 85 K=2,NM
85         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        YA0=CDABS(CBY0)
        LB=0
        CG0=CBY0
        CG1=CBY1
        DO 90 K=2,NM
           CYK=2.0D0*(K-1.0D0)/Z*CG1-CG0
           IF (CDABS(CYK).GT.1.0D+290) GO TO 90            
           YAK=CDABS(CYK)
           YA1=CDABS(CG0)
           IF (YAK.LT.YA0.AND.YAK.LT.YA1) LB=K
           CBY(K)=CYK
           CG0=CG1
           CG1=CYK
90      CONTINUE
        IF (LB.LE.4.OR.DIMAG(Z).EQ.0.0D0) GO TO 125
95      IF (LB.EQ.LB0) GO TO 125
        CH2=(1.0D0,0.0D0)
        CH1=(0.0D0,0.0D0)
        LB0=LB
        DO 100 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.NM) CBJ(LB+1)=2.0D0*LB/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CBY0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CBY0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CBY1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CBY1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,NM-1
           CYLK=2.0D0*K/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,NM
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CONTINUE
        DO 130 K=2,NM
130        CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END


        SUBROUTINE CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
!
!       ===========================================================
!       Purpose: Compute complex Bessel functions J0(z), J1(z)
!                Y0(z), Y1(z), and their derivatives
!       Input :  z --- Complex argument
!       Output:  CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
!       ===========================================================
!
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CDJ0=(0.0D0,0.0D0)
           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR/CBJ0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR/CBJ1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP/CS).LT.1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,&
                 -.5725014209747314D+00,.6074042001273483D+01,&
                 -.1100171402692467D+03,.3038090510922384D+04,&
                 -.1188384262567832D+06,.6252951493434797D+07,&
                 -.4259392165047669D+09,.3646840080706556D+11,&
                 -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,&
                  .1727727502584457D+01,-.2438052969955606D+02,&
                  .5513358961220206D+03,-.1825775547429318D+05,&
                  .8328593040162893D+06,-.5006958953198893D+08,&
                  .3836255180230433D+10,-.3649010818849833D+12,&
                  .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,&
                  .6765925884246826D+00,-.6883914268109947D+01,&
                  .1215978918765359D+03,-.3302272294480852D+04,&
                  .1276412726461746D+06,-.6656367718817688D+07,&
                  .4502786003050393D+09,-.3833857520742790D+11,&
                  .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,&
                  -.1993531733751297D+01,.2724882731126854D+02,&
                  -.6038440767050702D+03,.1971837591223663D+05,&
                  -.8902978767070678D+06,.5310411010968522D+08,&
                  -.4043620325107754D+10,.3827011346598605D+12,&
                  -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDJ0=-CBJ1
        CDJ1=CBJ0-1.0D0/Z*CBJ1
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
        RETURN
        END

        
        SUBROUTINE CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
!#< g   CJYNB   Evaluate a sequence of Bessel functions of the first and 
!       second kinds and their derivatives with integer orders and 
!       complex arguments (method 2).
!
!       =======================================================
!       Purpose: Compute Bessel functions Jn(z), Yn(z) and
!                their derivatives for a complex argument
!       Input :  z --- Complex argument of Jn(z) and Yn(z)
!                n --- Order of Jn(z) and Yn(z)
!                      ( n = 0,1,˙˙˙, n <= 250 )
!       Output:  CBJ(n) --- Jn(z)
!                CDJ(n) --- Jn'(z)
!                CBY(n) --- Yn(z)
!                CDY(n) --- Yn'(z)
!                NM --- Highest order computed
!       Routines called:
!                MSTA1 and MSTA2 to calculate the starting
!                point for backward recurrence
!       ======================================================= !#>
!#< b       Eaxmple: z = 4.0 + i 2.0
!
!     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
!    -------------------------------------------------------------------
!     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
!     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
!     2   .93050039D+00  -.77959350D+00  -.72363400D+00  -.72836666D+00 
!
!     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
!   --------------------------------------------------------------------
!     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
!     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
!     2   .79074211D+00   .86842120D+00   .78932897D+00  -.70142425D+00!#>

        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N),&
                  A(4),B(4),A1(4),B1(4)
        EL=0.5772156649015329D0
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        Y0=DABS(DIMAG(Z))
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        IF (A0.LE.300.D0.OR.N.GT.INT(0.25*A0)) THEN
           IF (N.EQ.0) NM=1
           M=MSTA1(A0,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(A0,NM,15)
           ENDIF
           CBS=(0.0D0,0.0D0)
           CSU=(0.0D0,0.0D0)
           CSV=(0.0D0,0.0D0)
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 IF (Y0.LE.1.0D0) THEN
                    CBS=CBS+2.0D0*CF
                 ELSE
                    CBS=CBS+(-1)**(K/2)*2.0D0*CF
                 ENDIF
                 CSU=CSU+(-1)**(K/2)*CF/K
              ELSE IF (K.GT.1) THEN
                 CSV=CSV+(-1)**(K/2)*K/(K*K-1.0D0)*CF
              ENDIF
              CF2=CF1
15            CF1=CF
           IF (Y0.LE.1.0D0) THEN
              CS0=CBS+CF
           ELSE
              CS0=(CBS+CF)/CDCOS(Z)
           ENDIF
           DO 20 K=0,NM
20            CBJ(K)=CBJ(K)/CS0
           CE=CDLOG(Z/2.0D0)+EL
           CBY(0)=R2P*(CE*CBJ(0)-4.0D0*CSU/CS0)
           CBY(1)=R2P*(-CBJ(0)/Z+(CE-1.0D0)*CBJ(1)-4.0D0*CSV/CS0)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,&
                  -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,&
                   .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,&
                   .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,&
                   -.1993531733751297D+01,.2724882731126854D+02/
           CT1=Z-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 25 K=1,4
25            CP0=CP0+A(K)*Z**(-2*K)
           CQ0=-0.125D0/Z
           DO 30 K=1,4
30            CQ0=CQ0+B(K)*Z**(-2*K-1)
           CU=CDSQRT(R2P/Z)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CBJ(0)=CBJ0
           CBY(0)=CBY0
           CT2=Z-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 35 K=1,4
35            CP1=CP1+A1(K)*Z**(-2*K)
           CQ1=0.375D0/Z
           DO 40 K=1,4
40            CQ1=CQ1+B1(K)*Z**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
           CBJ(1)=CBJ1
           CBY(1)=CBY1
           DO 45 K=2,NM
              CBJK=2.0D0*(K-1.0D0)/Z*CBJ1-CBJ0
              CBJ(K)=CBJK
              CBJ0=CBJ1
45            CBJ1=CBJK
        ENDIF
        CDJ(0)=-CBJ(1)
        DO 50 K=1,NM
50         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        IF (CDABS(CBJ(0)).GT.1.0D0) THEN
           CBY(1)=(CBJ(1)*CBY(0)-2.0D0/(PI*Z))/CBJ(0)
        ENDIF
        DO 55 K=2,NM
           IF (CDABS(CBJ(K-1)).GE.CDABS(CBJ(K-2))) THEN
              CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           ELSE
              CYY=(CBJ(K)*CBY(K-2)-4.0D0*(K-1.0D0)/(PI*Z*Z))/CBJ(K-2)
           ENDIF
           CBY(K)=CYY
55      CONTINUE
        CDY(0)=-CBY(1)
        DO 60 K=1,NM
60         CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END        
        
        SUBROUTINE CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
!#< g   CJYVB   Evaluate a sequence of Bessel functions of the first and 
!       second kinds and their derivatives with arbitrary real orders 
!       and complex arguments (method 2).
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
!       =========================================================== !#>
!#< b       Example:
!                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
!
!     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
!    -------------------------------------------------------------------
!     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
!     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
!     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
!     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
!     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
!     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00
!
!     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
!    -------------------------------------------------------------------
!     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
!     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
!     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
!     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
!     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
!     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
!       ============================================================= !#>

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
              M1=INT(X-1)
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
           NN=INT(N1-(N1-N0)/(1.0D0-F0/F1))              
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
           NN=INT(N1-(N1-N0)/(1.0D0-F0/F1))
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

      

      module hank
      contains
      
      SUBROUTINE HANKELS(Z,H0,H1)
!#define comparar
!#ifdef comparar
!     use specfun
!     use glovars,only:UI
!     integer :: NM,K
!     complex*16, dimension(0:10) :: CBJ,CDJ,CBY,CDY
!#endif
!     Z = COMPLEX ARGUMENT
!
!     COMPUTE SECOND KIND HANKEL FUNCTIONS H0 AND H1
!
      COMPLEX*16 :: Z,H0,H1,C,A,E,E2,ZH,P
      real*8 :: X,Y,R,PHI,J,AR
      X=REAL(Z)
      Y=AIMAG(Z)
      R=SQRT(X*X+Y*Y)
      PHI=ATAN2(Y,X)
      IF(R.LE.10.0)GO TO 20
      J=2.0*R
      C=(0.0,0.1250)/Z
      K=2
      P=C*C
      A=4.5*P
      P=7.5*P
      H0=1.0+C+A
      H1=1.0-3.0*C-P
10    I=4*K
      K=K+1
      DI=I
      DK=K
      A=A*C*(DI+1.0/DK)
      P=P*C*(DI-3.0/DK)
      H0=H0+A
      H1=H1-P
      AR=ABS(REAL(P))+ABS(AIMAG(P))
      IF(AR.GT.1.E-16.AND.K.LT.J)GO TO 10
      AR=0.785398163397448-X-PHI/2.0
      E=0.0
      IF(Y.GT.-160.0) E=0.7978845608028650/SQRT(R)*EXP(Y)*CMPLX(COS(AR),SIN(AR),8)
!     IF(X.EQ.0.0)E=CMPLX(0.0,AIMAG(E))
      IF(abs(X) .lt. 0.00001)E=CMPLX(0.0,AIMAG(E),8)
      H0=H0*E
      H1=H1*E*(0.0,1.0)
      GO TO 23
20    ZH=Z/2.0
      C=-ZH*ZH
      E=CMPLX(0.0,0.3183098861837910,8)
      E2=E*2.0
      A=1.0-E2*(0.5772156649015330+LOG(R/2.0))+PHI*0.636619772367582
      P=1.0
      K=1
      H0=A
      H1=A+E*(1.0-1.0/C)
25    A=A+E2/K
      P=P*C
      H0=H0+A*P
      K=K+1
      P=P/(K*K)
      H1=H1+(A*K+E)*P
      IF(ABS(REAL(P))+ABS(AIMAG(P)).GT.1.E-16) GO TO 25
      H1=H1*ZH
!     IF(X.NE.0.0)GO TO 23
      IF(abs(X) .gt. 0.00001) GO TO 23
      H0=CMPLX(0.0,AIMAG(H0),8)
      H1=CMPLX(REAL(H1),0.0,8)

23    K=K

!#ifdef comparar
!     print*,""
!     print*,"Z=",Z
!     print*,"H_0^(2)", H0
!     print*,"H_1^(2)", H1
!     call CJYNB(2,Z,NM,CBJ,CDJ,CBY,CDY)
!       WRITE(*,*)
!       WRITE(*,*)'   n     Re[Jn(z)]       Im[Jn(z)]',&
!                 '       Re[Yn(z)]      Im[Yn(z)]'
!       WRITE(*,*)' -------------------------------------',&
!                 '-------------------------------'
!       DO K=0,NM,1
!          WRITE(*,'(1X,I4,4D16.8)')K,CBJ(K),CBY(K)
!       end do
!       WRITE(*,*)
!       WRITE(*,*)'   n     Re[H_n^(2)(z)]       Im[H_n^(2)(z)]'
!       WRITE(*,*)' -------------------------------------',&
!                 '-------------------------------'
!       DO K=0,NM,1
!          WRITE(*,'(1X,I4,2D26.16)')K,CBJ(K)-UI*CBY(K)
!       end do
!     stop "HANK"
!#endif      
      RETURN
      END SUBROUTINE HANKELS
      end module hank
      
                 
      

      module debugStuff
      contains
      
      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      if (outpf .eq. 6) then
      do i = 1,m 
        do j = 1,n 
          write(outpf,'(A,E9.2,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E9.2,A)',advance='no') AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      else !a un archivo
      do i = 1,m
        do j = 1,n
          write(outpf,'(EN26.9,1X)',advance='no') REAL(MAT(i,j))
          write(outpf,'(EN26.9,3X)',advance='no') AIMAG(MAT(i,j))
        end do
        write(outpf,'(A)',advance='yes')''
      end do
      end if
      end subroutine
      !
      subroutine scripToMatlabMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n, outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=32), intent(in) :: name
      
      write(outpf,'(A,A)') trim(name), " = ["
      do i = 1,m
        write(outpf,'(a)',advance='no') "["
        do j = 1,n
          write(outpf,'(a,EN26.9,a)',advance='no') "(",REAL(MAT(i,j)),")+("
          write(outpf,'(EN26.9,a,2X)',advance='no') AIMAG(MAT(i,j)),")*1i"
        end do
        write(outpf,'(A)',advance='yes')'];'
      end do
      write(outpf,'(a)') "];"
      
      end subroutine scripToMatlabMNmatrixZ
      
      subroutine showMNmatrixZabs(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(E15.5,2x)',advance='no') ABS(MAT(i,j))
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      
      !
      subroutine showMNmatrixR(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n ,outpf
      real, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(F10.5,2x)',advance='no') MAT(i,j)
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine
      !
      subroutine showMNmatrixI(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n, outpf
      integer, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') ""
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(I0,3x)',advance='no') MAT(i,j)
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      write(outpf,'(A)') ""
      end subroutine
      
      end module
      
      !another set of functions in case needed
      

      module fitting
      contains
      
      function splitatY(surf_poly,degree,Y,aIN,bIN) !SINGLE PRES, REAL
      ! there is only one intersection.
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      real, intent(in) :: Y,aIN,bIN
      real, allocatable, dimension(:) :: surf_poly0
      integer, intent(in) :: degree
      integer :: i
      real :: a,b,splitatY,af,bf,temp,tmp2,tmp, & 
                 cx,d,cf,errorTol,s,sf
      logical :: mflag = .true.
      errorTol = real(0.000001,8)
      
      allocate(surf_poly0(degree+1))
      a = aIN
      b = bIN
      
      surf_poly0 = surf_poly
      surf_poly0(1) = surf_poly0(1) - Y
      
      !encontramos el cero de surf_poly0 entre a y b
      af = polyVal(surf_poly0,degree,a)
      bf = polyVal(surf_poly0,degree,b)
      ! usamos el método de Brent.
      if (af*bf >= real(0,8)) then
        if (af < bf ) then
          splitatY = af
        else
          splitatY = bf
        end if
        return
      else
        
        if (abs(af) < abs(bf)) then
          temp = b
          b = a
          a = temp
          temp = bf
          bf = af
          af = temp
        end if
        cx = a
        cf = af
        mflag = .true.
        i = 0
        d = 0.
        do while( (abs(bf) > errortol) .and. (abs(a-b) > errorTol ) )
!       do while((.not.(bf == real(0,8) )) .and. (abs(a-b) > errorTol ))
!          print*,"go_",i
          if ((abs(af-cf) > errorTol) .and. (abs(bf-cf)>errorTol)) then
!          if ((af /= cf) .and. (bf /= cf)) then
          ! interpolación cuadrática inversa
            s = a * bf * cf / (af-bf) / (af-cf) + & 
            b*af*cf/(bf-af)/(bf-cf)+cx*af*bf/(cf-af)/(cf-bf)
          else
          ! regla de la secante
            s = b - bf * (b-a)/(bf-af)
          end if
          tmp2 = (3.0*a + b)/4.0
          if ( (.not.(((s > tmp2) .and. (s < b)) .or. & 
          ((s < tmp2) .and. (s > b))) ) .or. &
          (mflag .and. ((abs(s-b)) .ge. (abs(b-cx)/2.0 ))) .or. &
          ((.not. (mflag)) .and. ((abs(s-b)) .ge. (abs(cx-d)/2.0 ))) ) then
            s = (a+b) / 2.0
            mflag = .true.
          else
            if ((mflag .and. (abs(b-cx)< errorTol)) .or. &
            ((.not. (mflag)) .and. (abs(cx-d) < errorTol))) then
              s = (a+b) / 2.0
              mflag = .true.
            else
              mflag = .false.
            end if
          end if
           sf = polyVal(surf_poly0,degree,s)
           d = cx
           cx = b
           cf = bf
!          if (af * sf < real(0,8)) then
           if (af * sf < errorTol) then
             b = s
             bf = sf
           else
             a = s
             af = sf
           end if
           if (abs(af) < abs(bf)) then
             tmp = a
             a = b
             b = tmp
             tmp = af
             af = bf
             bf = tmp
           end if
           i = i + 1
           if (i> 1000) then
             !error
             b = real(0.123456789)
             exit
           end if
        end do
      splitatY = b
      end if
      end function splitatY
      
      
      ! function evaluation
      function polyVal(surf_poly,degree,X) ! SINGLE, REAL
      implicit none
!     integer, parameter           :: dp = selected_real_kind(15, 307)
      real, dimension(:),intent(in) :: surf_poly
      integer, intent(in) :: degree 
      real, intent(in) :: X
      real :: polyVal
      integer :: i
      
      !surfo_poly are the polynomial coefficients: A0 A1 A2...An
      polyVal = surf_poly(1) !A0
      DO i = 1,degree
      polyVal = polyVal + X**i * surf_poly(i+1)
      end do
      
      end function
      
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
      ! http://rosettacode.org/wiki/Polynomial_regression#Fortran
      function polyfit(vx,vy,Ln,d,verbose,outpf) ! DOUBLE, REAL
      implicit none
      integer, intent(in)               :: verbose,outpf
      integer, intent(in)               :: Ln, d
      real*8, dimension(d+1)              :: polyfit
      real*8, dimension(:), intent(in)    :: vx, vy
      real*8, dimension(:,:), allocatable :: X
      real*8, dimension(:,:), allocatable :: XT
      real*8, dimension(:,:), allocatable :: XTX
      integer :: i, j
      integer :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      real*8, dimension(:), allocatable :: work
      
      n = d+1
      lda = n
      lwork = n
 
      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, Ln))
      allocate(X(Ln, n))
      allocate(XTX(n, n))
      
      if (verbose >= 1) then
        write(outpf,*)"fitting curve.."
      end if !
      if (verbose >= 4) then
       write(outpf,'(a)')"begin curve fit with:"
       write(outpf,'(a,I4)')"vx: of size:",size(vx)
       write (outpf, '(F9.4)') vx
       write(outpf,'(a,I4)')"vy: of size:",size(vy)
       write (outpf, '(F9.4)') vy
      end if
      
      ! prepare the matrix
      do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
      end do
      if (verbose >= 4) then
       write(outpf,'(a,I4,a,I4,a)') & 
                              "X: size, (",size(X,1),',',size(X,2),')'
       write(outpf,*) X
      end if
      
      XT  = transpose(X)
      XTX = matmul(XT, X)
      if (verbose >= 4) then
       write(outpf,'(a,I4,a,I4,a)') & 
                         "XTX: size, (",size(XTX,1),',',size(XTX,2),')'
       write(outpf,*)XTX
      end if
      
      ! calls to LAPACK subs DGETRF and DGETRI
      ! factorizacion LU
      call DGETRF(n, n, XTX, & 
                  lda, ipiv, info)
      if (verbose >= 4) then
       write(outpf,'(a)')"DGETRF (XTX):"
       write (outpf, '(E12.3)') XTX
      end if
      !
      if ( info /= 0 ) then
       write(outpf,*) "problem DGETRF =",info
       stop 1
      end if
      ! inversa de la matriz 
      call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
      if (verbose >= 4) then
       write(outpf,'(a)')"DGETRI (XTX):"
       write (outpf, '(F9.4)') XTX
      end if
      
      if ( info /= 0 ) then
       write(outpf,'(a)') "problem DGETRI =",info
       stop 1
      end if
 
      polyfit = matmul( matmul(XTX, XT), vy)
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      if (verbose >= 4) then
       write(outpf,'(a)')'polyfit='
       write (outpf, '(E12.3)') polyfit
      end if!
      if (verbose .ge. 4) then
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
       write(outpf,'(a)', ADVANCE = "NO") "0 ="
       do i=1,d+1
          write(outpf,'(ES12.3,a,I0)', ADVANCE = "NO") polyfit(i),"x^",i-1
       end do
      end if
      end function polyfit
      
      function Zpolyfit(vx,vy,Ln,d) ! complex*16, dimension(d+1) 
      use glovars, only : verbose, outpf => PrintNum
      implicit none
      real*8,     dimension(:), intent(in)    :: vx
      complex*16, dimension(:), intent(in)    :: vy
      integer, intent(in)                     :: d, Ln
      complex*16, dimension(d+1)              :: zpolyfit
      
      complex*16, dimension(:,:), allocatable :: X
      complex*16, dimension(:,:), allocatable :: XT
      complex*16, dimension(:,:), allocatable :: XTX
      integer :: i, j
      integer :: n, lda, lwork
      integer :: info
      integer, dimension(:), allocatable :: ipiv
      complex*16, dimension(:), allocatable :: work
      
      n = d+1
      lda = n
      lwork = n
 
      allocate(ipiv(n))
      allocate(work(lwork))
      allocate(XT(n, Ln))
      allocate(X(Ln, n))
      allocate(XTX(n, n))
      
      if (verbose >= 4) then
        write(outpf,*)"fitting curve.."
      end if !
      
      ! prepare the matrix
      do i = 0, d
       do j = 1, Ln
          X(j, i+1) = vx(j)**i
       end do
      end do
      XT  = transpose(X)
      XTX = matmul(XT, X)
      
      ! calls to LAPACK subs DGETRF and DGETRI
      ! factorizacion LU
      call ZGETRF(n, n, XTX, lda, ipiv, info)
      if ( info /= 0 ) then
       
       print*,Ln
       print*,vx
       write(6,*) "Zpolyfit :problem ZGETRF =",info
       stop 
      end if
      ! inversa de la matriz 
      call ZGETRI(n, XTX, lda, ipiv, work, lwork, info)
      if ( info /= 0 ) then
       write(6,'(a)') "Zpolyfit :problem ZGETRI =",info
       stop 
      end if
 
      zpolyfit = matmul( matmul(XTX, XT), vy(1:Ln))
      deallocate(ipiv)
      deallocate(work)
      deallocate(X)
      deallocate(XT)
      deallocate(XTX)
      if (verbose .ge. 4) then
      !fit a polynomial Anx^n + ... + A2x^2 + A1x + A0 = 0
      !of degree 'd' to data: (vx,vy)
      ! coefficients are orthered from the lower power: A0 A1 A2 .. AN
       write(outpf,'(a)', ADVANCE = "NO") "0 ="
       do i=1,d+1
          write(outpf,'(a,ES12.3,a,ES12.3,a,I0)', ADVANCE = "NO") & 
          'polyfit= ',real(zpolyfit(i))," i",aimag(zpolyfit(i)),"x^",i-1
       end do
      end if
      end function Zpolyfit
      
      function Zevapol(po,d,x)
      implicit none
      integer, intent(in)                     :: d
      complex*16, dimension(d+1)              :: po
      real*8 :: x
      complex*16 :: Zevapol
      integer :: i
      Zevapol = po(1)
      do i =2,d+1
        Zevapol = Zevapol + po(i)* x**(i-1)
      end do
      end function Zevapol
      
      !normal vectors to surf_poly at the nXI points (x,y) 
      function normalto(x,y,nXI,surf_poly,degree,verbose,outpf)
      implicit none
      integer,intent(in)      :: verbose,outpf
      integer, intent(in)     :: degree 
!     integer, parameter      :: dp = selected_real_kind(15, 307)
      real, intent(in), dimension(degree+1) :: surf_poly !surface 
      integer :: nXI !number of points
      real, intent(in), dimension(nXI) :: x,y !points coordinates
      
      integer :: i,j
      real*8, allocatable :: normalto(:,:) !function result value
      real*8, parameter      :: a = real(1,8) ! normal poing up parameter
!     real*8, parameter      :: tolerance = real(0.001,8) !to tell if its vertical
      real*8, dimension(degree) :: fprime !surf_poly derivative coeficients
      real*8 :: fprimeX, x0, mag 
      
      ALLOCATE (normalto(nXI,2))
      if(verbose >= 1)then
       write(outpf,'(a)',advance='no') "  ...getting the normal vectors"
      end if
      ! the normal line to curve f(x) = An x^n + ... + A1 x + A0 
      ! is y = y1 - 1 / (f'(x1)) (x - x1) 
      ! if we want the normal vectors pointing up (bigger 'y')
      ! we force   y = y1 + a   and we find x in terms of x1
      
      ! f'(x) coefficients
      do i = degree,1,-1
!      write(outpf,*)'i=',i
       fprime(i)= real(i,8) * surf_poly(i+1) 
      end do
      if (verbose >= 4 ) then
       write(outpf,'(a)')'surf_poly=A0,A1,...,AN '
       write(outpf,'(F9.4)') surf_poly
       write(outpf,'(a)')'f_prime_(x) A0 A1 ... An '
       write(outpf,'(ES9.4/)') fprime
      end if
      
      do i = 1,nXI !for every point
      !the derivative at (xi)
       fprimeX = real(0,8)
       do j = size(fprime),2,-1
         fprimeX = fprimeX + fprime(j) * x(i)**(j-1)
       end do
       fprimeX = fprimeX + fprime(1)
       x0 = x(i) - a * fprimeX
       
       !normalizar y poner como vector
       mag = sqrt((x(i)-x0)**2 + a**2)
       
       normalto(i,1)= (x0-x(i))/mag
       normalto(i,2)= a/mag
       
       if (verbose >= 4) then
         write(outpf,'(A,f7.2,A,f7.2,A,/A,f6.1,A,f6.1,/A,f6.2,A,f6.2)') &
         "(",x(i),",",y(i),")","x0= ",x0,"  y0= ",y(i)+a, &
         "nx=",normalto(i,1),"  ny=",normalto(i,2)
       end if
       
       
       
      end do
      
      if(verbose >= 1)then
       write(outpf,'(a)',advance='yes') " ... done"
      end if
      
      ! if by incrementing  y = y1 + a the value becomes infinite, then
      ! it is a vertical surface.                                     TODO
      
      
      end function normalto
      
      
!     subroutine makeTaperFuncs_cutF_cutK!(dir)
!     use resultvars, only: Hf!,Hk 
!     use waveNumVars, only : NFREC,DFREC!,NMAX,DK
!     use dislin
!     implicit none
!     integer, intent(in) :: dir
!     integer, parameter :: NpolF = 19
!     integer, parameter :: NpolK = 20
!     real, parameter :: cutoff_fracFmax = 0.92
!     real, parameter :: cutoff_fracKmax = 0.8
!     integer :: i
!     
!     if(.not. allocated(Hf)) allocate(Hf(NFREC+1))
!     if(.not. allocated(Hk)) allocate(Hk(NMAX*2))
!        
!        Hf = (/ (i,i=0,nfrec) /) * DFREC
!        Hf = 1.0/sqrt(1.0+(Hf/(cutoff_fracFmax*NFREC*DFREC))**(2*NpolF)) 
         
         
!!        Hf = 1.0_8
!        
!        Hk(1:nmax+1) = (/ (cmplx(i,0,8),i=0,nmax) /) * DK
!        Hk(nmax+2:2*nmax) = (/ (cmplx(i,0,8),i=nmax-1,1,-1) /) * DK 
!        HK = 1.0/sqrt(1.0+(Hk/(cutoff_fracKmax*NMAX*DK))**(2*NpolK))
!        
!        if (dir .eq. 2) then
!          Hk = 1.0_8/1.135_8  !SH
!        else
!!          Hk = 1.0_8/3.385_8 ! dk = 0.00343750
! !         Hk = 1.0_8/2.3_8    !PSV ! 256  Kmax = 14.07 nK = 2048 dk = 0.006875
!  !        Hk = 1.0_8/1.155_8  !PSV ! 512  Kmax = 28.15 nK = 2048 dk = 0.01375
!   !       Hk = 1.0_8/0.595_8  !PSV ! 1024 Kmax = 56.29 nK = 2048 dk = 0.0275
!          HK = HK * 1.0_8/(6100.275482093669_8 * dk * dk & 
!                      - 0292.363636363636 * dk &
!                      + 0004.021666666667)
!        end if
         
!!        CALL SETFIL('tapperK.pdf')
! !       call qplot((/(1.0*i,i=1,nmax*2)/),real(Hk,4),2*nmax)
!  !      stop
!         CALL SETFIL('tapperF.pdf')
!         call qplot((/(1.0*i,i=1,NFREC+1)/),real(Hf,4),NFREC+1)
!     !   stop
!     end subroutine makeTaperFuncs_cutF_cutK  
      
      subroutine spline3p(x,y,n01,n12,x0,y0,x1,y1,x2,y2)
      ! interpolación con spline de 3 puntos
      !    0     1         2    :  nodos
      !    1 2 3 4 5 6 7 8 9    :  vector interpolado
      !    <---> <--------->
      !     n01      n12 
      implicit none
      integer, intent(in) :: n01,n12 !puntos entre nodos,incluidos los 3 nodos
      real*8, intent(inout), dimension(n01+n12) :: x
      complex*16, intent(inout), dimension(n01+n12) :: y
      complex*16, intent(in) :: y0,y1,y2
      real*8, intent(in) :: x0,x1,x2
      integer :: i
      complex*16, dimension(3,3) :: A
      complex*16, dimension(3) :: B
      complex*16 :: a1,b1,a2,b2,t
      integer, dimension(3) :: ipiv
      complex*16, dimension(3) :: work
      integer :: info
      integer :: lwork
      A = 0
      A(1,1) = 2 / (x1-x0)
      A(1,2) = 1 / (x1-x0)
      A(2,1) = 1 / (x1-x0)
      A(2,2) = 2*(1/(x1-x0) + 1/(x2-x1))
      A(2,3) = 1 / (x2-x1)
      A(3,2) = 1 / (x2-x1)
      A(3,3) = 2 / (x2-x1)
      B(1) = 3*(y1-y0) / (x1-x0)**2
      B(2) = 3*((y1-y0)/(x1-x0)**2 + (y2-y1)/(x2-x1)**2)
      B(3) = 3*(y2-y1)/ (x2-x1)**2
      ! resolver y obtener ko, k1, k2
      lwork = 3*3
      call zgetrf(3,3,A,3,ipiv,info)
      call zgetri(3,A,3,ipiv,work,lwork,info)
      if(info .ne. 0) stop "Problem at inverse of SPLINE3P matrix "
      B = matmul(A,B)!     print*,B
      a1 = B(1)*(x1-x0) - (y1-y0)!; print*,a1
      b1 = -B(2)*(x1-x0) + (y1-y0)!; print*,b1
      a2 = B(2)*(x2-x1) - (y2-y1)!; print*,a2
      b2 = -B(3)*(x2-x1) + (y2-y1)!; print*,b2
      y(1) = y0
      do i = 2,n01
      t = (x(i)-x0) / (x1-x0)
      ! en puntos intermedios a x0 y x1
      y(i) = (1-t)*y0 + t * y1 + t * (1- t)* (a1*(1-t)+b1*t)
      end do!
      y(n01+1) = y1
      do i = n01+2,n01+n12-1
      t = (x(i)-x1) / (x2-x1)
      ! en puntos intermedios a x1 y x2
      y(i) = (1-t)*y1 + t * y2 + t * (1- t)* (a2*(1-t)+b2*t)
      end do
      y(n01+n12) = y2
      end subroutine spline3p
      
      subroutine splineIn(x,y,n01,n12,x0,y0,x1,y1,x2,y2)
      use debugstuff
      ! interpolación con spline de 3 puntos
      !    0     1         2    :  nodos
      !    1 2 3 4 5 6 7 8 9    :  vector interpolado
      !    <---> <--------->
      !     n01      n12 
      implicit none
      integer, intent(in) :: n01,n12 !puntos entre nodos,incluidos los 3 nodos
      integer, intent(inout), dimension(n01+n12) :: x
      integer, intent(inout), dimension(n01+n12) :: y
      integer, intent(in) :: y0,y1,y2
      integer, intent(in) :: x0,x1,x2
      integer :: i
      real, dimension(3,3) :: A
      real, dimension(3) :: B
      real :: a1,b1,a2,b2,t
      integer, dimension(3) :: ipiv
      real, dimension(3) :: work
      integer :: info
      integer :: lwork
!     print*,"spline3pIn"
!     print*,n01,n12
!     print*,x0,y0
!     print*,x1,y1
!     print*,x2,y2
      A(1,1) = 2. / (x1-x0)
      A(1,2) = 1. / (x1-x0)
      A(1,3) = 0.
      A(2,1) = 1. / (x1-x0)
      A(2,2) = 2.*(1./(x1-x0) + 1./(x2-x1))
      A(2,3) = 1. / (x2-x1)
      A(3,1) = 0.
      A(3,2) = 1. / (x2-x1)
      A(3,3) = 2. / (x2-x1)
      B(1) = 3.*(y1-y0) / (x1-x0)**2.
      B(2) = 3.*(1.*(y1-y0)/(x1-x0)**2. + 1.*(y2-y1)/(x2-x1)**2.)
      B(3) = 3.*(y2-y1)/ (x2-x1)**2.
      
      ! resolver y obtener ko, k1, k2
      lwork = 3*3
!     call showMNmatrixR(3,3,A,"  A  ",6)
!     call showMNmatrixR(3,1,B,"  B  ",6)
      call sgetrf(3,3,A,3,ipiv,info)
!     if(info .ne. 0) stop "Problem at LU of SPLINE3P matrix "
      call sgetri(3,A,3,ipiv,work,lwork,info)
!     if(info .ne. 0) stop "Problem at inverse of SPLINE3P matrix "
      B = matmul(A,B)!     print*,B
      a1 = B(1)*(x1-x0) - 1.*(y1-y0)!; print*,a1
      b1 = -B(2)*(x1-x0) + 1.*(y1-y0)!; print*,b1
      a2 = B(2)*(x2-x1) - 1.*(y2-y1)!; print*,a2
      b2 = -B(3)*(x2-x1) + 1.*(y2-y1)!; print*,b2
      y(1) = y0
      do i = 2,n01
      t = 1.*(x(i)-x0) / (x1-x0)
      ! en puntos intermedios a x0 y x1
      y(i) = int((1.-t)*y0 + t * 1.*y1 + t * (1.- t)* (a1*(1.-t)+1.*b1*t))
      end do!
      y(n01+1) = y1
      do i = n01+2,n01+n12-1
      t = 1.*(x(i)-x1) / (x2-x1)
      ! en puntos intermedios a x1 y x2
      y(i) = int((1.-t)*y1 + t * 1.*y2 + t * (1.- t)* (a2*(1.-t)+1.*b2*t))
      end do
      y(n01+n12) = y2
      end subroutine splineIn
      
      end module 
      module ploteo10pesos
      contains
      
      subroutine plotSpectrum(y_in,Df,full_n,n,titleN,xAx,yAx,logflag,W,H,maxfrec)
      ! (Uo,DFREC,size(Uo),size(Uo)/2.0,titleN,xAx,yAx,logflag,1200,800,maxfrec)
      use DISLIN
      implicit none
      real, intent(in)                              :: Df,maxfrec
      integer, intent(in)                           :: full_n,n,H,W
      character(LEN=100), intent(in)                 :: xAx
      character(LEN=100), intent(in)                 :: yAx
      character(LEN=9)                             :: logflag
      character(LEN=100)                            :: titleN
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
!     
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
      call color ('RED') 
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
      
      subroutine plotXYabs(y_in,Dt,n,titleN,xAx,yAx,W,H)
      ! (Uo,Dt,size(Uo),'FIGURE_NAME.pdf','time[sec]','amplitude',1200,800) 
      USE DISLIN
      implicit none
      real, intent(in)                              :: Dt
      integer, intent(in)                           :: n,H,W
      character(LEN=9), intent(in)                 :: xAx
      character(LEN=100)                            :: titleN,yAx
      COMPLEX*16, DIMENSION(n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      
      real, dimension(n) :: x
      real maxY,minY,xstep,ystep
      integer :: i
!     integer :: Modo,nx,ny
      character(LEN=100) :: dumb
      CHARACTER(LEN=30) :: CBUF
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
!     
      
!     allocate(x(n))
!     print*,size(y)
      DO i = 1,n
        x(i) = Dt*(i-1)
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
!       write(6,*) x(i), y(i)
      END DO
      
!     minY=MINVAL(real(y(:)),1)
!     write(6,*)"MinReal Val= ",minY
!     maxY=MAXVAL(real(y(:)),1)
!     write(6,*)"MaxReal Val= ",maxY
!     minYc=MINVAL(aimag(y(:)),1)
!     write(6,*)"MinComplex Val= ",minYc
!     maxYc=MAXVAL(aimag(y(:)),1)
!     write(6,*)"MaxComplex Val= ",maxYc
      minY =MINval(abs(y(:)))
      maxY =MAXval(abs(y(:)))
      
      
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
!     print*,"disini"
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL TEXMOD ('ON') ! latex!!
      CALL HWFONT()
      CALL axspos (int(370,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W,4) ,int(H,4)) !size of the axis system.
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(5,4) ,'XY') 
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
      call graf(real(x(1),4), & 
                real(x(n),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
      
!     call title 
      call color ('RED') 
      call curve(real(x,4) ,real(abs(y),4) ,int(n,4))
      call color ('FORE') 
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
      call legini(CBUF,int(1,4),int(20,4))
 !     nx = nxposn(x(n)*n + x(n)*n / 20.)
 !     ny = nyposn(minY + (maxY-minY)*0.7)
 !     print*,nx
 !     print*,ny
      call legpos(int(1600,4),int(320,4))
      write(dumb,'(a)') 'abs(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(1,4))

      call disfin()      
      
!     print*,'plotted ',trim(titleN)
      end subroutine plotXYabs
      
      subroutine plotFK(thisFK,x,z,tt,xAx,yAx,outpf, imecS ,imecE,onlythisJ,JJ) 
      use DISLIN
      use waveNumVars, only : NFREC,NMAX,DFREC,DK
      use glovars
      implicit none
      integer ,intent(in) :: outpf, imecS ,imecE
      
      !                     ,--- 0...frec
      !                     | ,--- -k...+k
      !                     | | ,--- iMec: 1:5
      !                     | | | 
      complex*16, dimension(NFREC+1,NMAX,3),intent(in) :: thisFK
      real, intent(in) :: x,z
      character(LEN=100),intent(in) :: tt
      character(LEN=9), intent(in) :: xAx,yAx
      
      real, dimension(NMAX)             :: vVert
      real, dimension(NFREC)            :: vHorz
      real, dimension(NFREC,NMAX)       :: Mab,Mre,Mim
      character(LEN=10), dimension(3)   :: nombre
      character(LEN=100)                :: titulo,t2
      integer :: i,ik,iMec!,Sf
!     real :: k
      real :: minX,minY,maxX,maxY,xstep,ystep,miV,maV,Vstep,mama!,x_i,z_i
      real, dimension(41)   :: ZLVRAY
      real, parameter :: p = 27. ! sharpness parameter
      logical :: onlythisJ
      integer :: JJ
      
      if(verbose>=2)write(outpf,'(a,a,a)') "will plot ",trim(tt),"..."
            
      Mab=0;Mre=0;Mim=0
      nombre(1)= '_w.png'
      nombre(2)= '_u.png'
      nombre(3)= '_v.png'
      
      do i=1,NFREC
         vHorz(i) = 0 + (i-1) * real(DFREC,4)  ! Hz
      end do
      !
      do ik=1,nmax
         vVert(ik) = 0 + (ik-1) * real(DK,4)
      end do
      
      minY = 0.
      maxY = vVert(size(vVert))
      ystep = maxY / 10.
      
      minX = vHorz(1)
      maxX = vHorz(size(vHorz))
      xstep = maxX / 5.
      
      do iMec = imecS,imecE   
      write(titulo,'(a,a,I0,a,I0,a,I0,a,I0,a,a)') trim(tt),'[', &
      int(x),'.',abs(int((x-int(x))*10)),';', & 
      int(z),'.',abs(int((z-int(z))*10)),']', & 
      trim(nombre(iMec))
!     print*,titulo,thisFK(1:NFREC,1:NMAX,iMec);cycle
       Mre = real(thisFK(1:NFREC,1:NMAX,iMec),4)
       Mim = real(aimag(thisFK(1:NFREC,1:NMAX,iMec)),4)
       Mab = real(log(1. + exp(p)*abs(thisFK(1:NFREC,1:NMAX,iMec))) / & 
           log(exp(p)+1.),4)
       Mab = Mab / maxval(Mab)
      
      if (onlythisJ .eqv. .true.) then
      CALL METAFL('PNG')
      call filmod('DELETE')
      CALL SETPAG('DA4L')
      write(t2,'(a,a)') "0_Abs_",trim(titulo)
      print*,trim(t2),'abs',sum(abs(thisFK(JJ,1:nmax,iMec)))*dK
      CALL SETFIL(trim(t2))
      call qplot(vVert(1:nmax),real(abs(thisFK(JJ,1:nmax,iMec)),4), nmax)
      
      write(t2,'(a,a)') "0_Rea_",trim(titulo)
      print*,trim(t2),'rea',sum(real(thisFK(JJ,1:nmax,iMec)))*dK
      CALL SETFIL(trim(t2))
      call qplot(vVert(1:nmax),real(real(thisFK(JJ,1:nmax,iMec)),4), nmax)
      
      write(t2,'(a,a)') "0_Ima_",trim(titulo)
      print*,trim(t2),'img',sum(aimag(thisFK(JJ,1:nmax,iMec)))*dK
      CALL SETFIL(trim(t2))
      call qplot(vVert(1:nmax),real(aimag(thisFK(JJ,1:nmax,iMec)),4), nmax)
      cycle !iMec
      end if
      
      ! shaded contour plot
      CALL METAFL('PNG') !'PDF'
      CALL SETFIL(trim(titulo))
      if(verbose>=2) write(outpf,'(a)')trim(titulo)
      call filmod('DELETE') ! para sobreescribir el archivo
!     CALL PAGE (2000, 1100)
!     Sf = 2
      CALL PAGE (int(5500,4) , int(2400,4))
      
      call imgfmt('RGB')
      call winsiz(int(1800,4),int(800,4))
      CALL SCRMOD('REVERS') !fondo blanco
      CALL DISINI()
      call errmod ("all", "off")
      CALL HEIGHT(int(30,4))
      ! REAL
      miV = minval(Mre)*0.1 !;print *, imec,"min",miV
      maV = maxval(Mre)*0.1 !;print *, imec,"max",maV
      if (isnan(mav)) then 
      print*,"ruined real",tt
      maV=1.0;miV=0.0
      print*, Mre(30,:)
      end if 
      
      mama = max(abs(miV),abs(maV))
      miV = -mama
      maV = mama
      Vstep = (maV-miV)/40.0
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      CALL axspos (int(360,4) ,int(2200,4)) !the position of an axis system. Lower left corner
      call axslen (int(1000,4),int(2000,4)) !size of the axis system.
      call labdig(int(2,4),'XY') !number of decimal places for labels
      call labdig(int(1,4),'Z') !number of decimal places for labels
      call labels('EXP','Z')
      call name(trim(xAx),'X') 
      call name(trim(yAx),'Y')
      CALL SETVLT ('SPEC')
     
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(minY,4),real(maxY,4),real(minY,4),real(ystep,4), & 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
                 
      CALL CONSHD(real(vHorz,4), int(size(vHorz),4), & 
                  real(vVert,4), int(size(vVert),4), & 
                  real(Mre,4), real(ZLVRAY,4),int(41,4))
      CALL ENDGRF
      ! imag
      miV = minval(Mim)*0.1!;print *, imec,"min",miV
      maV = maxval(Mim)*0.1!;print *, imec,"max",maV
      if (isnan(mav)) then 
      print*,"ruined imag",tt;maV=1.0;miV=0.0
      print*, Mim(30,:)
      end if 
      
      mama = max(abs(miV),abs(maV))
      miV = -mama
      maV = mama
      Vstep = (maV-miV)/40.0
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      CALL axspos (int(2150,4) ,int(2200,4)) !the position of an axis system. Lower left corner
      call axslen (int(1000,4),int(2000,4)) !size of the axis system.
      call labdig(int(2,4),'XY') !number of decimal places for labels
      call labdig(int(1,4),'Z') !number of decimal places for labels
      call labels('EXP','Z')
      call name(trim(xAx),'X') 
      call name(trim(" "),'Y')
      CALL SETVLT ('SPEC')
     
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(minY,4),real(maxY,4),real(minY,4),real(ystep,4), & 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
                 
      CALL CONSHD(real(vHorz,4), int(size(vHorz),4), & 
                  real(vVert,4), int(size(vVert),4), & 
                  real(Mim,4), real(ZLVRAY,4),int(41,4))
      CALL ENDGRF
      ! abs
      miV = minval(Mab) !;print *, imec,"min",miV
      maV = maxval(Mab) !;print *, imec,"max",maV
      if (isnan(mav)) then 
      print*,"ruined abs ",tt
      CALL DISFIN()
      cycle
      end if 
      Vstep = (maV-miV)/40.0
      
      DO i = 1,41
        ZLVRAY(i) = miV + Vstep * (i-1)
      end do
      CALL axspos (int(4000,4) ,int(2200,4)) !the position of an axis system. Lower left corner
      call axslen (int(1000,4),int(2000,4)) !size of the axis system.
      call labdig(int(2,4),'XY') !number of decimal places for labels
      call labdig(int(1,4),'Z') !number of decimal places for labels
      call labels('EXP','Z')
      call name(trim(xAx),'X') 
      call name(trim(" "),'Y')
      CALL SETVLT ('TEMP')
     
      call graf3(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(minY,4),real(maxY,4),real(minY,4),real(ystep,4), & 
                 real(miV,4),real(maV,4),real(miV,4),real(Vstep*4.0,4)) 
                 
      CALL CONSHD(real(vHorz,4), int(size(vHorz),4), & 
                  real(vVert,4), int(size(vVert),4), & 
                  real(Mab,4), real(ZLVRAY,4),int(41,4))
      
      CALL ENDGRF
      
      CALL DISFIN()
      
      end do !iMec
!     end do !iP
      end subroutine plotFK
      
      subroutine drawPHI(titleN,comp)
      use DISLIN
      use resultVars, only : BP => BouPoints,nbpts, trac0vec ,n_top_sub,n_con_sub, n_val_sub
      use geometryvars, only : nXI,Xcoord_ER,normXI, midPoint
      use waveNumVars, only : smallestWL
      use meshVars, only : MeshMaxX, MeshMaxZ, MeshMinX, MeshMinZ,nmarkZ,nmarkX, MeshVecLen
      implicit none
      character(LEN=100) :: titleN
      integer :: comp
      real :: maxY,minY,maxX,minX,xstep,zstep, encuadre
      integer :: i,J,tam,initial,final,Jinitial,Jfinal
      real, dimension(5,2) :: rec
      real, dimension(nBpts) :: phi
      real :: maxPhi
      integer*4 :: lentitle,sc
      character(LEN=60) :: CTIT,mr,mi
      sc=5
       
!     minX = minval(real(Xcoord_ER(:,1,:),4))
!     maxX = maxval(real(Xcoord_ER(:,1,:),4))
!     xstep = real(((maxX-minX) / 5.0 ))
!     minX = minX-xstep
!     maxX = maxX+xstep
!     
!     minY = minval(real(Xcoord_ER(:,2,:),4))
!     maxY = maxval(real(Xcoord_ER(:,2,:),4))
!     zstep = real(((maxY-minY) / 10. ))
!     minY = minY-zstep*2
!     maxY = maxY+zstep*2
!     
!     encuadre = (maxY-minY)/(maxX-minX)
      
      maxX = MeshMaxX
      minX = MeshMinX
      maxY = MeshMaxZ
      minY = MeshMinZ
      
      xstep = real(((maxX-minX) / nmarkX ))
      zstep = real(((maxY-minY) / nmarkZ ))
      encuadre = (maxY-minY)/(maxX-minX)

!     maxX = max(MeshMaxX,maxval(real(Xcoord_ER(:,1,:),4)))
!     minX = min(MeshMinX,minval(real(Xcoord_ER(:,1,:),4)))
!     maxY = max(MeshMaxZ,maxval(real(Xcoord_ER(:,2,:),4)))
!     minY = min(MeshMinZ,minval(real(Xcoord_ER(:,2,:),4)))
      
      
      maxphi = 1.0_8
      if (comp .gt. 0) then 
        tam = n_top_sub+n_con_sub
        initial = 0
        final = n_top_sub+n_con_sub
        Jinitial = 1
        Jfinal = n_top_sub+n_con_sub
      else if (comp .lt. 0) then 
        tam = n_con_sub+n_val_sub
        initial = n_top_sub+n_con_sub
        final = n_top_sub+2*n_con_sub+n_val_sub
        Jinitial = n_top_sub+1
        Jfinal = n_top_sub+n_con_sub+n_val_sub
      end if
      
      
      ! Dislin plotting routines 
      CALL METAFL('PDF') ! define intended display  XWIN,PS,EPS,PDF
      CALL SETFIL(trim(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(1200*sc,4),int(1100*sc,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize 
      call errmod ("all", "off")
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL TEXMOD ('ON') ! latex!!
      call incmrk (int(1,4))
      CALL HWFONT()
      CALL axspos (int(250*sc,4) ,int(950*sc,4)) !the position of an axis system. Lower left corner
      call axslen (int(800*sc,4) ,int(800*sc*encuadre,4)) !size of the axis system. 
      call labdig(int(1,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(10,4) ,'X') 
      call ticks (int(5,4) ,'Y') 
      CALL AXENDS ('ENDS', 'XY')
      !            low X   left Y   upp X   right Y
      call setgrf("NAME", "NAME", "NONE", "LINE")
      CALL SETVLT ('SPEC')
      ! increment for labels 
      
      ! xini xfin yfirstvalue xstep ymin ymax yminvalueTicks yTicks
      call graf(real(minX,4),real(maxX,4),real(minX,4),& 
                real(xstep,4),real(maxY,4),real(minY,4),& 
                real(maxY,4),real(-zstep,4))
      phi = 0.0
      if (comp .eq. 1) phi(1:tam) = real(real(trac0vec(2*initial+1:2*final:2)),4)
      if (comp .eq. 2) phi(1:tam) = real(real(trac0vec(2*initial+2:2*final:2)),4)
      if (comp .eq. 3) phi(1:tam) = real(real(trac0vec(initial+1:final)),4)
      if (comp .eq. -1) phi(1:tam) = real(real(trac0vec(2*initial+1:2*final:2)),4)
      if (comp .eq. -2) phi(1:tam) = real(real(trac0vec(2*initial+2:2*final:2)),4)
      if (comp .eq. -3) phi(1:tam) = real(real(trac0vec(initial:final)),4)
      maxPhi = maxval(abs(phi(1:tam)))
      write(mr,'(EN18.2)') maxPhi
      phi = phi / maxphi
      
      call color('RED')                                                    !
!     CALL SHDPAT (int(16,4)) ! filled shading                              ! 
      call mypat(45,5,3,1)
      i = 1
      do j=Jinitial,Jfinal
      rec(1,1) = real(BP(j)%bord_A%x ,4)
      rec(1,2) = real(BP(j)%bord_A%z ,4)
      
      rec(2,1) = real(BP(j)%bord_A%x + (phi(i))*BP(j)%normal%x * MeshVecLen ,4)
      rec(2,2) = real(BP(j)%bord_A%z + (phi(i))*BP(j)%normal%z * MeshVecLen ,4)
      
      rec(3,1) = real(BP(j)%bord_B%x + (phi(i))*BP(j)%normal%x * MeshVecLen ,4)
      rec(3,2) = real(BP(j)%bord_B%z + (phi(i))*BP(j)%normal%z * MeshVecLen ,4)
      
      rec(4,1) = real(BP(j)%bord_B%x ,4)
      rec(4,2) = real(BP(j)%bord_B%z ,4)
      
      rec(5,1) = real(BP(j)%bord_A%x ,4)
      rec(5,2) = real(BP(j)%bord_A%z ,4)
      i = i + 1
      CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      end do !j                                                                !
      phi = 0.0
      
      if (comp .eq. 1) phi(1:tam) = real(aimag(trac0vec(2*initial+1:2*final:2)),4)
      if (comp .eq. 2) phi(1:tam) = real(aimag(trac0vec(2*initial+2:2*final:2)),4)
      if (comp .eq. 3) phi(1:tam) = real(aimag(trac0vec(initial+1:final)),4)
      if (comp .eq. -1) phi(1:tam) = real(aimag(trac0vec(2*initial+1:2*final:2)),4)
      if (comp .eq. -2) phi(1:tam) = real(aimag(trac0vec(2*initial+2:2*final:2)),4)
      if (comp .eq. -3) phi(1:tam) = real(aimag(trac0vec(initial:final)),4)
      maxPhi = maxval(abs(phi(1:tam)))
      write(mi,'(EN18.2)') maxphi
      phi = phi / maxphi
      
      call color('BLUE')                                                    !
!     CALL SHDPAT (int(16,4)) ! filled shading                              ! 
      call mypat(45,5,3,1)
      i = 1
      do j=Jinitial,Jfinal
      rec(1,1) = real(BP(j)%bord_A%x ,4)
      rec(1,2) = real(BP(j)%bord_A%z ,4)
      
      rec(2,1) = real(BP(j)%bord_A%x + (phi(i))*BP(j)%normal%x * MeshVecLen ,4)
      rec(2,2) = real(BP(j)%bord_A%z + (phi(i))*BP(j)%normal%z * MeshVecLen ,4)
      
      rec(3,1) = real(BP(j)%bord_B%x + (phi(i))*BP(j)%normal%x * MeshVecLen ,4)
      rec(3,2) = real(BP(j)%bord_B%z + (phi(i))*BP(j)%normal%z * MeshVecLen ,4)
      
      rec(4,1) = real(BP(j)%bord_B%x ,4)
      rec(4,2) = real(BP(j)%bord_B%z ,4)
      
      rec(5,1) = real(BP(j)%bord_A%x ,4)
      rec(5,2) = real(BP(j)%bord_A%z ,4)
      i = i + 1
      CALL RLAREA(rec(:,1),rec(:,2),int(5,4))
      end do !j                                                                !
            
!     ! dibujar topografia original                                   !
!     call color ('FORE')                                             !
!     call marker(int(-1,4)) ! sin marcadores                         !
!     call PENWID(real(0.5,4))
!     call curve(real(Xcoord(:,1),4),real(Xcoord(:,2),4),int(nXI,4))  !
      
      !normales -------------------------------------------------------------
      call color('GREEN')                                                   !
      CALL HSYMBL(int(7,4)) !size of symbols                                ! 
      do j=1,nXI                                                            !
      CALL RLVEC (real(midPoint(j,1),4), real(midPoint(j,2),4), &           !
              real(midPoint(j,1)+normXI(j,1)* MeshVecLen*0.5,4), &                !
              real(midPoint(j,2)+normXI(j,2)* MeshVecLen*0.5,4), int(1001,4))     !
      end do                                                                !
      
      ! puntos centrales y gaussianos ------------------------------------
      call incmrk(-1) ! -1 : only symbols                                !
!     CALL HSYMBL(int(1,4)) !size of symbols                             !
!     call color('BACK')                                                 !
!     do j=1,nBpts                                                       !
!     !             cuadritos rellenos                                   !
!     CALL RLSYMB (16, real(BP(j)%bord_A%x,4), real(BP(j)%bord_A%z,4))   !
!     end do                                                             !
!     CALL HSYMBL(int(0,4)) !size of symbols                             !
!     do j=1,nBpts                                                       !
!     call color('FORE')                                                 !
      !             circulitos negros                                    !
!     CALL RLSYMB (17, real(BP(j)%center%x,4), real(BP(j)%center%z,4))   !
!     call color('ORANGE')                                               !
!     call marker(int(4,4)) !tachecitos                                  !
!     call curve(real(BP(j)%Gq_xXx_coords(:,1),4), &                     !
!                real(BP(j)%Gq_xXx_coords(:,2),4), &                     !
!                int(size(BP(J)%Gq_xXx_coords(:,1)),4))                  !
!     end do !j 
      
      ! longitud de onda
      call color ('FORE') 
      CALL HSYMBL(int(6,4)) !size of symbols      !
      CALL RLVEC (real((maxX+minX)/2- smallestWL/2,4), & 
          real(maxY-(maxY-minY)*0.03,4), &           !
                  real((maxX+minX)/2+ smallestWL/2,4), &     !
          real(maxY-(maxY-minY)*0.03,4), int(1122,4)) 
     
!     write(CTIT,'(a)') '$\lambda_{min}$'
!     lentitle = NLMESS(CTIT)
!     CALL HEIGHT (int(16*sc,4))
!     CALL MESSAG(CTIT,int((615*sc),4),int(900*sc,4))
      
      write(CTIT,'(a,a)') 'maxReal = ', trim(mr)
      lentitle = NLMESS(CTIT)
      CALL HEIGHT (int(12*sc,4))
      CALL MESSAG(CTIT,int((25*sc),4),int(30*sc,4))
      write(CTIT,'(a,a)') 'maxImag = ', trim(mi)
      lentitle = NLMESS(CTIT)
      CALL HEIGHT (int(12*sc,4))
      CALL MESSAG(CTIT,int((25*sc),4),int(50*sc,4))
      call disfin()
      end subroutine drawPHI
      end module ploteo10pesos

