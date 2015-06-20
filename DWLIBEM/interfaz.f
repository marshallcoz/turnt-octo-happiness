        subroutine globalmatrix_PSV(this_A,this_An,k,cOME_i)
          complex*16, intent(inout), dimension(:,:),pointer :: this_A,this_An
          real*8,     intent(in),pointer     :: k
          complex*16, intent(in),pointer     :: cOME_i 
        end subroutine globalmatrix_PSV
        
        subroutine gloMat_PSV(this_A,k,ik)
          complex*16,    intent(inout), dimension(:,:),pointer :: this_A
          real*8,     intent(in),pointer     :: k
          integer :: ik
!         complex*16, intent(in),pointer     :: cOME_i 
        end subroutine gloMat_PSV
        
        subroutine globalmatrix_SH(this_A,k,cOME_i)
          complex*16,    intent(inout), dimension(:,:),pointer :: this_A
          real*8,     intent(in),pointer     :: k
          complex*16, intent(in),pointer     :: cOME_i 
        end subroutine globalmatrix_SH
        
        subroutine inverseA(A,ipiv,work,n)
          integer, intent(in) :: n
          complex*16, dimension(:,:), intent(inout),pointer :: A
          integer, dimension(:), intent(inout),pointer :: ipiv
          complex*16, dimension(:),intent(inout),pointer :: work
        end subroutine inverseA
        
        subroutine subdivideTopo(J,FREC, onlythisJ ,minbet,bet,nbpts, BouPoints)
          use resultvars, only : Punto
          integer, intent(in) :: J
          real*8, intent(in) :: frec
          logical :: onlythisJ
          integer :: nbpts
          real*8 :: minbet
          real*8, dimension(:) :: BET
          type (Punto), dimension(:), allocatable,target ::  BouPoints
        end subroutine subdivideTopo
        
        subroutine preparePointerTable(pota,firstTime,smallestWL)
          integer, allocatable, dimension(:,:) :: pota
          logical,intent(in) :: firstTime
          real*8,intent(in) :: smallestWL
        end subroutine preparePointerTable
        
        subroutine reffField_by_(ipXi,dir_j,cOME)
          integer, intent(in) :: ipXi,dir_j 
          complex*16, intent(in),target  :: cOME
        end subroutine reffField_by_
        
        subroutine asociar(PX, i_Fuente ,itabla_z, itabla_x)
         use resultVars, only : Punto
         type(Punto), pointer :: PX
         integer, intent(in) :: itabla_x, itabla_z, i_Fuente
        end subroutine asociar
        
        subroutine PSVvectorB_force(i_zF,this_B,tam,pXi,direction,cOME,k,ik)
          use resultvars, only : Punto
          use soilVars, only : N
          integer, intent(in) :: i_zF,tam,ik
          complex*16, intent(inout), dimension(tam) :: this_B
          integer,    intent(in)    :: direction
          real*8,     intent(in)    :: k
          complex*16, intent(in)    :: cOME
          type(Punto),intent(in),target    :: pXi
        end subroutine PSVvectorB_force
        
        subroutine SHvectorB_force(i_zF,this_B,tam,pXi,cOME,k)
          use resultvars, only : Punto
          use soilVars, only : N
          integer, intent(in) :: i_zF,tam
          complex*16, intent(inout), dimension(tam) :: this_B
          real*8,     intent(in)    :: k
          complex*16, intent(in)    :: cOME
          type(Punto),intent(in),target    :: pXi
        end subroutine SHvectorB_force
        
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
        
        subroutine fill_termindep(auxK,come,mecS,mecE,p_x,pXi,dir_j)
         use resultvars, only:Punto
      type(Punto), pointer :: p_X,pXi
      integer :: mecS,mecE
      integer, intent(in) :: dir_j
      complex*16, dimension(1,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      end subroutine fill_termindep
      
      subroutine fill_ibemMat(i_zF,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto
      type(Punto), pointer :: pXi,p_X
      integer :: mecS,mecE,J,i_zF
      integer, intent(in) :: dir_j
      complex*16, dimension(1,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      end subroutine fill_ibemMat
      
      subroutine fill_diffbyStrata(i_zf,J,auxK,come,mecS,mecE,p_x,pXi,dir_j)
      use resultvars, only:Punto
      use waveNumVars, only : NMAX 
      type(Punto), pointer :: pXi,p_X
      integer :: i_zf,mecS,mecE,J
      integer, intent(in) :: dir_j
      complex*16, dimension(2*nmax,mecS:mecE), target :: auxK
      complex*16, intent(in)  :: cOME
      end subroutine fill_diffbyStrata
      

