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
        