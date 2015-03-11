C       Matriz Pre invertida para:     
C      un estrato sobre semiespacio:
C---------------------------------------

      subroutine psv1eHS(A0,cOME,k,Z,ALFA,BETA,AMU)
      IMPLICIT COMPLEX*16 (t) 
      integer N ! 1 para este caso
      PARAMETER (N=1)   
      complex*16 A0(4*N+2,4*N+2) !adjoint(A)
      complex*16 B0 ! 1/det(A)
C     variables de entrada : 
      complex*16 cOME
      real*8 k,Z(N+1)
      complex*16 ALFA(N+1),BETA(N+1),AMU(N+1)
C     auxiliares: 
      complex*16 h(N),egah(N),enuh(N)
      complex*16 ga(N+1),nu(N+1)
      complex*16  UI
      integer e
      complex*16 amue1,ga1,nu1,xi1,k2ga1,k2nu1,egah1,enuh1
      complex*16 amue2,ga2,nu2,xi2,k2ga2,k2nu2
      UI = (0.0,1.0)
      
      do e=1,N+1
        ga(e) = sqrt(cOME**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
        nu(e) = sqrt(cOME**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
        if(aimag(ga(e)).gt.0.0)ga(e) = conjg(ga(e))
        if(aimag(nu(e)).gt.0.0)nu(e) = conjg(nu(e))
      end do
C     para dentro del estrato
      do e=1,N
       h(e) = Z(e+1)-Z(e)
       egah(e) = exp(-UI * ga(e) * h(e))
       enuh(e) = exp(-UI * nu(e) * h(e))
      end do
      
      e = 1
      amue1 = amu(e)
      ga1 = ga(e)
      nu1 = nu(e)
      xi1 = k**2.0 - nu(e)**2.0 
      k2ga1 = 2.0_8*k*ga(e)
      k2nu1 = 2.0_8*k*nu(e)
      egah1 = egah(e)
      enuh1 = enuh(e)

C     ultimo medio      
      e = 2 
      amue2 = amu(e)
      ga2 = ga(e)
      nu2 = nu(e)
      xi2 = k**2.0 - nu(e)**2.0 
      k2ga2 = 2.0_8*k*ga(e)
      k2nu2 = 2.0_8*k*nu(e)
      
      t2 = amue1**2
      t3 = xi1**2
      t4 = k**2
      t5 = enuh1**2
      t6 = amue2**2
      t7 = xi2**2
      t8 = k2nu1**2
      t9 = k2ga1**2
      t10 = egah1**2
      t11 = amue2*t2*t3*t4*xi2*2.0D0
      t12 = amue2*ga1*k*k2nu2*t2*t3
      t13 = amue2*ga2*k*k2nu2*t2*t3
      t14 = amue2*k*k2ga2*nu1*t2*t3
      t15 = amue2*k*k2ga2*nu2*t2*t3
      t16 = amue1*egah1*enuh1*k2ga1*k2nu1*t2*t4*xi1*2.0D0
      t17 = amue2*ga1*k2ga2*k2nu1*nu2*t2*xi1
      t18 = amue2*ga2*k2ga1*k2nu2*nu1*t2*xi1
      t19 = amue2*ga1*k*k2nu1*t2*xi1*xi2
      t20 = amue2*ga2*k*k2nu1*t2*xi1*xi2
      t21 = amue2*k*k2ga1*nu1*t2*xi1*xi2
      t22 = amue2*k*k2ga1*nu2*t2*xi1*xi2
      t23 = amue1*egah1*enuh1*ga2*k2ga1*k2nu1*nu2*t2*xi1*2.0D0
      t24 = amue2*k*t3*t7
      t25 = amue2*k*k2ga1*k2nu1*t7
      t26 = amue2*k*k2ga2*k2nu2*t3
      t27 = amue2*k*k2ga1*k2ga2*k2nu1*k2nu2
      t28 = amue1*egah1*enuh1*k*k2ga1*k2nu1*xi1*xi2*2.0D0
      t29 = amue1*t3*t4*xi1
      t30 = amue1*ga2*nu2*t3*xi1
      t31 = amue1*k2ga1*k2nu1*t4*xi1
      t32 = amue1*ga2*k2ga1*k2nu1*nu2*xi1
      t33 = amue2*enuh1*k2ga1*k2nu1*t4*xi2
      t34 = amue2*enuh1*k2ga1*k2nu1*t4*t10*xi2
      t35 = amue1*egah1*k2ga1*t4*t8
      t36 = amue2*egah1*k2ga1*k2nu1*k2nu2*t4
      t37 = amue2*egah1*k2ga1*k2nu1*k2nu2*t4*t5
      t38 = t3**2
      t39 = amue1*t2*t3*t4*xi1
      t40 = amue1*t2*t3*t4*t5*xi1
      t41 = amue1*k2ga1*k2nu1*t2*t4*xi1
      t42 = amue1*k2ga1*k2nu1*t2*t4*t10*xi1
      t43 = amue1*k2ga1*t2*t4*t8
      t44 = amue1*k2nu1*t2*t3*t4
      t45 = amue1*k2nu1*t2*t3*t4*t5
      t46 = amue2*k2ga1*k2nu1*k2nu2*t2*t4
      t47 = amue2*k2ga1*k2nu1*k2nu2*t2*t4*t5
      t48 = amue1*enuh1*k2nu1*t4*t9
      t49 = amue2*enuh1*k2ga1*k2ga2*k2nu1*t4
      t50 = amue2*enuh1*k2ga1*k2ga2*k2nu1*t4*t10
      t51 = amue1*egah1*k2ga1*k2nu1*t4*xi1
      t52 = amue2*egah1*k2ga1*k2nu1*t4*xi2
      t53 = amue1*egah1*k2ga1*k2nu1*t4*t5*xi1
      t54 = amue2*egah1*k2ga1*k2nu1*t4*t5*xi2
      t55 = amue2*k*t3*xi1*xi2
      t56 = amue1*k*t10*t38
      t57 = amue1*k*t5*t38
      t58 = amue1*k*t8*t9*t10
      t59 = amue1*k*t5*t8*t9
      t60 = amue2*k*t3*t5*t10*xi1*xi2
      t61 = amue2*k*k2ga1*k2nu1*xi1*xi2
      t62 = amue1*egah1*enuh1*k*k2ga1*k2nu1*t3*8.0D0
      t63 = amue2*k*k2ga1*k2nu1*t10*xi1*xi2
      t64 = amue2*k*k2ga1*k2nu1*t5*xi1*xi2
      t65 = amue2*k*k2ga1*k2nu1*t5*t10*xi1*xi2
      t66 = amue1*k2nu1*t2*t4*t9
      t67 = amue1*k2ga1*t2*t3*t4
      t68 = amue1*k2ga1*t2*t3*t4*t10
      t69 = amue2*k2ga1*k2ga2*k2nu1*t2*t4
      t70 = amue2*k2ga1*k2ga2*k2nu1*t2*t4*t10
      t71 = amue1*t2*t3*t4*t5*t10*xi1
      t72 = amue2*t2*t3*t4*t10*xi2
      t73 = amue2*t2*t3*t4*t5*xi2
      t74 = amue1*k2ga1*k2nu1*t2*t4*t5*xi1
      t75 = amue1*k2ga1*k2nu1*t2*t4*t5*t10*xi1
      t76 = amue2*ga1*k2ga1*k2nu1*nu1*t2*t10*xi2
      t77 = amue2*ga1*k2ga1*k2nu1*nu1*t2*t5*xi2
      t78 = amue2*egah1*enuh1*ga1*k*k2nu1*t2*xi1*xi2*4.0D0
      t79 = amue2*egah1*enuh1*k*k2ga1*nu1*t2*xi1*xi2*4.0D0
      A0(1,1) = t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t40-
     +amue1*t2*t3*t4*xi1-amue1*t4*t6*t7*xi1-amue2*t2*t3*t4*t5*xi2*2.0D0+
     +amue1*t4*t5*t6*t7*xi1-amue1*ga2*nu2*t2*t3*xi1-amue2*ga1*nu2*t2*t3*
     +xi2-amue2*ga2*nu1*t2*t3*xi2-amue1*ga1*nu1*t6*t7*xi1-amue1*k2ga1*k2
     +nu1*t2*t4*xi1-amue2*k2ga1*k2nu2*t2*t4*xi1-amue2*k2ga2*k2nu1*t2*t4*
     +xi1-amue1*k2ga2*k2nu2*t4*t6*xi1-amue1*ga2*k2ga1*k2nu1*nu2*t2*xi1-a
     +mue1*ga1*k2ga2*k2nu2*nu1*t6*xi1-amue2*ga1*k*k2nu2*t2*t3*t5-amue2*g
     +a2*k*k2nu2*t2*t3*t5+amue2*k*k2ga2*nu1*t2*t3*t5-amue2*k*k2ga2*nu2*t
     +2*t3*t5+amue1*ga2*nu2*t2*t3*t5*xi1+amue2*ga1*nu2*t2*t3*t5*xi2-amue
     +2*ga2*nu1*t2*t3*t5*xi2-amue1*ga1*nu1*t5*t6*t7*xi1-amue1*k2ga1*k2nu
     +1*t2*t4*t5*xi1+amue2*k2ga1*k2nu2*t2*t4*t5*xi1-amue2*k2ga2*k2nu1*t2
     +*t4*t5*xi1+amue1*k2ga2*k2nu2*t4*t5*t6*xi1+amue1*egah1*enuh1*k*k2ga
     +1*nu1*t6*t7*2.0D0-amue2*egah1*enuh1*k2ga1*k2nu1*t2*t4*xi2*2.0D0-am
     +ue1*ga2*k2ga1*k2nu1*nu2*t2*t5*xi1+amue2*ga1*k2ga2*k2nu1*nu2*t2*t5*
     +xi1+amue2*ga2*k2ga1*k2nu2*nu1*t2*t5*xi1-amue1*ga1*k2ga2*k2nu2*nu1*
     +t5*t6*xi1+amue2*ga1*k*k2nu1*t2*t5*xi1*xi2+amue2*ga2*k*k2nu1*t2*t5*
     +xi1*xi2+amue2*k*k2ga1*nu1*t2*t5*xi1*xi2-amue2*k*k2ga1*nu2*t2*t5*xi
     +1*xi2-amue2*egah1*enuh1*k*k2ga1*k2ga2*k2nu1*nu2*t2*2.0D0+amue1*ega
     +h1*enuh1*k*k2ga1*k2ga2*k2nu2*nu1*t6*2.0D0-amue2*egah1*enuh1*ga2*k2
     +ga1*k2nu2*nu1*t2*xi1*2.0D0-amue2*egah1*enuh1*k*k2ga1*nu1*t2*xi1*xi
     +2*2.0D0
      A0(1,2) = t43+t44+t45+t46+t47+amue2*k2ga2*t2*t4*t8+amue1*k2nu1*t4*
     +t6*t7-amue2*k2nu1*t2*t4*xi1*xi2*2.0D0+amue1*ga2*k2ga1*nu2*t2*t8-am
     +ue2*ga1*k2ga2*nu2*t2*t8+amue1*ga2*k2nu1*nu2*t2*t3+amue1*ga1*k2nu1*
     +nu1*t6*t7+amue1*k2ga2*k2nu1*k2nu2*t4*t6-amue2*ga1*k*t2*t8*xi2-amue
     +2*ga2*k*t2*t8*xi2-amue1*k2ga1*t2*t4*t5*t8-amue2*k2ga2*t2*t4*t5*t8+
     +amue1*k2nu1*t4*t5*t6*t7-amue2*ga2*k2ga1*k2nu1*k2nu2*nu1*t2+amue1*g
     +a1*k2ga2*k2nu1*k2nu2*nu1*t6-amue1*egah1*enuh1*k2nu1*t2*t3*t4*2.0D0
     +-amue2*ga1*k*k2nu1*k2nu2*t2*xi1-amue2*ga2*k*k2nu1*k2nu2*t2*xi1-amu
     +e2*k*k2ga1*k2nu1*nu1*t2*xi2-amue2*k*k2ga2*k2nu1*nu1*t2*xi1-amue2*k
     +*k2ga1*k2nu1*nu2*t2*xi2-amue2*k*k2ga2*k2nu1*nu2*t2*xi1-amue1*ga2*k
     +2ga1*nu2*t2*t5*t8+amue2*ga1*k2ga2*nu2*t2*t5*t8+amue1*ga2*k2nu1*nu2
     +*t2*t3*t5-amue1*ga1*k2nu1*nu1*t5*t6*t7+amue1*k2ga2*k2nu1*k2nu2*t4*
     +t5*t6+amue2*ga1*k2nu1*nu2*t2*xi1*xi2+amue2*ga2*k2nu1*nu1*t2*xi1*xi
     +2+amue2*ga1*k*t2*t5*t8*xi2+amue2*ga2*k*t2*t5*t8*xi2-amue2*k2nu1*t2
     +*t4*t5*xi1*xi2*2.0D0-amue1*egah1*enuh1*ga2*k2nu1*nu2*t2*t3*2.0D0+a
     +mue2*egah1*enuh1*ga2*k2nu2*nu1*t2*t3*2.0D0+amue2*ga2*k2ga1*k2nu1*k
     +2nu2*nu1*t2*t5-amue1*ga1*k2ga2*k2nu1*k2nu2*nu1*t5*t6+amue2*egah1*e
     +nuh1*k*nu1*t2*t3*xi2*2.0D0-amue1*egah1*enuh1*k*nu1*t6*t7*xi1*2.0D0
     +-amue2*ga1*k*k2nu1*k2nu2*t2*t5*xi1-amue2*ga2*k*k2nu1*k2nu2*t2*t5*x
     +i1+amue2*egah1*enuh1*k2nu1*t2*t4*xi1*xi2*2.0D0+amue2*k*k2ga1*k2nu1
     +*nu1*t2*t5*xi2+amue2*k*k2ga2*k2nu1*nu1*t2*t5*xi1-amue2*k*k2ga1*k2n
     +u1*nu2*t2*t5*xi2-amue2*k*k2ga2*k2nu1*nu2*t2*t5*xi1+amue2*ga1*k2nu1
     +*nu2*t2*t5*xi1*xi2-amue2*ga2*k2nu1*nu1*t2*t5*xi1*xi2+amue2*egah1*e
     +nuh1*k*k2ga2*k2nu1*nu2*t2*xi1*2.0D0-amue1*egah1*enuh1*k*k2ga2*k2nu
     +2*nu1*t6*xi1*2.0D0
      A0(1,3) = amue2*t2*(amue2*egah1*nu1*t3*t7+amue1*egah1*k2ga1*k2ga2*
     +nu2*t8-amue1*egah1*k2ga2*k2nu1*nu2*t3+amue2*egah1*k2ga2*k2nu2*nu1*
     +t3-amue2*egah1*k2ga1*k2nu1*nu1*t7+amue1*enuh1*k2ga2*k2nu1*nu2*t3*2
     +.0D0+amue1*egah1*k*k2ga1*t8*xi2-amue1*egah1*k*k2nu1*t3*xi2-amue1*e
     +gah1*k*k2nu2*t3*xi1+amue1*enuh1*k*k2nu1*t3*xi2*2.0D0-amue2*enuh1*k
     +*k2nu1*t7*xi1*2.0D0+amue2*egah1*nu1*t3*t5*t7+amue1*egah1*nu2*t3*xi
     +1*xi2-amue2*egah1*k2ga1*k2ga2*k2nu1*k2nu2*nu1+amue1*egah1*k*k2ga1*
     +k2nu1*k2nu2*xi1-amue1*enuh1*k*k2ga1*k2nu1*k2nu2*xi1*2.0D0-amue2*en
     +uh1*k*k2ga2*k2nu1*k2nu2*xi1*2.0D0-amue1*egah1*k2ga1*k2ga2*nu2*t5*t
     +8-amue1*egah1*k2ga2*k2nu1*nu2*t3*t5+amue2*egah1*k2ga2*k2nu2*nu1*t3
     +*t5+amue2*egah1*k2ga1*k2nu1*nu1*t5*t7-amue1*egah1*k2ga1*k2nu1*nu2*
     +xi1*xi2+amue1*enuh1*k2ga1*k2nu1*nu2*xi1*xi2*2.0D0-amue1*egah1*k*k2
     +ga1*t5*t8*xi2-amue1*egah1*k*k2nu1*t3*t5*xi2+amue1*egah1*k*k2nu2*t3
     +*t5*xi1-amue1*egah1*nu2*t3*t5*xi1*xi2+amue2*egah1*k2ga1*k2ga2*k2nu
     +1*k2nu2*nu1*t5+amue1*egah1*k*k2ga1*k2nu1*k2nu2*t5*xi1-amue1*egah1*
     +k2ga1*k2nu1*nu2*t5*xi1*xi2)*(0.0D0,1.0D0)
      A0(1,4) = amue2*t2*(amue2*egah1*k*t3*t7-amue1*egah1*k*k2ga1*k2ga2*
     +t8+amue1*egah1*k*k2ga2*k2nu1*t3+amue2*egah1*k*k2ga2*k2nu2*t3-amue2
     +*egah1*k*k2ga1*k2nu1*t7-amue1*enuh1*k*k2ga2*k2nu1*t3*2.0D0+amue1*e
     +gah1*ga2*k2ga1*t8*xi2-amue1*egah1*ga2*k2nu1*t3*xi2-amue1*egah1*ga2
     +*k2nu2*t3*xi1+amue1*enuh1*ga2*k2nu1*t3*xi2*2.0D0+amue2*enuh1*ga1*k
     +2nu1*t7*xi1*2.0D0-amue2*egah1*k*t3*t5*t7-amue1*egah1*k*t3*xi1*xi2-
     +amue2*egah1*k*k2ga1*k2ga2*k2nu1*k2nu2+amue1*egah1*ga2*k2ga1*k2nu1*
     +k2nu2*xi1-amue1*enuh1*ga2*k2ga1*k2nu1*k2nu2*xi1*2.0D0+amue2*enuh1*
     +ga1*k2ga2*k2nu1*k2nu2*xi1*2.0D0+amue1*egah1*k*k2ga1*k2ga2*t5*t8+am
     +ue1*egah1*k*k2ga2*k2nu1*t3*t5-amue2*egah1*k*k2ga2*k2nu2*t3*t5-amue
     +2*egah1*k*k2ga1*k2nu1*t5*t7+amue1*egah1*k*k2ga1*k2nu1*xi1*xi2-amue
     +1*enuh1*k*k2ga1*k2nu1*xi1*xi2*2.0D0-amue1*egah1*ga2*k2ga1*t5*t8*xi
     +2-amue1*egah1*ga2*k2nu1*t3*t5*xi2+amue1*egah1*ga2*k2nu2*t3*t5*xi1+
     +amue1*egah1*k*t3*t5*xi1*xi2-amue2*egah1*k*k2ga1*k2ga2*k2nu1*k2nu2*
     +t5+amue1*egah1*ga2*k2ga1*k2nu1*k2nu2*t5*xi1+amue1*egah1*k*k2ga1*k2
     +nu1*t5*xi1*xi2)*(0.0D0,-1.0D0)
      A0(1,5) = -t2*(t52+t54+amue1*egah1*t3*t4*xi1-amue2*egah1*t3*t4*xi2
     +-amue2*egah1*k*k2ga2*nu1*t3-amue2*egah1*k*k2ga2*nu2*t3+amue1*egah1
     +*ga2*nu2*t3*xi1+amue2*egah1*ga2*nu1*t3*xi2-amue1*egah1*k2ga1*k2nu1
     +*t4*xi1+amue1*enuh1*k2ga1*k2nu1*t4*xi1*2.0D0+amue2*enuh1*k2ga2*k2n
     +u1*t4*xi1*2.0D0-amue1*egah1*t3*t4*t5*xi1+amue2*egah1*t3*t4*t5*xi2+
     +amue2*egah1*k*k2ga1*k2ga2*k2nu1*nu1+amue2*egah1*k*k2ga1*k2ga2*k2nu
     +1*nu2-amue1*egah1*ga2*k2ga1*k2nu1*nu2*xi1-amue2*egah1*ga2*k2ga1*k2
     +nu1*nu1*xi2+amue1*enuh1*ga2*k2ga1*k2nu1*nu2*xi1*2.0D0-amue2*enuh1*
     +ga1*k2ga2*k2nu1*nu2*xi1*2.0D0-amue2*egah1*k*k2ga2*nu1*t3*t5+amue2*
     +egah1*k*k2ga2*nu2*t3*t5-amue2*enuh1*ga1*k*k2nu1*xi1*xi2*2.0D0-amue
     +2*enuh1*ga2*k*k2nu1*xi1*xi2*2.0D0-amue1*egah1*ga2*nu2*t3*t5*xi1+am
     +ue2*egah1*ga2*nu1*t3*t5*xi2-amue1*egah1*k2ga1*k2nu1*t4*t5*xi1-amue
     +2*egah1*k*k2ga1*k2ga2*k2nu1*nu1*t5+amue2*egah1*k*k2ga1*k2ga2*k2nu1
     +*nu2*t5-amue1*egah1*ga2*k2ga1*k2nu1*nu2*t5*xi1+amue2*egah1*ga2*k2g
     +a1*k2nu1*nu1*t5*xi2)
      A0(1,6) = t2*(t35+t36+t37-amue1*egah1*k2nu1*t3*t4-amue2*egah1*k2nu
     +2*t3*t4+amue1*enuh1*k2nu1*t3*t4*2.0D0+amue1*egah1*ga2*k2ga1*nu2*t8
     +-amue1*egah1*ga2*k2nu1*nu2*t3+amue2*egah1*ga2*k2nu2*nu1*t3+amue1*e
     +nuh1*ga2*k2nu1*nu2*t3*2.0D0+amue2*egah1*k*nu1*t3*xi2+amue2*egah1*k
     +*nu2*t3*xi2-amue1*egah1*k2ga1*t4*t5*t8-amue1*egah1*k2nu1*t3*t4*t5+
     +amue2*egah1*k2nu2*t3*t4*t5-amue2*enuh1*k2nu1*t4*xi1*xi2*2.0D0-amue
     +2*egah1*ga2*k2ga1*k2nu1*k2nu2*nu1-amue2*enuh1*ga1*k*k2nu1*k2nu2*xi
     +1*2.0D0-amue2*enuh1*ga2*k*k2nu1*k2nu2*xi1*2.0D0-amue2*egah1*k*k2ga
     +1*k2nu1*nu1*xi2-amue2*egah1*k*k2ga1*k2nu1*nu2*xi2-amue1*egah1*ga2*
     +k2ga1*nu2*t5*t8-amue1*egah1*ga2*k2nu1*nu2*t3*t5+amue2*egah1*ga2*k2
     +nu2*nu1*t3*t5+amue2*enuh1*ga1*k2nu1*nu2*xi1*xi2*2.0D0+amue2*egah1*
     +k*nu1*t3*t5*xi2-amue2*egah1*k*nu2*t3*t5*xi2+amue2*egah1*ga2*k2ga1*
     +k2nu1*k2nu2*nu1*t5+amue2*egah1*k*k2ga1*k2nu1*nu1*t5*xi2-amue2*egah
     +1*k*k2ga1*k2nu1*nu2*t5*xi2)
      A0(2,1) = t66+t67+t68+t69+t70+amue1*k2ga1*t4*t6*t7+amue2*k2nu2*t2*
     +t4*t9-amue2*k2ga1*t2*t4*xi1*xi2*2.0D0+amue1*ga2*k2ga1*nu2*t2*t3+am
     +ue1*ga1*k2ga1*nu1*t6*t7+amue1*ga2*k2nu1*nu2*t2*t9-amue2*ga2*k2nu2*
     +nu1*t2*t9+amue1*k2ga1*k2ga2*k2nu2*t4*t6-amue2*k*nu1*t2*t9*xi2-amue
     +2*k*nu2*t2*t9*xi2+amue1*k2ga1*t4*t6*t7*t10-amue1*k2nu1*t2*t4*t9*t1
     +0-amue2*k2nu2*t2*t4*t9*t10-amue2*ga1*k2ga1*k2ga2*k2nu1*nu2*t2+amue
     +1*ga1*k2ga1*k2ga2*k2nu2*nu1*t6-amue1*egah1*enuh1*k2ga1*t2*t3*t4*2.
     +0D0-amue2*ga1*k*k2ga1*k2nu1*t2*xi2-amue2*ga1*k*k2ga1*k2nu2*t2*xi1-
     +amue2*ga2*k*k2ga1*k2nu1*t2*xi2-amue2*ga2*k*k2ga1*k2nu2*t2*xi1-amue
     +2*k*k2ga1*k2ga2*nu1*t2*xi1-amue2*k*k2ga1*k2ga2*nu2*t2*xi1+amue1*ga
     +2*k2ga1*nu2*t2*t3*t10-amue1*ga1*k2ga1*nu1*t6*t7*t10-amue1*ga2*k2nu
     +1*nu2*t2*t9*t10+amue2*ga2*k2nu2*nu1*t2*t9*t10+amue1*k2ga1*k2ga2*k2
     +nu2*t4*t6*t10+amue2*ga1*k2ga1*nu2*t2*xi1*xi2+amue2*ga2*k2ga1*nu1*t
     +2*xi1*xi2+amue2*k*nu1*t2*t9*t10*xi2+amue2*k*nu2*t2*t9*t10*xi2-amue
     +2*k2ga1*t2*t4*t10*xi1*xi2*2.0D0-amue1*egah1*enuh1*ga2*k2ga1*nu2*t2
     +*t3*2.0D0+amue2*egah1*enuh1*ga1*k2ga2*nu2*t2*t3*2.0D0+amue2*egah1*
     +enuh1*ga1*k*t2*t3*xi2*2.0D0-amue1*egah1*enuh1*ga1*k*t6*t7*xi1*2.0D
     +0+amue2*ga1*k2ga1*k2ga2*k2nu1*nu2*t2*t10-amue1*ga1*k2ga1*k2ga2*k2n
     +u2*nu1*t6*t10+amue2*ga1*k*k2ga1*k2nu1*t2*t10*xi2+amue2*ga1*k*k2ga1
     +*k2nu2*t2*t10*xi1-amue2*ga2*k*k2ga1*k2nu1*t2*t10*xi2-amue2*ga2*k*k
     +2ga1*k2nu2*t2*t10*xi1+amue2*egah1*enuh1*k2ga1*t2*t4*xi1*xi2*2.0D0-
     +amue2*k*k2ga1*k2ga2*nu1*t2*t10*xi1-amue2*k*k2ga1*k2ga2*nu2*t2*t10*
     +xi1-amue2*ga1*k2ga1*nu2*t2*t10*xi1*xi2+amue2*ga2*k2ga1*nu1*t2*t10*
     +xi1*xi2+amue2*egah1*enuh1*ga2*k*k2ga1*k2nu2*t2*xi1*2.0D0-amue1*ega
     +h1*enuh1*ga1*k*k2ga2*k2nu2*t6*xi1*2.0D0
      A0(2,2) = -t11-t12-t13-t14-t15-t16-t17-t18-t19-t20-t21-t22-t23+t39
     ++t41+t42+amue1*t4*t6*t7*xi1-amue1*t2*t3*t4*t10*xi1+amue2*t2*t3*t4*
     +t10*xi2*2.0D0-amue1*t4*t6*t7*t10*xi1+amue1*ga2*nu2*t2*t3*xi1+amue2
     +*ga1*nu2*t2*t3*xi2+amue2*ga2*nu1*t2*t3*xi2+amue1*ga1*nu1*t6*t7*xi1
     ++amue2*k2ga1*k2nu2*t2*t4*xi1+amue2*k2ga2*k2nu1*t2*t4*xi1+amue1*k2g
     +a2*k2nu2*t4*t6*xi1+amue1*ga2*k2ga1*k2nu1*nu2*t2*xi1+amue1*ga1*k2ga
     +2*k2nu2*nu1*t6*xi1-amue2*ga1*k*k2nu2*t2*t3*t10+amue2*ga2*k*k2nu2*t
     +2*t3*t10+amue2*k*k2ga2*nu1*t2*t3*t10+amue2*k*k2ga2*nu2*t2*t3*t10-a
     +mue1*ga2*nu2*t2*t3*t10*xi1+amue2*ga1*nu2*t2*t3*t10*xi2-amue2*ga2*n
     +u1*t2*t3*t10*xi2+amue1*ga1*nu1*t6*t7*t10*xi1+amue2*k2ga1*k2nu2*t2*
     +t4*t10*xi1-amue2*k2ga2*k2nu1*t2*t4*t10*xi1-amue1*k2ga2*k2nu2*t4*t6
     +*t10*xi1-amue1*egah1*enuh1*ga1*k*k2nu1*t6*t7*2.0D0+amue2*egah1*enu
     +h1*k2ga1*k2nu1*t2*t4*xi2*2.0D0+amue1*ga2*k2ga1*k2nu1*nu2*t2*t10*xi
     +1-amue2*ga1*k2ga2*k2nu1*nu2*t2*t10*xi1-amue2*ga2*k2ga1*k2nu2*nu1*t
     +2*t10*xi1+amue1*ga1*k2ga2*k2nu2*nu1*t6*t10*xi1-amue2*ga1*k*k2nu1*t
     +2*t10*xi1*xi2+amue2*ga2*k*k2nu1*t2*t10*xi1*xi2-amue2*k*k2ga1*nu1*t
     +2*t10*xi1*xi2-amue2*k*k2ga1*nu2*t2*t10*xi1*xi2+amue2*egah1*enuh1*g
     +a2*k*k2ga1*k2nu1*k2nu2*t2*2.0D0-amue1*egah1*enuh1*ga1*k*k2ga2*k2nu
     +1*k2nu2*t6*2.0D0+amue2*egah1*enuh1*ga1*k2ga2*k2nu1*nu2*t2*xi1*2.0D
     +0+amue2*egah1*enuh1*ga1*k*k2nu1*t2*xi1*xi2*2.0D0
      A0(2,3) = amue2*t2*(amue2*enuh1*k*t3*t7-amue1*egah1*k*k2ga1*k2nu2*
     +t3*2.0D0+amue1*enuh1*k*k2ga1*k2nu2*t3+amue2*enuh1*k*k2ga2*k2nu2*t3
     +-amue2*enuh1*k*k2ga1*k2nu1*t7-amue1*enuh1*k*k2nu1*k2nu2*t9+amue1*e
     +gah1*k2ga1*nu2*t3*xi2*2.0D0+amue2*egah1*k2ga1*nu1*t7*xi1*2.0D0-amu
     +e1*enuh1*k2ga1*nu2*t3*xi2-amue1*enuh1*k2ga2*nu2*t3*xi1+amue1*enuh1
     +*k2nu1*nu2*t9*xi2-amue2*enuh1*k*t3*t7*t10-amue1*enuh1*k*t3*xi1*xi2
     +-amue2*enuh1*k*k2ga1*k2ga2*k2nu1*k2nu2-amue1*egah1*k2ga1*k2ga2*k2n
     +u1*nu2*xi1*2.0D0+amue2*egah1*k2ga1*k2ga2*k2nu2*nu1*xi1*2.0D0+amue1
     +*enuh1*k2ga1*k2ga2*k2nu1*nu2*xi1+amue1*enuh1*k*k2ga1*k2nu2*t3*t10-
     +amue2*enuh1*k*k2ga2*k2nu2*t3*t10-amue2*enuh1*k*k2ga1*k2nu1*t7*t10+
     +amue1*enuh1*k*k2nu1*k2nu2*t9*t10-amue1*egah1*k*k2ga1*k2nu1*xi1*xi2
     +*2.0D0+amue1*enuh1*k*k2ga1*k2nu1*xi1*xi2-amue1*enuh1*k2ga1*nu2*t3*
     +t10*xi2+amue1*enuh1*k2ga2*nu2*t3*t10*xi1-amue1*enuh1*k2nu1*nu2*t9*
     +t10*xi2+amue1*enuh1*k*t3*t10*xi1*xi2-amue2*enuh1*k*k2ga1*k2ga2*k2n
     +u1*k2nu2*t10+amue1*enuh1*k2ga1*k2ga2*k2nu1*nu2*t10*xi1+amue1*enuh1
     +*k*k2ga1*k2nu1*t10*xi1*xi2)*(0.0D0,-1.0D0)
      A0(2,4) = amue2*t2*(amue2*enuh1*ga1*t3*t7+amue1*egah1*ga2*k2ga1*k2
     +nu2*t3*2.0D0-amue1*enuh1*ga2*k2ga1*k2nu2*t3+amue2*enuh1*ga1*k2ga2*
     +k2nu2*t3-amue2*enuh1*ga1*k2ga1*k2nu1*t7+amue1*enuh1*ga2*k2nu1*k2nu
     +2*t9+amue1*egah1*k*k2ga1*t3*xi2*2.0D0-amue2*egah1*k*k2ga1*t7*xi1*2
     +.0D0-amue1*enuh1*k*k2ga1*t3*xi2-amue1*enuh1*k*k2ga2*t3*xi1+amue1*e
     +nuh1*k*k2nu1*t9*xi2+amue2*enuh1*ga1*t3*t7*t10+amue1*enuh1*ga2*t3*x
     +i1*xi2-amue2*enuh1*ga1*k2ga1*k2ga2*k2nu1*k2nu2-amue1*egah1*k*k2ga1
     +*k2ga2*k2nu1*xi1*2.0D0-amue2*egah1*k*k2ga1*k2ga2*k2nu2*xi1*2.0D0+a
     +mue1*enuh1*k*k2ga1*k2ga2*k2nu1*xi1-amue1*enuh1*ga2*k2ga1*k2nu2*t3*
     +t10+amue2*enuh1*ga1*k2ga2*k2nu2*t3*t10+amue2*enuh1*ga1*k2ga1*k2nu1
     +*t7*t10-amue1*enuh1*ga2*k2nu1*k2nu2*t9*t10+amue1*egah1*ga2*k2ga1*k
     +2nu1*xi1*xi2*2.0D0-amue1*enuh1*ga2*k2ga1*k2nu1*xi1*xi2-amue1*enuh1
     +*k*k2ga1*t3*t10*xi2+amue1*enuh1*k*k2ga2*t3*t10*xi1-amue1*enuh1*k*k
     +2nu1*t9*t10*xi2-amue1*enuh1*ga2*t3*t10*xi1*xi2+amue2*enuh1*ga1*k2g
     +a1*k2ga2*k2nu1*k2nu2*t10+amue1*enuh1*k*k2ga1*k2ga2*k2nu1*t10*xi1-a
     +mue1*enuh1*ga2*k2ga1*k2nu1*t10*xi1*xi2)*(0.0D0,-1.0D0)
      A0(2,5) = t2*(t48+t49+t50+amue1*egah1*k2ga1*t3*t4*2.0D0-amue1*enuh
     +1*k2ga1*t3*t4-amue2*enuh1*k2ga2*t3*t4+amue1*egah1*ga2*k2ga1*nu2*t3
     +*2.0D0-amue1*enuh1*ga2*k2ga1*nu2*t3+amue2*enuh1*ga1*k2ga2*nu2*t3+a
     +mue1*enuh1*ga2*k2nu1*nu2*t9+amue2*enuh1*ga1*k*t3*xi2+amue2*enuh1*g
     +a2*k*t3*xi2-amue1*enuh1*k2ga1*t3*t4*t10+amue2*enuh1*k2ga2*t3*t4*t1
     +0-amue1*enuh1*k2nu1*t4*t9*t10-amue2*egah1*k2ga1*t4*xi1*xi2*2.0D0-a
     +mue2*enuh1*ga1*k2ga1*k2ga2*k2nu1*nu2-amue2*enuh1*ga1*k*k2ga1*k2nu1
     +*xi2-amue2*enuh1*ga2*k*k2ga1*k2nu1*xi2-amue2*egah1*k*k2ga1*k2ga2*n
     +u1*xi1*2.0D0-amue2*egah1*k*k2ga1*k2ga2*nu2*xi1*2.0D0-amue1*enuh1*g
     +a2*k2ga1*nu2*t3*t10+amue2*enuh1*ga1*k2ga2*nu2*t3*t10-amue1*enuh1*g
     +a2*k2nu1*nu2*t9*t10+amue2*egah1*ga2*k2ga1*nu1*xi1*xi2*2.0D0+amue2*
     +enuh1*ga1*k*t3*t10*xi2-amue2*enuh1*ga2*k*t3*t10*xi2+amue2*enuh1*ga
     +1*k2ga1*k2ga2*k2nu1*nu2*t10+amue2*enuh1*ga1*k*k2ga1*k2nu1*t10*xi2-
     +amue2*enuh1*ga2*k*k2ga1*k2nu1*t10*xi2)
      A0(2,6) = t2*(t33+t34+amue1*enuh1*t3*t4*xi1-amue2*enuh1*t3*t4*xi2-
     +amue2*enuh1*ga1*k*k2nu2*t3-amue2*enuh1*ga2*k*k2nu2*t3+amue1*enuh1*
     +ga2*nu2*t3*xi1+amue2*enuh1*ga1*nu2*t3*xi2+amue1*egah1*k2ga1*k2nu1*
     +t4*xi1*2.0D0+amue2*egah1*k2ga1*k2nu2*t4*xi1*2.0D0-amue1*enuh1*k2ga
     +1*k2nu1*t4*xi1-amue1*enuh1*t3*t4*t10*xi1+amue2*enuh1*t3*t4*t10*xi2
     ++amue2*enuh1*ga1*k*k2ga1*k2nu1*k2nu2+amue2*enuh1*ga2*k*k2ga1*k2nu1
     +*k2nu2+amue1*egah1*ga2*k2ga1*k2nu1*nu2*xi1*2.0D0-amue2*egah1*ga2*k
     +2ga1*k2nu2*nu1*xi1*2.0D0-amue1*enuh1*ga2*k2ga1*k2nu1*nu2*xi1-amue2
     +*enuh1*ga1*k2ga1*k2nu1*nu2*xi2-amue2*enuh1*ga1*k*k2nu2*t3*t10+amue
     +2*enuh1*ga2*k*k2nu2*t3*t10-amue2*egah1*k*k2ga1*nu1*xi1*xi2*2.0D0-a
     +mue2*egah1*k*k2ga1*nu2*xi1*xi2*2.0D0-amue1*enuh1*ga2*nu2*t3*t10*xi
     +1+amue2*enuh1*ga1*nu2*t3*t10*xi2-amue1*enuh1*k2ga1*k2nu1*t4*t10*xi
     +1-amue2*enuh1*ga1*k*k2ga1*k2nu1*k2nu2*t10+amue2*enuh1*ga2*k*k2ga1*
     +k2nu1*k2nu2*t10-amue1*enuh1*ga2*k2ga1*k2nu1*nu2*t10*xi1+amue2*enuh
     +1*ga1*k2ga1*k2nu1*nu2*t10*xi2)
      A0(3,1) = amue1*egah1*t2*t3*t4*xi1-amue2*egah1*t2*t3*t4*xi2*2.0D0+
     +amue1*egah1*t4*t6*t7*xi1+amue2*egah1*ga1*k*k2nu2*t2*t3-amue2*egah1
     +*ga2*k*k2nu2*t2*t3-amue2*egah1*k*k2ga2*nu1*t2*t3-amue2*egah1*k*k2g
     +a2*nu2*t2*t3+amue1*enuh1*k*k2ga1*nu1*t6*t7*2.0D0+amue1*egah1*ga2*n
     +u2*t2*t3*xi1-amue2*egah1*ga1*nu2*t2*t3*xi2+amue2*egah1*ga2*nu1*t2*
     +t3*xi2-amue1*egah1*ga1*nu1*t6*t7*xi1-amue1*egah1*k2ga1*k2nu1*t2*t4
     +*xi1-amue2*egah1*k2ga1*k2nu2*t2*t4*xi1+amue2*egah1*k2ga2*k2nu1*t2*
     +t4*xi1+amue1*egah1*k2ga2*k2nu2*t4*t6*xi1+amue1*enuh1*k2ga1*k2nu1*t
     +2*t4*xi1*2.0D0-amue2*enuh1*k2ga1*k2nu1*t2*t4*xi2*2.0D0-amue1*egah1
     +*t2*t3*t4*t5*xi1+amue2*egah1*t2*t3*t4*t5*xi2*2.0D0-amue1*egah1*t4*
     +t5*t6*t7*xi1-amue2*enuh1*k*k2ga1*k2ga2*k2nu1*nu2*t2*2.0D0+amue1*en
     +uh1*k*k2ga1*k2ga2*k2nu2*nu1*t6*2.0D0-amue1*egah1*ga2*k2ga1*k2nu1*n
     +u2*t2*xi1+amue2*egah1*ga1*k2ga2*k2nu1*nu2*t2*xi1+amue2*egah1*ga2*k
     +2ga1*k2nu2*nu1*t2*xi1-amue1*egah1*ga1*k2ga2*k2nu2*nu1*t6*xi1+amue1
     +*enuh1*ga2*k2ga1*k2nu1*nu2*t2*xi1*2.0D0-amue2*enuh1*ga2*k2ga1*k2nu
     +2*nu1*t2*xi1*2.0D0-amue2*egah1*ga1*k*k2nu2*t2*t3*t5+amue2*egah1*ga
     +2*k*k2nu2*t2*t3*t5-amue2*egah1*k*k2ga2*nu1*t2*t3*t5+amue2*egah1*k*
     +k2ga2*nu2*t2*t3*t5+amue2*egah1*ga1*k*k2nu1*t2*xi1*xi2-amue2*egah1*
     +ga2*k*k2nu1*t2*xi1*xi2+amue2*egah1*k*k2ga1*nu1*t2*xi1*xi2+amue2*eg
     +ah1*k*k2ga1*nu2*t2*xi1*xi2-amue2*enuh1*k*k2ga1*nu1*t2*xi1*xi2*2.0D
     +0-amue1*egah1*ga2*nu2*t2*t3*t5*xi1+amue2*egah1*ga1*nu2*t2*t3*t5*xi
     +2+amue2*egah1*ga2*nu1*t2*t3*t5*xi2-amue1*egah1*ga1*nu1*t5*t6*t7*xi
     +1-amue1*egah1*k2ga1*k2nu1*t2*t4*t5*xi1+amue2*egah1*k2ga1*k2nu2*t2*
     +t4*t5*xi1+amue2*egah1*k2ga2*k2nu1*t2*t4*t5*xi1-amue1*egah1*k2ga2*k
     +2nu2*t4*t5*t6*xi1-amue1*egah1*ga2*k2ga1*k2nu1*nu2*t2*t5*xi1+amue2*
     +egah1*ga1*k2ga2*k2nu1*nu2*t2*t5*xi1+amue2*egah1*ga2*k2ga1*k2nu2*nu
     +1*t2*t5*xi1-amue1*egah1*ga1*k2ga2*k2nu2*nu1*t5*t6*xi1+amue2*egah1*
     +ga1*k*k2nu1*t2*t5*xi1*xi2-amue2*egah1*ga2*k*k2nu1*t2*t5*xi1*xi2+am
     +ue2*egah1*k*k2ga1*nu1*t2*t5*xi1*xi2-amue2*egah1*k*k2ga1*nu2*t2*t5*
     +xi1*xi2
      A0(3,2) = amue1*egah1*k2ga1*t2*t4*t8-amue2*egah1*k2ga2*t2*t4*t8-am
     +ue1*egah1*k2nu1*t2*t3*t4-amue1*egah1*k2nu1*t4*t6*t7+amue1*enuh1*k2
     +nu1*t2*t3*t4*2.0D0+amue1*egah1*ga2*k2ga1*nu2*t2*t8-amue2*egah1*ga1
     +*k2ga2*nu2*t2*t8-amue1*egah1*ga2*k2nu1*nu2*t2*t3+amue1*egah1*ga1*k
     +2nu1*nu1*t6*t7+amue1*enuh1*ga2*k2nu1*nu2*t2*t3*2.0D0-amue2*enuh1*g
     +a2*k2nu2*nu1*t2*t3*2.0D0+amue2*egah1*k2ga1*k2nu1*k2nu2*t2*t4-amue1
     +*egah1*k2ga2*k2nu1*k2nu2*t4*t6-amue2*egah1*ga1*k*t2*t8*xi2+amue2*e
     +gah1*ga2*k*t2*t8*xi2-amue2*enuh1*k*nu1*t2*t3*xi2*2.0D0+amue1*enuh1
     +*k*nu1*t6*t7*xi1*2.0D0-amue1*egah1*k2ga1*t2*t4*t5*t8+amue2*egah1*k
     +2ga2*t2*t4*t5*t8-amue1*egah1*k2nu1*t2*t3*t4*t5-amue1*egah1*k2nu1*t
     +4*t5*t6*t7+amue2*egah1*k2nu1*t2*t4*xi1*xi2*2.0D0-amue2*enuh1*k2nu1
     +*t2*t4*xi1*xi2*2.0D0-amue2*egah1*ga2*k2ga1*k2nu1*k2nu2*nu1*t2+amue
     +1*egah1*ga1*k2ga2*k2nu1*k2nu2*nu1*t6-amue2*egah1*ga1*k*k2nu1*k2nu2
     +*t2*xi1+amue2*egah1*ga2*k*k2nu1*k2nu2*t2*xi1-amue2*egah1*k*k2ga1*k
     +2nu1*nu1*t2*xi2+amue2*egah1*k*k2ga2*k2nu1*nu1*t2*xi1-amue2*egah1*k
     +*k2ga1*k2nu1*nu2*t2*xi2+amue2*egah1*k*k2ga2*k2nu1*nu2*t2*xi1-amue2
     +*enuh1*k*k2ga2*k2nu1*nu2*t2*xi1*2.0D0+amue1*enuh1*k*k2ga2*k2nu2*nu
     +1*t6*xi1*2.0D0-amue1*egah1*ga2*k2ga1*nu2*t2*t5*t8+amue2*egah1*ga1*
     +k2ga2*nu2*t2*t5*t8-amue1*egah1*ga2*k2nu1*nu2*t2*t3*t5-amue1*egah1*
     +ga1*k2nu1*nu1*t5*t6*t7+amue2*egah1*k2ga1*k2nu1*k2nu2*t2*t4*t5-amue
     +1*egah1*k2ga2*k2nu1*k2nu2*t4*t5*t6+amue2*egah1*ga1*k2nu1*nu2*t2*xi
     +1*xi2-amue2*egah1*ga2*k2nu1*nu1*t2*xi1*xi2+amue2*egah1*ga1*k*t2*t5
     +*t8*xi2-amue2*egah1*ga2*k*t2*t5*t8*xi2+amue2*egah1*k2nu1*t2*t4*t5*
     +xi1*xi2*2.0D0+amue2*egah1*ga2*k2ga1*k2nu1*k2nu2*nu1*t2*t5-amue1*eg
     +ah1*ga1*k2ga2*k2nu1*k2nu2*nu1*t5*t6-amue2*egah1*ga1*k*k2nu1*k2nu2*
     +t2*t5*xi1+amue2*egah1*ga2*k*k2nu1*k2nu2*t2*t5*xi1+amue2*egah1*k*k2
     +ga1*k2nu1*nu1*t2*t5*xi2-amue2*egah1*k*k2ga2*k2nu1*nu1*t2*t5*xi1-am
     +ue2*egah1*k*k2ga1*k2nu1*nu2*t2*t5*xi2+amue2*egah1*k*k2ga2*k2nu1*nu
     +2*t2*t5*xi1+amue2*egah1*ga1*k2nu1*nu2*t2*t5*xi1*xi2+amue2*egah1*ga
     +2*k2nu1*nu1*t2*t5*xi1*xi2
      A0(3,3) = amue2*t2*(amue2*nu1*t3*t7-amue1*k2ga1*k2ga2*nu2*t8-amue1
     +*k2ga2*k2nu1*nu2*t3+amue2*k2ga2*k2nu2*nu1*t3+amue2*k2ga1*k2nu1*nu1
     +*t7-amue1*k*k2ga1*t8*xi2-amue1*k*k2nu1*t3*xi2-amue1*k*k2nu2*t3*xi1
     ++amue2*nu1*t3*t5*t7+amue1*nu2*t3*xi1*xi2-amue1*nu2*t3*t5*xi1*xi2+a
     +mue2*k2ga1*k2ga2*k2nu1*k2nu2*nu1-amue1*k*k2ga1*k2nu1*k2nu2*xi1+amu
     +e1*k2ga1*k2ga2*nu2*t5*t8-amue1*k2ga2*k2nu1*nu2*t3*t5+amue2*k2ga2*k
     +2nu2*nu1*t3*t5-amue2*k2ga1*k2nu1*nu1*t5*t7+amue1*k2ga1*k2nu1*nu2*x
     +i1*xi2+amue1*k*k2ga1*t5*t8*xi2-amue1*k*k2nu1*t3*t5*xi2+amue1*k*k2n
     +u2*t3*t5*xi1+amue1*egah1*enuh1*k2ga2*k2nu1*nu2*t3*2.0D0+amue1*egah
     +1*enuh1*k*k2nu1*t3*xi2*2.0D0-amue2*egah1*enuh1*k*k2nu1*t7*xi1*2.0D
     +0-amue2*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t5-amue1*k*k2ga1*k2nu1*k2nu2*t
     +5*xi1+amue1*k2ga1*k2nu1*nu2*t5*xi1*xi2+amue1*egah1*enuh1*k*k2ga1*k
     +2nu1*k2nu2*xi1*2.0D0-amue2*egah1*enuh1*k*k2ga2*k2nu1*k2nu2*xi1*2.0
     +D0-amue1*egah1*enuh1*k2ga1*k2nu1*nu2*xi1*xi2*2.0D0)*(0.0D0,-1.0D0)
      A0(3,4) = amue2*t2*(t24+t25+t26+t27+t28+amue1*k*k2ga1*k2ga2*t8+amu
     +e1*k*k2ga2*k2nu1*t3-amue1*ga2*k2ga1*t8*xi2-amue1*ga2*k2nu1*t3*xi2-
     +amue1*ga2*k2nu2*t3*xi1-amue2*k*t3*t5*t7-amue1*k*t3*xi1*xi2+amue1*k
     +*t3*t5*xi1*xi2-amue1*ga2*k2ga1*k2nu1*k2nu2*xi1-amue1*k*k2ga1*k2ga2
     +*t5*t8+amue1*k*k2ga2*k2nu1*t3*t5-amue2*k*k2ga2*k2nu2*t3*t5+amue2*k
     +*k2ga1*k2nu1*t5*t7-amue1*k*k2ga1*k2nu1*xi1*xi2+amue1*ga2*k2ga1*t5*
     +t8*xi2-amue1*ga2*k2nu1*t3*t5*xi2+amue1*ga2*k2nu2*t3*t5*xi1-amue1*e
     +gah1*enuh1*k*k2ga2*k2nu1*t3*2.0D0+amue1*egah1*enuh1*ga2*k2nu1*t3*x
     +i2*2.0D0-amue2*egah1*enuh1*ga1*k2nu1*t7*xi1*2.0D0+amue2*k*k2ga1*k2
     +ga2*k2nu1*k2nu2*t5-amue1*ga2*k2ga1*k2nu1*k2nu2*t5*xi1-amue1*k*k2ga
     +1*k2nu1*t5*xi1*xi2+amue1*egah1*enuh1*ga2*k2ga1*k2nu1*k2nu2*xi1*2.0
     +D0-amue2*egah1*enuh1*ga1*k2ga2*k2nu1*k2nu2*xi1*2.0D0)*(0.0D0,1.0D0
     +)
      A0(3,5) = t2*(t29+t30+t31+t32-amue2*t3*t4*xi2-amue2*k*k2ga2*nu1*t3
     +-amue2*k*k2ga2*nu2*t3+amue2*ga2*nu1*t3*xi2-amue2*k2ga1*k2nu1*t4*xi
     +2-amue1*t3*t4*t5*xi1+amue2*t3*t4*t5*xi2-amue2*k*k2ga1*k2ga2*k2nu1*
     +nu1-amue2*k*k2ga1*k2ga2*k2nu1*nu2+amue2*ga2*k2ga1*k2nu1*nu1*xi2-am
     +ue2*k*k2ga2*nu1*t3*t5+amue2*k*k2ga2*nu2*t3*t5-amue1*ga2*nu2*t3*t5*
     +xi1+amue2*ga2*nu1*t3*t5*xi2+amue1*k2ga1*k2nu1*t4*t5*xi1-amue2*k2ga
     +1*k2nu1*t4*t5*xi2-amue1*egah1*enuh1*k2ga1*k2nu1*t4*xi1*2.0D0+amue2
     +*egah1*enuh1*k2ga2*k2nu1*t4*xi1*2.0D0+amue2*k*k2ga1*k2ga2*k2nu1*nu
     +1*t5-amue2*k*k2ga1*k2ga2*k2nu1*nu2*t5+amue1*ga2*k2ga1*k2nu1*nu2*t5
     +*xi1-amue2*ga2*k2ga1*k2nu1*nu1*t5*xi2-amue1*egah1*enuh1*ga2*k2ga1*
     +k2nu1*nu2*xi1*2.0D0+amue2*egah1*enuh1*ga1*k2ga2*k2nu1*nu2*xi1*2.0D
     +0+amue2*egah1*enuh1*ga1*k*k2nu1*xi1*xi2*2.0D0-amue2*egah1*enuh1*ga
     +2*k*k2nu1*xi1*xi2*2.0D0)
      A0(3,6) = t2*(amue1*k2ga1*t4*t8+amue1*k2nu1*t3*t4+amue2*k2nu2*t3*t
     +4+amue1*ga2*k2ga1*nu2*t8+amue1*ga2*k2nu1*nu2*t3-amue2*ga2*k2nu2*nu
     +1*t3+amue2*k2ga1*k2nu1*k2nu2*t4-amue2*k*nu1*t3*xi2-amue2*k*nu2*t3*
     +xi2-amue1*k2ga1*t4*t5*t8+amue1*k2nu1*t3*t4*t5-amue2*k2nu2*t3*t4*t5
     +-amue2*ga2*k2ga1*k2nu1*k2nu2*nu1-amue1*egah1*enuh1*k2nu1*t3*t4*2.0
     +D0-amue2*k*k2ga1*k2nu1*nu1*xi2-amue2*k*k2ga1*k2nu1*nu2*xi2-amue1*g
     +a2*k2ga1*nu2*t5*t8+amue1*ga2*k2nu1*nu2*t3*t5-amue2*ga2*k2nu2*nu1*t
     +3*t5+amue2*k2ga1*k2nu1*k2nu2*t4*t5-amue2*k*nu1*t3*t5*xi2+amue2*k*n
     +u2*t3*t5*xi2-amue1*egah1*enuh1*ga2*k2nu1*nu2*t3*2.0D0+amue2*ga2*k2
     +ga1*k2nu1*k2nu2*nu1*t5+amue2*egah1*enuh1*k2nu1*t4*xi1*xi2*2.0D0+am
     +ue2*k*k2ga1*k2nu1*nu1*t5*xi2-amue2*k*k2ga1*k2nu1*nu2*t5*xi2-amue2*
     +egah1*enuh1*ga1*k*k2nu1*k2nu2*xi1*2.0D0+amue2*egah1*enuh1*ga2*k*k2
     +nu1*k2nu2*xi1*2.0D0+amue2*egah1*enuh1*ga1*k2nu1*nu2*xi1*xi2*2.0D0)
      A0(4,1) = amue1*egah1*k2ga1*t2*t3*t4*2.0D0-amue1*enuh1*k2ga1*t2*t3
     +*t4-amue1*enuh1*k2ga1*t4*t6*t7+amue1*enuh1*k2nu1*t2*t4*t9-amue2*en
     +uh1*k2nu2*t2*t4*t9+amue1*egah1*ga2*k2ga1*nu2*t2*t3*2.0D0-amue2*ega
     +h1*ga1*k2ga2*nu2*t2*t3*2.0D0-amue1*enuh1*ga2*k2ga1*nu2*t2*t3+amue1
     +*enuh1*ga1*k2ga1*nu1*t6*t7+amue1*enuh1*ga2*k2nu1*nu2*t2*t9-amue2*e
     +nuh1*ga2*k2nu2*nu1*t2*t9+amue2*enuh1*k2ga1*k2ga2*k2nu1*t2*t4-amue1
     +*enuh1*k2ga1*k2ga2*k2nu2*t4*t6-amue2*egah1*ga1*k*t2*t3*xi2*2.0D0+a
     +mue1*egah1*ga1*k*t6*t7*xi1*2.0D0-amue2*enuh1*k*nu1*t2*t9*xi2+amue2
     +*enuh1*k*nu2*t2*t9*xi2-amue1*enuh1*k2ga1*t2*t3*t4*t10-amue1*enuh1*
     +k2ga1*t4*t6*t7*t10-amue1*enuh1*k2nu1*t2*t4*t9*t10+amue2*enuh1*k2nu
     +2*t2*t4*t9*t10-amue2*egah1*k2ga1*t2*t4*xi1*xi2*2.0D0+amue2*enuh1*k
     +2ga1*t2*t4*xi1*xi2*2.0D0-amue2*enuh1*ga1*k2ga1*k2ga2*k2nu1*nu2*t2+
     +amue1*enuh1*ga1*k2ga1*k2ga2*k2nu2*nu1*t6-amue2*egah1*ga2*k*k2ga1*k
     +2nu2*t2*xi1*2.0D0+amue1*egah1*ga1*k*k2ga2*k2nu2*t6*xi1*2.0D0-amue2
     +*enuh1*ga1*k*k2ga1*k2nu1*t2*xi2+amue2*enuh1*ga1*k*k2ga1*k2nu2*t2*x
     +i1-amue2*enuh1*ga2*k*k2ga1*k2nu1*t2*xi2+amue2*enuh1*ga2*k*k2ga1*k2
     +nu2*t2*xi1-amue2*enuh1*k*k2ga1*k2ga2*nu1*t2*xi1+amue2*enuh1*k*k2ga
     +1*k2ga2*nu2*t2*xi1-amue1*enuh1*ga2*k2ga1*nu2*t2*t3*t10-amue1*enuh1
     +*ga1*k2ga1*nu1*t6*t7*t10-amue1*enuh1*ga2*k2nu1*nu2*t2*t9*t10+amue2
     +*enuh1*ga2*k2nu2*nu1*t2*t9*t10+amue2*enuh1*k2ga1*k2ga2*k2nu1*t2*t4
     +*t10-amue1*enuh1*k2ga1*k2ga2*k2nu2*t4*t6*t10-amue2*enuh1*ga1*k2ga1
     +*nu2*t2*xi1*xi2+amue2*enuh1*ga2*k2ga1*nu1*t2*xi1*xi2+amue2*enuh1*k
     +*nu1*t2*t9*t10*xi2-amue2*enuh1*k*nu2*t2*t9*t10*xi2+amue2*enuh1*k2g
     +a1*t2*t4*t10*xi1*xi2*2.0D0+amue2*enuh1*ga1*k2ga1*k2ga2*k2nu1*nu2*t
     +2*t10-amue1*enuh1*ga1*k2ga1*k2ga2*k2nu2*nu1*t6*t10+amue2*enuh1*ga1
     +*k*k2ga1*k2nu1*t2*t10*xi2-amue2*enuh1*ga1*k*k2ga1*k2nu2*t2*t10*xi1
     +-amue2*enuh1*ga2*k*k2ga1*k2nu1*t2*t10*xi2+amue2*enuh1*ga2*k*k2ga1*
     +k2nu2*t2*t10*xi1-amue2*enuh1*k*k2ga1*k2ga2*nu1*t2*t10*xi1+amue2*en
     +uh1*k*k2ga1*k2ga2*nu2*t2*t10*xi1+amue2*enuh1*ga1*k2ga1*nu2*t2*t10*
     +xi1*xi2+amue2*enuh1*ga2*k2ga1*nu1*t2*t10*xi1*xi2
      A0(4,2) = -amue1*enuh1*t2*t3*t4*xi1+amue2*enuh1*t2*t3*t4*xi2*2.0D0
     +-amue1*enuh1*t4*t6*t7*xi1-amue1*egah1*ga1*k*k2nu1*t6*t7*2.0D0+amue
     +2*enuh1*ga1*k*k2nu2*t2*t3+amue2*enuh1*ga2*k*k2nu2*t2*t3-amue2*enuh
     +1*k*k2ga2*nu1*t2*t3+amue2*enuh1*k*k2ga2*nu2*t2*t3-amue1*enuh1*ga2*
     +nu2*t2*t3*xi1-amue2*enuh1*ga1*nu2*t2*t3*xi2+amue2*enuh1*ga2*nu1*t2
     +*t3*xi2+amue1*enuh1*ga1*nu1*t6*t7*xi1-amue1*egah1*k2ga1*k2nu1*t2*t
     +4*xi1*2.0D0+amue2*egah1*k2ga1*k2nu1*t2*t4*xi2*2.0D0+amue1*enuh1*k2
     +ga1*k2nu1*t2*t4*xi1-amue2*enuh1*k2ga1*k2nu2*t2*t4*xi1+amue2*enuh1*
     +k2ga2*k2nu1*t2*t4*xi1-amue1*enuh1*k2ga2*k2nu2*t4*t6*xi1+amue1*enuh
     +1*t2*t3*t4*t10*xi1-amue2*enuh1*t2*t3*t4*t10*xi2*2.0D0+amue1*enuh1*
     +t4*t6*t7*t10*xi1+amue2*egah1*ga2*k*k2ga1*k2nu1*k2nu2*t2*2.0D0-amue
     +1*egah1*ga1*k*k2ga2*k2nu1*k2nu2*t6*2.0D0-amue1*egah1*ga2*k2ga1*k2n
     +u1*nu2*t2*xi1*2.0D0+amue2*egah1*ga1*k2ga2*k2nu1*nu2*t2*xi1*2.0D0+a
     +mue1*enuh1*ga2*k2ga1*k2nu1*nu2*t2*xi1-amue2*enuh1*ga1*k2ga2*k2nu1*
     +nu2*t2*xi1-amue2*enuh1*ga2*k2ga1*k2nu2*nu1*t2*xi1+amue1*enuh1*ga1*
     +k2ga2*k2nu2*nu1*t6*xi1+amue2*enuh1*ga1*k*k2nu2*t2*t3*t10-amue2*enu
     +h1*ga2*k*k2nu2*t2*t3*t10+amue2*enuh1*k*k2ga2*nu1*t2*t3*t10-amue2*e
     +nuh1*k*k2ga2*nu2*t2*t3*t10+amue2*egah1*ga1*k*k2nu1*t2*xi1*xi2*2.0D
     +0-amue2*enuh1*ga1*k*k2nu1*t2*xi1*xi2-amue2*enuh1*ga2*k*k2nu1*t2*xi
     +1*xi2-amue2*enuh1*k*k2ga1*nu1*t2*xi1*xi2+amue2*enuh1*k*k2ga1*nu2*t
     +2*xi1*xi2+amue1*enuh1*ga2*nu2*t2*t3*t10*xi1-amue2*enuh1*ga1*nu2*t2
     +*t3*t10*xi2-amue2*enuh1*ga2*nu1*t2*t3*t10*xi2+amue1*enuh1*ga1*nu1*
     +t6*t7*t10*xi1+amue1*enuh1*k2ga1*k2nu1*t2*t4*t10*xi1-amue2*enuh1*k2
     +ga1*k2nu2*t2*t4*t10*xi1-amue2*enuh1*k2ga2*k2nu1*t2*t4*t10*xi1+amue
     +1*enuh1*k2ga2*k2nu2*t4*t6*t10*xi1+amue1*enuh1*ga2*k2ga1*k2nu1*nu2*
     +t2*t10*xi1-amue2*enuh1*ga1*k2ga2*k2nu1*nu2*t2*t10*xi1-amue2*enuh1*
     +ga2*k2ga1*k2nu2*nu1*t2*t10*xi1+amue1*enuh1*ga1*k2ga2*k2nu2*nu1*t6*
     +t10*xi1-amue2*enuh1*ga1*k*k2nu1*t2*t10*xi1*xi2+amue2*enuh1*ga2*k*k
     +2nu1*t2*t10*xi1*xi2-amue2*enuh1*k*k2ga1*nu1*t2*t10*xi1*xi2+amue2*e
     +nuh1*k*k2ga1*nu2*t2*t10*xi1*xi2
      A0(4,3) = amue2*t2*(t24+t25+t26+t27+t28+amue1*k*k2ga1*k2nu2*t3+amu
     +e1*k*k2nu1*k2nu2*t9-amue1*k2ga1*nu2*t3*xi2-amue1*k2ga2*nu2*t3*xi1-
     +amue1*k2nu1*nu2*t9*xi2-amue2*k*t3*t7*t10-amue1*k*t3*xi1*xi2+amue1*
     +k*t3*t10*xi1*xi2-amue1*k2ga1*k2ga2*k2nu1*nu2*xi1+amue1*k*k2ga1*k2n
     +u2*t3*t10-amue2*k*k2ga2*k2nu2*t3*t10+amue2*k*k2ga1*k2nu1*t7*t10-am
     +ue1*k*k2nu1*k2nu2*t9*t10-amue1*k*k2ga1*k2nu1*xi1*xi2-amue1*k2ga1*n
     +u2*t3*t10*xi2+amue1*k2ga2*nu2*t3*t10*xi1+amue1*k2nu1*nu2*t9*t10*xi
     +2-amue1*egah1*enuh1*k*k2ga1*k2nu2*t3*2.0D0+amue2*k*k2ga1*k2ga2*k2n
     +u1*k2nu2*t10+amue1*egah1*enuh1*k2ga1*nu2*t3*xi2*2.0D0-amue2*egah1*
     +enuh1*k2ga1*nu1*t7*xi1*2.0D0-amue1*k2ga1*k2ga2*k2nu1*nu2*t10*xi1-a
     +mue1*k*k2ga1*k2nu1*t10*xi1*xi2+amue1*egah1*enuh1*k2ga1*k2ga2*k2nu1
     +*nu2*xi1*2.0D0-amue2*egah1*enuh1*k2ga1*k2ga2*k2nu2*nu1*xi1*2.0D0)*
     +(0.0D0,1.0D0)
      A0(4,4) = amue2*t2*(amue2*ga1*t3*t7-amue1*ga2*k2ga1*k2nu2*t3+amue2
     +*ga1*k2ga2*k2nu2*t3+amue2*ga1*k2ga1*k2nu1*t7-amue1*ga2*k2nu1*k2nu2
     +*t9-amue1*k*k2ga1*t3*xi2-amue1*k*k2ga2*t3*xi1-amue1*k*k2nu1*t9*xi2
     ++amue2*ga1*t3*t7*t10+amue1*ga2*t3*xi1*xi2-amue1*ga2*t3*t10*xi1*xi2
     ++amue2*ga1*k2ga1*k2ga2*k2nu1*k2nu2-amue1*k*k2ga1*k2ga2*k2nu1*xi1-a
     +mue1*ga2*k2ga1*k2nu2*t3*t10+amue2*ga1*k2ga2*k2nu2*t3*t10-amue2*ga1
     +*k2ga1*k2nu1*t7*t10+amue1*ga2*k2nu1*k2nu2*t9*t10+amue1*ga2*k2ga1*k
     +2nu1*xi1*xi2-amue1*k*k2ga1*t3*t10*xi2+amue1*k*k2ga2*t3*t10*xi1+amu
     +e1*k*k2nu1*t9*t10*xi2+amue1*egah1*enuh1*ga2*k2ga1*k2nu2*t3*2.0D0-a
     +mue2*ga1*k2ga1*k2ga2*k2nu1*k2nu2*t10+amue1*egah1*enuh1*k*k2ga1*t3*
     +xi2*2.0D0-amue2*egah1*enuh1*k*k2ga1*t7*xi1*2.0D0-amue1*k*k2ga1*k2g
     +a2*k2nu1*t10*xi1+amue1*ga2*k2ga1*k2nu1*t10*xi1*xi2+amue1*egah1*enu
     +h1*k*k2ga1*k2ga2*k2nu1*xi1*2.0D0-amue2*egah1*enuh1*k*k2ga1*k2ga2*k
     +2nu2*xi1*2.0D0-amue1*egah1*enuh1*ga2*k2ga1*k2nu1*xi1*xi2*2.0D0)*(0
     +.0D0,1.0D0)
      A0(4,5) = t2*(amue1*k2ga1*t3*t4+amue2*k2ga2*t3*t4+amue1*k2nu1*t4*t
     +9+amue1*ga2*k2ga1*nu2*t3-amue2*ga1*k2ga2*nu2*t3+amue1*ga2*k2nu1*nu
     +2*t9+amue2*k2ga1*k2ga2*k2nu1*t4-amue2*ga1*k*t3*xi2-amue2*ga2*k*t3*
     +xi2+amue1*k2ga1*t3*t4*t10-amue2*k2ga2*t3*t4*t10-amue1*k2nu1*t4*t9*
     +t10-amue2*ga1*k2ga1*k2ga2*k2nu1*nu2-amue1*egah1*enuh1*k2ga1*t3*t4*
     +2.0D0-amue2*ga1*k*k2ga1*k2nu1*xi2-amue2*ga2*k*k2ga1*k2nu1*xi2+amue
     +1*ga2*k2ga1*nu2*t3*t10-amue2*ga1*k2ga2*nu2*t3*t10-amue1*ga2*k2nu1*
     +nu2*t9*t10+amue2*k2ga1*k2ga2*k2nu1*t4*t10-amue2*ga1*k*t3*t10*xi2+a
     +mue2*ga2*k*t3*t10*xi2-amue1*egah1*enuh1*ga2*k2ga1*nu2*t3*2.0D0+amu
     +e2*ga1*k2ga1*k2ga2*k2nu1*nu2*t10+amue2*ga1*k*k2ga1*k2nu1*t10*xi2-a
     +mue2*ga2*k*k2ga1*k2nu1*t10*xi2+amue2*egah1*enuh1*k2ga1*t4*xi1*xi2*
     +2.0D0-amue2*egah1*enuh1*k*k2ga1*k2ga2*nu1*xi1*2.0D0+amue2*egah1*en
     +uh1*k*k2ga1*k2ga2*nu2*xi1*2.0D0+amue2*egah1*enuh1*ga2*k2ga1*nu1*xi
     +1*xi2*2.0D0)
      A0(4,6) = -t2*(t29+t30+t31+t32-amue2*t3*t4*xi2-amue2*ga1*k*k2nu2*t
     +3-amue2*ga2*k*k2nu2*t3+amue2*ga1*nu2*t3*xi2-amue2*k2ga1*k2nu1*t4*x
     +i2-amue1*t3*t4*t10*xi1+amue2*t3*t4*t10*xi2-amue2*ga1*k*k2ga1*k2nu1
     +*k2nu2-amue2*ga2*k*k2ga1*k2nu1*k2nu2+amue2*ga1*k2ga1*k2nu1*nu2*xi2
     +-amue2*ga1*k*k2nu2*t3*t10+amue2*ga2*k*k2nu2*t3*t10-amue1*ga2*nu2*t
     +3*t10*xi1+amue2*ga1*nu2*t3*t10*xi2+amue1*k2ga1*k2nu1*t4*t10*xi1-am
     +ue2*k2ga1*k2nu1*t4*t10*xi2+amue2*ga1*k*k2ga1*k2nu1*k2nu2*t10-amue2
     +*ga2*k*k2ga1*k2nu1*k2nu2*t10-amue1*egah1*enuh1*k2ga1*k2nu1*t4*xi1*
     +2.0D0+amue2*egah1*enuh1*k2ga1*k2nu2*t4*xi1*2.0D0+amue1*ga2*k2ga1*k
     +2nu1*nu2*t10*xi1-amue2*ga1*k2ga1*k2nu1*nu2*t10*xi2-amue1*egah1*enu
     +h1*ga2*k2ga1*k2nu1*nu2*xi1*2.0D0+amue2*egah1*enuh1*ga2*k2ga1*k2nu2
     +*nu1*xi1*2.0D0+amue2*egah1*enuh1*k*k2ga1*nu1*xi1*xi2*2.0D0-amue2*e
     +gah1*enuh1*k*k2ga1*nu2*xi1*xi2*2.0D0)
      A0(5,1) = t2*(t33+t34+t51+t53-amue2*egah1*ga1*k*k2nu2*t3-amue1*ega
     +h1*k*k2ga1*nu1*t3-amue1*egah1*k*k2ga1*nu2*t3+amue1*enuh1*k*k2ga1*n
     +u1*t3-amue1*enuh1*k*k2nu1*nu2*t9+amue2*enuh1*k*k2nu2*nu1*t9+amue1*
     +egah1*ga1*nu2*t3*xi1+amue2*egah1*ga1*nu1*t3*xi2+amue2*egah1*k2ga1*
     +k2nu2*t4*xi1-amue1*enuh1*k2ga1*k2nu1*t4*xi1+amue1*enuh1*ga1*k2ga1*
     +k2nu1*nu2*xi1-amue2*enuh1*ga1*k2ga1*k2nu2*nu1*xi1+amue2*egah1*ga1*
     +k*k2nu2*t3*t5-amue1*egah1*k*k2ga1*nu1*t3*t5+amue1*egah1*k*k2ga1*nu
     +2*t3*t5+amue1*enuh1*k*k2ga1*nu1*t3*t10+amue1*enuh1*k*k2nu1*nu2*t9*
     +t10-amue2*enuh1*k*k2nu2*nu1*t9*t10-amue2*egah1*ga1*k*k2nu1*xi1*xi2
     +-amue2*enuh1*k*k2ga1*nu1*xi1*xi2-amue1*egah1*ga1*nu2*t3*t5*xi1+amu
     +e2*egah1*ga1*nu1*t3*t5*xi2-amue2*egah1*k2ga1*k2nu2*t4*t5*xi1-amue1
     +*enuh1*k2ga1*k2nu1*t4*t10*xi1-amue1*enuh1*ga1*k2ga1*k2nu1*nu2*t10*
     +xi1+amue2*enuh1*ga1*k2ga1*k2nu2*nu1*t10*xi1-amue2*egah1*ga1*k*k2nu
     +1*t5*xi1*xi2-amue2*enuh1*k*k2ga1*nu1*t10*xi1*xi2)*(-2.0D0)
      A0(5,2) = t2*(t35+t36+t37+amue1*enuh1*k2nu1*t3*t4+amue1*egah1*ga1*
     +k2nu1*nu2*t3-amue1*enuh1*ga1*k2nu1*nu2*t3+amue2*enuh1*ga1*k2nu2*nu
     +1*t3-amue2*egah1*ga1*k*t8*xi2-amue1*enuh1*k*nu1*t3*xi1+amue2*enuh1
     +*k*nu1*t3*xi2-amue1*egah1*k2ga1*t4*t5*t8-amue1*enuh1*k2nu1*t3*t4*t
     +10-amue2*enuh1*k2nu1*t4*xi1*xi2-amue2*egah1*ga1*k*k2nu1*k2nu2*xi1-
     +amue1*egah1*k*k2ga1*k2nu1*nu1*xi1-amue1*egah1*k*k2ga1*k2nu1*nu2*xi
     +1+amue1*enuh1*k*k2ga1*k2nu1*nu2*xi1-amue2*enuh1*k*k2ga1*k2nu2*nu1*
     +xi1+amue1*egah1*ga1*k2nu1*nu2*t3*t5-amue1*enuh1*ga1*k2nu1*nu2*t3*t
     +10+amue2*enuh1*ga1*k2nu2*nu1*t3*t10+amue2*egah1*ga1*k2nu1*nu1*xi1*
     +xi2+amue2*egah1*ga1*k*t5*t8*xi2+amue1*enuh1*k*nu1*t3*t10*xi1-amue2
     +*enuh1*k*nu1*t3*t10*xi2+amue2*enuh1*k2nu1*t4*t10*xi1*xi2-amue2*ega
     +h1*ga1*k*k2nu1*k2nu2*t5*xi1+amue1*egah1*k*k2ga1*k2nu1*nu1*t5*xi1-a
     +mue1*egah1*k*k2ga1*k2nu1*nu2*t5*xi1+amue1*enuh1*k*k2ga1*k2nu1*nu2*
     +t10*xi1-amue2*enuh1*k*k2ga1*k2nu2*nu1*t10*xi1-amue2*egah1*ga1*k2nu
     +1*nu1*t5*xi1*xi2)*2.0D0
      A0(5,3) = amue1*t2*(-amue1*nu2*t38-amue1*nu2*t8*t9+amue1*nu2*t5*t3
     +8+amue1*nu2*t10*t38-amue1*k2ga1*k2nu1*nu2*t3*2.0D0+amue2*k2ga1*k2n
     +u2*nu1*t3+amue2*k2nu1*k2nu2*nu1*t9+amue2*k*k2ga1*t8*xi2+amue2*k*k2
     +nu1*t3*xi2+amue2*k*k2nu2*t3*xi1+amue1*nu2*t5*t8*t9+amue1*nu2*t8*t9
     +*t10-amue1*nu2*t5*t10*t38-amue2*nu1*t3*xi1*xi2-amue2*nu1*t3*t5*xi1
     +*xi2+amue2*nu1*t3*t10*xi1*xi2+amue2*k*k2ga1*k2nu1*k2nu2*xi1-amue1*
     +k2ga1*k2nu1*nu2*t3*t5*2.0D0+amue2*k2ga1*k2nu2*nu1*t3*t5-amue1*k2ga
     +1*k2nu1*nu2*t3*t10*2.0D0+amue2*k2ga1*k2nu2*nu1*t3*t10-amue2*k2nu1*
     +k2nu2*nu1*t5*t9-amue2*k2nu1*k2nu2*nu1*t9*t10-amue2*k2ga1*k2nu1*nu1
     +*xi1*xi2-amue2*k*k2ga1*t5*t8*xi2+amue2*k*k2ga1*t8*t10*xi2+amue2*k*
     +k2nu1*t3*t5*xi2-amue2*k*k2nu2*t3*t5*xi1-amue2*k*k2nu1*t3*t10*xi2-a
     +mue2*k*k2nu2*t3*t10*xi1-amue1*nu2*t5*t8*t9*t10+amue1*egah1*enuh1*k
     +2ga1*k2nu1*nu2*t3*8.0D0-amue2*egah1*enuh1*k2ga1*k2nu2*nu1*t3*4.0D0
     ++amue2*k*k2ga1*k2nu1*k2nu2*t5*xi1+amue2*k*k2ga1*k2nu1*k2nu2*t10*xi
     +1-amue1*k2ga1*k2nu1*nu2*t3*t5*t10*2.0D0+amue2*k2ga1*k2nu2*nu1*t3*t
     +5*t10+amue2*k2nu1*k2nu2*nu1*t5*t9*t10+amue2*k2ga1*k2nu1*nu1*t5*xi1
     +*xi2-amue2*k2ga1*k2nu1*nu1*t10*xi1*xi2-amue2*k*k2ga1*t5*t8*t10*xi2
     +-amue2*k*k2nu1*t3*t5*t10*xi2+amue2*k*k2nu2*t3*t5*t10*xi1+amue2*nu1
     +*t3*t5*t10*xi1*xi2-amue2*egah1*enuh1*k*k2ga1*k2nu1*k2nu2*xi1*4.0D0
     ++amue2*k*k2ga1*k2nu1*k2nu2*t5*t10*xi1+amue2*k2ga1*k2nu1*nu1*t5*t10
     +*xi1*xi2)*(0.0D0,1.0D0)
      A0(5,4) = amue1*t2*(t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65-am
     +ue1*k*t38-amue1*k*t8*t9-amue1*k*k2ga1*k2nu1*t3*2.0D0-amue2*k*k2ga1
     +*k2nu2*t3-amue2*k*k2nu1*k2nu2*t9+amue2*ga1*k2ga1*t8*xi2+amue2*ga1*
     +k2nu1*t3*xi2+amue2*ga1*k2nu2*t3*xi1-amue1*k*t5*t10*t38-amue2*k*t3*
     +t5*xi1*xi2-amue2*k*t3*t10*xi1*xi2+amue2*ga1*k2ga1*k2nu1*k2nu2*xi1-
     +amue1*k*k2ga1*k2nu1*t3*t5*2.0D0+amue2*k*k2ga1*k2nu2*t3*t5-amue1*k*
     +k2ga1*k2nu1*t3*t10*2.0D0-amue2*k*k2ga1*k2nu2*t3*t10-amue2*k*k2nu1*
     +k2nu2*t5*t9+amue2*k*k2nu1*k2nu2*t9*t10-amue2*ga1*k2ga1*t5*t8*xi2-a
     +mue2*ga1*k2ga1*t8*t10*xi2+amue2*ga1*k2nu1*t3*t5*xi2-amue2*ga1*k2nu
     +2*t3*t5*xi1+amue2*ga1*k2nu1*t3*t10*xi2+amue2*ga1*k2nu2*t3*t10*xi1-
     +amue1*k*t5*t8*t9*t10-amue2*egah1*enuh1*ga1*k2nu1*t3*xi2*4.0D0+amue
     +2*ga1*k2ga1*k2nu1*k2nu2*t5*xi1-amue2*ga1*k2ga1*k2nu1*k2nu2*t10*xi1
     +-amue1*k*k2ga1*k2nu1*t3*t5*t10*2.0D0+amue2*k*k2ga1*k2nu2*t3*t5*t10
     ++amue2*k*k2nu1*k2nu2*t5*t9*t10+amue2*ga1*k2ga1*t5*t8*t10*xi2+amue2
     +*ga1*k2nu1*t3*t5*t10*xi2-amue2*ga1*k2nu2*t3*t5*t10*xi1-amue2*egah1
     +*enuh1*k*k2ga1*k2nu1*xi1*xi2*4.0D0-amue2*ga1*k2ga1*k2nu1*k2nu2*t5*
     +t10*xi1)*(0.0D0,1.0D0)
      A0(5,5) = t39-t40+t41+t42+t71+t72+t73+t74+t75+t76+t77+t78+t79-amue
     +2*t2*t3*t4*xi2-amue1*t2*t3*t4*t10*xi1+amue1*k*k2ga1*nu1*t2*t3+amue
     +1*k*k2ga1*nu2*t2*t3+amue1*k*k2nu1*nu1*t2*t9+amue1*k*k2nu1*nu2*t2*t
     +9-amue1*ga1*nu2*t2*t3*xi1-amue2*ga1*nu1*t2*t3*xi2-amue2*k2ga1*k2nu
     +1*t2*t4*xi2-amue1*ga1*k2ga1*k2nu1*nu2*t2*xi1-amue2*ga1*k2ga1*k2nu1
     +*nu1*t2*xi2+amue1*k*k2ga1*nu1*t2*t3*t5-amue1*k*k2ga1*nu2*t2*t3*t5+
     +amue1*k*k2ga1*nu1*t2*t3*t10+amue1*k*k2ga1*nu2*t2*t3*t10-amue1*k*k2
     +nu1*nu1*t2*t5*t9+amue1*k*k2nu1*nu2*t2*t5*t9-amue1*k*k2nu1*nu1*t2*t
     +9*t10-amue1*k*k2nu1*nu2*t2*t9*t10+amue1*ga1*nu2*t2*t3*t5*xi1-amue2
     +*ga1*nu1*t2*t3*t5*xi2-amue1*ga1*nu2*t2*t3*t10*xi1-amue2*ga1*nu1*t2
     +*t3*t10*xi2-amue2*k2ga1*k2nu1*t2*t4*t5*xi2-amue2*k2ga1*k2nu1*t2*t4
     +*t10*xi2-amue2*t2*t3*t4*t5*t10*xi2-amue1*egah1*enuh1*k*k2ga1*nu1*t
     +2*t3*4.0D0-amue1*egah1*enuh1*k2ga1*k2nu1*t2*t4*xi1*4.0D0-amue1*ga1
     +*k2ga1*k2nu1*nu2*t2*t5*xi1+amue1*ga1*k2ga1*k2nu1*nu2*t2*t10*xi1+am
     +ue1*k*k2ga1*nu1*t2*t3*t5*t10-amue1*k*k2ga1*nu2*t2*t3*t5*t10+amue1*
     +k*k2nu1*nu1*t2*t5*t9*t10-amue1*k*k2nu1*nu2*t2*t5*t9*t10+amue1*ga1*
     +nu2*t2*t3*t5*t10*xi1-amue2*ga1*nu1*t2*t3*t5*t10*xi2-amue2*k2ga1*k2
     +nu1*t2*t4*t5*t10*xi2+amue1*ga1*k2ga1*k2nu1*nu2*t2*t5*t10*xi1-amue2
     +*ga1*k2ga1*k2nu1*nu1*t2*t5*t10*xi2
      A0(5,6) = t43+t44+t45+t46+t47+amue2*k2nu2*t2*t3*t4-amue1*ga1*k2ga1
     +*nu2*t2*t8-amue1*ga1*k2nu1*nu2*t2*t3+amue2*ga1*k2nu2*nu1*t2*t3-amu
     +e1*k*nu1*t2*t3*xi1-amue1*k*nu2*t2*t3*xi1-amue1*k2ga1*t2*t4*t5*t8+a
     +mue1*k2ga1*t2*t4*t8*t10-amue2*k2nu2*t2*t3*t4*t5-amue1*k2nu1*t2*t3*
     +t4*t10-amue2*k2nu2*t2*t3*t4*t10+amue2*ga1*k2ga1*k2nu1*k2nu2*nu1*t2
     +-amue1*k*k2ga1*k2nu1*nu1*t2*xi1-amue1*k*k2ga1*k2nu1*nu2*t2*xi1+amu
     +e1*ga1*k2ga1*nu2*t2*t5*t8+amue1*ga1*k2ga1*nu2*t2*t8*t10-amue1*ga1*
     +k2nu1*nu2*t2*t3*t5+amue2*ga1*k2nu2*nu1*t2*t3*t5-amue1*ga1*k2nu1*nu
     +2*t2*t3*t10+amue2*ga1*k2nu2*nu1*t2*t3*t10+amue2*k2ga1*k2nu1*k2nu2*
     +t2*t4*t10-amue1*k*nu1*t2*t3*t5*xi1+amue1*k*nu2*t2*t3*t5*xi1+amue1*
     +k*nu1*t2*t3*t10*xi1+amue1*k*nu2*t2*t3*t10*xi1-amue1*k2ga1*t2*t4*t5
     +*t8*t10-amue1*k2nu1*t2*t3*t4*t5*t10+amue2*k2nu2*t2*t3*t4*t5*t10+am
     +ue1*egah1*enuh1*ga1*k2nu1*nu2*t2*t3*4.0D0-amue2*ga1*k2ga1*k2nu1*k2
     +nu2*nu1*t2*t5-amue2*ga1*k2ga1*k2nu1*k2nu2*nu1*t2*t10+amue1*k*k2ga1
     +*k2nu1*nu1*t2*t5*xi1-amue1*k*k2ga1*k2nu1*nu2*t2*t5*xi1-amue1*k*k2g
     +a1*k2nu1*nu1*t2*t10*xi1-amue1*k*k2ga1*k2nu1*nu2*t2*t10*xi1-amue1*g
     +a1*k2ga1*nu2*t2*t5*t8*t10-amue1*ga1*k2nu1*nu2*t2*t3*t5*t10+amue2*g
     +a1*k2nu2*nu1*t2*t3*t5*t10+amue2*k2ga1*k2nu1*k2nu2*t2*t4*t5*t10+amu
     +e1*k*nu1*t2*t3*t5*t10*xi1-amue1*k*nu2*t2*t3*t5*t10*xi1-amue2*egah1
     +*enuh1*ga1*k*k2nu1*k2nu2*t2*xi1*4.0D0+amue1*egah1*enuh1*k*k2ga1*k2
     +nu1*nu2*t2*xi1*4.0D0-amue2*egah1*enuh1*k*k2ga1*k2nu2*nu1*t2*xi1*4.
     +0D0+amue2*ga1*k2ga1*k2nu1*k2nu2*nu1*t2*t5*t10+amue1*k*k2ga1*k2nu1*
     +nu1*t2*t5*t10*xi1-amue1*k*k2ga1*k2nu1*nu2*t2*t5*t10*xi1
      A0(6,1) = t2*(t48+t49+t50+amue1*egah1*k2ga1*t3*t4-amue1*egah1*ga2*
     +k2ga1*nu1*t3+amue2*egah1*ga1*k2ga2*nu1*t3+amue1*enuh1*ga2*k2ga1*nu
     +1*t3-amue1*egah1*ga1*k*t3*xi1+amue2*egah1*ga1*k*t3*xi2-amue2*enuh1
     +*k*nu1*t9*xi2-amue1*egah1*k2ga1*t3*t4*t5-amue1*enuh1*k2nu1*t4*t9*t
     +10-amue2*egah1*k2ga1*t4*xi1*xi2+amue1*egah1*ga2*k*k2ga1*k2nu1*xi1-
     +amue2*egah1*ga1*k*k2ga2*k2nu1*xi1-amue1*enuh1*ga1*k*k2ga1*k2nu1*xi
     +1-amue1*enuh1*ga2*k*k2ga1*k2nu1*xi1-amue2*enuh1*k*k2ga1*k2ga2*nu1*
     +xi1-amue1*egah1*ga2*k2ga1*nu1*t3*t5+amue2*egah1*ga1*k2ga2*nu1*t3*t
     +5+amue1*enuh1*ga2*k2ga1*nu1*t3*t10+amue2*enuh1*ga1*k2ga1*nu1*xi1*x
     +i2+amue1*egah1*ga1*k*t3*t5*xi1-amue2*egah1*ga1*k*t3*t5*xi2+amue2*e
     +nuh1*k*nu1*t9*t10*xi2+amue2*egah1*k2ga1*t4*t5*xi1*xi2+amue1*egah1*
     +ga2*k*k2ga1*k2nu1*t5*xi1-amue2*egah1*ga1*k*k2ga2*k2nu1*t5*xi1+amue
     +1*enuh1*ga1*k*k2ga1*k2nu1*t10*xi1-amue1*enuh1*ga2*k*k2ga1*k2nu1*t1
     +0*xi1-amue2*enuh1*k*k2ga1*k2ga2*nu1*t10*xi1-amue2*enuh1*ga1*k2ga1*
     +nu1*t10*xi1*xi2)*2.0D0
      A0(6,2) = t2*(t51-t52+t53-t54+amue1*egah1*ga2*k*k2ga1*t8-amue2*ega
     +h1*ga1*k*k2ga2*t8-amue1*egah1*ga1*k*k2nu1*t3+amue1*enuh1*ga1*k*k2n
     +u1*t3+amue1*enuh1*ga2*k*k2nu1*t3+amue2*enuh1*k*k2ga2*nu1*t3-amue1*
     +enuh1*ga2*nu1*t3*xi1-amue2*enuh1*ga1*nu1*t3*xi2-amue1*enuh1*k2ga1*
     +k2nu1*t4*xi1-amue2*enuh1*k2ga2*k2nu1*t4*xi1-amue1*egah1*ga2*k2ga1*
     +k2nu1*nu1*xi1+amue2*egah1*ga1*k2ga2*k2nu1*nu1*xi1-amue1*egah1*ga2*
     +k*k2ga1*t5*t8+amue2*egah1*ga1*k*k2ga2*t5*t8-amue1*egah1*ga1*k*k2nu
     +1*t3*t5+amue1*enuh1*ga1*k*k2nu1*t3*t10-amue1*enuh1*ga2*k*k2nu1*t3*
     +t10-amue2*enuh1*k*k2ga2*nu1*t3*t10+amue2*egah1*ga1*k*k2nu1*xi1*xi2
     ++amue2*enuh1*k*k2ga1*nu1*xi1*xi2+amue1*enuh1*ga2*nu1*t3*t10*xi1-am
     +ue2*enuh1*ga1*nu1*t3*t10*xi2-amue1*enuh1*k2ga1*k2nu1*t4*t10*xi1+am
     +ue2*enuh1*k2ga2*k2nu1*t4*t10*xi1+amue1*egah1*ga2*k2ga1*k2nu1*nu1*t
     +5*xi1-amue2*egah1*ga1*k2ga2*k2nu1*nu1*t5*xi1+amue2*egah1*ga1*k*k2n
     +u1*t5*xi1*xi2+amue2*enuh1*k*k2ga1*nu1*t10*xi1*xi2)*(-2.0D0)
      A0(6,3) = amue1*t2*(t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65-am
     +ue1*k*t38-amue1*k*t8*t9-amue2*k*k2ga1*k2ga2*t8-amue1*k*k2ga1*k2nu1
     +*t3*2.0D0-amue2*k*k2ga2*k2nu1*t3+amue2*k2ga1*nu1*t3*xi2+amue2*k2ga
     +2*nu1*t3*xi1+amue2*k2nu1*nu1*t9*xi2-amue1*k*t5*t10*t38-amue2*k*t3*
     +t5*xi1*xi2-amue2*k*t3*t10*xi1*xi2+amue2*k2ga1*k2ga2*k2nu1*nu1*xi1+
     +amue2*k*k2ga1*k2ga2*t5*t8-amue2*k*k2ga1*k2ga2*t8*t10-amue1*k*k2ga1
     +*k2nu1*t3*t5*2.0D0-amue2*k*k2ga2*k2nu1*t3*t5-amue1*k*k2ga1*k2nu1*t
     +3*t10*2.0D0+amue2*k*k2ga2*k2nu1*t3*t10+amue2*k2ga1*nu1*t3*t5*xi2+a
     +mue2*k2ga2*nu1*t3*t5*xi1+amue2*k2ga1*nu1*t3*t10*xi2-amue2*k2ga2*nu
     +1*t3*t10*xi1-amue2*k2nu1*nu1*t5*t9*xi2-amue2*k2nu1*nu1*t9*t10*xi2-
     +amue1*k*t5*t8*t9*t10-amue2*egah1*enuh1*k2ga1*nu1*t3*xi2*4.0D0-amue
     +2*k2ga1*k2ga2*k2nu1*nu1*t5*xi1+amue2*k2ga1*k2ga2*k2nu1*nu1*t10*xi1
     ++amue2*k*k2ga1*k2ga2*t5*t8*t10-amue1*k*k2ga1*k2nu1*t3*t5*t10*2.0D0
     ++amue2*k*k2ga2*k2nu1*t3*t5*t10+amue2*k2ga1*nu1*t3*t5*t10*xi2-amue2
     +*k2ga2*nu1*t3*t5*t10*xi1+amue2*k2nu1*nu1*t5*t9*t10*xi2-amue2*egah1
     +*enuh1*k*k2ga1*k2nu1*xi1*xi2*4.0D0-amue2*k2ga1*k2ga2*k2nu1*nu1*t5*
     +t10*xi1)*(0.0D0,1.0D0)
      A0(6,4) = amue1*t2*(-amue1*ga2*t38-amue1*ga2*t8*t9+amue1*ga2*t5*t3
     +8+amue1*ga2*t10*t38+amue2*ga1*k2ga1*k2ga2*t8-amue1*ga2*k2ga1*k2nu1
     +*t3*2.0D0+amue2*ga1*k2ga2*k2nu1*t3+amue2*k*k2ga1*t3*xi2+amue2*k*k2
     +ga2*t3*xi1+amue2*k*k2nu1*t9*xi2+amue1*ga2*t5*t8*t9+amue1*ga2*t8*t9
     +*t10-amue1*ga2*t5*t10*t38-amue2*ga1*t3*xi1*xi2+amue2*ga1*t3*t5*xi1
     +*xi2-amue2*ga1*t3*t10*xi1*xi2+amue2*k*k2ga1*k2ga2*k2nu1*xi1-amue2*
     +ga1*k2ga1*k2ga2*t5*t8-amue2*ga1*k2ga1*k2ga2*t8*t10-amue1*ga2*k2ga1
     +*k2nu1*t3*t5*2.0D0+amue2*ga1*k2ga2*k2nu1*t3*t5-amue1*ga2*k2ga1*k2n
     +u1*t3*t10*2.0D0+amue2*ga1*k2ga2*k2nu1*t3*t10-amue2*ga1*k2ga1*k2nu1
     +*xi1*xi2-amue2*k*k2ga1*t3*t5*xi2-amue2*k*k2ga2*t3*t5*xi1+amue2*k*k
     +2ga1*t3*t10*xi2-amue2*k*k2ga2*t3*t10*xi1+amue2*k*k2nu1*t5*t9*xi2-a
     +mue2*k*k2nu1*t9*t10*xi2-amue1*ga2*t5*t8*t9*t10+amue1*egah1*enuh1*g
     +a2*k2ga1*k2nu1*t3*8.0D0-amue2*egah1*enuh1*ga1*k2ga2*k2nu1*t3*4.0D0
     ++amue2*k*k2ga1*k2ga2*k2nu1*t5*xi1+amue2*k*k2ga1*k2ga2*k2nu1*t10*xi
     +1+amue2*ga1*k2ga1*k2ga2*t5*t8*t10-amue1*ga2*k2ga1*k2nu1*t3*t5*t10*
     +2.0D0+amue2*ga1*k2ga2*k2nu1*t3*t5*t10-amue2*ga1*k2ga1*k2nu1*t5*xi1
     +*xi2+amue2*ga1*k2ga1*k2nu1*t10*xi1*xi2-amue2*k*k2ga1*t3*t5*t10*xi2
     ++amue2*k*k2ga2*t3*t5*t10*xi1-amue2*k*k2nu1*t5*t9*t10*xi2+amue2*ga1
     +*t3*t5*t10*xi1*xi2-amue2*egah1*enuh1*k*k2ga1*k2ga2*k2nu1*xi1*4.0D0
     ++amue2*k*k2ga1*k2ga2*k2nu1*t5*t10*xi1+amue2*ga1*k2ga1*k2nu1*t5*t10
     +*xi1*xi2)*(0.0D0,-1.0D0)
      A0(6,5) = t66+t67+t68+t69+t70+amue2*k2ga2*t2*t3*t4-amue1*ga2*k2ga1
     +*nu1*t2*t3+amue2*ga1*k2ga2*nu1*t2*t3-amue1*ga2*k2nu1*nu1*t2*t9-amu
     +e1*ga1*k*t2*t3*xi1-amue1*ga2*k*t2*t3*xi1-amue1*k2ga1*t2*t3*t4*t5-a
     +mue2*k2ga2*t2*t3*t4*t5-amue2*k2ga2*t2*t3*t4*t10+amue1*k2nu1*t2*t4*
     +t5*t9-amue1*k2nu1*t2*t4*t9*t10+amue2*ga1*k2ga1*k2ga2*k2nu1*nu1*t2-
     +amue1*ga1*k*k2ga1*k2nu1*t2*xi1-amue1*ga2*k*k2ga1*k2nu1*t2*xi1-amue
     +1*ga2*k2ga1*nu1*t2*t3*t5+amue2*ga1*k2ga2*nu1*t2*t3*t5-amue1*ga2*k2
     +ga1*nu1*t2*t3*t10+amue2*ga1*k2ga2*nu1*t2*t3*t10+amue1*ga2*k2nu1*nu
     +1*t2*t5*t9+amue1*ga2*k2nu1*nu1*t2*t9*t10+amue2*k2ga1*k2ga2*k2nu1*t
     +2*t4*t5+amue1*ga1*k*t2*t3*t5*xi1+amue1*ga2*k*t2*t3*t5*xi1-amue1*ga
     +1*k*t2*t3*t10*xi1+amue1*ga2*k*t2*t3*t10*xi1-amue1*k2ga1*t2*t3*t4*t
     +5*t10+amue2*k2ga2*t2*t3*t4*t5*t10-amue1*k2nu1*t2*t4*t5*t9*t10+amue
     +1*egah1*enuh1*ga2*k2ga1*nu1*t2*t3*4.0D0-amue2*ga1*k2ga1*k2ga2*k2nu
     +1*nu1*t2*t5-amue2*ga1*k2ga1*k2ga2*k2nu1*nu1*t2*t10-amue1*ga1*k*k2g
     +a1*k2nu1*t2*t5*xi1-amue1*ga2*k*k2ga1*k2nu1*t2*t5*xi1+amue1*ga1*k*k
     +2ga1*k2nu1*t2*t10*xi1-amue1*ga2*k*k2ga1*k2nu1*t2*t10*xi1-amue1*ga2
     +*k2ga1*nu1*t2*t3*t5*t10+amue2*ga1*k2ga2*nu1*t2*t3*t5*t10-amue1*ga2
     +*k2nu1*nu1*t2*t5*t9*t10+amue2*k2ga1*k2ga2*k2nu1*t2*t4*t5*t10+amue1
     +*ga1*k*t2*t3*t5*t10*xi1-amue1*ga2*k*t2*t3*t5*t10*xi1+amue1*egah1*e
     +nuh1*ga2*k*k2ga1*k2nu1*t2*xi1*4.0D0-amue2*egah1*enuh1*ga1*k*k2ga2*
     +k2nu1*t2*xi1*4.0D0-amue2*egah1*enuh1*k*k2ga1*k2ga2*nu1*t2*xi1*4.0D
     +0+amue2*ga1*k2ga1*k2ga2*k2nu1*nu1*t2*t5*t10+amue1*ga1*k*k2ga1*k2nu
     +1*t2*t5*t10*xi1-amue1*ga2*k*k2ga1*k2nu1*t2*t5*t10*xi1
      A0(6,6) = -t39+t40-t41-t42-t71-t72-t73-t74-t75-t76-t77-t78-t79+amu
     +e2*t2*t3*t4*xi2+amue1*t2*t3*t4*t10*xi1-amue1*ga1*k*k2ga1*t2*t8-amu
     +e1*ga2*k*k2ga1*t2*t8-amue1*ga1*k*k2nu1*t2*t3-amue1*ga2*k*k2nu1*t2*
     +t3+amue1*ga2*nu1*t2*t3*xi1+amue2*ga1*nu1*t2*t3*xi2+amue2*k2ga1*k2n
     +u1*t2*t4*xi2+amue1*ga2*k2ga1*k2nu1*nu1*t2*xi1+amue2*ga1*k2ga1*k2nu
     +1*nu1*t2*xi2+amue1*ga1*k*k2ga1*t2*t5*t8+amue1*ga2*k*k2ga1*t2*t5*t8
     ++amue1*ga1*k*k2ga1*t2*t8*t10-amue1*ga2*k*k2ga1*t2*t8*t10-amue1*ga1
     +*k*k2nu1*t2*t3*t5-amue1*ga2*k*k2nu1*t2*t3*t5-amue1*ga1*k*k2nu1*t2*
     +t3*t10+amue1*ga2*k*k2nu1*t2*t3*t10+amue1*ga2*nu1*t2*t3*t5*xi1+amue
     +2*ga1*nu1*t2*t3*t5*xi2-amue1*ga2*nu1*t2*t3*t10*xi1+amue2*ga1*nu1*t
     +2*t3*t10*xi2+amue2*k2ga1*k2nu1*t2*t4*t5*xi2+amue2*k2ga1*k2nu1*t2*t
     +4*t10*xi2+amue2*t2*t3*t4*t5*t10*xi2+amue1*egah1*enuh1*ga1*k*k2nu1*
     +t2*t3*4.0D0+amue1*egah1*enuh1*k2ga1*k2nu1*t2*t4*xi1*4.0D0-amue1*ga
     +2*k2ga1*k2nu1*nu1*t2*t5*xi1+amue1*ga2*k2ga1*k2nu1*nu1*t2*t10*xi1-a
     +mue1*ga1*k*k2ga1*t2*t5*t8*t10+amue1*ga2*k*k2ga1*t2*t5*t8*t10-amue1
     +*ga1*k*k2nu1*t2*t3*t5*t10+amue1*ga2*k*k2nu1*t2*t3*t5*t10-amue1*ga2
     +*nu1*t2*t3*t5*t10*xi1+amue2*ga1*nu1*t2*t3*t5*t10*xi2+amue2*k2ga1*k
     +2nu1*t2*t4*t5*t10*xi2-amue1*ga2*k2ga1*k2nu1*nu1*t2*t5*t10*xi1+amue
     +2*ga1*k2ga1*k2nu1*nu1*t2*t5*t10*xi2
      
      t2 = amue1**2
      t3 = xi1**2
      t4 = t2**2
      t5 = k**2
      t6 = t3**2
      t7 = egah1**2
      t8 = k2ga1**2
      t9 = k2nu1**2
      t10 = enuh1**2
      t11 = amue2**2
      t12 = xi2**2
      B0 = -1.0D0/(-t4*t5*t6-ga2*nu2*t4*t6+t4*t5*t6*t7+t4*t5*t6*t10-t4*t
     +5*t8*t9+ga2*nu2*t4*t6*t7+ga2*nu2*t4*t6*t10-ga2*nu2*t4*t8*t9-k2ga1*
     +k2nu1*t3*t4*t5*2.0D0-t4*t5*t6*t7*t10-t2*t3*t5*t11*t12+t4*t5*t7*t8*
     +t9+t4*t5*t8*t9*t10-ga1*nu1*t2*t3*t11*t12-ga2*nu2*t4*t6*t7*t10+ga2*
     +nu2*t4*t7*t8*t9+ga2*nu2*t4*t8*t9*t10-k2ga1*k2nu1*t3*t4*t5*t7*2.0D0
     +-k2ga1*k2nu1*t3*t4*t5*t10*2.0D0-k2ga2*k2nu2*t2*t3*t5*t11-k2ga1*k2n
     +u1*t2*t5*t11*t12+t2*t3*t5*t7*t11*t12+t2*t3*t5*t10*t11*t12-t4*t5*t7
     +*t8*t9*t10-ga2*k2ga1*k2nu1*nu2*t3*t4*2.0D0-amue1*amue2*k2ga1*k2ga2
     +*t2*t5*t9-amue1*amue2*k2ga1*k2nu2*t2*t3*t5-amue1*amue2*k2ga2*k2nu1
     +*t2*t3*t5-amue1*amue2*k2nu1*k2nu2*t2*t5*t8+egah1*enuh1*k2ga1*k2nu1
     +*t3*t4*t5*8.0D0-ga2*k2ga1*k2nu1*nu2*t3*t4*t7*2.0D0-ga1*k2ga2*k2nu2
     +*nu1*t2*t3*t11-ga2*k2ga1*k2nu1*nu2*t3*t4*t10*2.0D0-ga1*k2ga1*k2nu1
     +*nu1*t2*t11*t12-k2ga1*k2ga2*k2nu1*k2nu2*t2*t5*t11+amue1*amue2*t2*t
     +3*t5*xi1*xi2*2.0D0-ga1*nu1*t2*t3*t7*t11*t12-ga1*nu1*t2*t3*t10*t11*
     +t12-ga2*nu2*t4*t7*t8*t9*t10-k2ga1*k2nu1*t3*t4*t5*t7*t10*2.0D0+k2ga
     +2*k2nu2*t2*t3*t5*t7*t11+k2ga2*k2nu2*t2*t3*t5*t10*t11-k2ga1*k2nu1*t
     +2*t5*t7*t11*t12-k2ga1*k2nu1*t2*t5*t10*t11*t12-t2*t3*t5*t7*t10*t11*
     +t12+amue1*amue2*ga1*k2ga1*k2ga2*nu2*t2*t9+amue1*amue2*ga1*k2ga2*k2
     +nu1*nu2*t2*t3+amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t3+amue1*amue2*ga
     +2*k2nu1*k2nu2*nu1*t2*t8+egah1*enuh1*ga2*k2ga1*k2nu1*nu2*t3*t4*8.0D
     +0+amue1*amue2*ga1*k*k2ga1*t2*t9*xi2+amue1*amue2*ga2*k*k2ga1*t2*t9*
     +xi2+amue1*amue2*ga1*k*k2nu1*t2*t3*xi2+amue1*amue2*ga1*k*k2nu2*t2*t
     +3*xi1+amue1*amue2*ga2*k*k2nu1*t2*t3*xi2+amue1*amue2*ga2*k*k2nu2*t2
     +*t3*xi1+amue1*amue2*k*k2ga1*nu1*t2*t3*xi2+amue1*amue2*k*k2ga2*nu1*
     +t2*t3*xi1+amue1*amue2*k*k2ga1*nu2*t2*t3*xi2+amue1*amue2*k*k2ga2*nu
     +2*t2*t3*xi1+amue1*amue2*k*k2nu1*nu1*t2*t8*xi2+amue1*amue2*k*k2nu1*
     +nu2*t2*t8*xi2-amue1*amue2*k2ga1*k2ga2*t2*t5*t7*t9+amue1*amue2*k2ga
     +1*k2ga2*t2*t5*t9*t10-amue1*amue2*k2ga1*k2nu2*t2*t3*t5*t7+amue1*amu
     +e2*k2ga2*k2nu1*t2*t3*t5*t7+amue1*amue2*k2ga1*k2nu2*t2*t3*t5*t10-am
     +ue1*amue2*k2ga2*k2nu1*t2*t3*t5*t10+amue1*amue2*k2nu1*k2nu2*t2*t5*t
     +7*t8-amue1*amue2*k2nu1*k2nu2*t2*t5*t8*t10-ga1*k2ga1*k2ga2*k2nu1*k2
     +nu2*nu1*t2*t11-amue1*amue2*ga1*nu2*t2*t3*xi1*xi2-amue1*amue2*ga2*n
     +u1*t2*t3*xi1*xi2+amue1*amue2*k2ga1*k2nu1*t2*t5*xi1*xi2*2.0D0-ga1*k
     +2ga2*k2nu2*nu1*t2*t3*t7*t11-ga2*k2ga1*k2nu1*nu2*t3*t4*t7*t10*2.0D0
     +-ga1*k2ga2*k2nu2*nu1*t2*t3*t10*t11+ga1*k2ga1*k2nu1*nu1*t2*t7*t11*t
     +12+ga1*k2ga1*k2nu1*nu1*t2*t10*t11*t12-k2ga1*k2ga2*k2nu1*k2nu2*t2*t
     +5*t7*t11-k2ga1*k2ga2*k2nu1*k2nu2*t2*t5*t10*t11-amue1*amue2*t2*t3*t
     +5*t7*xi1*xi2*2.0D0-amue1*amue2*t2*t3*t5*t10*xi1*xi2*2.0D0-ga1*nu1*
     +t2*t3*t7*t10*t11*t12-k2ga2*k2nu2*t2*t3*t5*t7*t10*t11-k2ga1*k2nu1*t
     +2*t5*t7*t10*t11*t12+amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t2*xi1+amu
     +e1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t2*xi1+amue1*amue2*k*k2ga1*k2ga2*
     +k2nu1*nu1*t2*xi1+amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t2*xi1-amue1*
     +amue2*ga1*k2ga1*k2ga2*nu2*t2*t7*t9-amue1*amue2*ga1*k2ga1*k2ga2*nu2
     +*t2*t9*t10+amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t3*t7+amue1*amue2*ga
     +2*k2ga1*k2nu2*nu1*t2*t3*t7+amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t3*t
     +10+amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t3*t10-amue1*amue2*ga2*k2nu1
     +*k2nu2*nu1*t2*t7*t8-amue1*amue2*ga2*k2nu1*k2nu2*nu1*t2*t8*t10-amue
     +1*amue2*ga1*k2ga1*k2nu1*nu2*t2*xi1*xi2-amue1*amue2*ga2*k2ga1*k2nu1
     +*nu1*t2*xi1*xi2-amue1*amue2*ga1*k*k2ga1*t2*t7*t9*xi2+amue1*amue2*g
     +a2*k*k2ga1*t2*t7*t9*xi2-amue1*amue2*ga1*k*k2ga1*t2*t9*t10*xi2-amue
     +1*amue2*ga2*k*k2ga1*t2*t9*t10*xi2+amue1*amue2*ga1*k*k2nu1*t2*t3*t7
     +*xi2+amue1*amue2*ga1*k*k2nu2*t2*t3*t7*xi1-amue1*amue2*ga2*k*k2nu1*
     +t2*t3*t7*xi2-amue1*amue2*ga2*k*k2nu2*t2*t3*t7*xi1+amue1*amue2*ga1*
     +k*k2nu1*t2*t3*t10*xi2-amue1*amue2*ga1*k*k2nu2*t2*t3*t10*xi1+amue1*
     +amue2*ga2*k*k2nu1*t2*t3*t10*xi2-amue1*amue2*ga2*k*k2nu2*t2*t3*t10*
     +xi1+amue1*amue2*k*k2ga1*nu1*t2*t3*t7*xi2-amue1*amue2*k*k2ga2*nu1*t
     +2*t3*t7*xi1+amue1*amue2*k*k2ga1*nu2*t2*t3*t7*xi2-amue1*amue2*k*k2g
     +a2*nu2*t2*t3*t7*xi1+amue1*amue2*k*k2ga1*nu1*t2*t3*t10*xi2+amue1*am
     +ue2*k*k2ga2*nu1*t2*t3*t10*xi1-amue1*amue2*k*k2ga1*nu2*t2*t3*t10*xi
     +2-amue1*amue2*k*k2ga2*nu2*t2*t3*t10*xi1-amue1*amue2*k*k2nu1*nu1*t2
     +*t7*t8*xi2-amue1*amue2*k*k2nu1*nu2*t2*t7*t8*xi2-amue1*amue2*k*k2nu
     +1*nu1*t2*t8*t10*xi2+amue1*amue2*k*k2nu1*nu2*t2*t8*t10*xi2+egah1*en
     +uh1*ga1*k*k2nu1*t2*t11*t12*xi1*4.0D0+amue1*amue2*k2ga1*k2ga2*t2*t5
     +*t7*t9*t10+amue1*amue2*k2ga1*k2nu2*t2*t3*t5*t7*t10+amue1*amue2*k2g
     +a2*k2nu1*t2*t3*t5*t7*t10+amue1*amue2*k2nu1*k2nu2*t2*t5*t7*t8*t10+g
     +a1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t2*t7*t11+ga1*k2ga1*k2ga2*k2nu1*k2n
     +u2*nu1*t2*t10*t11+egah1*enuh1*k*k2ga1*nu1*t2*t11*t12*xi1*4.0D0-amu
     +e1*amue2*ga1*nu2*t2*t3*t7*xi1*xi2+amue1*amue2*ga2*nu1*t2*t3*t7*xi1
     +*xi2+amue1*amue2*ga1*nu2*t2*t3*t10*xi1*xi2-amue1*amue2*ga2*nu1*t2*
     +t3*t10*xi1*xi2+amue1*amue2*k2ga1*k2nu1*t2*t5*t7*xi1*xi2*2.0D0+amue
     +1*amue2*k2ga1*k2nu1*t2*t5*t10*xi1*xi2*2.0D0-ga1*k2ga2*k2nu2*nu1*t2
     +*t3*t7*t10*t11-ga1*k2ga1*k2nu1*nu1*t2*t7*t10*t11*t12-k2ga1*k2ga2*k
     +2nu1*k2nu2*t2*t5*t7*t10*t11+amue1*amue2*t2*t3*t5*t7*t10*xi1*xi2*2.
     +0D0-amue1*amue2*egah1*enuh1*ga1*k2ga2*k2nu1*nu2*t2*t3*4.0D0-amue1*
     +amue2*egah1*enuh1*ga2*k2ga1*k2nu2*nu1*t2*t3*4.0D0-amue1*amue2*egah
     +1*enuh1*ga1*k*k2nu1*t2*t3*xi2*4.0D0-amue1*amue2*egah1*enuh1*k*k2ga
     +1*nu1*t2*t3*xi2*4.0D0-amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t2*t7*xi
     +1+amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t2*t7*xi1+amue1*amue2*ga1*k*
     +k2ga1*k2nu1*k2nu2*t2*t10*xi1+amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t
     +2*t10*xi1-amue1*amue2*egah1*enuh1*k2ga1*k2nu1*t2*t5*xi1*xi2*8.0D0+
     +amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*t2*t7*xi1+amue1*amue2*k*k2ga1*
     +k2ga2*k2nu1*nu2*t2*t7*xi1-amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*t2*t
     +10*xi1+amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t2*t10*xi1+egah1*enuh1*
     +ga1*k*k2ga2*k2nu1*k2nu2*t2*t11*xi1*4.0D0+amue1*amue2*ga1*k2ga1*k2g
     +a2*nu2*t2*t7*t9*t10+amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t3*t7*t10+a
     +mue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t3*t7*t10+amue1*amue2*ga2*k2nu1*
     +k2nu2*nu1*t2*t7*t8*t10+egah1*enuh1*k*k2ga1*k2ga2*k2nu2*nu1*t2*t11*
     +xi1*4.0D0+amue1*amue2*ga1*k2ga1*k2nu1*nu2*t2*t7*xi1*xi2-amue1*amue
     +2*ga2*k2ga1*k2nu1*nu1*t2*t7*xi1*xi2-amue1*amue2*ga1*k2ga1*k2nu1*nu
     +2*t2*t10*xi1*xi2+amue1*amue2*ga2*k2ga1*k2nu1*nu1*t2*t10*xi1*xi2+am
     +ue1*amue2*ga1*k*k2ga1*t2*t7*t9*t10*xi2-amue1*amue2*ga2*k*k2ga1*t2*
     +t7*t9*t10*xi2+amue1*amue2*ga1*k*k2nu1*t2*t3*t7*t10*xi2-amue1*amue2
     +*ga1*k*k2nu2*t2*t3*t7*t10*xi1-amue1*amue2*ga2*k*k2nu1*t2*t3*t7*t10
     +*xi2+amue1*amue2*ga2*k*k2nu2*t2*t3*t7*t10*xi1+amue1*amue2*k*k2ga1*
     +nu1*t2*t3*t7*t10*xi2-amue1*amue2*k*k2ga2*nu1*t2*t3*t7*t10*xi1-amue
     +1*amue2*k*k2ga1*nu2*t2*t3*t7*t10*xi2+amue1*amue2*k*k2ga2*nu2*t2*t3
     +*t7*t10*xi1+amue1*amue2*k*k2nu1*nu1*t2*t7*t8*t10*xi2-amue1*amue2*k
     +*k2nu1*nu2*t2*t7*t8*t10*xi2-ga1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t2*t7*
     +t10*t11+amue1*amue2*ga1*nu2*t2*t3*t7*t10*xi1*xi2+amue1*amue2*ga2*n
     +u1*t2*t3*t7*t10*xi1*xi2+amue1*amue2*k2ga1*k2nu1*t2*t5*t7*t10*xi1*x
     +i2*2.0D0-amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t2*t7*t10*xi1+amue1*a
     +mue2*ga2*k*k2ga1*k2nu1*k2nu2*t2*t7*t10*xi1-amue1*amue2*k*k2ga1*k2g
     +a2*k2nu1*nu1*t2*t7*t10*xi1+amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t2*
     +t7*t10*xi1+amue1*amue2*ga1*k2ga1*k2nu1*nu2*t2*t7*t10*xi1*xi2+amue1
     +*amue2*ga2*k2ga1*k2nu1*nu1*t2*t7*t10*xi1*xi2-amue1*amue2*egah1*enu
     +h1*ga2*k*k2ga1*k2nu1*k2nu2*t2*xi1*4.0D0-amue1*amue2*egah1*enuh1*k*
     +k2ga1*k2ga2*k2nu1*nu2*t2*xi1*4.0D0)
      
      A0 = A0 * B0
      
      end 

      subroutine psv1eHSfullInverse(A0,cOME,k,Z,ALFA,BETA,AMU)
      IMPLICIT COMPLEX*16 (t) 
      integer N ! 1 para este caso
      PARAMETER (N=1)   
      complex*16 A0(4*N+2,4*N+2)
C     variables de entrada : 
      complex*16 cOME
      real*8 k,Z(N+1)
      complex*16 ALFA(N+1),BETA(N+1),AMU(N+1)
C     auxiliares: 
      complex*16 h(N),egah(N),enuh(N)
      complex*16 ga(N+1),nu(N+1)
      complex*16  UI
      integer e
      complex*16 amue1,ga1,nu1,xi1,k2ga1,k2nu1,egah1,enuh1
      complex*16 amue2,ga2,nu2,xi2,k2ga2,k2nu2
      UI = (0.0,1.0)
      
      do e=1,N+1
        ga(e) = sqrt(cOME**2.0_8/ALFA(e)**2.0_8 - k**2.0_8)
        nu(e) = sqrt(cOME**2.0_8/BETA(e)**2.0_8 - k**2.0_8)
        if(aimag(ga(e)).gt.0.0)ga(e) = conjg(ga(e))
        if(aimag(nu(e)).gt.0.0)nu(e) = conjg(nu(e))
      end do
C     para dentro del estrato
      do e=1,N
       h(e) = Z(e+1)-Z(e)
       egah(e) = exp(-UI * ga(e) * h(e))
       enuh(e) = exp(-UI * nu(e) * h(e))
      end do
      
      e = 1
      amue1 = amu(e)
      ga1 = ga(e)
      nu1 = nu(e)
      xi1 = k**2.0 - nu(e)**2.0 
      k2ga1 = 2.0_8*k*ga(e)
      k2nu1 = 2.0_8*k*nu(e)
      egah1 = egah(e)
      enuh1 = enuh(e)

C     ultimo medio      
      e = 2 
      amue2 = amu(e)
      ga2 = ga(e)
      nu2 = nu(e)
      xi2 = k**2.0 - nu(e)**2.0 
      k2ga2 = 2.0_8*k*ga(e)
      k2nu2 = 2.0_8*k*nu(e)
      
      t2 = xi1**2
      t3 = k**2
      t4 = enuh1**2
      t5 = ga2*nu2
      t6 = amue2**2
      t7 = xi2**2
      t8 = amue1**2
      t9 = t2**2
      t10 = egah1**2
      t11 = k2ga1**2
      t12 = k2nu1**2
      t13 = 1.0D0/amue1
      t14 = t3*t8*t9*t10
      t15 = t3*t4*t8*t9
      t16 = amue1*amue2*t2*t3*xi1*xi2*2.0D0
      t17 = t3*t8*t10*t11*t12
      t18 = t3*t4*t8*t11*t12
      t19 = t2*t3*t6*t7*t10
      t20 = t2*t3*t4*t6*t7
      t21 = ga2*nu2*t8*t9*t10
      t22 = ga2*nu2*t4*t8*t9
      t23 = amue1*amue2*ga1*k*k2nu2*t2*xi1
      t24 = amue1*amue2*ga2*k*k2nu2*t2*xi1
      t25 = amue1*amue2*k*k2ga2*nu1*t2*xi1
      t26 = amue1*amue2*k*k2ga2*nu2*t2*xi1
      t27 = ga2*nu2*t8*t10*t11*t12
      t28 = ga2*nu2*t4*t8*t11*t12
      t29 = k2ga2*k2nu2*t2*t3*t6*t10
      t30 = k2ga2*k2nu2*t2*t3*t4*t6
      t31 = amue1*amue2*k*k2ga2*nu1*t2*t4*xi1
      t32 = amue1*amue2*ga2*nu1*t2*t10*xi1*xi2
      t33 = amue1*amue2*ga1*nu2*t2*t4*xi1*xi2
      t34 = amue1*amue2*k2nu1*k2nu2*t3*t10*t11
      t35 = amue1*amue2*k2ga1*k2ga2*t3*t4*t12
      t36 = amue1*amue2*k2ga2*k2nu1*t2*t3*t10
      t37 = amue1*amue2*k2ga1*k2nu2*t2*t3*t4
      t38 = egah1*enuh1*k2ga1*k2nu1*t2*t3*t8*8.0D0
      t39 = ga1*k2ga1*k2nu1*nu1*t6*t7*t10
      t40 = ga1*k2ga1*k2nu1*nu1*t4*t6*t7
      t41 = amue1*amue2*ga1*k2ga1*k2ga2*nu2*t12
      t42 = amue1*amue2*ga2*k2nu1*k2nu2*nu1*t11
      t43 = amue1*amue2*ga1*k*k2ga1*t12*xi2
      t44 = amue1*amue2*ga2*k*k2ga1*t12*xi2
      t45 = amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2
      t46 = amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2
      t47 = amue1*amue2*k*k2nu1*nu1*t11*xi2
      t48 = amue1*amue2*k*k2nu1*nu2*t11*xi2
      t49 = amue1*amue2*ga1*k*k2nu1*t2*xi2
      t50 = amue1*amue2*ga2*k*k2nu1*t2*xi2
      t51 = amue1*amue2*k2ga1*k2nu1*t3*xi1*xi2*2.0D0
      t52 = amue1*amue2*k*k2ga1*nu1*t2*xi2
      t53 = amue1*amue2*k*k2ga1*nu2*t2*xi2
      t54 = amue1*amue2*t2*t3*t4*t10*xi1*xi2*2.0D0
      t55 = amue1*amue2*ga1*k*k2nu2*t2*t10*xi1
      t56 = amue1*amue2*k2ga1*k2ga2*t3*t4*t10*t12
      t57 = amue1*amue2*k2nu1*k2nu2*t3*t4*t10*t11
      t58 = amue1*amue2*k2ga1*k2nu2*t2*t3*t4*t10
      t59 = amue1*amue2*k2ga2*k2nu1*t2*t3*t4*t10
      t60 = ga1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t6*t10
      t61 = ga1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t4*t6
      t62 = amue1*amue2*ga2*k*k2ga1*t10*t12*xi2
      t63 = amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t10
      t64 = amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t10
      t65 = amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t4
      t66 = amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t4
      t67 = amue1*amue2*k*k2nu1*nu2*t4*t11*xi2
      t68 = egah1*enuh1*ga2*k2ga1*k2nu1*nu2*t2*t8*8.0D0
      t69 = amue1*amue2*ga1*k*k2nu1*t2*t10*xi2
      t70 = amue1*amue2*ga1*k*k2nu1*t2*t4*xi2
      t71 = amue1*amue2*ga2*k*k2nu1*t2*t4*xi2
      t72 = amue1*amue2*k2ga1*k2nu1*t3*t10*xi1*xi2*2.0D0
      t73 = amue1*amue2*k2ga1*k2nu1*t3*t4*xi1*xi2*2.0D0
      t74 = egah1*enuh1*ga1*k*k2nu1*t6*t7*xi1*4.0D0
      t75 = amue1*amue2*k*k2ga1*nu1*t2*t10*xi2
      t76 = amue1*amue2*k*k2ga1*nu2*t2*t10*xi2
      t77 = amue1*amue2*k*k2ga1*nu1*t2*t4*xi2
      t78 = egah1*enuh1*k*k2ga1*nu1*t6*t7*xi1*4.0D0
      t79 = amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*xi1
      t80 = amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*xi1
      t81 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*xi1
      t82 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*xi1
      t83 = amue1*amue2*ga2*k*k2nu2*t2*t4*t10*xi1
      t84 = amue1*amue2*k*k2ga2*nu2*t2*t4*t10*xi1
      t85 = amue1*amue2*ga1*nu2*t2*t4*t10*xi1*xi2
      t86 = amue1*amue2*ga2*nu1*t2*t4*t10*xi1*xi2
      t87 = amue1*amue2*ga1*k2ga1*k2ga2*nu2*t4*t10*t12
      t88 = amue1*amue2*ga2*k2nu1*k2nu2*nu1*t4*t10*t11
      t89 = amue1*amue2*ga1*k*k2ga1*t4*t10*t12*xi2
      t90 = amue1*amue2*ga1*k2ga2*k2nu1*nu2*t2*t4*t10
      t91 = amue1*amue2*ga2*k2ga1*k2nu2*nu1*t2*t4*t10
      t92 = amue1*amue2*k*k2nu1*nu1*t4*t10*t11*xi2
      t93 = amue1*amue2*ga1*k*k2nu1*t2*t4*t10*xi2
      t94 = amue1*amue2*k2ga1*k2nu1*t3*t4*t10*xi1*xi2*2.0D0
      t95 = amue1*amue2*k*k2ga1*nu1*t2*t4*t10*xi2
      t96 = amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t10*xi1
      t97 = amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t4*xi1
      t98 = amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t4*xi1
      t99 = egah1*enuh1*ga1*k*k2ga2*k2nu1*k2nu2*t6*xi1*4.0D0
      t100 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*t10*xi1
      t101 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t10*xi1
      t102 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t4*xi1
      t103 = egah1*enuh1*k*k2ga1*k2ga2*k2nu2*nu1*t6*xi1*4.0D0
      t104 = amue1*amue2*ga1*k2ga1*k2nu1*nu2*t10*xi1*xi2
      t105 = amue1*amue2*ga2*k2ga1*k2nu1*nu1*t4*xi1*xi2
      t106 = amue1*amue2*ga2*k*k2ga1*k2nu1*k2nu2*t4*t10*xi1
      t107 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu2*t4*t10*xi1
      t108 = amue1*amue2*ga1*k2ga1*k2nu1*nu2*t4*t10*xi1*xi2
      t109 = amue1*amue2*ga2*k2ga1*k2nu1*nu1*t4*t10*xi1*xi2
      t112 = t3*t8*t9
      t113 = t3*t8*t11*t12
      t114 = t2*t3*t6*t7
      t115 = ga2*nu2*t8*t9
      t116 = t3*t4*t8*t9*t10
      t117 = ga2*nu2*t8*t11*t12
      t118 = k2ga1*k2nu1*t2*t3*t8*2.0D0
      t119 = k2ga1*k2nu1*t3*t6*t7
      t120 = k2ga2*k2nu2*t2*t3*t6
      t121 = ga1*nu1*t2*t6*t7
      t122 = amue1*amue2*t2*t3*t10*xi1*xi2*2.0D0
      t123 = amue1*amue2*t2*t3*t4*xi1*xi2*2.0D0
      t124 = t3*t4*t8*t10*t11*t12
      t125 = t2*t3*t4*t6*t7*t10
      t126 = amue1*amue2*ga1*nu2*t2*xi1*xi2
      t127 = amue1*amue2*ga2*nu1*t2*xi1*xi2
      t128 = ga2*nu2*t4*t8*t9*t10
      t129 = k2ga1*k2nu1*t2*t3*t8*t10*2.0D0
      t130 = k2ga1*k2nu1*t3*t6*t7*t10
      t131 = k2ga1*k2nu1*t2*t3*t4*t8*2.0D0
      t132 = k2ga1*k2nu1*t3*t4*t6*t7
      t133 = ga1*nu1*t2*t6*t7*t10
      t134 = ga1*nu1*t2*t4*t6*t7
      t135 = amue1*amue2*k2ga1*k2ga2*t3*t12
      t136 = amue1*amue2*k2nu1*k2nu2*t3*t11
      t137 = k2ga1*k2ga2*k2nu1*k2nu2*t3*t6
      t138 = amue1*amue2*k2ga1*k2nu2*t2*t3
      t139 = amue1*amue2*k2ga2*k2nu1*t2*t3
      t140 = ga2*k2ga1*k2nu1*nu2*t2*t8*2.0D0
      t141 = ga1*k2ga1*k2nu1*nu1*t6*t7
      t142 = ga1*k2ga2*k2nu2*nu1*t2*t6
      t143 = amue1*amue2*k*k2ga2*nu1*t2*t10*xi1
      t144 = amue1*amue2*k*k2ga2*nu2*t2*t10*xi1
      t145 = amue1*amue2*k*k2ga2*nu2*t2*t4*xi1
      t146 = amue1*amue2*ga1*nu2*t2*t10*xi1*xi2
      t147 = amue1*amue2*ga2*nu1*t2*t4*xi1*xi2
      t148 = ga2*nu2*t4*t8*t10*t11*t12
      t149 = k2ga1*k2nu1*t2*t3*t4*t8*t10*2.0D0
      t150 = k2ga1*k2nu1*t3*t4*t6*t7*t10
      t151 = k2ga2*k2nu2*t2*t3*t4*t6*t10
      t152 = ga1*nu1*t2*t4*t6*t7*t10
      t153 = amue1*amue2*k2ga1*k2ga2*t3*t10*t12
      t154 = amue1*amue2*k2nu1*k2nu2*t3*t4*t11
      t155 = k2ga1*k2ga2*k2nu1*k2nu2*t3*t6*t10
      t156 = k2ga1*k2ga2*k2nu1*k2nu2*t3*t4*t6
      t157 = amue1*amue2*k2ga1*k2nu2*t2*t3*t10
      t158 = amue1*amue2*k2ga2*k2nu1*t2*t3*t4
      t159 = ga2*k2ga1*k2nu1*nu2*t2*t8*t10*2.0D0
      t160 = ga1*k2ga2*k2nu2*nu1*t2*t6*t10
      t161 = ga2*k2ga1*k2nu1*nu2*t2*t4*t8*2.0D0
      t162 = ga1*k2ga2*k2nu2*nu1*t2*t4*t6
      t163 = ga1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t6
      t164 = amue1*amue2*ga2*k*k2nu2*t2*t10*xi1
      t165 = amue1*amue2*ga1*k*k2nu2*t2*t4*xi1
      t166 = amue1*amue2*ga2*k*k2nu2*t2*t4*xi1
      t167 = k2ga1*k2ga2*k2nu1*k2nu2*t3*t4*t6*t10
      t168 = ga2*k2ga1*k2nu1*nu2*t2*t4*t8*t10*2.0D0
      t169 = ga1*k2ga1*k2nu1*nu1*t4*t6*t7*t10
      t170 = ga1*k2ga2*k2nu2*nu1*t2*t4*t6*t10
      t171 = amue1*amue2*ga1*k2ga1*k2ga2*nu2*t10*t12
      t172 = amue1*amue2*ga2*k2nu1*k2nu2*nu1*t10*t11
      t173 = amue1*amue2*ga1*k2ga1*k2ga2*nu2*t4*t12
      t174 = amue1*amue2*ga2*k2nu1*k2nu2*nu1*t4*t11
      t175 = amue1*amue2*ga1*k*k2ga1*t10*t12*xi2
      t176 = amue1*amue2*ga1*k*k2ga1*t4*t12*xi2
      t177 = amue1*amue2*ga2*k*k2ga1*t4*t12*xi2
      t178 = amue1*amue2*k*k2nu1*nu1*t10*t11*xi2
      t179 = amue1*amue2*k*k2nu1*nu2*t10*t11*xi2
      t180 = amue1*amue2*k*k2nu1*nu1*t4*t11*xi2
      t181 = amue1*amue2*ga2*k*k2nu1*t2*t10*xi2
      t182 = amue1*amue2*k*k2ga1*nu2*t2*t4*xi2
      t183 = amue1*amue2*ga1*k2ga1*k2nu1*nu2*xi1*xi2
      t184 = amue1*amue2*ga2*k2ga1*k2nu1*nu1*xi1*xi2
      t185 = amue1*amue2*ga1*k*k2nu2*t2*t4*t10*xi1
      t186 = amue1*amue2*k*k2ga2*nu1*t2*t4*t10*xi1
      t187 = ga1*k2ga1*k2ga2*k2nu1*k2nu2*nu1*t4*t6*t10
      t188 = amue1*amue2*ga2*k*k2ga1*t4*t10*t12*xi2
      t189 = amue1*amue2*k*k2nu1*nu2*t4*t10*t11*xi2
      t190 = amue1*amue2*ga2*k*k2nu1*t2*t4*t10*xi2
      t191 = amue1*amue2*k*k2ga1*nu2*t2*t4*t10*xi2
      t192 = amue1*amue2*egah1*enuh1*ga1*k2ga2*k2nu1*nu2*t2*4.0D0
      t193 = amue1*amue2*egah1*enuh1*ga2*k2ga1*k2nu2*nu1*t2*4.0D0
      t194 = amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t10*xi1
      t195 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*t4*xi1
      t196 = amue1*amue2*egah1*enuh1*ga1*k*k2nu1*t2*xi2*4.0D0
      t197 = amue1*amue2*egah1*enuh1*k2ga1*k2nu1*t3*xi1*xi2*8.0D0
      t198 = amue1*amue2*egah1*enuh1*k*k2ga1*nu1*t2*xi2*4.0D0
      t199 = amue1*amue2*ga2*k2ga1*k2nu1*nu1*t10*xi1*xi2
      t200 = amue1*amue2*ga1*k2ga1*k2nu1*nu2*t4*xi1*xi2
      t201 = amue1*amue2*ga1*k*k2ga1*k2nu1*k2nu2*t4*t10*xi1
      t202 = amue1*amue2*k*k2ga1*k2ga2*k2nu1*nu1*t4*t10*xi1
      t203 = amue1*amue2*egah1*enuh1*ga2*k*k2ga1*k2nu1*k2nu2*xi1*4.0D0
      t204 = amue1*amue2*egah1*enuh1*k*k2ga1*k2ga2*k2nu1*nu2*xi1*4.0D0
      t110 = t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28
     ++t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t
     +45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61
     ++t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t
     +78+t79+t80+t81+t82+t83+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94
     ++t95+t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108+
     +t109-t112-t113-t114-t115-t116-t117-t118-t119-t120-t121-t122-t123-t
     +124-t125-t126-t127-t128-t129-t130-t131-t132-t133-t134-t135-t136-t1
     +37-t138-t139-t140-t141-t142-t143-t144-t145-t146-t147-t148-t149-t15
     +0-t151-t152-t153-t154-t155-t156-t157-t158-t159-t160-t161-t162-t163
     +-t164-t165-t166-t167-t168-t169-t170-t171-t172-t173-t174-t175-t176-
     +t177-t178-t179-t180-t181-t182-t183-t184-t185-t186-t187-t188-t189-t
     +190-t191-t192-t193-t194-t195-t196-t197-t198-t199-t200-t201-t202-t2
     +03-t204
      t111 = 1.0D0/t110
      t205 = t3*t4
      t206 = ga2*nu2*t4
      t207 = amue2*ga2*k2nu2*nu1
      t208 = amue2*k*nu1*xi2
      t209 = amue2*k*nu2*xi2
      t210 = amue2*t3*xi2*2.0D0
      t211 = amue2*k*k2ga2*nu1
      t212 = amue2*k*k2ga2*nu2
      t213 = amue2*ga1*k*k2nu2
      t214 = amue2*ga2*k*k2nu2
      t215 = amue2*ga1*k2ga2*nu2
      t216 = amue2*ga1*k*xi2
      t217 = amue2*ga2*k*xi2
      t218 = t3*t6*t7
      t219 = k2ga2*k2nu2*t3*t6
      t220 = ga1*nu1*t6*t7
      t221 = ga1*k2ga2*k2nu2*nu1*t6
      t222 = k2ga1*t3
      t223 = ga2*k2ga1*nu2
      t224 = k2nu1*t3
      t225 = ga2*k2nu1*nu2
      t226 = t3*xi2
      t227 = t3*t10
      t228 = ga2*nu2*t10
      t229 = amue2*k*k2ga2*nu1*t4
      t230 = amue2*ga1*nu2*t4*xi2
      t231 = amue2*enuh1*ga2*k2nu2*nu1*2.0D0
      t232 = amue2*enuh1*k*nu1*xi2*2.0D0
      t233 = t231+t232
      t234 = amue2*k2nu2*t3*t4
      t235 = amue2*ga2*k2nu2*nu1*t4
      t236 = amue2*k*nu1*t4*xi2
      t237 = enuh1*k2nu1*t3*2.0D0
      t238 = enuh1*ga2*k2nu1*nu2*2.0D0
      t239 = t237+t238
      t240 = k2nu1*t3*t4
      t241 = ga2*k2nu1*nu2*t4
      t242 = t224+t225+t240+t241
      t243 = t3+t5-t205-t206
      t244 = enuh1*k*nu1*t6*t7*2.0D0
      t245 = enuh1*k*k2ga2*k2nu2*nu1*t6*2.0D0
      t246 = t244+t245
      t247 = ga1*nu1*t4*t6*t7
      t248 = ga1*k2ga2*k2nu2*nu1*t4*t6
      t249 = amue2*enuh1*k*k2nu1*t7*(0.0D0,2.0D0)
      t250 = amue2*enuh1*k*k2ga2*k2nu1*k2nu2*(0.0D0,2.0D0)
      t251 = t249+t250
      t252 = amue2*t251
      t253 = k2ga1*t3*xi2
      t254 = k*k2ga1*k2ga2*nu1
      t255 = k*k2ga1*k2ga2*nu2
      t256 = k2ga1*t3*t4*xi2
      t257 = k*k2ga1*k2ga2*nu2*t4
      t258 = ga2*k2ga1*nu1*t4*xi2
      t259 = t253+t254+t255+t256+t257+t258-ga2*k2ga1*nu1*xi2-k*k2ga1*k2g
     +a2*nu1*t4
      t260 = amue2*t259
      t261 = k2ga1*t3*t4
      t262 = ga2*k2ga1*nu2*t4
      t263 = k*k2ga2*nu1
      t264 = k*k2ga2*nu2
      t265 = k*k2ga2*nu1*t4
      t266 = t226+t263+t264+t265-ga2*nu1*xi2-t3*t4*xi2-k*k2ga2*nu2*t4-ga
     +2*nu1*t4*xi2
      t267 = amue1*t2*t243*xi1
      t268 = t3*t12
      t269 = ga2*nu2*t12
      t353 = t3*t4*t12
      t270 = t268+t269-t353-ga2*nu2*t4*t12
      t271 = amue1*t270
      t272 = k2nu1*k2nu2*t3
      t273 = k2nu1*k2nu2*t3*t4
      t274 = ga2*k2nu1*k2nu2*nu1*t4
      t275 = k*k2nu1*nu1*t4*xi2
      t276 = t272+t273+t274+t275-ga2*k2nu1*k2nu2*nu1-k*k2nu1*nu1*xi2-k*k
     +2nu1*nu2*xi2-k*k2nu1*nu2*t4*xi2
      t277 = amue2*t276
      t278 = t271+t277
      t279 = k2ga1*t278
      t280 = k2nu2*t3*t4
      t281 = ga2*k2nu2*nu1
      t282 = k*nu1*xi2
      t283 = k*nu2*xi2
      t284 = ga2*k2nu2*nu1*t4
      t285 = k*nu1*t4*xi2
      t352 = k2nu2*t3
      t286 = t280+t281+t282+t283+t284+t285-t352-k*nu2*t4*xi2
      t287 = amue2*egah1*t3*xi2*2.0D0
      t288 = amue2*egah1*ga2*k*k2nu2*2.0D0
      t289 = t287+t288
      t290 = egah1*t3*2.0D0
      t291 = egah1*ga2*nu2*2.0D0
      t292 = egah1*ga1*k*t6*t7*2.0D0
      t293 = egah1*ga1*k*k2ga2*k2nu2*t6*2.0D0
      t294 = t292+t293
      t295 = amue2*enuh1*t3*xi2*2.0D0
      t296 = amue2*egah1*ga1*k2ga2*nu2*2.0D0
      t297 = amue2*egah1*ga1*k*xi2*2.0D0
      t298 = enuh1*t3
      t299 = enuh1*t3*t10
      t300 = enuh1*ga2*nu2
      t301 = enuh1*ga2*nu2*t10
      t302 = amue2*t7*(0.0D0,1.0D0)
      t303 = amue2*k2ga2*k2nu2*(0.0D0,1.0D0)
      t304 = amue1*xi2*(0.0D0,1.0D0)
      t305 = amue2*egah1*k*k2ga1*t7*(0.0D0,2.0D0)
      t306 = amue2*egah1*k*k2ga1*k2ga2*k2nu2*(0.0D0,2.0D0)
      t307 = t305+t306
      t308 = amue2*t307
      t309 = egah1*k*k2ga1*k2ga2*(0.0D0,2.0D0)
      t310 = egah1*ga2*k2ga1*xi2*(0.0D0,-2.0D0)
      t311 = t3*t11
      t312 = ga2*nu2*t11
      t387 = t3*t10*t11
      t313 = t311+t312-t387-ga2*nu2*t10*t11
      t314 = amue1*t313
      t315 = k2ga1*k2ga2*t3
      t316 = k2ga1*k2ga2*t3*t10
      t317 = ga1*k2ga1*k2ga2*nu2*t10
      t318 = ga1*k*k2ga1*t10*xi2
      t319 = t315+t316+t317+t318-ga1*k2ga1*k2ga2*nu2-ga1*k*k2ga1*xi2-ga2
     +*k*k2ga1*xi2-ga2*k*k2ga1*t10*xi2
      t320 = amue2*t319
      t321 = t314+t320
      t322 = k2nu1*t321
      t323 = k2ga2*t3*t10
      t324 = ga1*k2ga2*nu2
      t325 = ga1*k*xi2
      t326 = ga2*k*xi2
      t327 = ga1*k2ga2*nu2*t10
      t328 = ga1*k*t10*xi2
      t386 = k2ga2*t3
      t329 = t323+t324+t325+t326+t327+t328-t386-ga2*k*t10*xi2
      t330 = amue2*t329
      t331 = k2ga1*t3*t10
      t332 = ga2*k2ga1*nu2*t10
      t333 = k2nu1*t3*xi2
      t334 = ga1*k*k2nu1*k2nu2
      t335 = ga2*k*k2nu1*k2nu2
      t336 = k2nu1*t3*t10*xi2
      t337 = ga2*k*k2nu1*k2nu2*t10
      t338 = ga1*k2nu1*nu2*t10*xi2
      t339 = t333+t334+t335+t336+t337+t338-ga1*k2nu1*nu2*xi2-ga1*k*k2nu1
     +*k2nu2*t10
      t340 = amue2*t339
      t341 = k2nu1*t3*t10
      t342 = ga2*k2nu1*nu2*t10
      t343 = ga1*k*k2nu2
      t344 = ga2*k*k2nu2
      t345 = ga1*k*k2nu2*t10
      t346 = t226+t343+t344+t345-ga1*nu2*xi2-t3*t10*xi2-ga2*k*k2nu2*t10-
     +ga1*nu2*t10*xi2
      t347 = t3+t5-t227-t228
      t348 = amue1*t2*t347*xi1
      t349 = enuh1*k2ga1*t3*2.0D0
      t350 = k2ga1*t3*t4*t10
      t351 = t3*t4*t10
      t354 = amue1*egah1*k2ga1*t3*2.0D0
      t355 = amue1*egah1*k2ga1*t3*t4*2.0D0
      t356 = egah1*ga1*k*2.0D0
      t357 = egah1*ga1*k*t4*2.0D0
      t358 = egah1*k2ga1*t3*2.0D0
      t359 = amue1*enuh1*k2nu1*t3*2.0D0
      t360 = egah1*k2nu1*t3*2.0D0
      t361 = enuh1*k*nu1*2.0D0
      t362 = amue1*t9*(0.0D0,1.0D0)
      t363 = amue1*t11*t12*(0.0D0,1.0D0)
      t364 = amue1*t9*t10*(0.0D0,-1.0D0)
      t365 = amue1*t4*t9*(0.0D0,-1.0D0)
      t366 = amue1*t10*t11*t12*(0.0D0,-1.0D0)
      t367 = amue1*t4*t11*t12*(0.0D0,-1.0D0)
      t368 = amue1*t4*t9*t10*(0.0D0,1.0D0)
      t369 = amue1*k2ga1*k2nu1*t2*(0.0D0,2.0D0)
      t370 = amue1*k2ga1*k2nu1*t2*t10*(0.0D0,2.0D0)
      t371 = amue1*k2ga1*k2nu1*t2*t4*(0.0D0,2.0D0)
      t372 = amue1*t4*t10*t11*t12*(0.0D0,1.0D0)
      t373 = amue1*egah1*enuh1*k2ga1*k2nu1*t2*(0.0D0,-8.0D0)
      t374 = amue1*k2ga1*k2nu1*t2*t4*t10*(0.0D0,2.0D0)
      t375 = t362+t363+t364+t365+t366+t367+t368+t369+t370+t371+t372+t373
     ++t374
      t376 = amue1*t375
      t377 = t2*xi1*xi2*(0.0D0,1.0D0)
      t378 = t2*t10*xi1*xi2*(0.0D0,-1.0D0)
      t379 = t2*t4*xi1*xi2*(0.0D0,-1.0D0)
      t380 = k2ga1*k2nu1*xi1*xi2*(0.0D0,1.0D0)
      t381 = t2*t4*t10*xi1*xi2*(0.0D0,1.0D0)
      t382 = k2ga1*k2nu1*t10*xi1*xi2*(0.0D0,1.0D0)
      t383 = k2ga1*k2nu1*t4*xi1*xi2*(0.0D0,1.0D0)
      t384 = k2ga1*k2nu1*t4*t10*xi1*xi2*(0.0D0,1.0D0)
      t385 = egah1*enuh1*k2ga1*k2nu1*xi1*xi2*(0.0D0,-4.0D0)
      t388 = egah1*enuh1*ga1*k*k2nu1*4.0D0
      t389 = k2nu1*t3*t4*t10
      t390 = egah1*enuh1*k*k2ga1*nu1*4.0D0
      t391 = t388+t390
      t392 = ga1*nu1
      t393 = ga1*nu1*t10
      t394 = ga1*nu1*t4
      t395 = ga1*nu1*t4*t10
      t396 = t3-t205-t227+t351+t392+t393+t394+t395
      t397 = t2*t396
      A0(1,1) = t13*t111*(t8*(t2*xi1*(t3+t5-t3*t4-ga2*nu2*t4)+k2ga1*k2nu
     +1*xi1*(t3+t5+t205+t206-egah1*(enuh1*t3*2.0D0+enuh1*ga2*nu2*2.0D0))
     +)-amue1*(k2nu1*(xi1*(t215+t216+t217-amue2*k2ga2*t3-amue2*k2ga2*t3*
     +t4+amue2*ga1*k2ga2*nu2*t4+amue2*ga1*k*t4*xi2+amue2*ga2*k*t4*xi2)-e
     +gah1*k2ga1*(t295+amue2*enuh1*k*k2ga2*nu2*2.0D0))+t2*(t210+t211+t21
     +2+t213+t214+t229+t230-amue2*ga1*nu2*xi2-amue2*ga2*nu1*xi2-amue2*t3
     +*t4*xi2*2.0D0-amue2*ga1*k*k2nu2*t4-amue2*ga2*k*k2nu2*t4-amue2*k*k2
     +ga2*nu2*t4-amue2*ga2*nu1*t4*xi2)+k2ga1*xi1*(t207+t208+t209+t234+t2
     +35+t236-egah1*t233-amue2*k2nu2*t3-amue2*k*nu2*t4*xi2))+xi1*(t218+t
     +219+t220+t221+t247+t248-t3*t4*t6*t7-k2ga2*k2nu2*t3*t4*t6)-egah1*k2
     +ga1*t246)
      A0(1,2) = -t13*t111*(k2ga1*t3*t8*t12+k2nu1*t2*t3*t8+k2nu1*t3*t6*t7
     ++amue1*amue2*k2ga2*t3*t12+ga2*k2ga1*nu2*t8*t12+ga2*k2nu1*nu2*t2*t8
     ++ga1*k2nu1*nu1*t6*t7+k2ga2*k2nu1*k2nu2*t3*t6-k2ga1*t3*t4*t8*t12+k2
     +nu1*t2*t3*t4*t8+k2nu1*t3*t4*t6*t7-amue1*amue2*ga1*k2ga2*nu2*t12+am
     +ue1*amue2*k2ga1*k2nu1*k2nu2*t3-amue1*amue2*ga1*k*t12*xi2-amue1*amu
     +e2*ga2*k*t12*xi2-amue1*amue2*k2ga2*t3*t4*t12+ga1*k2ga2*k2nu1*k2nu2
     +*nu1*t6-egah1*enuh1*k2nu1*t2*t3*t8*2.0D0-amue1*amue2*k2nu1*t3*xi1*
     +xi2*2.0D0-ga2*k2ga1*nu2*t4*t8*t12+ga2*k2nu1*nu2*t2*t4*t8-ga1*k2nu1
     +*nu1*t4*t6*t7+k2ga2*k2nu1*k2nu2*t3*t4*t6-amue1*amue2*ga2*k2ga1*k2n
     +u1*k2nu2*nu1-amue1*amue2*ga1*k*k2nu1*k2nu2*xi1-amue1*amue2*ga2*k*k
     +2nu1*k2nu2*xi1-amue1*amue2*k*k2ga1*k2nu1*nu1*xi2-amue1*amue2*k*k2g
     +a2*k2nu1*nu1*xi1-amue1*amue2*k*k2ga1*k2nu1*nu2*xi2-amue1*amue2*k*k
     +2ga2*k2nu1*nu2*xi1+amue1*amue2*ga1*k2ga2*nu2*t4*t12+amue1*amue2*k2
     +ga1*k2nu1*k2nu2*t3*t4-egah1*enuh1*ga2*k2nu1*nu2*t2*t8*2.0D0+amue1*
     +amue2*ga1*k2nu1*nu2*xi1*xi2+amue1*amue2*ga2*k2nu1*nu1*xi1*xi2+amue
     +1*amue2*ga1*k*t4*t12*xi2+amue1*amue2*ga2*k*t4*t12*xi2-ga1*k2ga2*k2
     +nu1*k2nu2*nu1*t4*t6-egah1*enuh1*k*nu1*t6*t7*xi1*2.0D0-amue1*amue2*
     +k2nu1*t3*t4*xi1*xi2*2.0D0+amue1*amue2*egah1*enuh1*ga2*k2nu2*nu1*t2
     +*2.0D0+amue1*amue2*ga2*k2ga1*k2nu1*k2nu2*nu1*t4+amue1*amue2*egah1*
     +enuh1*k*nu1*t2*xi2*2.0D0-amue1*amue2*ga1*k*k2nu1*k2nu2*t4*xi1-amue
     +1*amue2*ga2*k*k2nu1*k2nu2*t4*xi1+amue1*amue2*egah1*enuh1*k2nu1*t3*
     +xi1*xi2*2.0D0+amue1*amue2*k*k2ga1*k2nu1*nu1*t4*xi2+amue1*amue2*k*k
     +2ga2*k2nu1*nu1*t4*xi1-amue1*amue2*k*k2ga1*k2nu1*nu2*t4*xi2-amue1*a
     +mue2*k*k2ga2*k2nu1*nu2*t4*xi1-egah1*enuh1*k*k2ga2*k2nu2*nu1*t6*xi1
     +*2.0D0+amue1*amue2*ga1*k2nu1*nu2*t4*xi1*xi2-amue1*amue2*ga2*k2nu1*
     +nu1*t4*xi1*xi2+amue1*amue2*egah1*enuh1*k*k2ga2*k2nu1*nu2*xi1*2.0D0
     +)
      A0(1,3) = t111*(k2ga1*(amue2*(amue2*egah1*k2nu1*nu1*t7*(0.0D0,1.0D
     +0)+amue2*egah1*k2ga2*k2nu1*k2nu2*nu1*(0.0D0,1.0D0)-amue2*egah1*k2n
     +u1*nu1*t4*t7*(0.0D0,1.0D0)-amue2*egah1*k2ga2*k2nu1*k2nu2*nu1*t4*(0
     +.0D0,1.0D0))-amue1*amue2*(egah1*k2ga2*nu2*t12*(0.0D0,1.0D0)+egah1*
     +k*t12*xi2*(0.0D0,1.0D0)-egah1*k2ga2*nu2*t4*t12*(0.0D0,1.0D0)-egah1
     +*k*t4*t12*xi2*(0.0D0,1.0D0)))+xi1*(t252-amue1*amue2*k2ga1*(egah1*k
     +*k2nu1*k2nu2*(0.0D0,1.0D0)-enuh1*k*k2nu1*k2nu2*(0.0D0,2.0D0)-egah1
     +*k2nu1*nu2*xi2*(0.0D0,1.0D0)+enuh1*k2nu1*nu2*xi2*(0.0D0,2.0D0)+ega
     +h1*k*k2nu1*k2nu2*t4*(0.0D0,1.0D0)-egah1*k2nu1*nu2*t4*xi2*(0.0D0,1.
     +0D0)))-t2*(amue2*(amue2*egah1*nu1*t7*(0.0D0,1.0D0)+amue2*egah1*k2g
     +a2*k2nu2*nu1*(0.0D0,1.0D0)+amue2*egah1*nu1*t4*t7*(0.0D0,1.0D0)+amu
     +e2*egah1*k2ga2*k2nu2*nu1*t4*(0.0D0,1.0D0))-amue1*amue2*(egah1*k2ga
     +2*k2nu1*nu2*(0.0D0,1.0D0)-enuh1*k2ga2*k2nu1*nu2*(0.0D0,2.0D0)+egah
     +1*k*k2nu1*xi2*(0.0D0,1.0D0)-enuh1*k*k2nu1*xi2*(0.0D0,2.0D0)+egah1*
     +k2ga2*k2nu1*nu2*t4*(0.0D0,1.0D0)+egah1*k*k2nu1*t4*xi2*(0.0D0,1.0D0
     +)))+amue1*amue2*t2*xi1*(egah1*k*k2nu2*(0.0D0,1.0D0)-egah1*nu2*xi2*
     +(0.0D0,1.0D0)-egah1*k*k2nu2*t4*(0.0D0,1.0D0)+egah1*nu2*t4*xi2*(0.0
     +D0,1.0D0)))
      A0(1,4) = -t111*(amue2*(amue2*egah1*k*t2*t7*(0.0D0,-1.0D0)-amue2*e
     +gah1*k*k2ga2*k2nu2*t2*(0.0D0,1.0D0)+amue2*egah1*k*k2ga1*k2nu1*t7*(
     +0.0D0,1.0D0)-amue2*enuh1*ga1*k2nu1*t7*xi1*(0.0D0,2.0D0)+amue2*egah
     +1*k*t2*t4*t7*(0.0D0,1.0D0)+amue2*egah1*k*k2ga1*k2ga2*k2nu1*k2nu2*(
     +0.0D0,1.0D0)-amue2*enuh1*ga1*k2ga2*k2nu1*k2nu2*xi1*(0.0D0,2.0D0)+a
     +mue2*egah1*k*k2ga2*k2nu2*t2*t4*(0.0D0,1.0D0)+amue2*egah1*k*k2ga1*k
     +2nu1*t4*t7*(0.0D0,1.0D0)+amue2*egah1*k*k2ga1*k2ga2*k2nu1*k2nu2*t4*
     +(0.0D0,1.0D0))-amue1*amue2*(egah1*k*k2ga1*k2ga2*t12*(0.0D0,-1.0D0)
     ++egah1*k*k2ga2*k2nu1*t2*(0.0D0,1.0D0)-enuh1*k*k2ga2*k2nu1*t2*(0.0D
     +0,2.0D0)+egah1*ga2*k2ga1*t12*xi2*(0.0D0,1.0D0)-egah1*ga2*k2nu1*t2*
     +xi2*(0.0D0,1.0D0)-egah1*ga2*k2nu2*t2*xi1*(0.0D0,1.0D0)+enuh1*ga2*k
     +2nu1*t2*xi2*(0.0D0,2.0D0)-egah1*k*t2*xi1*xi2*(0.0D0,1.0D0)+egah1*k
     +*t2*t4*xi1*xi2*(0.0D0,1.0D0)+egah1*ga2*k2ga1*k2nu1*k2nu2*xi1*(0.0D
     +0,1.0D0)-enuh1*ga2*k2ga1*k2nu1*k2nu2*xi1*(0.0D0,2.0D0)+egah1*k*k2g
     +a1*k2ga2*t4*t12*(0.0D0,1.0D0)+egah1*k*k2ga2*k2nu1*t2*t4*(0.0D0,1.0
     +D0)+egah1*k*k2ga1*k2nu1*xi1*xi2*(0.0D0,1.0D0)-enuh1*k*k2ga1*k2nu1*
     +xi1*xi2*(0.0D0,2.0D0)-egah1*ga2*k2ga1*t4*t12*xi2*(0.0D0,1.0D0)-ega
     +h1*ga2*k2nu1*t2*t4*xi2*(0.0D0,1.0D0)+egah1*ga2*k2nu2*t2*t4*xi1*(0.
     +0D0,1.0D0)+egah1*ga2*k2ga1*k2nu1*k2nu2*t4*xi1*(0.0D0,1.0D0)+egah1*
     +k*k2ga1*k2nu1*t4*xi1*xi2*(0.0D0,1.0D0)))
      A0(1,5) = t111*(egah1*(t267+k2nu1*(t260-amue1*xi1*(t222+t223+t261+
     +t262))-amue2*t2*t266)+k2nu1*xi1*(amue1*(t349+enuh1*ga2*k2ga1*nu2*2
     +.0D0)-amue2*(enuh1*k2ga2*t3*(-2.0D0)+enuh1*ga1*k2ga2*nu2*2.0D0+
     +enuh1*ga1*k*xi2*2.0D0+enuh1*ga2*k*xi2*2.0D0)))
      A0(1,6) = -t111*(egah1*(t279-t2*(amue1*t242-amue2*t286))-amue2*xi1
     +*(enuh1*k2nu1*t3*xi2*2.0D0+enuh1*ga1*k*k2nu1*k2nu2*2.0D0+enuh1*ga2
     +*k*k2nu1*k2nu2*2.0D0-enuh1*ga1*k2nu1*nu2*xi2*2.0D0)+amue1*t2*t239)
      A0(2,1) = -t13*t111*(k2ga1*t2*t3*t8+k2ga1*t3*t6*t7+k2nu1*t3*t8*t11
     ++amue1*amue2*k2nu2*t3*t11+ga2*k2ga1*nu2*t2*t8+ga1*k2ga1*nu1*t6*t7+
     +ga2*k2nu1*nu2*t8*t11+k2ga1*k2ga2*k2nu2*t3*t6+k2ga1*t2*t3*t8*t10+k2
     +ga1*t3*t6*t7*t10-k2nu1*t3*t8*t10*t11-amue1*amue2*ga2*k2nu2*nu1*t11
     ++amue1*amue2*k2ga1*k2ga2*k2nu1*t3-amue1*amue2*k*nu1*t11*xi2-amue1*
     +amue2*k*nu2*t11*xi2-amue1*amue2*k2nu2*t3*t10*t11+ga1*k2ga1*k2ga2*k
     +2nu2*nu1*t6-egah1*enuh1*k2ga1*t2*t3*t8*2.0D0-amue1*amue2*k2ga1*t3*
     +xi1*xi2*2.0D0+ga2*k2ga1*nu2*t2*t8*t10-ga1*k2ga1*nu1*t6*t7*t10-ga2*
     +k2nu1*nu2*t8*t10*t11+k2ga1*k2ga2*k2nu2*t3*t6*t10-amue1*amue2*ga1*k
     +2ga1*k2ga2*k2nu1*nu2-amue1*amue2*ga1*k*k2ga1*k2nu1*xi2-amue1*amue2
     +*ga1*k*k2ga1*k2nu2*xi1-amue1*amue2*ga2*k*k2ga1*k2nu1*xi2-amue1*amu
     +e2*ga2*k*k2ga1*k2nu2*xi1-amue1*amue2*k*k2ga1*k2ga2*nu1*xi1-amue1*a
     +mue2*k*k2ga1*k2ga2*nu2*xi1+amue1*amue2*ga2*k2nu2*nu1*t10*t11+amue1
     +*amue2*k2ga1*k2ga2*k2nu1*t3*t10-egah1*enuh1*ga2*k2ga1*nu2*t2*t8*2.
     +0D0+amue1*amue2*ga1*k2ga1*nu2*xi1*xi2+amue1*amue2*ga2*k2ga1*nu1*xi
     +1*xi2+amue1*amue2*k*nu1*t10*t11*xi2+amue1*amue2*k*nu2*t10*t11*xi2-
     +egah1*enuh1*ga1*k*t6*t7*xi1*2.0D0-ga1*k2ga1*k2ga2*k2nu2*nu1*t6*t10
     +-amue1*amue2*k2ga1*t3*t10*xi1*xi2*2.0D0+amue1*amue2*egah1*enuh1*ga
     +1*k2ga2*nu2*t2*2.0D0+amue1*amue2*egah1*enuh1*ga1*k*t2*xi2*2.0D0+am
     +ue1*amue2*ga1*k2ga1*k2ga2*k2nu1*nu2*t10+amue1*amue2*ga1*k*k2ga1*k2
     +nu1*t10*xi2+amue1*amue2*ga1*k*k2ga1*k2nu2*t10*xi1-amue1*amue2*ga2*
     +k*k2ga1*k2nu1*t10*xi2-amue1*amue2*ga2*k*k2ga1*k2nu2*t10*xi1+amue1*
     +amue2*egah1*enuh1*k2ga1*t3*xi1*xi2*2.0D0-amue1*amue2*k*k2ga1*k2ga2
     +*nu1*t10*xi1-amue1*amue2*k*k2ga1*k2ga2*nu2*t10*xi1-egah1*enuh1*ga1
     +*k*k2ga2*k2nu2*t6*xi1*2.0D0-amue1*amue2*ga1*k2ga1*nu2*t10*xi1*xi2+
     +amue1*amue2*ga2*k2ga1*nu1*t10*xi1*xi2+amue1*amue2*egah1*enuh1*ga2*
     +k*k2ga1*k2nu2*xi1*2.0D0)
      A0(2,2) = -t13*t111*(xi1*(t218+t219+t220+t221-t3*t6*t7*t10+ga1*nu1
     +*t6*t7*t10-k2ga2*k2nu2*t3*t6*t10+ga1*k2ga2*k2nu2*nu1*t6*t10)-amue1
     +*(t2*(t210+t211+t212+t213+t214-amue2*ga1*nu2*xi2-amue2*ga2*nu1*xi2
     +-amue2*t3*t10*xi2*2.0D0+amue2*ga1*k*k2nu2*t10-amue2*ga2*k*k2nu2*t1
     +0-amue2*k*k2ga2*nu1*t10-amue2*k*k2ga2*nu2*t10-amue2*ga1*nu2*t10*xi
     +2+amue2*ga2*nu1*t10*xi2)+k2ga1*(xi1*(t207+t208+t209-amue2*k2nu2*t3
     +-amue2*k2nu2*t3*t10+amue2*ga2*k2nu2*nu1*t10+amue2*k*nu1*t10*xi2+am
     +ue2*k*nu2*t10*xi2)-enuh1*k2nu1*t289)+k2nu1*xi1*(t215+t216+t217-enu
     +h1*(t296+t297)-amue2*k2ga2*t3+amue2*k2ga2*t3*t10+amue2*ga1*k2ga2*n
     +u2*t10+amue2*ga1*k*t10*xi2-amue2*ga2*k*t10*xi2))+t8*(t2*xi1*(t3+t5
     +-t3*t10-ga2*nu2*t10)+k2ga1*k2nu1*xi1*(t3+t5+t227+t228-enuh1*(t290+
     +t291)))-enuh1*k2nu1*t294)
      A0(2,3) = -t111*(amue2*(amue2*enuh1*k*t2*t7*(0.0D0,-1.0D0)-amue2*e
     +nuh1*k*k2ga2*k2nu2*t2*(0.0D0,1.0D0)+amue2*enuh1*k*k2ga1*k2nu1*t7*(
     +0.0D0,1.0D0)-amue2*egah1*k2ga1*nu1*t7*xi1*(0.0D0,2.0D0)+amue2*enuh
     +1*k*t2*t7*t10*(0.0D0,1.0D0)+amue2*enuh1*k*k2ga1*k2ga2*k2nu1*k2nu2*
     +(0.0D0,1.0D0)-amue2*egah1*k2ga1*k2ga2*k2nu2*nu1*xi1*(0.0D0,2.0D0)+
     +amue2*enuh1*k*k2ga2*k2nu2*t2*t10*(0.0D0,1.0D0)+amue2*enuh1*k*k2ga1
     +*k2nu1*t7*t10*(0.0D0,1.0D0)+amue2*enuh1*k*k2ga1*k2ga2*k2nu1*k2nu2*
     +t10*(0.0D0,1.0D0))-amue1*amue2*(egah1*k*k2ga1*k2nu2*t2*(0.0D0,-2.0
     +D0)+enuh1*k*k2ga1*k2nu2*t2*(0.0D0,1.0D0)-enuh1*k*k2nu1*k2nu2*t11*(
     +0.0D0,1.0D0)+egah1*k2ga1*nu2*t2*xi2*(0.0D0,2.0D0)-enuh1*k2ga1*nu2*
     +t2*xi2*(0.0D0,1.0D0)-enuh1*k2ga2*nu2*t2*xi1*(0.0D0,1.0D0)+enuh1*k2
     +nu1*nu2*t11*xi2*(0.0D0,1.0D0)-enuh1*k*t2*xi1*xi2*(0.0D0,1.0D0)+enu
     +h1*k*t2*t10*xi1*xi2*(0.0D0,1.0D0)-egah1*k2ga1*k2ga2*k2nu1*nu2*xi1*
     +(0.0D0,2.0D0)+enuh1*k2ga1*k2ga2*k2nu1*nu2*xi1*(0.0D0,1.0D0)+enuh1*
     +k*k2ga1*k2nu2*t2*t10*(0.0D0,1.0D0)+enuh1*k*k2nu1*k2nu2*t10*t11*(0.
     +0D0,1.0D0)-egah1*k*k2ga1*k2nu1*xi1*xi2*(0.0D0,2.0D0)+enuh1*k*k2ga1
     +*k2nu1*xi1*xi2*(0.0D0,1.0D0)-enuh1*k2ga1*nu2*t2*t10*xi2*(0.0D0,1.0
     +D0)+enuh1*k2ga2*nu2*t2*t10*xi1*(0.0D0,1.0D0)-enuh1*k2nu1*nu2*t10*t
     +11*xi2*(0.0D0,1.0D0)+enuh1*k2ga1*k2ga2*k2nu1*nu2*t10*xi1*(0.0D0,1.
     +0D0)+enuh1*k*k2ga1*k2nu1*t10*xi1*xi2*(0.0D0,1.0D0)))
      A0(2,4) = -t111*(k2nu1*(amue2*(amue2*enuh1*ga1*k2ga1*t7*(0.0D0,1.0
     +D0)+amue2*enuh1*ga1*k2ga1*k2ga2*k2nu2*(0.0D0,1.0D0)-amue2*enuh1*ga
     +1*k2ga1*t7*t10*(0.0D0,1.0D0)-amue2*enuh1*ga1*k2ga1*k2ga2*k2nu2*t10
     +*(0.0D0,1.0D0))-amue1*amue2*(enuh1*ga2*k2nu2*t11*(0.0D0,1.0D0)+enu
     +h1*k*t11*xi2*(0.0D0,1.0D0)-enuh1*ga2*k2nu2*t10*t11*(0.0D0,1.0D0)-e
     +nuh1*k*t10*t11*xi2*(0.0D0,1.0D0)))+xi1*(t308+amue1*amue2*k2nu1*(t3
     +09+t310-enuh1*k*k2ga1*k2ga2*(0.0D0,1.0D0)+enuh1*ga2*k2ga1*xi2*(0.0
     +D0,1.0D0)-enuh1*k*k2ga1*k2ga2*t10*(0.0D0,1.0D0)+enuh1*ga2*k2ga1*t1
     +0*xi2*(0.0D0,1.0D0)))-t2*(amue2*(amue2*enuh1*ga1*t7*(0.0D0,1.0D0)+
     +amue2*enuh1*ga1*k2ga2*k2nu2*(0.0D0,1.0D0)+amue2*enuh1*ga1*t7*t10*(
     +0.0D0,1.0D0)+amue2*enuh1*ga1*k2ga2*k2nu2*t10*(0.0D0,1.0D0))-amue1*
     +amue2*(egah1*ga2*k2ga1*k2nu2*(0.0D0,-2.0D0)+enuh1*ga2*k2ga1*k2nu2*
     +(0.0D0,1.0D0)-egah1*k*k2ga1*xi2*(0.0D0,2.0D0)+enuh1*k*k2ga1*xi2*(0
     +.0D0,1.0D0)+enuh1*ga2*k2ga1*k2nu2*t10*(0.0D0,1.0D0)+enuh1*k*k2ga1*
     +t10*xi2*(0.0D0,1.0D0)))+amue1*amue2*t2*xi1*(enuh1*k*k2ga2*(0.0D0,1
     +.0D0)-enuh1*ga2*xi2*(0.0D0,1.0D0)-enuh1*k*k2ga2*t10*(0.0D0,1.0D0)+
     +enuh1*ga2*t10*xi2*(0.0D0,1.0D0)))
      A0(2,5) = -t111*(enuh1*(t322+t2*(t330-amue1*(t222+t223+t331+t332))
     +)+amue1*t2*(t358+egah1*ga2*k2ga1*nu2*2.0D0)-amue2*xi1*(egah1*k2ga1
     +*t3*xi2*2.0D0+egah1*k*k2ga1*k2ga2*nu1*2.0D0+egah1*k*k2ga1*k2ga2*nu
     +2*2.0D0-egah1*ga2*k2ga1*nu1*xi2*2.0D0))
      A0(2,6) = -t111*(enuh1*(t348+k2ga1*(t340-amue1*xi1*(t224+t225+t341
     ++t342))-amue2*t2*t346)+k2ga1*xi1*(amue1*(t360+egah1*ga2*k2nu1*nu2*
     +2.0D0)-amue2*(egah1*k2nu2*t3*(-2.0D0)+egah1*ga2*k2nu2*nu1*2.0D0+
     +egah1*k*nu1*xi2*2.0D0+egah1*k*nu2*xi2*2.0D0)))
      A0(3,1) = -t13*t111*(k2ga1*t246+t8*(k2ga1*(t239*xi1-egah1*t242*xi1
     +)+egah1*t2*t243*xi1)-amue1*(egah1*(t2*(t210+t211+t212-t213+t214+t2
     +29-t230+amue2*ga1*nu2*xi2-amue2*ga2*nu1*xi2-amue2*t3*t4*xi2*2.0D0+
     +amue2*ga1*k*k2nu2*t4-amue2*ga2*k*k2nu2*t4-amue2*k*k2ga2*nu2*t4-amu
     +e2*ga2*nu1*t4*xi2)-xi1*(amue2*k2ga2*k2nu1*t3+amue2*ga1*k2ga2*k2nu1
     +*nu2+amue2*ga1*k*k2nu1*xi2-amue2*ga2*k*k2nu1*xi2+amue2*k2ga2*k2nu1
     +*t3*t4+amue2*ga1*k2ga2*k2nu1*nu2*t4+amue2*ga1*k*k2nu1*t4*xi2-amue2
     +*ga2*k*k2nu1*t4*xi2))+k2ga1*(t233*xi1-egah1*xi1*(t207+t208+t209+t2
     +34+t235+t236-amue2*k2nu2*t3-amue2*k*nu2*t4*xi2)+amue2*enuh1*k2nu1*
     +t3*xi2*2.0D0+amue2*enuh1*k*k2ga2*k2nu1*nu2*2.0D0))-egah1*xi1*(-t21
     +8-t219+t220+t221+t247+t248+t3*t4*t6*t7+k2ga2*k2nu2*t3*t4*t6))
      A0(3,2) = t13*t111*(-egah1*k2ga1*t3*t8*t12+egah1*k2nu1*t2*t3*t8+eg
     +ah1*k2nu1*t3*t6*t7-enuh1*k2nu1*t2*t3*t8*2.0D0+amue1*amue2*egah1*k2
     +ga2*t3*t12-egah1*ga2*k2ga1*nu2*t8*t12+egah1*ga2*k2nu1*nu2*t2*t8-eg
     +ah1*ga1*k2nu1*nu1*t6*t7-enuh1*ga2*k2nu1*nu2*t2*t8*2.0D0+egah1*k2ga
     +2*k2nu1*k2nu2*t3*t6-enuh1*k*nu1*t6*t7*xi1*2.0D0+egah1*k2ga1*t3*t4*
     +t8*t12+egah1*k2nu1*t2*t3*t4*t8+egah1*k2nu1*t3*t4*t6*t7+amue1*amue2
     +*egah1*ga1*k2ga2*nu2*t12+amue1*amue2*enuh1*ga2*k2nu2*nu1*t2*2.0D0-
     +amue1*amue2*egah1*k2ga1*k2nu1*k2nu2*t3+amue1*amue2*egah1*ga1*k*t12
     +*xi2-amue1*amue2*egah1*ga2*k*t12*xi2+amue1*amue2*enuh1*k*nu1*t2*xi
     +2*2.0D0-amue1*amue2*egah1*k2ga2*t3*t4*t12-egah1*ga1*k2ga2*k2nu1*k2
     +nu2*nu1*t6-amue1*amue2*egah1*k2nu1*t3*xi1*xi2*2.0D0+amue1*amue2*en
     +uh1*k2nu1*t3*xi1*xi2*2.0D0-enuh1*k*k2ga2*k2nu2*nu1*t6*xi1*2.0D0+eg
     +ah1*ga2*k2ga1*nu2*t4*t8*t12+egah1*ga2*k2nu1*nu2*t2*t4*t8+egah1*ga1
     +*k2nu1*nu1*t4*t6*t7+egah1*k2ga2*k2nu1*k2nu2*t3*t4*t6+amue1*amue2*e
     +gah1*ga2*k2ga1*k2nu1*k2nu2*nu1+amue1*amue2*egah1*ga1*k*k2nu1*k2nu2
     +*xi1-amue1*amue2*egah1*ga2*k*k2nu1*k2nu2*xi1+amue1*amue2*egah1*k*k
     +2ga1*k2nu1*nu1*xi2-amue1*amue2*egah1*k*k2ga2*k2nu1*nu1*xi1+amue1*a
     +mue2*egah1*k*k2ga1*k2nu1*nu2*xi2-amue1*amue2*egah1*k*k2ga2*k2nu1*n
     +u2*xi1+amue1*amue2*enuh1*k*k2ga2*k2nu1*nu2*xi1*2.0D0-amue1*amue2*e
     +gah1*ga1*k2ga2*nu2*t4*t12-amue1*amue2*egah1*k2ga1*k2nu1*k2nu2*t3*t
     +4-amue1*amue2*egah1*ga1*k2nu1*nu2*xi1*xi2+amue1*amue2*egah1*ga2*k2
     +nu1*nu1*xi1*xi2-amue1*amue2*egah1*ga1*k*t4*t12*xi2+amue1*amue2*ega
     +h1*ga2*k*t4*t12*xi2+egah1*ga1*k2ga2*k2nu1*k2nu2*nu1*t4*t6-amue1*am
     +ue2*egah1*k2nu1*t3*t4*xi1*xi2*2.0D0-amue1*amue2*egah1*ga2*k2ga1*k2
     +nu1*k2nu2*nu1*t4+amue1*amue2*egah1*ga1*k*k2nu1*k2nu2*t4*xi1-amue1*
     +amue2*egah1*ga2*k*k2nu1*k2nu2*t4*xi1-amue1*amue2*egah1*k*k2ga1*k2n
     +u1*nu1*t4*xi2+amue1*amue2*egah1*k*k2ga2*k2nu1*nu1*t4*xi1+amue1*amu
     +e2*egah1*k*k2ga1*k2nu1*nu2*t4*xi2-amue1*amue2*egah1*k*k2ga2*k2nu1*
     +nu2*t4*xi1-amue1*amue2*egah1*ga1*k2nu1*nu2*t4*xi1*xi2-amue1*amue2*
     +egah1*ga2*k2nu1*nu1*t4*xi1*xi2)
      A0(3,3) = -t111*(-k2ga1*(amue2*(amue2*k2nu1*nu1*t7*(0.0D0,1.0D0)+a
     +mue2*k2ga2*k2nu1*k2nu2*nu1*(0.0D0,1.0D0)-amue2*k2nu1*nu1*t4*t7*(0.
     +0D0,1.0D0)-amue2*k2ga2*k2nu1*k2nu2*nu1*t4*(0.0D0,1.0D0))-amue1*amu
     +e2*(k2ga2*nu2*t12*(0.0D0,1.0D0)+k*t12*xi2*(0.0D0,1.0D0)-k2ga2*nu2*
     +t4*t12*(0.0D0,1.0D0)-k*t4*t12*xi2*(0.0D0,1.0D0)))-t2*(amue2*(amue2
     +*nu1*t7*(0.0D0,1.0D0)+amue2*k2ga2*k2nu2*nu1*(0.0D0,1.0D0)+amue2*nu
     +1*t4*t7*(0.0D0,1.0D0)+amue2*k2ga2*k2nu2*nu1*t4*(0.0D0,1.0D0))-amue
     +1*amue2*(k2ga2*k2nu1*nu2*(0.0D0,1.0D0)+k*k2nu1*xi2*(0.0D0,1.0D0)+k
     +2ga2*k2nu1*nu2*t4*(0.0D0,1.0D0)+k*k2nu1*t4*xi2*(0.0D0,1.0D0)))+ega
     +h1*(xi1*(t252-amue1*amue2*k2ga1*(enuh1*k*k2nu1*k2nu2*(0.0D0,2.0D0)
     +-enuh1*k2nu1*nu2*xi2*(0.0D0,2.0D0)))-amue1*amue2*t2*(enuh1*k2ga2*k
     +2nu1*nu2*(0.0D0,2.0D0)+enuh1*k*k2nu1*xi2*(0.0D0,2.0D0)))+amue1*amu
     +e2*t2*xi1*(k*k2nu2*(0.0D0,1.0D0)-nu2*xi2*(0.0D0,1.0D0)-k*k2nu2*t4*
     +(0.0D0,1.0D0)+nu2*t4*xi2*(0.0D0,1.0D0))+amue1*amue2*k2ga1*xi1*(k*k
     +2nu1*k2nu2*(0.0D0,1.0D0)-k2nu1*nu2*xi2*(0.0D0,1.0D0)+k*k2nu1*k2nu2
     +*t4*(0.0D0,1.0D0)-k2nu1*nu2*t4*xi2*(0.0D0,1.0D0)))
      A0(3,4) = t111*(k2ga1*(amue2*(amue1*ga2*t12*xi2*(0.0D0,1.0D0)-amue
     +1*ga2*t4*t12*xi2*(0.0D0,1.0D0))-amue2*k*(amue1*k2ga2*t12*(0.0D0,1.
     +0D0)+amue2*k2nu1*t7*(0.0D0,1.0D0)+amue2*k2ga2*k2nu1*k2nu2*(0.0D0,1
     +.0D0)-amue1*k2ga2*t4*t12*(0.0D0,1.0D0)+amue2*k2nu1*t4*t7*(0.0D0,1.
     +0D0)+amue2*k2ga2*k2nu1*k2nu2*t4*(0.0D0,1.0D0)))+xi1*(amue2*(amue2*
     +egah1*enuh1*ga1*k2nu1*t7*(0.0D0,2.0D0)+amue2*egah1*enuh1*ga1*k2ga2
     +*k2nu1*k2nu2*(0.0D0,2.0D0))+k2ga1*(amue2*(amue1*ga2*k2nu1*k2nu2*(0
     +.0D0,1.0D0)+amue1*ga2*k2nu1*k2nu2*t4*(0.0D0,1.0D0)-amue1*egah1*enu
     +h1*ga2*k2nu1*k2nu2*(0.0D0,2.0D0))+amue2*k*(amue1*k2nu1*xi2*(0.0D0,
     +1.0D0)+amue1*k2nu1*t4*xi2*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2nu1*xi
     +2*(0.0D0,2.0D0))))+t2*(amue2*(amue1*ga2*k2nu1*xi2*(0.0D0,1.0D0)+am
     +ue1*ga2*k2nu1*t4*xi2*(0.0D0,1.0D0)-amue1*egah1*enuh1*ga2*k2nu1*xi2
     +*(0.0D0,2.0D0))-amue2*k*(t302+t303+amue1*k2ga2*k2nu1*(0.0D0,1.0D0)
     +-amue2*t4*t7*(0.0D0,1.0D0)+amue1*k2ga2*k2nu1*t4*(0.0D0,1.0D0)-amue
     +2*k2ga2*k2nu2*t4*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2ga2*k2nu1*(0.0D
     +0,2.0D0)))+t2*xi1*(amue2*(amue1*ga2*k2nu2*(0.0D0,1.0D0)-amue1*ga2*
     +k2nu2*t4*(0.0D0,1.0D0))+amue2*k*(t304-amue1*t4*xi2*(0.0D0,1.0D0)))
     +)
      A0(3,5) = t111*(-t267+k2nu1*(t260-xi1*(amue1*(t222+t223+t261+t262-
     +egah1*enuh1*k2ga1*t3*2.0D0-egah1*enuh1*ga2*k2ga1*nu2*2.0D0)+amue2*
     +(egah1*enuh1*k2ga2*t3*2.0D0+egah1*enuh1*ga1*k2ga2*nu2*2.0D0+egah1*
     +enuh1*ga1*k*xi2*2.0D0-egah1*enuh1*ga2*k*xi2*2.0D0)))+amue2*t2*t266
     +)
      A0(3,6) = -t111*(t279-t2*(amue2*t286-amue1*(t224+t225+t240+t241-eg
     +ah1*enuh1*k2nu1*t3*2.0D0-egah1*enuh1*ga2*k2nu1*nu2*2.0D0))+amue2*x
     +i1*(egah1*enuh1*k2nu1*t3*xi2*2.0D0-egah1*enuh1*ga1*k*k2nu1*k2nu2*2
     +.0D0+egah1*enuh1*ga2*k*k2nu1*k2nu2*2.0D0+egah1*enuh1*ga1*k2nu1*nu2
     +*xi2*2.0D0))
      A0(4,1) = t13*t111*(egah1*k2ga1*t2*t3*t8*(-2.0D0)+
     +enuh1*k2ga1*t2*t3*
     +t8+enuh1*k2ga1*t3*t6*t7-enuh1*k2nu1*t3*t8*t11+amue1*amue2*enuh1*k2
     +nu2*t3*t11-egah1*ga2*k2ga1*nu2*t2*t8*2.0D0+enuh1*ga2*k2ga1*nu2*t2*
     +t8-enuh1*ga1*k2ga1*nu1*t6*t7-enuh1*ga2*k2nu1*nu2*t8*t11+enuh1*k2ga
     +1*k2ga2*k2nu2*t3*t6-egah1*ga1*k*t6*t7*xi1*2.0D0+enuh1*k2ga1*t2*t3*
     +t8*t10+enuh1*k2ga1*t3*t6*t7*t10+enuh1*k2nu1*t3*t8*t10*t11+amue1*am
     +ue2*egah1*ga1*k2ga2*nu2*t2*2.0D0+amue1*amue2*enuh1*ga2*k2nu2*nu1*t
     +11-amue1*amue2*enuh1*k2ga1*k2ga2*k2nu1*t3+amue1*amue2*egah1*ga1*k*
     +t2*xi2*2.0D0+amue1*amue2*enuh1*k*nu1*t11*xi2-amue1*amue2*enuh1*k*n
     +u2*t11*xi2-amue1*amue2*enuh1*k2nu2*t3*t10*t11-enuh1*ga1*k2ga1*k2ga
     +2*k2nu2*nu1*t6+amue1*amue2*egah1*k2ga1*t3*xi1*xi2*2.0D0-amue1*amue
     +2*enuh1*k2ga1*t3*xi1*xi2*2.0D0-egah1*ga1*k*k2ga2*k2nu2*t6*xi1*2.0D
     +0+enuh1*ga2*k2ga1*nu2*t2*t8*t10+enuh1*ga1*k2ga1*nu1*t6*t7*t10+enuh
     +1*ga2*k2nu1*nu2*t8*t10*t11+enuh1*k2ga1*k2ga2*k2nu2*t3*t6*t10+amue1
     +*amue2*enuh1*ga1*k2ga1*k2ga2*k2nu1*nu2+amue1*amue2*egah1*ga2*k*k2g
     +a1*k2nu2*xi1*2.0D0+amue1*amue2*enuh1*ga1*k*k2ga1*k2nu1*xi2-amue1*a
     +mue2*enuh1*ga1*k*k2ga1*k2nu2*xi1+amue1*amue2*enuh1*ga2*k*k2ga1*k2n
     +u1*xi2-amue1*amue2*enuh1*ga2*k*k2ga1*k2nu2*xi1+amue1*amue2*enuh1*k
     +*k2ga1*k2ga2*nu1*xi1-amue1*amue2*enuh1*k*k2ga1*k2ga2*nu2*xi1-amue1
     +*amue2*enuh1*ga2*k2nu2*nu1*t10*t11-amue1*amue2*enuh1*k2ga1*k2ga2*k
     +2nu1*t3*t10+amue1*amue2*enuh1*ga1*k2ga1*nu2*xi1*xi2-amue1*amue2*en
     +uh1*ga2*k2ga1*nu1*xi1*xi2-amue1*amue2*enuh1*k*nu1*t10*t11*xi2+amue
     +1*amue2*enuh1*k*nu2*t10*t11*xi2+enuh1*ga1*k2ga1*k2ga2*k2nu2*nu1*t6
     +*t10-amue1*amue2*enuh1*k2ga1*t3*t10*xi1*xi2*2.0D0-amue1*amue2*enuh
     +1*ga1*k2ga1*k2ga2*k2nu1*nu2*t10-amue1*amue2*enuh1*ga1*k*k2ga1*k2nu
     +1*t10*xi2+amue1*amue2*enuh1*ga1*k*k2ga1*k2nu2*t10*xi1+amue1*amue2*
     +enuh1*ga2*k*k2ga1*k2nu1*t10*xi2-amue1*amue2*enuh1*ga2*k*k2ga1*k2nu
     +2*t10*xi1+amue1*amue2*enuh1*k*k2ga1*k2ga2*nu1*t10*xi1-amue1*amue2*
     +enuh1*k*k2ga1*k2ga2*nu2*t10*xi1-amue1*amue2*enuh1*ga1*k2ga1*nu2*t1
     +0*xi1*xi2-amue1*amue2*enuh1*ga2*k2ga1*nu1*t10*xi1*xi2)
      A0(4,2) = -t13*t111*(xi1*(-enuh1*t3*t6*t7+enuh1*ga1*nu1*t6*t7-enuh
     +1*k2ga2*k2nu2*t3*t6+enuh1*t3*t6*t7*t10+enuh1*ga1*k2ga2*k2nu2*nu1*t
     +6+enuh1*ga1*nu1*t6*t7*t10+enuh1*k2ga2*k2nu2*t3*t6*t10+enuh1*ga1*k2
     +ga2*k2nu2*nu1*t6*t10)-k2nu1*t294+amue1*(t2*(t295+amue2*enuh1*ga1*k
     +*k2nu2+amue2*enuh1*ga2*k*k2nu2-amue2*enuh1*k*k2ga2*nu1+amue2*enuh1
     +*k*k2ga2*nu2-amue2*enuh1*ga1*nu2*xi2+amue2*enuh1*ga2*nu1*xi2-amue2
     +*enuh1*t3*t10*xi2*2.0D0+amue2*enuh1*ga1*k*k2nu2*t10-amue2*enuh1*ga
     +2*k*k2nu2*t10+amue2*enuh1*k*k2ga2*nu1*t10-amue2*enuh1*k*k2ga2*nu2*
     +t10-amue2*enuh1*ga1*nu2*t10*xi2-amue2*enuh1*ga2*nu1*t10*xi2)-k2nu1
     +*xi1*(-t296-t297-amue2*enuh1*k2ga2*t3+amue2*enuh1*ga1*k2ga2*nu2+am
     +ue2*enuh1*ga1*k*xi2+amue2*enuh1*ga2*k*xi2+amue2*enuh1*k2ga2*t3*t10
     ++amue2*enuh1*ga1*k2ga2*nu2*t10+amue2*enuh1*ga1*k*t10*xi2-amue2*enu
     +h1*ga2*k*t10*xi2))+k2ga1*(amue1*(k2nu1*t289-xi1*(amue2*enuh1*k2nu2
     +*t3+amue2*enuh1*ga2*k2nu2*nu1+amue2*enuh1*k*nu1*xi2-amue2*enuh1*k*
     +nu2*xi2+amue2*enuh1*k2nu2*t3*t10+amue2*enuh1*ga2*k2nu2*nu1*t10+amu
     +e2*enuh1*k*nu1*t10*xi2-amue2*enuh1*k*nu2*t10*xi2))+k2nu1*t8*xi1*(-
     +t290-t291+t298+t299+t300+t301))-t2*t8*xi1*(t298-t299+t300-t301))
      A0(4,3) = t111*(k2nu1*(amue2*(amue1*nu2*t11*xi2*(0.0D0,1.0D0)-amue
     +1*nu2*t10*t11*xi2*(0.0D0,1.0D0))-amue2*k*(amue2*k2ga1*t7*(0.0D0,1.
     +0D0)+amue1*k2nu2*t11*(0.0D0,1.0D0)+amue2*k2ga1*k2ga2*k2nu2*(0.0D0,
     +1.0D0)+amue2*k2ga1*t7*t10*(0.0D0,1.0D0)-amue1*k2nu2*t10*t11*(0.0D0
     +,1.0D0)+amue2*k2ga1*k2ga2*k2nu2*t10*(0.0D0,1.0D0)))+xi1*(amue2*(am
     +ue2*egah1*enuh1*k2ga1*nu1*t7*(0.0D0,2.0D0)+amue2*egah1*enuh1*k2ga1
     +*k2ga2*k2nu2*nu1*(0.0D0,2.0D0))+k2nu1*(amue2*(amue1*k2ga1*k2ga2*nu
     +2*(0.0D0,1.0D0)+amue1*k2ga1*k2ga2*nu2*t10*(0.0D0,1.0D0)-amue1*egah
     +1*enuh1*k2ga1*k2ga2*nu2*(0.0D0,2.0D0))+amue2*k*(amue1*k2ga1*xi2*(0
     +.0D0,1.0D0)+amue1*k2ga1*t10*xi2*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2
     +ga1*xi2*(0.0D0,2.0D0))))+t2*(amue2*(amue1*k2ga1*nu2*xi2*(0.0D0,1.0
     +D0)+amue1*k2ga1*nu2*t10*xi2*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2ga1*
     +nu2*xi2*(0.0D0,2.0D0))-amue2*k*(t302+t303+amue1*k2ga1*k2nu2*(0.0D0
     +,1.0D0)-amue2*t7*t10*(0.0D0,1.0D0)+amue1*k2ga1*k2nu2*t10*(0.0D0,1.
     +0D0)-amue2*k2ga2*k2nu2*t10*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2ga1*k
     +2nu2*(0.0D0,2.0D0)))+t2*xi1*(amue2*(amue1*k2ga2*nu2*(0.0D0,1.0D0)-
     +amue1*k2ga2*nu2*t10*(0.0D0,1.0D0))+amue2*k*(t304-amue1*t10*xi2*(0.
     +0D0,1.0D0))))
      A0(4,4) = t111*(-k2nu1*(amue2*(amue2*ga1*k2ga1*t7*(0.0D0,1.0D0)+am
     +ue2*ga1*k2ga1*k2ga2*k2nu2*(0.0D0,1.0D0)-amue2*ga1*k2ga1*t7*t10*(0.
     +0D0,1.0D0)-amue2*ga1*k2ga1*k2ga2*k2nu2*t10*(0.0D0,1.0D0))-amue1*am
     +ue2*(ga2*k2nu2*t11*(0.0D0,1.0D0)+k*t11*xi2*(0.0D0,1.0D0)-ga2*k2nu2
     +*t10*t11*(0.0D0,1.0D0)-k*t10*t11*xi2*(0.0D0,1.0D0)))-t2*(amue2*(am
     +ue2*ga1*t7*(0.0D0,1.0D0)+amue2*ga1*t7*t10*(0.0D0,1.0D0)+amue2*ga1*
     +k2ga2*k2nu2*(0.0D0,1.0D0)+amue2*ga1*k2ga2*k2nu2*t10*(0.0D0,1.0D0))
     +-amue1*amue2*(ga2*k2ga1*k2nu2*(0.0D0,1.0D0)+k*k2ga1*xi2*(0.0D0,1.0
     +D0)+ga2*k2ga1*k2nu2*t10*(0.0D0,1.0D0)+k*k2ga1*t10*xi2*(0.0D0,1.0D0
     +)))+enuh1*(xi1*(t308-amue1*amue2*k2nu1*(t309+t310))-amue1*amue2*t2
     +*(egah1*ga2*k2ga1*k2nu2*(0.0D0,2.0D0)+egah1*k*k2ga1*xi2*(0.0D0,2.0
     +D0)))+amue1*amue2*t2*xi1*(k*k2ga2*(0.0D0,1.0D0)-ga2*xi2*(0.0D0,1.0
     +D0)-k*k2ga2*t10*(0.0D0,1.0D0)+ga2*t10*xi2*(0.0D0,1.0D0))+amue1*amu
     +e2*k2nu1*xi1*(k*k2ga1*k2ga2*(0.0D0,1.0D0)-ga2*k2ga1*xi2*(0.0D0,1.0
     +D0)+k*k2ga1*k2ga2*t10*(0.0D0,1.0D0)-ga2*k2ga1*t10*xi2*(0.0D0,1.0D0
     +)))
      A0(4,5) = -t111*(t322-t2*(t330-amue1*(t222+t223+t331+t332-egah1*en
     +uh1*k2ga1*t3*2.0D0-egah1*enuh1*ga2*k2ga1*nu2*2.0D0))+amue2*xi1*(eg
     +ah1*enuh1*k2ga1*t3*xi2*2.0D0-egah1*enuh1*k*k2ga1*k2ga2*nu1*2.0D0+e
     +gah1*enuh1*k*k2ga1*k2ga2*nu2*2.0D0+egah1*enuh1*ga2*k2ga1*nu1*xi2*2
     +.0D0))
      A0(4,6) = -t111*(-t348+k2ga1*(t340-xi1*(amue1*(t224+t225+t341+t342
     +-egah1*enuh1*k2nu1*t3*2.0D0-egah1*enuh1*ga2*k2nu1*nu2*2.0D0)+amue2
     +*(egah1*enuh1*k2nu2*t3*2.0D0+egah1*enuh1*ga2*k2nu2*nu1*2.0D0+egah1
     +*enuh1*k*nu1*xi2*2.0D0-egah1*enuh1*k*nu2*xi2*2.0D0)))+amue2*t2*t34
     +6)
      A0(5,1) = t111*(-t2*(amue1*egah1*k*k2ga1*nu1*2.0D0+amue1*egah1*k*k
     +2ga1*nu2*2.0D0-amue1*enuh1*k*k2ga1*nu1*2.0D0+amue1*egah1*k*k2ga1*n
     +u1*t4*2.0D0-amue1*egah1*k*k2ga1*nu2*t4*2.0D0-amue1*enuh1*k*k2ga1*n
     +u1*t10*2.0D0)+k2nu1*(xi1*(t354+t355-amue1*enuh1*k2ga1*t3*2.0D0+amu
     +e1*enuh1*ga1*k2ga1*nu2*2.0D0-amue1*enuh1*k2ga1*t3*t10*2.0D0-amue1*
     +enuh1*ga1*k2ga1*nu2*t10*2.0D0)-amue1*enuh1*k*nu2*t11*2.0D0+amue1*e
     +nuh1*k*nu2*t10*t11*2.0D0)+amue2*(k2nu1*(xi2*(t349+enuh1*k2ga1*t3*t
     +10*2.0D0)-xi1*xi2*(t356+t357))+t2*(xi2*(egah1*ga1*nu1*2.0D0+egah1*
     +ga1*nu1*t4*2.0D0)-egah1*ga1*k*k2nu2*2.0D0+egah1*ga1*k*k2nu2*t4*2.0
     +D0)-xi1*(xi2*(enuh1*k*k2ga1*nu1*2.0D0+enuh1*k*k2ga1*nu1*t10*2.0D0)
     +-egah1*k2ga1*k2nu2*t3*2.0D0+enuh1*ga1*k2ga1*k2nu2*nu1*2.0D0+egah1*
     +k2ga1*k2nu2*t3*t4*2.0D0-enuh1*ga1*k2ga1*k2nu2*nu1*t10*2.0D0)+enuh1
     +*k*k2nu2*nu1*t11*2.0D0-enuh1*k*k2nu2*nu1*t10*t11*2.0D0)+t2*xi1*(am
     +ue1*egah1*ga1*nu2*2.0D0-amue1*egah1*ga1*nu2*t4*2.0D0))
      A0(5,2) = -t111*(amue2*(t2*(xi2*(t361-enuh1*k*nu1*t10*2.0D0)+enuh1
     +*ga1*k2nu2*nu1*2.0D0+enuh1*ga1*k2nu2*nu1*t10*2.0D0)+k2ga1*(-xi1*(e
     +nuh1*k*k2nu2*nu1*2.0D0+enuh1*k*k2nu2*nu1*t10*2.0D0)+egah1*k2nu1*k2
     +nu2*t3*2.0D0+egah1*k2nu1*k2nu2*t3*t4*2.0D0)-xi2*(egah1*ga1*k*t12*2
     +.0D0-egah1*ga1*k*t4*t12*2.0D0)-xi1*(xi2*(t237-egah1*ga1*k2nu1*nu1*
     +2.0D0-enuh1*k2nu1*t3*t10*2.0D0+egah1*ga1*k2nu1*nu1*t4*2.0D0)+egah1
     +*ga1*k*k2nu1*k2nu2*2.0D0+egah1*ga1*k*k2nu1*k2nu2*t4*2.0D0))+t2*(t3
     +59+amue1*egah1*ga1*k2nu1*nu2*2.0D0-amue1*enuh1*ga1*k2nu1*nu2*2.0D0
     +-amue1*enuh1*k2nu1*t3*t10*2.0D0+amue1*egah1*ga1*k2nu1*nu2*t4*2.0D0
     +-amue1*enuh1*ga1*k2nu1*nu2*t10*2.0D0)-k2ga1*(xi1*(amue1*egah1*k*k2
     +nu1*nu1*2.0D0+amue1*egah1*k*k2nu1*nu2*2.0D0-amue1*enuh1*k*k2nu1*nu
     +2*2.0D0-amue1*egah1*k*k2nu1*nu1*t4*2.0D0+amue1*egah1*k*k2nu1*nu2*t
     +4*2.0D0-amue1*enuh1*k*k2nu1*nu2*t10*2.0D0)-amue1*egah1*t3*t12*2.0D
     +0+amue1*egah1*t3*t4*t12*2.0D0)-t2*xi1*(amue1*enuh1*k*nu1*2.0D0-amu
     +e1*enuh1*k*nu1*t10*2.0D0))
      A0(5,3) = t111*(amue1*(amue1*nu2*t9*(0.0D0,1.0D0)-amue1*nu2*t4*t9*
     +(0.0D0,1.0D0)-amue1*nu2*t9*t10*(0.0D0,1.0D0)+amue1*nu2*t11*t12*(0.
     +0D0,1.0D0)+amue1*k2ga1*k2nu1*nu2*t2*(0.0D0,2.0D0)+amue1*nu2*t4*t9*
     +t10*(0.0D0,1.0D0)-amue1*nu2*t4*t11*t12*(0.0D0,1.0D0)-amue1*nu2*t10
     +*t11*t12*(0.0D0,1.0D0)+amue1*k2ga1*k2nu1*nu2*t2*t4*(0.0D0,2.0D0)+a
     +mue1*k2ga1*k2nu1*nu2*t2*t10*(0.0D0,2.0D0)+amue1*nu2*t4*t10*t11*t12
     +*(0.0D0,1.0D0)-amue1*egah1*enuh1*k2ga1*k2nu1*nu2*t2*(0.0D0,8.0D0)+
     +amue1*k2ga1*k2nu1*nu2*t2*t4*t10*(0.0D0,2.0D0))-amue2*(amue1*(k2ga1
     +*k2nu2*nu1*t2*(0.0D0,1.0D0)+k2nu1*k2nu2*nu1*t11*(0.0D0,1.0D0)-nu1*
     +t2*xi1*xi2*(0.0D0,1.0D0)+k2ga1*k2nu2*nu1*t2*t4*(0.0D0,1.0D0)+k2ga1
     +*k2nu2*nu1*t2*t10*(0.0D0,1.0D0)-k2nu1*k2nu2*nu1*t4*t11*(0.0D0,1.0D
     +0)-k2nu1*k2nu2*nu1*t10*t11*(0.0D0,1.0D0)-k2ga1*k2nu1*nu1*xi1*xi2*(
     +0.0D0,1.0D0)-nu1*t2*t4*xi1*xi2*(0.0D0,1.0D0)+nu1*t2*t10*xi1*xi2*(0
     +.0D0,1.0D0)+k2ga1*k2nu1*nu1*t4*xi1*xi2*(0.0D0,1.0D0)-k2ga1*k2nu1*n
     +u1*t10*xi1*xi2*(0.0D0,1.0D0)+nu1*t2*t4*t10*xi1*xi2*(0.0D0,1.0D0)-e
     +gah1*enuh1*k2ga1*k2nu2*nu1*t2*(0.0D0,4.0D0)+k2ga1*k2nu2*nu1*t2*t4*
     +t10*(0.0D0,1.0D0)+k2nu1*k2nu2*nu1*t4*t10*t11*(0.0D0,1.0D0)+k2ga1*k
     +2nu1*nu1*t4*t10*xi1*xi2*(0.0D0,1.0D0))+amue1*k*(k2ga1*t12*xi2*(0.0
     +D0,1.0D0)+k2nu1*t2*xi2*(0.0D0,1.0D0)+k2nu2*t2*xi1*(0.0D0,1.0D0)+k2
     +ga1*k2nu1*k2nu2*xi1*(0.0D0,1.0D0)-k2ga1*t4*t12*xi2*(0.0D0,1.0D0)+k
     +2ga1*t10*t12*xi2*(0.0D0,1.0D0)+k2nu1*t2*t4*xi2*(0.0D0,1.0D0)-k2nu2
     +*t2*t4*xi1*(0.0D0,1.0D0)-k2nu1*t2*t10*xi2*(0.0D0,1.0D0)-k2nu2*t2*t
     +10*xi1*(0.0D0,1.0D0)+k2ga1*k2nu1*k2nu2*t4*xi1*(0.0D0,1.0D0)+k2ga1*
     +k2nu1*k2nu2*t10*xi1*(0.0D0,1.0D0)-k2ga1*t4*t10*t12*xi2*(0.0D0,1.0D
     +0)-k2nu1*t2*t4*t10*xi2*(0.0D0,1.0D0)+k2nu2*t2*t4*t10*xi1*(0.0D0,1.
     +0D0)+k2ga1*k2nu1*k2nu2*t4*t10*xi1*(0.0D0,1.0D0)-egah1*enuh1*k2ga1*
     +k2nu1*k2nu2*xi1*(0.0D0,4.0D0))))
      A0(5,4) = t111*(k*(t376-amue1*amue2*(t377+t378+t379+t380+t381+t382
     ++t383+t384+t385-k2ga1*k2nu2*t2*(0.0D0,1.0D0)-k2nu1*k2nu2*t11*(0.0D
     +0,1.0D0)+k2ga1*k2nu2*t2*t4*(0.0D0,1.0D0)-k2ga1*k2nu2*t2*t10*(0.0D0
     +,1.0D0)-k2nu1*k2nu2*t4*t11*(0.0D0,1.0D0)+k2nu1*k2nu2*t10*t11*(0.0D
     +0,1.0D0)+k2ga1*k2nu2*t2*t4*t10*(0.0D0,1.0D0)+k2nu1*k2nu2*t4*t10*t1
     +1*(0.0D0,1.0D0)))-amue1*amue2*(ga1*k2ga1*t12*xi2*(0.0D0,1.0D0)+ga1
     +*k2nu1*t2*xi2*(0.0D0,1.0D0)+ga1*k2nu2*t2*xi1*(0.0D0,1.0D0)+ga1*k2g
     +a1*k2nu1*k2nu2*xi1*(0.0D0,1.0D0)-ga1*k2ga1*t4*t12*xi2*(0.0D0,1.0D0
     +)-ga1*k2ga1*t10*t12*xi2*(0.0D0,1.0D0)+ga1*k2nu1*t2*t4*xi2*(0.0D0,1
     +.0D0)-ga1*k2nu2*t2*t4*xi1*(0.0D0,1.0D0)+ga1*k2nu1*t2*t10*xi2*(0.0D
     +0,1.0D0)+ga1*k2nu2*t2*t10*xi1*(0.0D0,1.0D0)+ga1*k2ga1*t4*t10*t12*x
     +i2*(0.0D0,1.0D0)+ga1*k2nu1*t2*t4*t10*xi2*(0.0D0,1.0D0)-ga1*k2nu2*t
     +2*t4*t10*xi1*(0.0D0,1.0D0)-egah1*enuh1*ga1*k2nu1*t2*xi2*(0.0D0,4.0
     +D0)+ga1*k2ga1*k2nu1*k2nu2*t4*xi1*(0.0D0,1.0D0)-ga1*k2ga1*k2nu1*k2n
     +u2*t10*xi1*(0.0D0,1.0D0)-ga1*k2ga1*k2nu1*k2nu2*t4*t10*xi1*(0.0D0,1
     +.0D0)))
      A0(5,5) = -t111*(amue1*(k2nu1*(k*nu1*t11+k*nu2*t11-k*nu1*t4*t11+k*
     +nu2*t4*t11-k*nu1*t10*t11-k*nu2*t10*t11+k*nu1*t4*t10*t11-k*nu2*t4*t
     +10*t11)+t2*(k*k2ga1*nu1+k*k2ga1*nu2+k*k2ga1*nu1*t4-k*k2ga1*nu2*t4+
     +k*k2ga1*nu1*t10+k*k2ga1*nu2*t10-egah1*enuh1*k*k2ga1*nu1*4.0D0+k*k2
     +ga1*nu1*t4*t10-k*k2ga1*nu2*t4*t10)+k2nu1*xi1*(t222+t261+t331+t350-
     +ga1*k2ga1*nu2-egah1*enuh1*k2ga1*t3*4.0D0-ga1*k2ga1*nu2*t4+ga1*k2ga
     +1*nu2*t10+ga1*k2ga1*nu2*t4*t10)+t2*xi1*(t3-t205-t227+t351-ga1*nu2+
     +ga1*nu2*t4-ga1*nu2*t10+ga1*nu2*t4*t10))-amue2*xi2*(t397+k2nu1*(t22
     +2+t261+t331+t350+ga1*k2ga1*nu1-ga1*k2ga1*nu1*t4-ga1*k2ga1*nu1*t10+
     +ga1*k2ga1*nu1*t4*t10)-t391*xi1))
      A0(5,6) = t111*(amue1*(t2*(-t224-t240+t341+t389+ga1*k2nu1*nu2+ga1*
     +k2nu1*nu2*t4+ga1*k2nu1*nu2*t10-egah1*enuh1*ga1*k2nu1*nu2*4.0D0+ga1
     +*k2nu1*nu2*t4*t10)-k2ga1*(t268-t353-ga1*nu2*t12+t3*t10*t12+ga1*nu2
     +*t4*t12+ga1*nu2*t10*t12-t3*t4*t10*t12-ga1*nu2*t4*t10*t12)+t2*xi1*(
     +k*nu1+k*nu2+k*nu1*t4-k*nu2*t4-k*nu1*t10-k*nu2*t10-k*nu1*t4*t10+k*n
     +u2*t4*t10)+k2ga1*xi1*(k*k2nu1*nu1+k*k2nu1*nu2-k*k2nu1*nu1*t4+k*k2n
     +u1*nu2*t4+k*k2nu1*nu1*t10+k*k2nu1*nu2*t10-egah1*enuh1*k*k2nu1*nu2*
     +4.0D0-k*k2nu1*nu1*t4*t10+k*k2nu1*nu2*t4*t10))-amue2*(t2*(-t280+t35
     +2+ga1*k2nu2*nu1-k2nu2*t3*t10+ga1*k2nu2*nu1*t4+ga1*k2nu2*nu1*t10+k2
     +nu2*t3*t4*t10+ga1*k2nu2*nu1*t4*t10)+k2ga1*(t272+t273+ga1*k2nu1*k2n
     +u2*nu1+k2nu1*k2nu2*t3*t10-ga1*k2nu1*k2nu2*nu1*t4-ga1*k2nu1*k2nu2*n
     +u1*t10+k2nu1*k2nu2*t3*t4*t10+ga1*k2nu1*k2nu2*nu1*t4*t10)-xi1*(egah
     +1*enuh1*ga1*k*k2nu1*k2nu2*4.0D0+egah1*enuh1*k*k2ga1*k2nu2*nu1*4.0D
     +0)))
      A0(6,1) = -t111*(amue2*(t2*(xi2*(t356-t357)+egah1*ga1*k2ga2*nu1*2.
     +0D0+egah1*ga1*k2ga2*nu1*t4*2.0D0)+k2nu1*(-xi1*(egah1*ga1*k*k2ga2*2
     +.0D0+egah1*ga1*k*k2ga2*t4*2.0D0)+enuh1*k2ga1*k2ga2*t3*2.0D0+enuh1*
     +k2ga1*k2ga2*t3*t10*2.0D0)-xi2*(enuh1*k*nu1*t11*2.0D0-enuh1*k*nu1*t
     +10*t11*2.0D0)-xi1*(xi2*(t358-enuh1*ga1*k2ga1*nu1*2.0D0-egah1*k2ga1
     +*t3*t4*2.0D0+enuh1*ga1*k2ga1*nu1*t10*2.0D0)+enuh1*k*k2ga1*k2ga2*nu
     +1*2.0D0+enuh1*k*k2ga1*k2ga2*nu1*t10*2.0D0))+k2nu1*(xi1*(amue1*egah
     +1*ga2*k*k2ga1*2.0D0-amue1*enuh1*ga1*k*k2ga1*2.0D0-amue1*enuh1*ga2*
     +k*k2ga1*2.0D0+amue1*egah1*ga2*k*k2ga1*t4*2.0D0+amue1*enuh1*ga1*k*k
     +2ga1*t10*2.0D0-amue1*enuh1*ga2*k*k2ga1*t10*2.0D0)+amue1*enuh1*t3*t
     +11*2.0D0-amue1*enuh1*t3*t10*t11*2.0D0)+t2*(t354-t355-amue1*egah1*g
     +a2*k2ga1*nu1*2.0D0+amue1*enuh1*ga2*k2ga1*nu1*2.0D0-amue1*egah1*ga2
     +*k2ga1*nu1*t4*2.0D0+amue1*enuh1*ga2*k2ga1*nu1*t10*2.0D0)-t2*xi1*(a
     +mue1*egah1*ga1*k*2.0D0-amue1*egah1*ga1*k*t4*2.0D0))
      A0(6,2) = -t111*(amue2*(k2ga1*(xi2*(t360+egah1*k2nu1*t3*t4*2.0D0)-
     +xi1*xi2*(t361+enuh1*k*nu1*t10*2.0D0))+t2*(xi2*(enuh1*ga1*nu1*2.0D0
     ++enuh1*ga1*nu1*t10*2.0D0)-enuh1*k*k2ga2*nu1*2.0D0+enuh1*k*k2ga2*nu
     +1*t10*2.0D0)-xi1*(xi2*(egah1*ga1*k*k2nu1*2.0D0+egah1*ga1*k*k2nu1*t
     +4*2.0D0)-enuh1*k2ga2*k2nu1*t3*2.0D0+egah1*ga1*k2ga2*k2nu1*nu1*2.0D
     +0+enuh1*k2ga2*k2nu1*t3*t10*2.0D0-egah1*ga1*k2ga2*k2nu1*nu1*t4*2.0D
     +0)+egah1*ga1*k*k2ga2*t12*2.0D0-egah1*ga1*k*k2ga2*t4*t12*2.0D0)+k2g
     +a1*(xi1*(t359-amue1*egah1*k2nu1*t3*2.0D0+amue1*egah1*ga2*k2nu1*nu1
     +*2.0D0-amue1*egah1*k2nu1*t3*t4*2.0D0+amue1*enuh1*k2nu1*t3*t10*2.0D
     +0-amue1*egah1*ga2*k2nu1*nu1*t4*2.0D0)-amue1*egah1*ga2*k*t12*2.0D0+
     +amue1*egah1*ga2*k*t4*t12*2.0D0)+t2*(amue1*egah1*ga1*k*k2nu1*2.0D0-
     +amue1*enuh1*ga1*k*k2nu1*2.0D0-amue1*enuh1*ga2*k*k2nu1*2.0D0+amue1*
     +egah1*ga1*k*k2nu1*t4*2.0D0-amue1*enuh1*ga1*k*k2nu1*t10*2.0D0+amue1
     +*enuh1*ga2*k*k2nu1*t10*2.0D0)+t2*xi1*(amue1*enuh1*ga2*nu1*2.0D0-am
     +ue1*enuh1*ga2*nu1*t10*2.0D0))
      A0(6,3) = t111*(k*(t376-amue1*amue2*(t377+t378+t379+t380+t381+t382
     ++t383+t384+t385-k2ga1*k2ga2*t12*(0.0D0,1.0D0)-k2ga2*k2nu1*t2*(0.0D
     +0,1.0D0)+k2ga1*k2ga2*t4*t12*(0.0D0,1.0D0)-k2ga1*k2ga2*t10*t12*(0.0
     +D0,1.0D0)-k2ga2*k2nu1*t2*t4*(0.0D0,1.0D0)+k2ga2*k2nu1*t2*t10*(0.0D
     +0,1.0D0)+k2ga1*k2ga2*t4*t10*t12*(0.0D0,1.0D0)+k2ga2*k2nu1*t2*t4*t1
     +0*(0.0D0,1.0D0)))-amue1*amue2*(k2ga1*nu1*t2*xi2*(0.0D0,1.0D0)+k2ga
     +2*nu1*t2*xi1*(0.0D0,1.0D0)+k2nu1*nu1*t11*xi2*(0.0D0,1.0D0)+k2ga1*k
     +2ga2*k2nu1*nu1*xi1*(0.0D0,1.0D0)+k2ga1*nu1*t2*t4*xi2*(0.0D0,1.0D0)
     ++k2ga2*nu1*t2*t4*xi1*(0.0D0,1.0D0)+k2ga1*nu1*t2*t10*xi2*(0.0D0,1.0
     +D0)-k2ga2*nu1*t2*t10*xi1*(0.0D0,1.0D0)-k2nu1*nu1*t4*t11*xi2*(0.0D0
     +,1.0D0)-k2nu1*nu1*t10*t11*xi2*(0.0D0,1.0D0)+k2ga1*nu1*t2*t4*t10*xi
     +2*(0.0D0,1.0D0)-k2ga2*nu1*t2*t4*t10*xi1*(0.0D0,1.0D0)+k2nu1*nu1*t4
     +*t10*t11*xi2*(0.0D0,1.0D0)-egah1*enuh1*k2ga1*nu1*t2*xi2*(0.0D0,4.0
     +D0)-k2ga1*k2ga2*k2nu1*nu1*t4*xi1*(0.0D0,1.0D0)+k2ga1*k2ga2*k2nu1*n
     +u1*t10*xi1*(0.0D0,1.0D0)-k2ga1*k2ga2*k2nu1*nu1*t4*t10*xi1*(0.0D0,1
     +.0D0)))
      A0(6,4) = -t111*(amue1*(amue1*ga2*t9*(0.0D0,1.0D0)-amue1*ga2*t4*t9
     +*(0.0D0,1.0D0)-amue1*ga2*t9*t10*(0.0D0,1.0D0)+amue1*ga2*t11*t12*(0
     +.0D0,1.0D0)+amue1*ga2*k2ga1*k2nu1*t2*(0.0D0,2.0D0)+amue1*ga2*t4*t9
     +*t10*(0.0D0,1.0D0)-amue1*ga2*t4*t11*t12*(0.0D0,1.0D0)-amue1*ga2*t1
     +0*t11*t12*(0.0D0,1.0D0)+amue1*ga2*k2ga1*k2nu1*t2*t4*(0.0D0,2.0D0)+
     +amue1*ga2*k2ga1*k2nu1*t2*t10*(0.0D0,2.0D0)+amue1*ga2*t4*t10*t11*t1
     +2*(0.0D0,1.0D0)-amue1*egah1*enuh1*ga2*k2ga1*k2nu1*t2*(0.0D0,8.0D0)
     ++amue1*ga2*k2ga1*k2nu1*t2*t4*t10*(0.0D0,2.0D0))-amue2*(amue1*(ga1*
     +k2ga1*k2ga2*t12*(0.0D0,1.0D0)+ga1*k2ga2*k2nu1*t2*(0.0D0,1.0D0)-ga1
     +*t2*xi1*xi2*(0.0D0,1.0D0)-ga1*k2ga1*k2ga2*t4*t12*(0.0D0,1.0D0)-ga1
     +*k2ga1*k2ga2*t10*t12*(0.0D0,1.0D0)+ga1*k2ga2*k2nu1*t2*t4*(0.0D0,1.
     +0D0)+ga1*k2ga2*k2nu1*t2*t10*(0.0D0,1.0D0)-ga1*k2ga1*k2nu1*xi1*xi2*
     +(0.0D0,1.0D0)+ga1*t2*t4*xi1*xi2*(0.0D0,1.0D0)-ga1*t2*t10*xi1*xi2*(
     +0.0D0,1.0D0)-ga1*k2ga1*k2nu1*t4*xi1*xi2*(0.0D0,1.0D0)+ga1*k2ga1*k2
     +nu1*t10*xi1*xi2*(0.0D0,1.0D0)+ga1*t2*t4*t10*xi1*xi2*(0.0D0,1.0D0)-
     +egah1*enuh1*ga1*k2ga2*k2nu1*t2*(0.0D0,4.0D0)+ga1*k2ga1*k2ga2*t4*t1
     +0*t12*(0.0D0,1.0D0)+ga1*k2ga2*k2nu1*t2*t4*t10*(0.0D0,1.0D0)+ga1*k2
     +ga1*k2nu1*t4*t10*xi1*xi2*(0.0D0,1.0D0))+amue1*k*(k2ga1*t2*xi2*(0.0
     +D0,1.0D0)+k2ga2*t2*xi1*(0.0D0,1.0D0)+k2nu1*t11*xi2*(0.0D0,1.0D0)+k
     +2ga1*k2ga2*k2nu1*xi1*(0.0D0,1.0D0)-k2ga1*t2*t4*xi2*(0.0D0,1.0D0)-k
     +2ga2*t2*t4*xi1*(0.0D0,1.0D0)+k2ga1*t2*t10*xi2*(0.0D0,1.0D0)-k2ga2*
     +t2*t10*xi1*(0.0D0,1.0D0)+k2nu1*t4*t11*xi2*(0.0D0,1.0D0)-k2nu1*t10*
     +t11*xi2*(0.0D0,1.0D0)+k2ga1*k2ga2*k2nu1*t4*xi1*(0.0D0,1.0D0)+k2ga1
     +*k2ga2*k2nu1*t10*xi1*(0.0D0,1.0D0)-k2ga1*t2*t4*t10*xi2*(0.0D0,1.0D
     +0)+k2ga2*t2*t4*t10*xi1*(0.0D0,1.0D0)-k2nu1*t4*t10*t11*xi2*(0.0D0,1
     +.0D0)+k2ga1*k2ga2*k2nu1*t4*t10*xi1*(0.0D0,1.0D0)-egah1*enuh1*k2ga1
     +*k2ga2*k2nu1*xi1*(0.0D0,4.0D0))))
      A0(6,5) = -t111*(amue2*(t2*(-t323+t386+ga1*k2ga2*nu1-k2ga2*t3*t4+g
     +a1*k2ga2*nu1*t4+ga1*k2ga2*nu1*t10+k2ga2*t3*t4*t10+ga1*k2ga2*nu1*t4
     +*t10)+k2nu1*(t315+t316+ga1*k2ga1*k2ga2*nu1+k2ga1*k2ga2*t3*t4-ga1*k
     +2ga1*k2ga2*nu1*t4-ga1*k2ga1*k2ga2*nu1*t10+k2ga1*k2ga2*t3*t4*t10+ga
     +1*k2ga1*k2ga2*nu1*t4*t10)-xi1*(egah1*enuh1*ga1*k*k2ga2*k2nu1*4.0D0
     ++egah1*enuh1*k*k2ga1*k2ga2*nu1*4.0D0))-amue1*(t2*(-t222+t261-t331+
     +t350+ga2*k2ga1*nu1+ga2*k2ga1*nu1*t4+ga2*k2ga1*nu1*t10-egah1*enuh1*
     +ga2*k2ga1*nu1*4.0D0+ga2*k2ga1*nu1*t4*t10)-k2nu1*(t311-t387-ga2*nu1
     +*t11+t3*t4*t11+ga2*nu1*t4*t11+ga2*nu1*t10*t11-t3*t4*t10*t11-ga2*nu
     +1*t4*t10*t11)+k2nu1*xi1*(ga1*k*k2ga1+ga2*k*k2ga1+ga1*k*k2ga1*t4+ga
     +2*k*k2ga1*t4-ga1*k*k2ga1*t10+ga2*k*k2ga1*t10-egah1*enuh1*ga2*k*k2g
     +a1*4.0D0-ga1*k*k2ga1*t4*t10+ga2*k*k2ga1*t4*t10)+t2*xi1*(ga1*k+ga2*
     +k-ga1*k*t4-ga2*k*t4+ga1*k*t10-ga2*k*t10-ga1*k*t4*t10+ga2*k*t4*t10)
     +))
      A0(6,6) = t111*(amue1*(t2*(-t388+ga1*k*k2nu1+ga2*k*k2nu1+ga1*k*k2n
     +u1*t4+ga2*k*k2nu1*t4+ga1*k*k2nu1*t10-ga2*k*k2nu1*t10+ga1*k*k2nu1*t
     +4*t10-ga2*k*k2nu1*t4*t10)+k2ga1*(ga1*k*t12+ga2*k*t12-ga1*k*t4*t12-
     +ga2*k*t4*t12-ga1*k*t10*t12+ga2*k*t10*t12+ga1*k*t4*t10*t12-ga2*k*t4
     +*t10*t12)+k2ga1*xi1*(t224+t240+t341+t389-ga2*k2nu1*nu1-egah1*enuh1
     +*k2nu1*t3*4.0D0+ga2*k2nu1*nu1*t4-ga2*k2nu1*nu1*t10+ga2*k2nu1*nu1*t
     +4*t10)+t2*xi1*(t3-t205-t227+t351-ga2*nu1-ga2*nu1*t4+ga2*nu1*t10+ga
     +2*nu1*t4*t10))-amue2*xi2*(t397+k2ga1*(t224+t240+t341+t389+ga1*k2nu
     +1*nu1-ga1*k2nu1*nu1*t4-ga1*k2nu1*nu1*t10+ga1*k2nu1*nu1*t4*t10)-t39
     +1*xi1))
      
      end 
