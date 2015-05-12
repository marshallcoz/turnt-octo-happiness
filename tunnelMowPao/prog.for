        PROGRAM MCJYVB
C
C       =============================================================
C       Purpose: This program computes Bessel functions Jv(z), Yv(z),
C                and their derivatives for a complex argument using
C                subroutine CJYVB
C       Input :  z --- Complex argument
C                v --- Order of Jv(z) and Yv(z)
C                      ( v = n+v0, 0 Û n Û 250, 0 Û v0 < 1 )
C       Output:  CBJ(n) --- Jn+v0(z)
C                CDJ(n) --- Jn+v0'(z)
C                CBY(n) --- Yn+v0(z)
C                CDY(n) --- Yn+v0'(z)
C       Example:
C                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
C
C     n     Re[Jv(z)]       Im[Jv(z)]      Re[Jv'(z)]      Im[Jv'(z)]
C    -------------------------------------------------------------------
C     0  -.13829878D+01  -.30855145D+00  -.18503756D+00   .13103689D+01
C     1   .82553327D-01  -.12848394D+01  -.12336901D+01   .45079506D-01
C     2   .10843924D+01  -.39871046D+00  -.33046401D+00  -.84574964D+00
C     3   .74348135D+00   .40665987D+00   .45318486D+00  -.42198992D+00
C     4   .17802266D+00   .44526939D+00   .39624497D+00   .97902890D-01
C     5  -.49008598D-01   .21085409D+00   .11784299D+00   .19422044D+00
C
C     n     Re[Yv(z)]      Im[Yv(z)]       Re[Yv'(z)]      Im[Yv'(z)]
C    -------------------------------------------------------------------
C     0   .34099851D+00  -.13440666D+01  -.13544477D+01  -.15470699D+00
C     1   .13323787D+01   .53735934D-01  -.21467271D-01  -.11807457D+01
C     2   .38393305D+00   .10174248D+01   .91581083D+00  -.33147794D+00
C     3  -.49924295D+00   .71669181D+00   .47786442D+00   .37321597D+00
C     4  -.57179578D+00   .27099289D+00  -.12111686D+00   .23405313D+00
C     5  -.25700924D+00   .24858555D+00  -.43023156D+00  -.13123662D+00
C       =============================================================
C
        USE SPECFUN
        IMPLICIT DOUBLE PRECISION (V,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250)
        WRITE(*,*)'  Please enter v, x and y ( z=x+iy )'
        READ(*,*)V,X,Y
        Z=CMPLX(X,Y)
        N=INT(V)
        V0=V-N
        WRITE(*,25)V0,X,Y
        IF (N.LE.8) THEN
           NS=1
        ELSE
           WRITE(*,*)'  Please enter order step Ns'
           READ(*,*)NS
        ENDIF
        CALL CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
        NM=INT(VM)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Jv(z)]       Im[Jv(z)]',
     &            '       Re[Jv''(z)]      Im[Jv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 10 K=0,NM,NS
10         WRITE(*,20) K,CBJ(K),CDJ(K)
        WRITE(*,*)
        WRITE(*,*)'  n       Re[Yv(z)]       Im[Yv(z)]',
     &            '       Re[Yv''(z)]      Im[Yv''(z)]'
        WRITE(*,*)' ----------------------------------',
     &            '-----------------------------------'
        DO 15 K=0,NM,NS
15         WRITE(*,20) K,CBY(K),CDY(K)
20      FORMAT(1X,I3,2X,4D16.8)
25      FORMAT(8X,'v = n+v0',',  v0 =',F5.2,',  z =',F7.2,' +',F7.2,'i')
        END