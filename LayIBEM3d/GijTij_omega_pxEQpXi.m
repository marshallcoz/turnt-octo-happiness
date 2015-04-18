function [GNIJ,T] = GijTij_omega_pxEQpXi(m,f,p_x,pXi)
% Greenex: Green function for displacement by analytic integration
% when the source and receiver are very near or at the same point.
T = eye(3)*0.5;
GNIJ = zeros(3);
%
VN = pXi.normal(1:3); %normal de la superficie
RNM = distancia(p_x,pXi); %distancia
R = pXi.radio; %radio del disco
if (RNM >= 0.1*R)
    EPS = RNM/R;
    G = (p_x.center(1:3) - pXi.center(1:3))/RNM; %cosenos directores
else
    EPS = 0;
    G = 0;
end
clear RNM 

CAKA = f.come/m.beta;
%CAQA = f.come/m.alfa;

BEA2  = m.bealf2;
BEA3  = m.bealf3;
BEA4  = m.bealf4;

% GREEN FUNCTION (WITH ANALYTIC INTEGRATION)
%CAQR=CAQA*R;
CAKR=CAKA*R;
%AQR2=CAQR*CAQR;
AKR2=CAKR*CAKR;
%AQR3=AQR2*CAQR;
AKR3=AKR2*CAKR;
if (EPS == 0)
F1=1-1i*(2.0+BEA3)*CAKR/6.0-1*(1.0+0.5*BEA4)*AKR2/9.0+1i*BEA2*BEA3*AKR3/24.0;
F1=F1/(pi*R);
F2=1*(1.0+BEA2)/2.0-1i*(2.0+BEA3)*CAKR/6.0-1*(2.0+BEA4)*AKR2/18.0+1i*AKR3/24.0;
F2=F2/(pi*R);
for I=1:3
    for J=1:3
        AUX=(F2-F1)*VN(I)*VN(J);
        if(I == J),AUX=AUX+F1+F2;end
        GNIJ(I,J)=AUX/4.0/m.amu;
    end
end

else
PG(1)=VN(2)*G(3)-VN(3)*G(2);
PG(2)=VN(3)*G(1)-VN(1)*G(3);
PG(3)=VN(1)*G(2)-VN(2)*G(1);
EPS2=EPS*EPS;
EPS4=EPS2*EPS2;
EPS6=EPS4*EPS2;
C1=1.0-EPS2*0.25-EPS4*3./64.-EPS6*5./192.;
C2=1.0;
C3=1.0+EPS2*0.75-EPS4*3./64.-EPS6/256.;
C4=1.0+EPS2*2.;
F2C=1*(1.+BEA2)/2.*C1-1i*(2.+BEA3)*CAKR/6.*C2-1*(2.+BEA4)*AKR2/18.*C3+1i*AKR3/24.*C4;
F2C=F2C/(pi*R);
D1=1.0-EPS2*0.125-EPS4/64.-EPS6*5./1024.;
D2=1.0+EPS2*0.5;
D3=1.0+EPS2*21./8.-EPS4*21./64.-EPS6*43./1024.;
D4=1.0+EPS2*4.;
E1=1.0-EPS2*0.375-EPS4*5./64.-EPS6*35./1024.;
E2=1.0-EPS2*0.5;
E3=1.0+EPS2*1.875-EPS4*57./64.-EPS6*55./1024.;
E4=1.0;
F1D=1*D1-1i*(2.+BEA3)*CAKR/6.*D2-1*(1.+.5*BEA4)*AKR2/9.*D3+1i*BEA2*BEA3*AKR3/24.*D4;
F1D=F1D/(pi*R);
F2D=1*(1.+BEA2)/2.*D1-1i*(2.+BEA3)*CAKR/6.*D2-1*(2.+BEA4)*AKR2/18.*D3+1i*AKR3/24.*D4;
F2D=F2D/(pi*R);

F1E=1*E1-1i*(2.+BEA3)*CAKR/6.*E2-1*(1.+.5*BEA4)*AKR2/9.*E3+1i*BEA2*BEA3*AKR3/24.*E4;
F1E=F1E/(pi*R);
F2E=1*(1.+BEA2)/2.*E1-1i*(2.+BEA3)*CAKR/6.*E2 -1*(2.+BEA4)*AKR2/18.*E3+1i*AKR3/24.*E4;
F2E=F2E/(pi*R);

for I=1:3
    for J=1:3
        AUX=(F1D-F2D)*G(I)*G(J)+(F1E-F2E)*PG(I)*PG(J);
        if(I == J), AUX=AUX+F2C*2.0;end
        GNIJ(I,J)=AUX/4.0/m.amu;
    end
end
end