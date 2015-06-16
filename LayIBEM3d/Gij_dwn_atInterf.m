%% funcion de Green 3d Espacio completo sin fase radial 
% (para el término independiente)
function [PSV,SH] = Gij_dwn_atInterf(m,f,z)
PSV = zeros(4,2,f.nmax+1); %Szz,Szr,Uz,Ur; Horz, Vert; ik
SH =  zeros(2,1,f.nmax+1); %    Szr,   Ur; Horz      ; ik
k = 0:f.dk:f.dk*f.nmax; k(1) = 0.01*f.dk;
k2 = k.^2;
ome = f.come;
alf = m.alfa;
bet = m.beta;
mu = m.amu;
b2_4piom2m = bet^2/(4*pi*ome^2*mu);
b2_4pio2 = bet^2/(4*pi*ome^2);
gam =((ome/alf)^2 - k2).^0.5;
nu = ((ome/bet)^2 - k2).^0.5;
[found] = find(imag(gam)>0);  gam(found) = -(gam(found)); clear found
[found] = find(imag(nu)>0);  nu(found) = -(nu(found)); clear found
nu2 = nu.^2;
k2_ga = k2./gam;
k2_nu2 = k2-nu2;

ega = exp(-1i*gam*abs(z));
enu = exp(-1i*nu*abs(z));
om2be2 = ome^2/bet^2;

%fuerza vertical (2) -------------------------------------
% esfuerzos:
cons1 = sign(z)*b2_4pio2;
PSV(1,2,:) = cons1*(k2_nu2.*ega-2*k2.*enu).*k; %szz
PSV(2,2,:) = -b2_4pio2*(2*gam.*ega+k2_nu2./nu.*enu).*k2;%szr
% desplazamientos:
PSV(3,2,:) = -b2_4piom2m*((gam.*ega + k2./nu.*enu).*k); %uz
PSV(4,2,:) = -sign(z)*b2_4piom2m*((ega - enu).*k2); %ur

%fuerza horizontal (1) -------------------------------------------
%PSV
% esfuerzos:
PSV(1,1,:) = -1i*b2_4pio2*(k2_nu2./gam.*ega + 2*nu.*enu).*k2; %szz
PSV(2,1,:) = -1i*cons1*(2*k2.*ega-k2_nu2.*enu).*k; %szr
% desplazamientos:
PSV(3,1,:) = sign(z)*1i*b2_4piom2m*(ega-enu).*k2; %uz
PSV(4,1,:) = -1i*b2_4piom2m*(k2_ga.*ega+nu.*enu).*k; % ur

%SH
SH(1,1,:) = -1i*cons1*om2be2*enu.*k; %szr
SH(2,1,:) = -1i*b2_4piom2m*om2be2./nu.*enu.*k; % ur
