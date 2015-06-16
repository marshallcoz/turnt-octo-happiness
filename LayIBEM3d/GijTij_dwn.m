%% funcion de Green 3d Espacio completo 
% vía integral en num de onda radial
function [G,T,FK] = GijTij_dwn(m,f,p_x,pXi)
%m = m_vars;f=f_vars;
%% se resuelve en cilíndricas y se transforma a cartesianas
G = zeros(3);
T = zeros(3);
[r,th,thp,z] = cilindricas(p_x,pXi);
k = 0:f.dk:f.dk*f.nmax; k(1) = 0.01*f.dk;
k2 = k.^2;
k3 = k.^3;
kr = k*r;
ome = f.come;
alf = m.alfa;
bet = m.beta;
mu = m.amu;
lam = m.lam;
b2_4piom2m = bet^2/(4*pi*ome^2*mu);
b2_4pio2 = bet^2/(4*pi*ome^2);
gam =((ome/alf)^2 - k2).^0.5;
nu = ((ome/bet)^2 - k2).^0.5;
[found] = find(imag(gam)>0);  gam(found) = -(gam(found)); clear found
[found] = find(imag(nu)>0);  nu(found) = -(nu(found)); clear found
% figure;plot(real(gam),'r');hold on;plot(imag(gam),'b')
% figure;plot(real(nu),'r');hold on;plot(imag(nu),'b')
nu2 = nu.^2;
ga2 = gam.^2;
k2_ga = k2./gam;
J0kr= besselj(0,kr);
J1kr= besselj(1,kr);
ega = exp(-1i*gam*abs(z));
enu = exp(-1i*nu*abs(z));
om2be2 = ome^2/bet^2;
%figure;plot(real(ega),'r');hold on;plot(imag(ega),'b')
%figure;plot(real(enu),'r');hold on;plot(imag(enu),'b')
%% FK
FK = zeros(1,f.nmax);
for ik = 1:f.nmax
   FK(ik) = ((k2_ga(ik)*ega(ik)+nu(ik)*enu(ik))*(J0kr(ik)-J1kr(ik)/kr(ik))+...
    om2be2/nu(ik)*enu(ik)*J1kr(ik)/kr(ik))*k(ik);
end

%%
%fuerza vertical j=3(no importa th) -------------------------------------
%desplazamientos
ur = -sign(z)*b2_4piom2m*(sum((ega - enu).*k2.*J1kr)*f.dk);
uz = -1i*b2_4piom2m*(sum((gam.*ega + k2./nu.*enu).*k.*J0kr)*f.dk);
ux = ur*cos(th);
uy = ur*sin(th);
G(:,3) = [ux uy uz];
%% esfuerzos:
% rr rt rz
% tr tt tz
% zr zt zz
cons1 = sign(z)*b2_4pio2;
szz = cons1*sum(((k2-nu2).*ega-2*k2.*enu).*k.*J0kr)*f.dk;
szr = 1i*b2_4pio2*sum((2*gam.*ega+(k2-nu2)./nu.*enu).*k2.*J1kr)*f.dk;
szt = 0;
srr = -cons1*sum(((nu2-2*ga2+k2).*k.*J0kr-(k2/(mu*r)).*J1kr).*ega -...
    (2*k3.*J0kr-2/r*k2.*J1kr).*enu)*f.dk;
stt = -cons1*sum((ome^2*(1/bet^2-2/alf^2)*k.*J0kr + 2/r*k2.*J1kr).*ega -...
    (2/r*k2.*J1kr).*enu)*f.dk;
srt = 0;
esf_cil = [[srr srt szr]; ...
           [srt stt szt]; ...
           [szr szt szz]]; % en polares
[esf_car] = esf_cil2car(esf_cil,th);
%las tracciones en x
T(:,3) = tracciones(esf_car,p_x.normal);

%fuerza horizontal j=1 ---------------------------------------------------
%% desplazamientos:
ur = -1i*b2_4piom2m*(sum(((k2_ga.*ega+nu.*enu).*(J0kr-J1kr./kr)+...
    om2be2./nu.*enu.*J1kr./kr).*k)*f.dk*cos(th));
uth = 1i*b2_4piom2m*(sum(((k2_ga.*ega+nu.*enu).*J1kr./kr+...
    om2be2./nu.*enu.*(J0kr-J1kr./kr)).*k)*f.dk*sin(th));
uz = -sign(z)*b2_4piom2m*(sum((ega-enu).*k2.*J1kr)*f.dk*cos(th));
ux = ur*cos(th)-uth*sin(th);
uy = ur*sin(th)+uth*cos(th);
G(:,1) = [ux uy uz];
%% esfuerzos:
szz = -1i*b2_4pio2*sum(((k2-nu2)./gam.*ega + 2*nu.*enu).*k2.*J1kr)*f.dk*cos(th);
szr = -cons1*sum(((2*k2.*ega-(k2-nu2).*enu).*(J0kr-J1kr./kr)+om2be2*enu.*J1kr./kr).*k)*f.dk*cos(th);
szt = cons1*sum(((2*k2.*ega-(k2-nu2).*enu).*J1kr./kr+om2be2*enu.*(J0kr-J1kr./kr)).*k)*f.dk*sin(th);
srr = -1i*b2_4pio2*sum((ega.*(J1kr.*(2*ga2-om2be2+4/r^2)-2/r*J0kr.*k).*k./gam + ...
    enu.*(J1kr.*(-4/r^2*k./nu-2*nu.*k)+2/r*J0kr.*k2./nu)).*k)*f.dk*cos(th);
stt = -1i*b2_4piom2m*sum((ega.*(J1kr.*(-ome^2/alf^2*k./gam)*lam - 2*mu/r^2*k2.*J1kr./gam +...
    2*mu/r*(k2./gam.*(J0kr-J1kr./kr))) + enu.*(2*mu/r*(-nu.*J1kr./kr-ome^2/bet^2./nu.*(J0kr-J1kr./kr)+ ...
    nu.*(J0kr-J1kr./kr)+ome^2/bet^2./nu.*J1kr./kr))).*k)*f.dk*cos(th);
dUth_r = 1i*b2_4piom2m*sum(((k2./gam.*ega+nu.*enu).*(J0kr/r-2*J1kr./k/r^2) + ...
    ome^2/bet^2./nu.*enu.*(J1kr.*(2./kr/r-k)-J0kr/r)).*k)*f.dk*sin(th);
srt = mu*(-1/r*ur*tan(th)+ dUth_r - uth/r);
esf_cil = [[srr srt szr]; ...
     [srt stt szt]; ...
     [szr szt szz]]; % en polares
[esf_car] = esf_cil2car(esf_cil,th);
%las tracciones en x
T(:,1) = tracciones(esf_car,p_x.normal);

% fuerza horizontal j=2 --------------------------------------------------
%desplazamientos:
ur = -1i*b2_4piom2m*(sum(((k2_ga.*ega+nu.*enu).*(J0kr-J1kr./kr)+...
    om2be2./nu.*enu.*J1kr./kr).*k)*f.dk*cos(thp));
uth = 1i*b2_4piom2m*(sum(((k2_ga.*ega+nu.*enu).*J1kr./kr+...
    om2be2./nu.*enu.*(J0kr-J1kr./kr)).*k)*f.dk*sin(thp));
uz = -sign(z)*b2_4piom2m*(sum((ega-enu).*k2.*J1kr)*f.dk*cos(thp));
ux = ur*cos(th)-uth*sin(th);
uy = ur*sin(th)+uth*cos(th);
G(:,2) = [ux uy uz];
%esfuerzos:
szz = -1i*b2_4pio2*sum(((k2-nu2)./gam.*ega + 2*nu.*enu).*k2.*J1kr)*f.dk*cos(thp);
szr = -cons1*sum(((2*k2.*ega-(k2-nu2).*enu).*(J0kr-J1kr./kr)+om2be2*enu.*J1kr./kr).*k)*f.dk*cos(thp);
szt = cons1*sum(((2*k2.*ega-(k2-nu2).*enu).*J1kr./kr+om2be2*enu.*(J0kr-J1kr./kr)).*k)*f.dk*sin(thp);
srr = -1i*b2_4pio2*sum((ega.*(J1kr.*(2*ga2-om2be2+4/r^2)-2/r*J0kr.*k).*k./gam + ...
    enu.*(J1kr.*(-4/r^2*k./nu-2*nu.*k)+2/r*J0kr.*k2./nu)).*k)*f.dk*cos(thp);
stt = -1i*b2_4piom2m*sum((ega.*(J1kr.*(-ome^2/alf^2*k./gam)*lam - 2*mu/r^2*k2.*J1kr./gam +...
    2*mu/r*(k2./gam.*(J0kr-J1kr./kr))) + enu.*(2*mu/r*(-nu.*J1kr./kr-ome^2/bet^2./nu.*(J0kr-J1kr./kr)+ ...
    nu.*(J0kr-J1kr./kr)+ome^2/bet^2./nu.*J1kr./kr))).*k)*f.dk*cos(thp);
dUth_r = 1i*b2_4piom2m*sum(((k2./gam.*ega+nu.*enu).*(J0kr/r-2*J1kr./k/r^2) + ...
    ome^2/bet^2./nu.*enu.*(J1kr.*(2./kr/r-k)-J0kr/r)).*k)*f.dk*sin(thp);
srt = mu*(-1/r*ur*tan(thp)+ dUth_r - uth/r);
esf_cil = [[srr srt szr]; ...
     [srt stt szt]; ...
     [szr szt szz]]; % en polares
[esf_car] = esf_cil2car(esf_cil,th);
%las tracciones en x
T(:,2) = tracciones(esf_car,p_x.normal);


