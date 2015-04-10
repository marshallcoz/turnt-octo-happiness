function [G,T] = Diff_x(cellMat,CoefpsvH,CoefpsvV,Coefsh,p_x,pXi,m,f,ops,dir)
% Tensor en p_X del campo difractado por
% la estratigrafía + la incidencia real.
% Se arroja la respuesta en frecuencia y
% en coordenadas cartesianas.
G = zeros(3);
T = zeros(3);
[e] = el_estrato_es(m,ops.N,p_x.center(3));
[eXi] = el_estrato_es(m,ops.N,pXi.center(3));

% Coordenada polar entre fuente y receptor
[r,th,thp,z] = cilindricas(p_x,pXi);

% en el estrato del receptor
sD0 = cellMat{1,1}.sD0{1,1}(:,:,:,e);  %Uz,Ur,Ut
sS0 = cellMat{1,1}.sS0{1,1}(:,:,:,e);
sD0sh = cellMat{1,1}.sD0sh{1,1}(:,:,:,e); %Ur,Ut
sS0sh = cellMat{1,1}.sS0sh{1,1}(:,:,:,e);
dosOcuatro=4;
unoOdos = 2;
if (e == ops.N+1)
    dosOcuatro = 2;
    unoOdos = 1;
end
if (dir == 3)
    Coef = CoefpsvV(4*(e-1)+1:4*(e-1)+dosOcuatro,1,:);
else
    Coef = CoefpsvH(4*(e-1)+1:4*(e-1)+dosOcuatro,1,:);
    CoefSH = Coefsh(2*(e-1)+1:2*(e-1)+unoOdos,1,:);
end
ome = f.come;
alf = m.alfa;
bet = m.beta;
mu = m.amu;
lam = m.lam;
omealf = ome/alf;
%omebet = ome/bet;
k = cellMat{1,1}.k{1,1};
kr = k*r;
k2 = k.^2;
J0kr= besselj(0,kr);
J1kr= besselj(1,kr);
J1_kr = J1kr./kr;
Jo_J1_kr = J0kr-J1_kr;
% Campo difractado por estratigrafía en receptores
%% Uz:
a = (sD0(1,1,:).*Coef(1,1,:) + ...
    sD0(1,2,:).*Coef(2,1,:));
if (e ~= ops.N+1)
    a = a + (sD0(1,3,:).*Coef(3,1,:) + ...
        sD0(1,4,:).*Coef(4,1,:));
end
if (dir == 3) % fuerza vertical
    uz = sum(a(:).*J0kr(:)/(-1i))*f.dk;
else % fuerza horizontal
    uz = sum(a(:).*J1kr(:)/(-1i))*f.dk*cos(th);
end

%% Ur,Uth:
a = (sD0(2,1,:).*Coef(1,1,:) + ...
    sD0(2,2,:).*Coef(2,1,:));
if (e ~= ops.N+1)
    a = a + (sD0(2,3,:).*Coef(3,1,:) + ...
        sD0(2,4,:).*Coef(4,1,:));
end
if (dir == 3) % fuerza vertical
    ur = sum(a(:).*J1kr(:))*f.dk;
    ux = ur*cos(th);
    uy = ur*sin(th);
    G(:,3) = [ux uy uz];
else %fuerza horizontal
    ur = sum(a(:).*Jo_J1_kr(:)); %PSV
    ut = sum(-a(:).*J1_kr(:));
    a = (sD0sh(1,1,:).*CoefSH(1)); %SH
    if (e ~= ops.N+1)
        a = a + (sD0sh(1,2,:).*CoefSH(2));
    end
    ur = ur + sum(a(:).*J1_kr(:));
    ut = ut + sum(-a(:).*Jo_J1_kr(:));
    ur = ur *f.dk*cos(th);
    ut = ut *f.dk*sin(th);
    
    ux = ur*cos(th)-ut*sin(th);
    uy = ur*sin(th)+ut*cos(th);
    G(:,dir) = [ux uy uz];
end

%% szz:
a = (sS0(1,1,:).*Coef(1,1,:) + ...
    sS0(1,2,:).*Coef(2,1,:));
if (e ~= ops.N+1)
    a = a + (sS0(1,3,:).*Coef(3,1,:) + ...
        sS0(1,4,:).*Coef(4,1,:));
end
if (dir == 3) % fuerza vertical
    szz = sum(a(:).*J0kr(:))*f.dk;
else %fuerza horizontal
    szz = sum(a(:).* J1kr(:))*f.dk*cos(th);
end
%% szr
a = (sS0(2,1,:).*Coef(1,1,:) + ...
    sS0(2,2,:).*Coef(2,1,:));
if (e ~= ops.N+1)
    a = a + (sS0(2,3,:).*Coef(3,1,:) + ...
        sS0(2,4,:).*Coef(4,1,:));
end
if (dir == 3) % fuerza vertical
    szr = sum(a(:).*J1_kr(:)/1i)*f.dk;
else %fuerza horizontal
    szr = sum(a(:).*Jo_J1_kr(:));
    a = (sS0sh(1,1,:).*CoefSH(1)); %SH
    if (e ~= ops.N+1)
        a = a + (sS0sh(1,2,:).*CoefSH(2));
    end
    szr = szr + sum(a(:).*J1_kr(:));
    szr = szr/1i *f.dk*cos(th);
end
%% szt
a = (sS0(3,1,:).*Coef(1,1,:) + ...
    sS0(3,2,:).*Coef(2,1,:));
if (e ~= ops.N+1)
    a = a + (sS0(3,3,:).*Coef(3,1,:) + ...
        sS0(3,4,:).*Coef(4,1,:));
end
if (dir == 3) % fuerza vertical
    szt = 0;
else %fuerza horizontal
    szt = sum(a(:).*J1_kr(:));
    a = (sS0sh(2,1,:).*CoefSH(1)); %SH
    if (e ~= ops.N+1)
        a = a + (sS0sh(2,2,:).*CoefSH(2));
    end
    szt = szt + sum(a(:).*Jo_J1_kr(:));
    szt = szt/1i *f.dk*sin(th);
end
%% fix coef A ikB A ikB -> A B A B
a = Coef(2,1,:);
Coef(2,1,:) = a(:)./(1i*k(:));
a = Coef(4,1,:);
Coef(4,1,:) = a(:)./(1i*k(:));
%% srr:

if (dir == 3) % fuerza vertical
    t1 = -lam*omealf^2*J0kr-2*mu*k2.*Jo_J1_kr;
    t2 = 2*mu*k2.*Jo_J1_kr;
    t4 = sS0(4,1,:).*Coef(1,1,:);
    t5 = sS0(4,2,:).*Coef(2,1,:);
    a = (t4(:).*t1(:) + ...
         t5(:).*t2(:));
    if (e ~= ops.N+1)
        t4 = sS0(4,3,:).*Coef(3,1,:);
        t5 = sS0(4,4,:).*Coef(4,1,:);
        a = a + (t4(:).*t1(:) + ...
                 t5(:).*t2(:));
    end
    srr = sum(a(:))*f.dk;
else %fuerza horizontal
    t2 = 2*mu*k.*(-J1kr.*k-J0kr/r+2*J1kr./(kr));
    t1 = -lam*omealf^2*J1kr+t2;
    t3 = 2*mu*k.*(J0kr/r-2*J1kr./(k*r*r));
    t4 = sS0(4,1,:).*Coef(1,1,:);
    t5 = sS0(4,2,:).*Coef(2,1,:);
    a = (t4(:).*t1(:) + ...
         t5(:).*t2(:));
    if (e ~= ops.N+1)
        t4 = sS0(4,3,:).*Coef(3,1,:);
        t5 = sS0(4,4,:).*Coef(4,1,:);
        a = a + (t4(:).*t1(:) + ...
                 t5(:).*t2(:));
    end
    srr = sum(a(:));
    t4 = sS0sh(3,1,:).*CoefSH(1);
    a = t4(:).*t3(:); %SH
    if (e ~= ops.N+1)
        t4 = sS0sh(3,2,:).*CoefSH(2);
        a = a + t4(:).*t3(:);
    end
    srr = (srr + sum(a(:)))*cos(th)*f.dk;
end
%% srt:
if (dir == 3) % fuerza vertical
    srt = 0;
else %fuerza horizontal
    t1 = 2*mu*J1kr./r^2;
    t3 = mu*(2*J0kr.*k/r+J1kr.*k2);
    t4 = sS0(5,1,:).*Coef(1,1,:);
    t5 = sS0(5,2,:).*Coef(2,1,:);
    a = (t4(:).*t1(:) + ...
         t5(:).*t1(:));
    if (e ~= ops.N+1)
        t4 = sS0(5,3,:).*Coef(3,1,:);
        t5 = sS0(5,4,:).*Coef(4,1,:);
        a = a + (t4(:).*t1(:) + ...
                 t5(:).*t1(:));
    end
    srt = sum(a(:));
    t4 = sS0sh(4,1,:).*CoefSH(1);
    a = t4(:).*t3(:); %SH
    if (e ~= ops.N+1)
        t4 = sS0sh(4,2,:).*CoefSH(2);
        a = a + t4(:).*t3(:);
    end
    srt = (srt + sum(a(:)))*sin(th)*f.dk;
end
%% stt:
if (dir == 3) % fuerza vertical
    t2 = 2*mu/r*k.*J1kr;
    t1 = -lam*omealf^2*J0kr-t2;
    t4 = sS0(6,1,:).*Coef(1,1,:);
    t5 = sS0(6,2,:).*Coef(2,1,:);
    a = (t4(:).*t1(:) + ...
         t5(:).*t2(:));
if (e ~= ops.N+1)
    t4 = sS0(6,3,:).*Coef(3,1,:);
    t5 = sS0(6,4,:).*Coef(4,1,:);
    a = a + (t4(:).*t1(:) + ...
             t5(:).*t2(:));
end
    stt = sum(a(:));
else %fuerza horizontal
    t2 = 2*mu*k/r.*(J0kr);
    t1 = -lam*omealf^2*J1kr+t2;
    t3 = 2*mu*J1kr/r^2;
    t4 = sS0(6,1,:).*Coef(1,1,:);
    t5 = sS0(6,2,:).*Coef(2,1,:);
    a = (t4(:).*t1(:) + ...
         t5(:).*t2(:));
    if (e ~= ops.N+1)
        t4 = sS0(6,3,:).*Coef(3,1,:);
        t5 = sS0(6,4,:).*Coef(4,1,:);
        a = a + (t4(:).*t1(:) + ...
                 t5(:).*t2(:));
    end
    stt = sum(a(:));
    t4 = sS0sh(5,1,:).*CoefSH(1);
    a = t4(:).*t3(:); %SH
    if (e ~= ops.N+1)
        t4 = sS0sh(5,2,:).*CoefSH(2);
        a = a + t4(:).*t3(:);
    end
    stt = (stt + sum(a(:)))*cos(th)*f.dk;
end


esf_cil = [[srr srt szr]; ...
           [srt stt szt]; ...
           [szr szt szz]]; % en polares
[esf_car] = esf_cil2car(esf_cil,th);
%las tracciones en x
T(:,dir) = tracciones(esf_car,p_x.normal);

%% adicionar campo incidente
if (e == eXi)
    %[G0,T0,nil] = GijTij_omega(m,f,p_x,pXi);
    %G = G + G0;
    %T = T + T0;
end


%% tensor completo de campo difractado (programar de algebras 6 elemntos faltantes)
%% fase horizontal e integración. (hacer como para el campo completo)
%% a cartesianas
%% + campo incidente
%% repartir

