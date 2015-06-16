function [Bpsv,B_sh] = Gij_B_dwn(m,f,ops,pXi)
Bpsv = zeros(4*ops.N+2,2,f.nmax+1); %Horz, Vert; ik
B_sh = zeros(2*ops.N+1,1,f.nmax+1); %Horz      ; ik
% las fuentes son los discos cargados

e_f = pXi.layer; %estrato de la fuente
% si el disco está sobre una interfaz de suelo
if (pXi.isoninterf)
    
else
    % interfaz de arriba *********************************
    zloc(1) = m(e_f).z - pXi.center(3); %Downward
    zloc(2) = m(e_f+1).z - pXi.center(3); %Upward
    [PSV,SH] = Gij_dwn_atInterf(m(e_f),f,zloc(1));
end
% acomodar en vector PSV: w u szz szx  SH: v sxx
if (e_f ~= 1)
    % fza horizontal
    Bpsv(1+4*(e_f-1)-2,1,:) = PSV(3,1,:); %  uz
    Bpsv(1+4*(e_f-1)-1,1,:) = PSV(4,1,:); %  ur
    B_sh(1+2*(e_f-1)-1,1,:) = SH(2,1,:);  %  ur
    % fza vertical
    Bpsv(1+4*(e_f-1)-2,2,:) = PSV(3,2,:); %  uz
    Bpsv(1+4*(e_f-1)-1,2,:) = PSV(4,2,:); %  ur
end
    % fza horizontal
    Bpsv(1+4*(e_f-1)  ,1,:) = PSV(1,1,:); % szz
    Bpsv(1+4*(e_f-1)+1,1,:) = PSV(2,1,:); % szr   ! delta
    B_sh(1+2*(e_f-1)  ,1,:) = SH(1,1,:);  % srr
    % fza vertical
    Bpsv(1+4*(e_f-1)  ,2,:) = PSV(1,2,:); % szz
    Bpsv(1+4*(e_f-1)+1,2,:) = PSV(2,2,:); % szr   ! delta
% interfaz de abajo **********************************
if (pXi.isoninterf)
    return
end
if (e_f == N+1)
    return
end
[PSV,SH] = Gij_dwn_atInterf(m(e_f),f,zloc(2));
    % fza horizontal
    Bpsv(1+4*(e_f-1)+2,1,:) = - PSV(3,1,:); %  uz
    Bpsv(1+4*(e_f-1)+3,1,:) = - PSV(4,1,:); %  ur
    B_sh(1+2*(e_f-1)+1,1,:) = - SH(2,1,:);  %  ur
    Bpsv(1+4*(e_f-1)+4,1,:) = - PSV(1,1,:); % szz
    Bpsv(1+4*(e_f-1)+5,1,:) = - PSV(2,1,:); % szr
    B_sh(1+2*(e_f-1)+2,1,:) = - SH(1,1,:);  % srr
    % fza vertical
    Bpsv(1+4*(e_f-1)+2,2,:) = - PSV(3,2,:); %  uz
    Bpsv(1+4*(e_f-1)+3,2,:) = - PSV(4,2,:); %  ur
    Bpsv(1+4*(e_f-1)+4,2,:) = - PSV(1,2,:); % szz
    Bpsv(1+4*(e_f-1)+5,2,:) = - PSV(2,2,:); % szr
