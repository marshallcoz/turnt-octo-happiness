% funcion de Green 3D semiespacio estratificado
function [Mpsv,M_sh,cellMat] = GijTij_HSestr_dwn(m,f,ops)
% Matriz global de continuidad P-SV y SH a cada k
k = 0:f.dk:f.dk*f.nmax; k(1) = 0.01*f.dk;
k2 = k.^2;
Z0 = k*0;
ome = f.come;

% matrices de coeficientes para P-SV y para SH
Mpsv = zeros(4*ops.N+2,4*ops.N+2,f.nmax+1);
M_sh = zeros(2*ops.N+1,2*ops.N+1,f.nmax+1);
ega = zeros(2,1,f.nmax+1); enu = zeros(2,1,f.nmax+1);
sD0 = zeros(3,4,f.nmax+1,ops.N+1); sS0 = zeros(6,4,f.nmax+1,ops.N+1);
sD = zeros(2,4,f.nmax+1);  sS = zeros(2,4,f.nmax+1);
sD0sh = zeros(2,2,f.nmax+1,ops.N+1); sS0sh = zeros(5,2,f.nmax+1,ops.N+1);
sDsh = zeros(1,2,f.nmax+1); sSsh = zeros(1,2,f.nmax+1);

iR= 0;iC= 0; ir= 0; ic= 0;
for e = 1:ops.N+1
    alf = m(e).alfa;
    bet = m(e).beta;
    amu = m(e).amu;
% %        ! Se debe cumplir que la parte imaginaria del n�mero de onda 
% %        ! vertical debe ser menor que cero. La parte imaginaria contri-
% %        ! buye ondas planas inhomog�neas con decaimiento expo-
% %        ! nencial a medida que z crece.
    gam =((ome/alf)^2 - k2).^0.5;
    nu = ((ome/bet)^2 - k2).^0.5;
    [found] = find(imag(gam)>0);  gam(found) = -(gam(found)); clear found
    [found] = find(imag(nu)>0);  nu(found) = -(nu(found)); clear found
    % figure;plot(real(gam),'r');hold on;plot(imag(gam),'b')
    % figure;plot(real(nu),'r');hold on;plot(imag(nu),'b')
    dKga = 2*k.*gam;
    dKnu = 2*k.*nu;
    xi = k.^2 - nu.^2;
    % P-SV
    sD0(1,1,:,e)=-gam; sD0(1,2,:,e)=-k; sD0(1,3,:,e)=gam; sD0(1,4,:,e)=-k;  %Uz
    sD0(2,1,:,e)=k;   sD0(2,2,:,e)=-nu; sD0(2,3,:,e)=k;  sD0(2,4,:,e)=nu;   %Ur
    sD0(3,1,:,e)=-k;   sD0(3,2,:,e)=nu; sD0(3,3,:,e)=-k;  sD0(3,4,:,e)=-nu; %Ut
    % SH
    sD0sh(1,1,:,e)= k; sD0sh(1,2,:,e)= -k; %Ur
    sD0sh(2,1,:,e)= -k; sD0sh(2,2,:,e)= -k; %Ut
    % P-SV
    %  subMatS0 = [[xi    -dKnu xi   dKnu];  szz
    %              [-dKga -xi   dKga -xi ]]; szr
    %        EPD                 ESD             EPU                ESU
    sS0(1,1,:,e)=xi;    sS0(1,2,:,e)=-dKnu; sS0(1,3,:,e)=xi;    sS0(1,4,:,e)=dKnu; %szz
    sS0(2,1,:,e)=-dKga; sS0(2,2,:,e)=-xi;   sS0(2,3,:,e)=dKga;  sS0(2,4,:,e)=-xi;  %szr
    sS0(3,1,:,e)= dKga; sS0(3,2,:,e)= xi;   sS0(3,3,:,e)=-dKga; sS0(3,4,:,e)=xi;   %szt
    sS0(1:3,:,:,e) = sS0(1:3,:,:,e) * amu;
    % (nota: los siguientes ocupan los coeficientes corregidos y ya est�n
    % corregidos)
    sS0(4,1,:,e)=1;   sS0(4,2,:,e)=-1i*nu; sS0(4,3,:,e)=1;   sS0(4,4,:,e)=1i*nu; % srr
    sS0(5,1,:,e)=1;   sS0(5,2,:,e)=-1i*nu; sS0(5,3,:,e)=1;   sS0(5,4,:,e)=1i*nu; % srt
    sS0(6,1,:,e)=1;   sS0(6,2,:,e)=-1i*nu; sS0(6,3,:,e)=1;   sS0(6,4,:,e)=1i*nu; % stt
    %sS0(4:6,:,:,e) = sS0(4:6,:,:,e) * 2*amu;
    % SH
    sS0sh(1,1,:,e)=  amu*nu.*k; sS0sh(1,2,:,e)= amu*nu.*k; %szr
    sS0sh(2,1,:,e)= -amu*nu.*k; sS0sh(2,2,:,e)=-amu*nu.*k; %szt
    % (nota idem)
    sS0sh(3,1,:,e)=1 ; sS0sh(3,2,:,e)=1; %srr
    sS0sh(4,1,:,e)=1 ; sS0sh(4,2,:,e)=1; %srt
    sS0sh(5,1,:,e)=1 ; sS0sh(5,2,:,e)=1; %stt
    
% ! la profundidad z de la frontera superior del estrato
% !         z_i = Z(e)   ! e=1  ->  z = z0 = 0
% !         z_f = Z(e+1) ! e=1  ->  z = z1 
    for bord = 0:1
        if ((e == 0) && (bord == 0)), continue,end
        if (e+bord > ops.N+1), break,end 
          % si 1+0;1+1;2+0;[2+1] > 2
        if (bord == 0) %---->---->---->---->---->---->
           if (e ~= ops.N+1) % (radiation condition lower HS)
               ega(1,1,:) = exp(-1i * gam * (m(e+1).z - m(e).z));
               enu(1,1,:) = exp(-1i * nu *  (m(e+1).z - m(e).z));
           else
              ega(1,1,:) = Z0;
              enu(1,1,:) = Z0;
           end
           ega(2,1,:)=ega(1,1,:);enu(2,1,:)=enu(1,1,:);
           % P-SV
           sD(:,1,:) = sD0(1:2,1,:,e);
           sD(:,2,:) = sD0(1:2,2,:,e);
           sD(:,3,:) = sD0(1:2,3,:,e) .* ega;
           sD(:,4,:) = sD0(1:2,4,:,e) .* enu;
           sS(:,1,:) = sS0(1:2,1,:,e);
           sS(:,2,:) = sS0(1:2,2,:,e);
           sS(:,3,:) = sS0(1:2,3,:,e) .* ega;
           sS(:,4,:) = sS0(1:2,4,:,e) .* enu;
           % SH
           sDsh(:,1,:) = sD0sh(1:1,1,:,e);
           sDsh(:,2,:) = sD0sh(1:1,2,:,e) .* enu(1,1,:);
           sSsh(:,1,:) = sS0sh(1:1,1,:,e);
           sSsh(:,2,:) = sS0sh(1:1,2,:,e) .* enu(1,1,:);
        else %bord .eq. 1 -----------------------------------
            if (e ~= 0) %(radiation condition upper HS)
                ega(1,1,:) = exp(-1i * gam * (m(e+1).z - m(e).z));
                enu(1,1,:) = exp(-1i * nu *  (m(e+1).z - m(e).z));      
            else
                ega(1,1,:) = Z0;
                enu(1,1,:) = Z0;
            end  %<----<----<----<----<----<----<----<----<--
            ega(2,1,:)=ega(1,1,:);enu(2,1,:)=enu(1,1,:);
            sD(:,1,:) = sD0(1:2,1,:,e).* ega;
            sD(:,2,:) = sD0(1:2,2,:,e).* enu;
            sD(:,3,:) = sD0(1:2,3,:,e);
            sD(:,4,:) = sD0(1:2,4,:,e);
            sS(:,1,:) = sS0(1:2,1,:,e).* ega;
            sS(:,2,:) = sS0(1:2,2,:,e).* enu;
            sS(:,3,:) = sS0(1:2,3,:,e);
            sS(:,4,:) = sS0(1:2,4,:,e);
           % SH
           sDsh(:,1,:) = sD0sh(1:1,1,:,e) .* enu(1,1,:);
           sDsh(:,2,:) = sD0sh(1:1,2,:,e);
           sSsh(:,1,:) = sS0sh(1:1,1,:,e) .* enu(1,1,:);
           sSsh(:,2,:) = sS0sh(1:1,2,:,e);
        end
        %ensamble de la macro columna i
        %evaluadas en el borde SUPERIOR del layer i
        if ((bord == 0) && (e ~= 0))  
            if (e ~= ops.N+1) %(radiation condition)
                if (e ~= 1) %(only stress bound.cond. in the surface
                    Mpsv( iR-1 : iR , iC+1 : iC+4 ,:) = -sD;
                    M_sh( ir   : ir , ic+1 : ic+2 ,:) = -sDsh;
                end
                Mpsv( iR+1 : iR+2 , iC+1 : iC+4,:) = -sS;
                M_sh( ir+1 : ir+1 , ic+1 : ic+2,:) = -sSsh;
            else
                if (e ~= 1) %(only stress bound.cond. in the surface
                    Mpsv( iR-1 : iR , iC+1 : iC+2,:) = -sD  (:,1:2,:);
                    M_sh( ir   : ir , ic+1 : ic+1,:) = -sDsh(:,1:1,:);
                end 
                Mpsv( iR+1 : iR+2 , iC+1 : iC+2,:) = -sS  (:,1:2,:);
                M_sh( ir+1 : ir+1 , ic+1 : ic+1,:) = -sSsh(:,1:1,:);
            end 
       end 
       %evaluadas en el borde INFERIOR del layer i
       if ((bord == 1) && (e ~= ops.N+1)) % cond de radiaci�n en HS
           if (e ~= 0) %(radiation condition upward)
               Mpsv( iR+3 : iR+4 , iC+1 : iC+4,:) = sD;
               Mpsv( iR+5 : iR+6 , iC+1 : iC+4,:) = sS;
               M_sh( ir+2 : ir+2 , ic+1 : ic+2,:) = sDsh;
               M_sh( ir+3 : ir+3 , ic+1 : ic+2,:) = sSsh;
           else
               Mpsv( iR+3 : iR+4 , iC+3 : iC+4,:) = sD  (:,3:4,:);
               Mpsv( iR+5 : iR+6 , iC+3 : iC+4,:) = sS  (:,3:4,:);
               M_sh( ir+2 : ir+2 , ic+2 : ic+2,:) = sDsh(:,2:2,:);
               M_sh( ir+3 : ir+3 , ic+2 : ic+2,:) = sSsh(:,2:2,:);
           end
       end
    end % bord
    iR= iR+4;
    iC= iC+4;
    ir= ir+2;
    ic= ic+2;
end % e

 % Guardar matrices para la posteridad
cellMat = cell(1);
cellMat{1,1}.sD0 = {sD0};
cellMat{1,1}.sS0 = {sS0};
cellMat{1,1}.sD0sh = {sD0sh};
cellMat{1,1}.sS0sh = {sS0sh};
cellMat{1,1}.k = {k};

% take the inverse of the global matrices
for ik=1:f.nmax+1
    Mpsv(:,:,ik) = inv(Mpsv(:,:,ik));
    M_sh(:,:,ik) = inv(M_sh(:,:,ik));
end
