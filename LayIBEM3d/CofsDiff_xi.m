function [ CoefpsvH,CoefpsvV,Coefsh ] = CofsDiff_xi(pXi,Mpsv,M_sh,m,ops,f,k)
% Hacer el vector de la fuente pXi y encontrar los coeficientes
% del campo difractado en la estratificación.
zXi = pXi.center(3);
[eXi,atInterf] = el_estrato_es(m,ops.N,zXi);
BpsvH = zeros(4*ops.N+2,f.nmax+1); %Horz
BpsvV = zeros(4*ops.N+2,f.nmax+1); %Vert
B_sh = zeros(2*ops.N+1,f.nmax+1);
if (atInterf) % si la fuente está en la interfaz (y estrato) eXi
    % caso P-SV
    % para fuerza vertical
        BpsvV(1+4*(eXi-1)  ,:) = -1/(2*pi)* k; %szz
    % para fuerza horizontal
        BpsvH(1+4*(eXi-1)+1,:) = 1i/(2*pi)* k; %szr
    % caso SH
    % para fuerza horizontal
        B_sh(2*eXi-1    ,:,1) =  1/(2*pi); %szt
else % la fuente está dentro del estrato
    % la fuente (en expansión en WN) se propaga a las interfaces adyacentes 
    % para la interfaz de arriba (prop hacia abajo)
    z_loc = m(eXi).z - zXi;
    [PSV,SH] = Gij_dwn_atInterf(m(eXi),f,k,z_loc);
     % interfaz de arriba
     % para fuerza horizontal 
     if(eXi ~= 1)
        BpsvH(1+4*(eXi-1)-2,:) = PSV(3,1,:);%Uz
        BpsvH(1+4*(eXi-1)-1,:) = PSV(4,1,:);%Ur
        B_sh(1+2*(eXi-1)-1,:) = SH(2,1,:); %Ur
     end
        BpsvH(1+4*(eXi-1)  ,:) = PSV(1,1,:);%szz
        BpsvH(1+4*(eXi-1)+1,:) = PSV(2,1,:);%szr
        B_sh(1+2*(eXi-1)   ,:) = SH(1,1,:); %szr
     % para fuerza vertical
     if(eXi ~= 1)
        BpsvV(1+4*(eXi-1)-2,:) = PSV(3,2,:);%Uz
        BpsvV(1+4*(eXi-1)-1,:) = PSV(4,2,:);%Ur
     end    
        BpsvV(1+4*(eXi-1)  ,:) = PSV(1,2,:);%szz
        BpsvV(1+4*(eXi-1)+1,:) = PSV(2,2,:);%szr
      
     if (eXi ~= ops.N+1)    
     % para la interfaz de abajo 
        z_loc = m(eXi+1).z - zXi;
        [PSV,SH] = Gij_dwn_atInterf(m(eXi),f,k,z_loc);
        % para fuerza horizontal
        BpsvH(1+4*(eXi-1)+2,:) = -PSV(1,1,:);%szz
        BpsvH(1+4*(eXi-1)+3,:) = -PSV(2,1,:);%szr
        BpsvH(1+4*(eXi-1)+4,:) = -PSV(3,1,:);%Uz
        BpsvH(1+4*(eXi-1)+5,:) = -PSV(4,1,:);%Ur
        B_sh(1+2*(eXi-1)+1,:) = -SH(2,1,:); %ur
        B_sh(1+2*(eXi-1)+2,:) = -SH(1,1,:); %szt
        % para fuerza vertical
        BpsvV(1+4*(eXi-1)+2,:) = -PSV(1,2,:);%szz
        BpsvV(1+4*(eXi-1)+3,:) = -PSV(2,2,:);%szr
        BpsvV(1+4*(eXi-1)+4,:) = -PSV(3,2,:);%Uz
        BpsvV(1+4*(eXi-1)+5,:) = -PSV(4,2,:);%Ur
     end
end % atInterf

%% Resolver el sistema Aa = B -> coeficientes
CoefpsvH = zeros(4*ops.N+2,1,f.nmax+1); %Horz
Coefsh = zeros(2*ops.N+1,1,f.nmax+1);
CoefpsvV = zeros(4*ops.N+2,1,f.nmax+1); %Vert
for ik=1:f.nmax+1
    CoefpsvH(:,1,ik) = Mpsv(:,:,ik)*BpsvH(:,ik); %Horz
    Coefsh(:,1,ik) = M_sh(:,:,ik)*B_sh(:,ik); 
    CoefpsvV(:,1,ik) = Mpsv(:,:,ik)*BpsvV(:,ik); %Vert
end
