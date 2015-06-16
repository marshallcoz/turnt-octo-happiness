function [receptor,Bou] = getGij_lay_psv_dwn(pota,receptor,Bou,...
                                                  Mpsv,M_sh,m,f,ops,mtrcs)
% obtner las tensores de Green entre segmentos:
% Bou.pt{i}.gG .gT
% y entre segmentos y receptores
% Bou.pt{i}.xG .xT
% en la frecuencia
% Allocate GreenFun variables
% now the green functions from each segment to each other segment
for i = 1:Bou.nBou
    Bou.pt{i}.gG = zeros(Bou.nBou,3,3); %will just valid each frequency
    Bou.pt{i}.gT = zeros(Bou.nBou,3,3); 
end
% now the green functions from each segment to each receiver
tam = size(receptor,1);
for i = 1:Bou.nBou
    Bou.pt{i}.xG = zeros(tam,3,3); % just valid each frequency
    Bou.pt{i}.xT = zeros(tam,3,3); 
end

nDeps = size(pota,1);
% para cada fuente virtual
for iDepF = 1:nDeps
    if (pota(iDepF,2) == 0)
        continue
    end
    for iFteSameZ = 1:pota(iDepF,2) %para cada fte a esa profundidad
    [pXi] = Bou.pt{-pota(iDepF,4+iFteSameZ)};
    
    % vector de terms indep para pXi:
    [Bpsv,B_sh] = Gij_B_dwn(m,f,ops,pXi);
    % coeficientes de ondas en cada estrato:
    for ik = 1:nmax+1
       Bpsv(:,1,ik) = Mpsv(:,:,ik) * Bpsv(:,1,ik); % horizontal
       Bpsv(:,2,ik) = Mpsv(:,:,ik) * Bpsv(:,2,ik); % vertical
       B_sh(:,1,ik) = M_sh(:,:,ik) * B_sh(:,1,ik); % horizontal
    end
    % campo difractado por estratigrafía:
    for iDepR = 1:nDeps %(todos los puntos en E son receptores)
        z = pota(iDepR,3);
        e = pota(iDepR,4);
        [mec] = campoDifractadoen(z,e,Bpsv,B_sh,m,f,ops.N,mtrcs);
        for iSameDep = 1:pota(iDepR,1)+pota(iDepR,2)
            if (pota(iDepR,4+iSameDep) > 0)
                p_X = receptor{pota(iDepR,4+iSameDep)};
            elseif (pota(iDepR,4+iSameDep) < 0)
                p_X = Bou.pt{-pota(iDepR,4+iSameDep)};
            else % 0
                break
            end
            % agregar fase horizontal r,th
            
            % integral en num de onda radial
            
            % agregar campo directo
            
        end
    end
    
    end %iFteSameZ
end %iDep