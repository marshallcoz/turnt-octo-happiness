function [Bou,Bmat,Bvec,med] = initBmat_allCont(res,Bou)
%Si toda la matriz es de continuidad
Bmat = zeros(3*2*Bou.nBou);
% 3 : componentes de fuerza y desplazamiento
% 2 : tracciones y desplazamientos
% n : fuentes/puntos de colocación
Bvec = zeros(3*2*Bou.nBou,3);
med = 3*Bou.nBou;

% now the green functions from each segment to each station
for i = 1:Bou.nBou
    Bou.pt{i}.gG = zeros(res.nrecep,3,3); %just valid each frequency
    Bou.pt{i}.gT = zeros(res.nrecep,3,3); 
end
