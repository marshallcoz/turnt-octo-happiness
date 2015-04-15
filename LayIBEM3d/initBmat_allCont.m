function [Bmat,Bvec,med] = initBmat_allCont(Bou)
%Si toda la matriz es de continuidad
Bmat = zeros(3*2*Bou.nBou);
% 3 : componentes de fuerza y desplazamiento
% 2 : tracciones y desplazamientos
% n : fuentes/puntos de colocación
Bvec = zeros(3*2*Bou.nBou,3);
med = 3*Bou.nBou;

