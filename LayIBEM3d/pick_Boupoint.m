function [p_x] = pick_Boupoint(i,Bou)
% Tomar uno de los puntos de colocación
p_x = Bou.pt{i};
% p_x.center(1:3) = Bou.pt{i}.center;
% p_x.normal(1:3) = Bou.pt{i}.normal;
% p_x.radio = Bou.pt{i}.radio;