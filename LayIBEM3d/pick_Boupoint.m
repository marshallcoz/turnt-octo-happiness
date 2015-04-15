function [p_x] = pick_Boupoint(i,Bou)
% Tomar uno de los puntos de colocación
p_x.center(1:3) = Bou.pt{i}.center;
p_x.normal(1:3) = Bou.pt{i}.normal;
p_x.radio = Bou.pt{i}.radio;
% p_x.greenG = zeros(3,3,f_vars.NFREC); 
% p_x.greenT = zeros(3,3,f_vars.NFREC); 
% p_x.sism = zeros(3,3,f_vars.ntiempo);