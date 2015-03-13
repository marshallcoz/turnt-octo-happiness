function [r] = distancia(p_x,pXi)
r = sqrt(sum((p_x.center(1:3)-pXi.center(1:3)).^2));