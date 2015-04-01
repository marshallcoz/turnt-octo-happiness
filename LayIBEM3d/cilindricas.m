function [r,th,thp,z] = cilindricas(p_x,pXi)
r = sqrt(sum((p_x.center(1:2)-pXi.center(1:2)).^2));
% desde el eje x
th = atan((p_x.center(2)-pXi.center(2))/(p_x.center(1)-pXi.center(1)));
% desde el eje y
thp = th-pi/2;
z = p_x.center(3)-pXi.center(3);