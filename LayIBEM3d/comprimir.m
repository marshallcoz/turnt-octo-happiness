function [s] = comprimir(s,mx)
p = 15; %8
s =  log((1. + exp(p)*abs(s))) / (log(exp(p)+1.));
%disp(s)
s = s / mx;