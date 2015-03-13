function [t] = tracciones(s,n)
% t_i = \sigma_{ij} n_j
t = zeros(3,1);
for i=1:3
    t(i) = dot(s(i,1:3)',n(1:3));
end