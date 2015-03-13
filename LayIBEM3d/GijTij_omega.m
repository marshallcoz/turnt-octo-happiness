%% funcion de Green 3d Espacio completo analitica
function [G,T,nil] = GijTij_omega(m,f,p_x,pXi)
nil = zeros(1,f.nmax);
come = f.come;
r = distancia(p_x,pXi);
gamma = (p_x.center - pXi.center)./r;

n = p_x.normal;
k = come/m.beta;
q = come/m.alfa;

ba = m.beta/m.alfa;
ba2 = ba^2;
qr = q * r;
kr = k * r;
kr2 = kr^2;
f1 = ba2*(1-2i/(qr)-2/(qr)^2)*exp(-1i*qr) + ...
    (2i/(kr)+2/(kr)^2)*exp(-1i*kr);
f2 = ba2*(1i/(qr)+1/(qr)^2)*exp(-1i*qr) + ...
    (1-1i/(kr)-1/(kr)^2)*exp(-1i*kr);

G = zeros(3);
for i=1:3
    for j=1:3
        G(i,j) = (f2*Kron(i,j) + (f1 - f2)*gamma(i)*gamma(j))/(4*pi*m.amu*r);
    end
end

A=zeros(2,3);B=A;C=A;D=A;
A(1,:) = [0 0 -1i];
A(2,:) = [-1i*ba 1i*(2*ba^3-ba) 0];
B(1,:) = [4 -2 -3];
B(2,:) = [-4*ba2-1 4*ba2-1 2*ba2];
C(1,:) = [-12i 6i 6i];
C(2,:) = [12i*ba -6i*ba -6i*ba];
D(1,:) = [-12 6 6];
D(2,:) = [12 -6 -6];

gj = zeros(3,1);
for j=1:3
    gj(j) = (kr*A(1,j) + B(1,j) + C(1,j)/kr + D(1,j)/kr2)*exp(-1i*kr) + ...
        (kr*A(2,j) + B(2,j) + C(2,j)/kr + D(2,j)/kr2)*exp(-1i*qr);
end
T = zeros(3);
for i=1:3
    for j=1:3
        T(i,j) = ((gj(1)-gj(2)-2*gj(3))*gamma(i)*gamma(j)*dot(gamma,n) + ...
            gj(3)*gamma(i)*n(j) + gj(2)*gamma(j)*n(i) + gj(3)*dot(gamma,n)*Kron(i,j))...
            /(4*pi*r*r);
    end
end

