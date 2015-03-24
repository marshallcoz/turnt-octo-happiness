%% funcion de Green 3D semiespacio estratificado
function [G,T,FK] = GijTij_HSestr_dwn(m,f,ops,p_x,pXi)
%m = m_vars;f=f_vars;
%% FK
FK = zeros(1,f.nmax);
%% se resuelve en cilíndricas y se transforma a cartesianas
G = zeros(3);
T = zeros(3);
[r,th,thp,z] = cilindricas(p_x,pXi);
k = 0:f.dk:f.dk*f.nmax; k(1) = 0.01*f.dk;


%% matriz de coeficientes polarizaciones P-SV
Mpsv = zeros(4*ops.N+2,4*ops.N+2);

