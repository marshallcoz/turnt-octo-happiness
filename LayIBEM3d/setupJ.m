function [m,f] = setupJ(J,m,f,ops)
f.frec = f.DFREC*max([J-1 0.01]);
f.ome = 2*pi*f.frec;
f.come = f.ome;

if f.dwn==true %amortiguamiento de las fuentes periódicas/y los polos
    f.omei = -1*pi/f.TW;
    f.come = f.ome + 1i * f.omei;
else
    f.omei = 0;
end
f.come = f.come * (1 - 1i/(2*f.Q)); %amortiguamiento histerético
disp(['J' num2str(J) ' : ' num2str(f.come)])
nm = size(m,2);
for im = 1:nm
if ops.UseAzimi
    m(im).alfa = m(im).alfa0*(1+1/pi/f.Q*log(f.come/2/pi))-1i/2/f.Q;
    m(im).beta = m(im).beta0*(1+1/pi/f.Q*log(f.come/2/pi))-1i/2/f.Q;
else
    m(im).alfa = m(im).alfa0;
    m(im).beta = m(im).beta0;
end
m(im).amu = m(im).rho * m(im).beta^2;
m(im).lambda = m(im).rho * m(im).alfa^2 - 2*m(im).amu;
m(im).bealf = m(im).beta / m(im).alfa;
m(im).bealf2 = (m(im).bealf)^2;
m(im).bealf3 = (m(im).bealf)^3;
m(im).bealf4 = (m(im).bealf)^4;
m(im).anu = (m(im).bealf2 - 0.5)/(-1 + m(im).bealf2);
m(im).lam = m(im).rho * m(im).alfa^2 - 2* m(im).amu;
end

