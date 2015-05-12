%%
clear
%z = 100.02142048715710 -5.0000002374872565E-004i;
%z = 299.47273931896569 -5.0000002374872565E-004i;
%z = 5.7736110083589434E+10 +0.0000i;
%z = 7.5574971632834371 -3.7787484971801066E-002i;
%z = 45 -0.3;
%z = 9.0689962445798342E-002 -4.5344980209359865E-004i;
z = 0.90689962445798342 -4.5344980209359865E-003i;
%z = 45 + 2i;
%n = -1:1:10;
n = 0:0.25:10;
%
J = besselj(n,z);
figure;
plot(n,real(J),'r-'); hold on
plot(n,imag(J),'b-')
%
Y = bessely(n(),z);
figure;
plot(n,real(Y),'r-'); hold on
plot(n,imag(Y),'b-')
%%
h1 = besselh(n,1,z);
figure;
plot(n,real(h1),'r-'); hold on
plot(n,imag(h1),'b-')
%%
h2 = besselh(n,2,z);
figure;
plot(n,real(h2),'r-'); hold on
plot(n,imag(h2),'b-')