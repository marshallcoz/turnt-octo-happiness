function [m,f,ops,res] = setUp
% Variables de frecuencia 
f.DFREC = 0.2;
f.NFREC = 200;
f.TW = 5;
f.Q = 1000;
f.dwn = true;
f.ntiempo = 2^10; % cantidad de puntos en tiempo
f.nmax = 2^11;
multDk = 1.1;

% Medio E homogeneo
m(1).beta0 = 0.5; %m/s
m(1).alfa0 = 1; %m/s
m(1).rho = 1.0; % T/m3

%Medio R homogeneo
m(2).beta0 = 0.5; %m/s
m(2).alfa0 = 1; %m/s
m(2).rho = 1.0; % T/m3

ops.UseAzimi = false;
ops.sacarSismogramas = true;
ops.sacarFotoramas = true;
ops.tp = 0.3; %ricker (ancho)  / frec central= 1/tp 
ops.ts = 1.5; %ricker (centro)
ops.sigma = 10.0; %gaussiana (% de fmax)

%Receptores
clear res
res.Del = 0.1;
res.box = [[-1 1];...
           [-1 1];...
           [1 1]];
res.norm =[1 1 1];

f.dk = 2*pi*(f.DFREC * f.NFREC)/(min(m(1:2).beta0) * f.nmax) * multDk;
f.dx = pi / (f.nmax * f.dk);
f.L = 2*pi/f.dk;
f.dt = 1 / (f.ntiempo * f.DFREC);
f.tmax = f.dt * f.ntiempo; % = 1/f.DFREC

