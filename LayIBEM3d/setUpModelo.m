function [m,f,ops,res] = setUpModelo
% de la frecuencia y n�mero de onda
f.DFREC = 0.2;
f.NFREC = 200;
f.TW = 5;
f.Q = 1000;
f.dwn = true;
f.ntiempo = 2^10; % cantidad de puntos en tiempo
f.nk = 100;
f.nmax = 2^8;
f.Dk_R = 2.5;

% Medios E estratificado
ops.N = 2; % Num de estratos
clear m
m(1).z = 0; 
m(1).beta0 = 0.5; % [m/s]
m(1).alfa0 = 1; % [m/s]
m(1).rho = 1.0; % [T/m3]

m(2).z = 2;
m(2).beta0 = 0.5; % [m/s]
m(2).alfa0 = 1; % [m/s]
m(2).rho = 1.0; % [T/m3]

% semiespacio
m(ops.N+1).z = 3;
m(ops.N+1).beta0 = 0.5; % [m/s]
m(ops.N+1).alfa0 = 1; % [m/s]
m(ops.N+1).rho = 1.0; % [T/m3]

%Medio R homogeneo
m(ops.N+2).beta0 = 0.5; % [m/s]
m(ops.N+2).alfa0 = 1; % [m/s]
m(ops.N+2).rho = 1.0; % [T/m3]

ops.UseAzimi = false;
ops.sacarSismogramas = true;
ops.sacarFotoramas = false;
ops.tp = 0.3; %ricker (ancho)  / frec central= 1/tp 
ops.ts = 1.5; %ricker (centro)

%Receptores -----------------------------------------
clear res
res.Del = 0.1; %espacio entre receptores [m]
res.box = [[1 1];...
           [1 1];...
           [1 1]];
res.norm =[1 1 1];

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% otras variables :: 
f.dk = 2*pi*(f.DFREC * f.NFREC)/(min(m(1:2).beta0) * f.nmax) * f.Dk_R;
f.dx = pi / (f.nmax * f.dk);
f.L = 2*pi/f.dk;
f.dt = 1 / (f.ntiempo * f.DFREC);
f.tmax = f.dt * f.ntiempo; % = 1/f.DFREC

