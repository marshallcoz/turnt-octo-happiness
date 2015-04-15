function [m,f,ops,res] = setUpModelo
% de la frecuencia y número de onda
f.DFREC = 0.05;
f.NFREC = 250;
f.TW = 5;
f.Q = 1000;
f.dwn = true;
f.ntiempo = 2^10; % cantidad de puntos en tiempo
f.nk = 100;
f.nmax = 2^8;
f.Dk_R = 2.5;



% Medios E estratificado
ops.N = 2; % Num de estratos
clear m                          % E
m(1).z = 0; 
m(1).beta0 = 3; % [m/s]
m(1).alfa0 = 6; % [m/s]
m(1).rho = 1.0; % [T/m3]

m(2).z = 2;                      % R
m(2).beta0 = 3; % [m/s]
m(2).alfa0 = 6; % [m/s]
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
ops.sigma = 10.0; %gaussiana (% de fmax)

%Receptores -----------------------------------------
clear res
% BoundaryFile
%res.BouFile = 'Ico5Unit.txt';
%res.BouFile = 'cuboUni.txt';
res.BouFile = 'cuboU24.txt';

res.Del = 0.5; %espacio entre receptores [m]
res.box = [[-1 1];...
           [2 2];...
           [0 0]];
res.norm =[1 1 1];

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% otras variables :: 
f.dk = 2*pi*(f.DFREC * f.NFREC)/(min(m(1:2).beta0) * f.nmax) * f.Dk_R;
f.dx = pi / (f.nmax * f.dk);
f.L = 2*pi/f.dk;
f.dt = 1 / (f.ntiempo * f.DFREC);
f.tmax = f.dt * f.ntiempo; % = 1/f.DFREC

