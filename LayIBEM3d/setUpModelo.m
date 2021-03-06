function [m,f,ops,res,p0] = setUpModelo
% de la frecuencia y n�mero de onda
f.DFREC = 0.05;
f.NFREC = 200;
f.TW = 5;
f.Q = 1000;
f.dwn = true;
f.ntiempo = max(2*f.NFREC,2^8); % cantidad de puntos en tiempo
f.nk = 200;
f.nmax = 2^8;
f.Dk_R = 2.5;

% Medios E estratificado
ops.N = 2; % Num de estratos
clear m                          % E
m(1).z = 0; 
m(1).beta0 = 1; % [m/s]
m(1).alfa0 = 2; % [m/s]
m(1).rho = 1.0; % [T/m3]

m(2).z = 0.5;                      % R
m(2).beta0 = 3; % [m/s]
m(2).alfa0 = 6; % [m/s]
m(2).rho = 1.0; % [T/m3]

% semiespacio
m(ops.N+1).z = 3;
m(ops.N+1).beta0 = 7; % [m/s]
m(ops.N+1).alfa0 = 14; % [m/s]
m(ops.N+1).rho = 1.0; % [T/m3]

%Medio R homogeneo
m(ops.N+2).beta0 = 3; % [m/s]
m(ops.N+2).alfa0 = 6; % [m/s]
m(ops.N+2).rho = 1.0; % [T/m3]
plotProfile(m,ops.N);

ops.UseAzimi = false;
ops.sacarSismogramas = false;
ops.sacarFotoramas = true;
ops.sacarPhiplot = false;
ops.tp = 0.3; %ricker (ancho)  / frec central= 1/tp 
ops.ts = 1.5; %ricker (centro)
ops.t0 = -1.5; %tiempo inicial
ops.sigma = 10.0; %gaussiana (% de fmax)

%Receptores -----------------------------------------
clear res
% BoundaryFile
%res.BouFile = 'Ico5Unit.txt';
% res.BouFile = 'icoshpere2.blend.txt';
% res.BouFileZ0 = 'icoshpere2.blendz0.0.txt';
%res.BouFile = 'icosphere3.blend.txt';
res.BouFileCont = 'mediaEsfera.blend.txt';
%res.BouFileZ0 = 'icosphere3.blendz0.0.txt';
res.BouFileZ0 = 'mediaEsfera.blendz0.0.txt';

res.Del = 0.2; %espacio entre receptores [m]
res.box = [[0 2];...
           [-2 2];...
           [0 0]];
res.norm =[1 1 1];

% fuente real :::::::::::::::::::::::::::::::
p0.center(1:3) =[3 -3 0];
p0.region = 1;

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% otras variables :: 
f.dk = 2*pi*(f.DFREC * f.NFREC)/(min(m(1:2).beta0) * f.nmax) * f.Dk_R;
f.dx = pi / (f.nmax * f.dk);
f.L = 2*pi/f.dk;
f.dt = 1 / (f.ntiempo * f.DFREC);
f.tmax = f.dt * f.ntiempo; % = 1/f.DFREC

