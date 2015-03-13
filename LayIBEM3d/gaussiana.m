% Ricker wavelet in frecuency
% En frecuencia se ha corregido que el tiempo inicial en ts
function [Uo] = gaussiana(f,sig)
s = sig/100 * f.ntiempo/2;
Uo = zeros(f.ntiempo,1);

Uo(1) = exp(-0.5*(0.01/s).^2);
Uo(2:f.ntiempo/2+1) = exp(-0.5*((1:f.ntiempo/2)/s).^2);
Uo(f.ntiempo/2+2:f.ntiempo) = conj(Uo(f.ntiempo/2:-1:2));

h_rick = figure('Name','Funcion de amplitud'); set(gcf,'Visible', 'off');
subplot(2,1,1);
plot((0:f.ntiempo-1).*f.DFREC,real(Uo),'r');hold on;
plot((0:f.ntiempo-1).*f.DFREC,imag(Uo),'b');
plot((0:f.ntiempo-1).*f.DFREC,abs(Uo),'k')
title('Espectro de pulso')
% xlim([0 f.ntiempo/2*f.DFREC]);
xlim([0 f.ntiempo*f.DFREC]);
xlabel('frecuencia en Hertz')

Ut = fft(Uo)/f.dt; %backward
subplot(2,1,2);
plot((0:f.ntiempo-1).*f.dt,real(Ut));
title('Ondícula de amplitud gaussiana')
xlabel('tiempo en segundos')

saveas(h_rick,'amplitud.jpg');
close(h_rick);
Uo = Uo';