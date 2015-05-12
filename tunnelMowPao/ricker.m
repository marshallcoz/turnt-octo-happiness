% Ricker wavelet in frecuency
% En frecuencia se ha corregido que el tiempo inicial en ts
function [rick] = ricker(f,ts,tp,t0)
a = pi*(-ts)/tp;
rick = zeros(f.ntiempo,1);
rick(1) = (a*a-0.5)*exp(-a*a);
for i=2:f.ntiempo/2+1
    a = pi*(f.dt*(i-1)-ts)/tp;
    a = a*a;
    rick(i) = (a-0.5)*exp(-a);
end
h_rick = figure('Name','Funcion de amplitud'); set(gcf,'Visible', 'off');
subplot(2,1,1);
plot((0:f.ntiempo-1).*f.dt,rick);
title('Ondícula de Ricker')
xlabel('tiempo en segundos')
rick = fft(rick)*f.dt; %forward

tvec = zeros(f.ntiempo,1);
tvec(1) = exp(-1i*(0.01)*( 2*pi*f.DFREC)*(t0));
tvec(2:f.ntiempo/2+1) = exp(-1i*(1:f.ntiempo/2)*( 2*pi*f.DFREC)*(t0));
tvec(f.ntiempo/2+2:f.ntiempo) = conj(tvec(f.ntiempo/2:-1:2));
rick = rick .* tvec;

subplot(2,1,2);
title('Espectro de pulso con t_0 = -ts')
plot((0:f.ntiempo-1).*f.DFREC,real(rick),'r');hold on;
plot((0:f.ntiempo-1).*f.DFREC,imag(rick),'b');
plot((0:f.ntiempo-1).*f.DFREC,abs(rick),'k')
[ymM] = get(gca,'ylim');
line([1 1]*f.DFREC*f.NFREC,[ymM(1) ymM(2)])
xlim([0 f.ntiempo/2*f.DFREC]);
xlabel('frecuencia en Hertz')

% tvec = zeros(f.ntiempo,1);
% tvec(1) = exp(-1i*(0.01)*( 2*pi*f.DFREC)*(ts));
% tvec(2:f.ntiempo/2+1) = exp(-1i*(1:f.ntiempo/2)*( 2*pi*f.DFREC)*(ts));
% tvec(f.ntiempo/2+2:f.ntiempo) = conj(tvec(f.ntiempo/2:-1:2));
% rick = rick .* tvec;
% rick = ifft(rick)/f.dt; %backward
% subplot(2,1,1);hold on
% plot((0:f.ntiempo-1).*f.dt,real(rick),'g--')
%chdir out
saveas(h_rick,'amplitud.jpg');
%chdir ..
close(h_rick);
rick = rick';