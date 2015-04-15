% La mitad de arriba del sismograma, con la frecuencia
function sismogramaA(iRecep,f_vars,s,comp,dirFza,center)
name = ['G' num2str(comp) num2str(dirFza) '_(' num2str(iRecep) ')'];
figure('Name',name);hold on; set(gcf,'Visible', 'off');
subplot(2,1,1)
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,real(s(comp,:)),'r');hold on;
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,imag(s(comp,:)),'b');
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,abs(s(comp,:)),'k')
xlim([0 f_vars.ntiempo/2*f_vars.DFREC]);
txcoor = ['x=(' num2str(center) ')'];
xlabel('Frecuencia en Hertz')
title(txcoor);