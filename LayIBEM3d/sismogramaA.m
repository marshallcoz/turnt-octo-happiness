% La mitad de arriba del sismograma, con la frecuencia
function sismogramaA(iRecep,f_vars,s,i,j,center)
    name = ['G' num2str(i) num2str(j) '_(' num2str(iRecep) ')'];
            figure('Name',name);hold on; set(gcf,'Visible', 'off');
            subplot(2,1,1)
            plot((0:f_vars.ntiempo-1).*f_vars.DFREC,real(s(i,:)),'r');hold on;
            plot((0:f_vars.ntiempo-1).*f_vars.DFREC,imag(s(i,:)),'b');
            plot((0:f_vars.ntiempo-1).*f_vars.DFREC,abs(s(i,:)),'k')
            xlim([0 f_vars.ntiempo/2*f_vars.DFREC]);
     txcoor = ['x=(' num2str(center) ')'];
     xlabel('Frecuencia en Hertz')
     title(txcoor);