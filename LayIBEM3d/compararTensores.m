i = 3;
j = 3;
name = ['T' num2str(i) num2str(j) 'anaVSdwn.jpg'];
som = p_x_frec.greenT(i,j,:); 
%sdm = p_x_dwn.greenT(i,j,:); 
sdm = p_x.greenT(i,j,:);
so = zeros(f_vars.NFREC,1);
sd = so;

for ij = 1:f_vars.NFREC
    so(ij) = som(ij);
    sd(ij) = sdm(ij);
end
cociente = real(so)./real(sd) + 1i*imag(so)./imag(sd);
disp(mean(cociente(170:200)))

han = figure;
subplot(5,1,1:4)
hold on;
plot((0:f_vars.NFREC-1).*f_vars.DFREC,real(so),'r');
plot((0:f_vars.NFREC-1).*f_vars.DFREC,imag(so),'b');
% plot((0:f_vars.NFREC-1).*f_vars.DFREC,abs(so),'k');
plot((0:f_vars.NFREC-1).*f_vars.DFREC,real(sd),'r--');
plot((0:f_vars.NFREC-1).*f_vars.DFREC,imag(sd),'b--');
% plot((0:f_vars.NFREC-1).*f_vars.DFREC,abs(sd),'k--');
%xlim([0 3])
xlabel('frecuencia en Hertz')
ylabel('amplitud del espectro')
title('Comparación de G_{ij} analítico (---) y con DWN (- -)')
subplot(5,1,5)
hold on;
plot((0:f_vars.NFREC-1).*f_vars.DFREC,real(so-sd),'r');
plot((0:f_vars.NFREC-1).*f_vars.DFREC,imag(so-sd),'b');
xlabel('frecuencia en Hertz')
ylabel('diferencia')
saveas(han,name);
close(han);
%%