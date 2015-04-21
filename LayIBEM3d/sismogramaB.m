% Parte de abajo del sismograma. Seña en tiempo
function sismogramaB(iRecep,f_vars,s,comp,dirFza)
tiem = (0:(f_vars.ntiempo-1))'.*f_vars.dt;
vec = zeros(f_vars.ntiempo,1);
name = ['G' num2str(comp) num2str(dirFza) '_(' num2str(iRecep) ')'];
h_figura = gcf;
hold on;
subplot(2,1,2)
for t = 1:f_vars.ntiempo
    vec(t) = s(t);
end
plot(tiem,vec,'k')
xlabel('tiempo en segundos')
cd out
cd Sis
%saveas(h_figura,[name '.jpg']);
hgsave(h_figura,[name '.fig']);
cd ..
cd ..
close(h_figura);