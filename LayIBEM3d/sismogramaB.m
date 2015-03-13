% Parte de abajo del sismograma. Seña en tiempo
function sismogramaB(iRecep,f_vars,s,i,j)
    tiem = (0:(f_vars.ntiempo-1))'.*f_vars.dt;
    vec = zeros(f_vars.ntiempo,1);
        name = ['G' num2str(i) num2str(j) '_(' num2str(iRecep) ')'];
        h_figura = gcf;
        hold on;
        subplot(2,1,2)
        for t = 1:f_vars.ntiempo
            vec(t) = s(t);
        end
        plot(tiem,vec,'k')
        xlabel('tiempo en segundos')
        saveas(h_figura,[name '.jpg']);
        close(h_figura);