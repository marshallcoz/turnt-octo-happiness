%% programa principal
%clear;clc;close('all','hidden')
cd '/Users/marshall/Documents/DOC/unit6/LayIBEM3d'
% del medio y la frecuencia
[m_vars,f_vars,ops,res] = setUp;

%fuente
pXi.center(1:3) =[0 0 0];
%[Uo] = ricker(f_vars,ops.ts,ops.tp);
[Uo] = gaussiana(f_vars,ops.sigma);
%
dirFza = 1;
pXi.normal(1:3) =[1 0 0];
mx =0;

% if f_vars.dwn==true
     Gij = @GijTij_dwn;
% else
%     Gij = @GijTij_omega;
% end

%receptores
[res] = initreceptores(res);
res.fotogramas = zeros(3,f_vars.ntiempo,res.nx,res.ny,res.nz);
FK = zeros(f_vars.nmax,f_vars.NFREC+1);
for iRecep = 1:res.nrecep
    [p_x] = pick_receptor(iRecep,res);
    p_x.greenG = zeros(3,3,f_vars.NFREC); 
    p_x.greenT = zeros(3,3,f_vars.NFREC); 
    p_x.sism = zeros(3,3,f_vars.ntiempo);
    
    for j=1:f_vars.NFREC+1
        [m_vars,f_vars] = setupJ(j,m_vars,f_vars,ops); 
        [G,T,FK(:,j)] = Gij(m_vars(1),f_vars,p_x,pXi);
        p_x.greenG(1:3,1:3,j) = G;
        p_x.greenT(1:3,1:3,j) = T;
    end
    % pasar al tiempo
    for j=dirFza %1:3 %cada direccion de la fuerza
        s = zeros(3,f_vars.ntiempo);
        s(1:3,1:f_vars.NFREC+1) =  p_x.greenG(1:3,j,1:f_vars.NFREC+1);
        s(1:3,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));
        for i=1:3
            s(i,:)=s(i,:).*Uo;
            if ops.sacarSismogramas,sismogramaA(iRecep,f_vars,s,i,j,p_x.center),end
            s(i,:) = real(ifft(s(i,:))/f_vars.dt);
            if f_vars.dwn==true, s(i,:) = s(i,:) .* exp(-f_vars.omei * f_vars.dt*(0:1:f_vars.ntiempo-1));end
            if ops.sacarSismogramas,sismogramaB(iRecep,f_vars,s(i,:),i,j),end
        end
        res.fotogramas(1:3,1:f_vars.ntiempo,p_x.p(1),p_x.p(2),p_x.p(3)) ... 
                   = s(1:3,1:f_vars.ntiempo);
        mx = max(mx,max(s));
    end   
end
mx = mean(mx);
%% graficar fotogramas
if ops.sacarFotoramas
disp('haciendo fotograms')
contadort=1;
for t=461:10:1024%161:10:1024%f_vars.ntiempo
    h_figura = figure('Name',['Gi1_' num2str(t)]);hold on; set(gcf,'Visible', 'off')
    quiver3(pXi.center(1),pXi.center(2),pXi.center(3),...
        pXi.normal(1),pXi.normal(2),pXi.normal(3),1,'r');
    for iRecep = 1:res.nrecep
        [p_x] = pick_receptor(iRecep,res);
        n = res.fotogramas(1:3,t,p_x.p(1),p_x.p(2),p_x.p(3));
        mag = comprimir(magnitud(n(1:3)),mx);
        if abs(mag) < 0.001/mx
            continue
        else
            quiver3(p_x.center(1),p_x.center(2),p_x.center(3),...
            n(1),n(2),n(3),mag,'k');
        end
    end
    view(-35,65)
    xlim(res.box(1,:));
    ylim(res.box(2,:));
    zlim([0 2]);
    xlabel('x');ylabel('y');zlabel('z')
    title(['t=' num2str(t*f_vars.dt) 'sec']);
    saveas(h_figura,['F' num2str(dirFza) '_' num2str(contadort) '.jpg']);
    contadort = contadort+1;
    close('all','hidden')
end
end
