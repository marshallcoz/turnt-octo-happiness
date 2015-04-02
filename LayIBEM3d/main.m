%% DWN+IBEM 3d region E: estratificado, R: homogeneo
clear;clc;close('all','hidden')
%cd '/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/LayIBEM3d'
% del medio y la frecuencia
[m_vars,f_vars,ops,res] = setUpModelo;

% fuente
pXi.center(1:3) =[0 0 0];
[Uo] = ricker(f_vars,ops.ts,ops.tp);
%[Uo] = gaussiana(f_vars,ops.sigma);
dirFza = 3; pXi.normal(1:3) = [0 0 0];
pXi.normal(dirFza) = 3; 

% receptores
[res] = initreceptores(res);
% 
res.fotogramas = zeros(3,f_vars.ntiempo,res.nx,res.ny,res.nz);
mx =0;
FK = zeros(f_vars.nmax,f_vars.NFREC+1);
%%
Gij = @GijTij_HSestr_dwn;
for iPx = 1:res.nrecep
    [p_x] = pick_receptor(iPx,res,f_vars);
    
    for j=1:f_vars.NFREC+1
        [m_vars,f_vars] = setupJ(j,m_vars,f_vars,ops); 
        
        [G,T,FK(:,j)] = Gij(m_vars(1),f_vars,ops,p_x,pXi);
        p_x.greenG(1:3,1:3,j) = G;
        p_x.greenT(1:3,1:3,j) = T;
    end
    % pasar al tiempo
    for j=dirFza %1:3 %cada direccion de la fuerza
        s = zeros(3,f_vars.ntiempo);
        s(1:3,1:f_vars.NFREC+1) =  p_x.greenG(1:3,j,1:f_vars.NFREC+1);
        s(1:3,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));
        for i=1:3 %componente de la respuesta
            s(i,:)=s(i,:).*Uo;
            if ops.sacarSismogramas,sismogramaA(iPx,f_vars,s,i,j,p_x.center),end
            s(i,:) = real(ifft(s(i,:))/f_vars.dt);
            if f_vars.dwn==true, s(i,:) = s(i,:) .* exp(-f_vars.omei * f_vars.dt*(0:1:f_vars.ntiempo-1));end
            if ops.sacarSismogramas,sismogramaB(iPx,f_vars,s(i,:),i,j),end
        end
        res.fotogramas(1:3,1:f_vars.ntiempo,p_x.p(1),p_x.p(2),p_x.p(3)) ... 
                   = s(1:3,1:f_vars.ntiempo);
        mx = max(mx,max(s));
    end   
end
mx = mean(mx);