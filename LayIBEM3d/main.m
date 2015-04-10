%% DWN+IBEM 3d region E: estratificado, R: estratificado
clear;clc;close('all','hidden')
%cd '/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/LayIBEM3d'
% del medio y la frecuencia
[m_vars,f_vars,ops,res] = setUpModelo;

% fuente real en E

%[Uo] = ricker(f_vars,ops.ts,ops.tp);
%[Uo] = gaussiana(f_vars,ops.sigma);

% receptores y puntos de colocación
[res] = initreceptores(res,f_vars);
[Bou] = initBoundary;

for J=1:f_vars.NFREC+1
[m_vars,f_vars] = setupJ(J,m_vars,f_vars,ops);

clear cellMat
% Las matrices glovales ya invertidas y en todo su esplendor
[Mpsv,M_sh,cellMat] = GijTij_HSestr_dwn(m_vars,f_vars,ops);

%*****************************************
%        Fuente real en región E         %
pXi.center(1:3) =[1.5 0 0];              %
dirFza = 3; pXi.normal(1:3) = [0 0 1];   %
%*****************************************

[ CoefpsvH,CoefpsvV,Coefsh ] = CofsDiff_xi(pXi,Mpsv,...
    M_sh,m_vars,ops,f_vars,cellMat{1,1}.k{1,1});
% Obtener campo incidente en receptores
for iPx = 1:res.nrecep
    [p_x] = pick_receptor(iPx,res);
    %if (p_x.region ~= 'E'), continue, end
    [G,T] = Diff_x(cellMat,CoefpsvH,CoefpsvV,Coefsh,...
                   p_x,pXi,m_vars,f_vars,ops,dirFza);
    [res] = pack_receptor(J,iPx,res,G,T);
end
% Obtener campo incidente en puntos de colocación (term indep IBEM)
%for iPx = 1:Bou.nBou
%    [p_x] = pick_Boupoint(1,Bou); % punto de colocación
%    
%end

end %frec loop
% 
%res.fotogramas = zeros(3,f_vars.ntiempo,res.nx,res.ny,res.nz);
%mx =0;
%FK = zeros(f_vars.nmax,f_vars.NFREC+1);


% 
% 
% 
% 
% for iPx = 1:res.nrecep
%     [p_x] = pick_receptor(iPx,res,f_vars);
%     
%     for j=1:f_vars.NFREC+1
%         [m_vars,f_vars] = setupJ(j,m_vars,f_vars,ops); 
%         
%         [G,T,FK(:,j)] = Gij(m_vars(1),f_vars,ops,p_x,pXi);
%         p_x.greenG(1:3,1:3,j) = G;
%         p_x.greenT(1:3,1:3,j) = T;
%     end
%     % pasar al tiempo
%     for j=dirFza %1:3 %cada direccion de la fuerza
%         s = zeros(3,f_vars.ntiempo);
%         s(1:3,1:f_vars.NFREC+1) =  p_x.greenG(1:3,j,1:f_vars.NFREC+1);
%         s(1:3,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));
%         for i=1:3 %componente de la respuesta
%             s(i,:)=s(i,:).*Uo;
%             if ops.sacarSismogramas,sismogramaA(iPx,f_vars,s,i,j,p_x.center),end
%             s(i,:) = real(ifft(s(i,:))/f_vars.dt);
%             if f_vars.dwn==true, s(i,:) = s(i,:) .* exp(-f_vars.omei * f_vars.dt*(0:1:f_vars.ntiempo-1));end
%             if ops.sacarSismogramas,sismogramaB(iPx,f_vars,s(i,:),i,j),end
%         end
%         res.fotogramas(1:3,1:f_vars.ntiempo,p_x.p(1),p_x.p(2),p_x.p(3)) ... 
%                    = s(1:3,1:f_vars.ntiempo);
%         mx = max(mx,max(s));
%     end   
% end
% mx = mean(mx);