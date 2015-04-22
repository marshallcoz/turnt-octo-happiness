%% programa principal IBEM INCLUSION HOMOGENEO
clear;clc;close('all','hidden')
cd '/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/LayIBEM3d'
%
rmdir('out','s')
mkdir out
mkdir out/phi
mkdir out/Sis
mkdir out/vid

% del medio y la frecuencia
[m_vars,f_vars,ops,res,p0] = setUpModelo;

[Uo] = ricker(f_vars,ops.ts,ops.tp,ops.t0);
%[Uo] = gaussiana(f_vars,ops.sigma);

% receptores y puntos de colocación
[res,receptor] = initreceptores(res,f_vars);

% if f_vars.dwn==true
%     Gij = @GijTij_dwn;
% else
Gij = @GijTij_omega;
GijX = @GijTij_omega_pxEQpXi;
% end
GE = zeros(3); TE = GE; G = GE; T = GE; mx =0;

[Bou] = initBoundary(res);
[xli,yli,zli] = plotGeom(1,Bou,res,receptor,p0);
phiVec = zeros(f_vars.NFREC,2*3*Bou.nBou,3);

% Frequency loop
jini = 1;
jfin = f_vars.NFREC;
%%
disp('incident field')
for j=jini:jfin
    [m,f] = setupJ(j,m_vars,f_vars,ops);
    for iPx = 1:res.nrecep %for each station
        [p_x] = receptor{iPx};%pick_receptor(iPx,res);
        if (p_x.region == p0.region)
            [G,T] = Gij(m(p0.region),f,p_x,p0);
            receptor{iPx}.greenT(:,:,j) = T(:,:);
            receptor{iPx}.greenG(:,:,j) = G(:,:);
        end
    end
end
%%
disp('diffracted field')
for j=jini:jfin
    [m,f] = setupJ(j,m_vars,f_vars,ops);
    
    [Bou,Bmat,Bvec,med] = initBmat_allCont(res,Bou);
    
    % fill the continuity conditions matrix
    j_m = 1:3;
    for iPxi = 1:Bou.nBou %for each BE (as a source) {column on E and column on R}
        [pXi] = pick_Boupoint(iPxi,Bou);
        DA = pi*pXi.radio^2;
        i_m = 1:3;
        for iPx = 1:Bou.nBou %for each BE (as a receiver at the boundary) {same region as pXi}
            [p_x] = pick_Boupoint(iPx,Bou);
            if (iPxi == iPx)
                % GreenEx y 0.5
                [GE,TE] = GijX(m(1),f,p_x,pXi); %(region E column)
                [GR,TR] = GijX(m(2),f,p_x,pXi); %(region R column)
                GE=GE*DA;GR=GR*DA;
                TR = -TR;
            else
                % full space Green function :
                [GE,TE] = Gij(m(1),f,p_x,pXi); %(region E column)
                [GR,TR] = Gij(m(2),f,p_x,pXi); %(region R column)
                GE=GE*DA;TE=TE*DA;GR=GR*DA;TR=TR*DA;
            end
            % fill the matrix
            Bmat(i_m,j_m) = TE;
            Bmat(i_m+med,j_m) = GE;
            Bmat(i_m,j_m+med) = -TR;
            Bmat(i_m+med,j_m+med) = -GR;
            i_m = i_m+3;
        end
        for iPx = 1:res.nrecep %for each station
            [p_x] = receptor{iPx};%pick_receptor(iPx,res);
            r = distancia(p_x,pXi);
            me = p_x.region; %= '1->E'
            if (r <= 0.7*pXi.radio)
                [G,T] = GijX(m(me),f,p_x,pXi);
                Bou.pt{iPxi}.gG(iPx,:,:) = G*DA;
                Bou.pt{iPxi}.gT(iPx,:,:) = T*DA;
            else
                % full space Green function:
                [G,T] = Gij(m(me),f,p_x,pXi);
                Bou.pt{iPxi}.gG(iPx,:,:) = G*DA;
                Bou.pt{iPxi}.gT(iPx,:,:) = T*DA;
            end
        end
        j_m = j_m+3;
    end %iPxi
    % fill independent term vector
    i_m = 1:3;
    for iPx = 1:Bou.nBou %for each BE (as a receiver at the boundary) {same region as pXi}
        [p_x] = pick_Boupoint(iPx,Bou);
        [G,T] = Gij(m(p0.region),f,p_x,p0); %(on region E)
        for dirFza=1:3 % 3 direcciones de la fuerza
            Bvec(i_m,dirFza) = -T(:,dirFza);
            Bvec(i_m+med,dirFza) = -G(:,dirFza);
        end
        i_m = i_m+3;
    end
    
    % solve diffracted field
    
    for dirFza=1:3 % 3 direcciones de la fuerza
        phiVec(j,:,dirFza) = linsolve(Bmat,Bvec(:,dirFza));
        if(ops.sacarPhiplot),plotphi(j,Bou,squeeze(phiVec(j,:,dirFza)),dirFza);end
        
        % diffracted/refracted field
        for iPx = 1:res.nrecep %for each station
            [p_x] = receptor{iPx};%pick_receptor(iPx,res);
            lmed = 0;
            if (p_x.region == 2),lmed = 1;end
            j_m = 1:3;
            for iPxi = 1:Bou.nBou %for each BE
                phirange = j_m+lmed*med;
                [pXi] = pick_Boupoint(iPxi,Bou);
                aux(1:3,1:3) = pXi.gG(iPx,:,:);
                receptor{iPx}.greenG(:,dirFza,j) = ...
                receptor{iPx}.greenG(:,dirFza,j) + ...
                    aux * phiVec(j,phirange,dirFza);
                aux(1:3,1:3) = pXi.gT(iPx,:,:);
                receptor{iPx}.greenT(:,dirFza,j) = ...
                receptor{iPx}.greenT(:,dirFza,j) + ...
                    aux * phiVec(phirange,dirFza);
                j_m = j_m+3;
            end
        end
    end
end %frequency loop

%% pasar al tiempo
disp('Pasar al tiempo')
for iPx = 1:res.nrecep %for each station
%     [p_xt] = pick_receptor(iPx,res);
    for dirFza=1:3 % dir fza
        s = zeros(3,f_vars.ntiempo);
        s(1:3,1:f_vars.NFREC) =  receptor{iPx}.greenG(1:3,dirFza,1:f_vars.NFREC);
        s(1:3,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));
        for comp=1:3 % componente
            thiS=s(comp,:);
            s(comp,:)=thiS(:).*Uo(:);
            if ops.sacarSismogramas,sismogramaA(iPx,f_vars,s,comp,dirFza,receptor{iPx}.center),end
            s(comp,:) = real(ifft(s(comp,:))/f_vars.dt);
            if ops.sacarSismogramas,sismogramaB(iPx,f_vars,s(comp,:),comp,dirFza),end
            if ops.sacarFotoramas
                res.fotogramas(dirFza,comp,1:f_vars.ntiempo,receptor{iPx}.p(1),receptor{iPx}.p(2),receptor{iPx}.p(3)) ...
                    = s(comp,1:f_vars.ntiempo);
                mx = max(mx,max(s));
            end
        end
        
    end
end



%% graficar fotogramas
if ops.sacarFotoramas
    disp('haciendo fotograms')
    mx = abs(mean(mx));
    for dirFza=1:3 % dir fza
        p0.normal(1:3) = 0; p0.normal(dirFza) = 1;
%         contadort=1;
        for t=30 :1:30%f_vars.ntiempo
            h_figura = figure('Name',['Foto_dirFza' num2str(dirFza) '_' num2str(t)]);hold on; set(gcf,'Visible', 'off')
            for i = 1:Bou.nBou
                plotPatch(Bou.pt{i}.vert,0.2);
            end
            quiver3(p0.center(1),p0.center(2),p0.center(3),...
                p0.normal(1),p0.normal(2),p0.normal(3),0.1*(xli(2)-xli(1)),'r','MarkerSize',25);
            for iPx = 1:res.nrecep
                %         [p_xt] = pick_receptor(iRecep,res);
                n = res.fotogramas(dirFza,1:3,t,receptor{iPx}.p(1),receptor{iPx}.p(2),receptor{iPx}.p(3));
                magN = magnitud(n(1:3));
                mag = real(comprimir(magN,mx));
%                 if (magN < 0.0001/mx)
%                     continue
%                 end
                n = n/magN;
%                 disp([n mag])
%                 if abs(mag) < 0.001/mx
%                     continue
%                 else
                    try
                    quiver3(receptor{iPx}.center(1),receptor{iPx}.center(2),receptor{iPx}.center(3),...
                        n(1),n(2),n(3),0.00050*mag,'k');
                    catch ER
                        
                    end
%                 end
            end
            view(45,36)
            light('Position',[xli(2)*2 yli(1)*2 zli(2)*5],'Style','infinite');
            axis equal
            xlim(xli)
            ylim(yli)
            zlim(zli)
            xlabel('x');ylabel('y');zlabel('z')
            title(['t=' num2str(t*f_vars.dt) 'sec']);
            set(gca,'Box','off')
            cd out
            cd vid
            saveas(h_figura,['F_fza' num2str(dirFza) '_' num2str(t) '.jpg']);
            cd ..
            cd ..
%             contadort = contadort+1;
            close('all','hidden')
        end
    end
end