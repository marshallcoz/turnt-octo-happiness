%% programa principal IBEM INCLUSION HOMOGENEO
%clear;clc;close('all','hidden')
%cd '/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/LayIBEM3d'
% del medio y la frecuencia
[m_vars,f_vars,ops,res] = setUpModelo;
mkdir out
mkdir out/phi
%fuente real
p0.center(1:3) =[0 -1.5 0];
p0.region = 1;
[Uo] = ricker(f_vars,ops.ts,ops.tp);
%[Uo] = gaussiana(f_vars,ops.sigma);

% receptores y puntos de colocación
[res] = initreceptores(res,f_vars);

% if f_vars.dwn==true
%     Gij = @GijTij_dwn;
% else
Gij = @GijTij_omega;
GijX = @GijTij_omega_pxEQpXi;
% end
GE = zeros(3); TE = GE; G = GE; T = GE; mx =0;
% Frequency loop
jini = 1;
jfin = f_vars.NFREC;
for j=jini:jfin
    [m,f] = setupJ(j,m_vars,f_vars,ops);
    % incident field
    regP0 = p0.region;
    mregP0 = m(regP0);
%     resoutT=zeros(res.nrecep,3,3);
%     resoutG=zeros(res.nrecep,3,3);
    for iPx = 1:res.nrecep %for each station
        [p_x] = pick_receptor(iPx,res);
        if (p_x.region == regP0)
            [G,T] = Gij(mregP0,f,p_x,p0);
%             resoutT(iPx,:,:) = T(:,:);
%             resoutG(iPx,:,:) = G(:,:);
            res.receptor{iPx}.greenT(:,:,j) = T(:,:);
            res.receptor{iPx}.greenG(:,:,j) = G(:,:);
        end
    end
%     for iPx = 1:res.nrecep
%         res.receptor{iPx}.greenT(:,:,j) =resoutT(iPx,:,:);
%         res.receptor{iPx}.greenG(:,:,j) =resoutG(iPx,:,:);
%     end
    
    [Bou] = initBoundary(res);
    if (j==jini),plotGeom(j,Bou,res,p0);end
    [Bmat,Bvec,med] = initBmat_allCont(Bou);
    
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
            [p_x] = pick_receptor(iPx,res);
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
    phiVec = zeros(2*med,3);
    for dirFza=1:3 % 3 direcciones de la fuerza
        phiVec(:,dirFza) = linsolve(Bmat,Bvec(:,dirFza));
        if(ops.sacarPhiplot),plotphi(j,Bou,phiVec,dirFza);end
        
        % diffracted/refracted field
        for iPx = 1:res.nrecep %for each station
            [p_x] = pick_receptor(iPx,res);
            lmed = 0;
            if (p_x.region == 2),lmed = 1;end
            j_m = 1:3;
            for iPxi = 1:Bou.nBou %for each BE
                phirange = j_m+lmed*med;
                [pXi] = pick_Boupoint(iPxi,Bou);
                aux(1:3,1:3) = pXi.gG(iPx,:,:);
                res.receptor{iPx}.greenG(:,dirFza,j) = ...
                res.receptor{iPx}.greenG(:,dirFza,j) + ...
                    aux * phiVec(phirange,dirFza);
                aux(1:3,1:3) = pXi.gT(iPx,:,:);
                res.receptor{iPx}.greenT(:,dirFza,j) = ...
                res.receptor{iPx}.greenT(:,dirFza,j) + ...
                    aux * phiVec(phirange,dirFza);
                j_m = j_m+3;
            end
        end
    end
end %frequency loop

%% pasar al tiempo
for iPx = 1:res.nrecep %for each station
    [p_x] = pick_receptor(iPx,res);
    
    for dirFza=1:3 % dir fza
        s = zeros(3,f_vars.ntiempo);
        s(1:3,1:f_vars.NFREC) =  p_x.greenG(1:3,dirFza,1:f_vars.NFREC);
        s(1:3,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));
        for comp=1:3 % componente
            thiS=s(comp,:);
            s(comp,:)=thiS(:).*Uo(:);
            if ops.sacarSismogramas,sismogramaA(iPx,f_vars,s,comp,dirFza,p_x.center),end
            s(comp,:) = real(ifft(s(comp,:))/f_vars.dt);
            if ops.sacarSismogramas,sismogramaB(iPx,f_vars,s(comp,:),comp,dirFza),end
        end
    end
%     if ops.sacarFotoramas
%         res.fotogramas(1:3,1:f_vars.ntiempo,p_x.p(1),p_x.p(2),p_x.p(3)) ...
%             = s(1:3,1:f_vars.ntiempo);
%         mx = max(mx,max(s));
%     end
end
% mx = mean(mx);


%% graficar fotogramas
% if ops.sacarFotoramas
% disp('haciendo fotograms')
% contadort=1;
% for t=461:10:1024%161:10:1024%f_vars.ntiempo
%     h_figura = figure('Name',['Gi1_' num2str(t)]);hold on; set(gcf,'Visible', 'off')
%     quiver3(pXi.center(1),pXi.center(2),pXi.center(3),...
%         pXi.normal(1),pXi.normal(2),pXi.normal(3),1,'r');
%     for iRecep = 1:res.nrecep
%         [p_x] = pick_receptor(iRecep,res);
%         n = res.fotogramas(1:3,t,p_x.p(1),p_x.p(2),p_x.p(3));
%         mag = comprimir(magnitud(n(1:3)),mx);
%         if abs(mag) < 0.001/mx
%             continue
%         else
%             quiver3(p_x.center(1),p_x.center(2),p_x.center(3),...
%             n(1),n(2),n(3),mag,'k');
%         end
%     end
%     view(-35,65)
%     xlim(res.box(1,:));
%     ylim(res.box(2,:));
%     zlim([0 2]);
%     xlabel('x');ylabel('y');zlabel('z')
%     title(['t=' num2str(t*f_vars.dt) 'sec']);
%     saveas(h_figura,['F' num2str(dirFza) '_' num2str(contadort) '.jpg']);
%     contadort = contadort+1;
%     close('all','hidden')
% end
% end