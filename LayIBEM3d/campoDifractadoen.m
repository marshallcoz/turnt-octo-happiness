function [mec] = campoDifractadoen(z,e,Bpsv,B_sh,m,f,N,mtrcs)
egammaN = zeros(f.nmax+1);
egammaP = egammaN;
enuN    = egammaN;
enuP    = egammaN;
coeffsPSV = zeros(4,f.nmax+1);
coeffs_SH = zeros(2,f.nmax+1);
alf = m(e).alfa;
bet = m(e).beta;
gam =((ome/alf)^2 - k2).^0.5;
nu = ((ome/bet)^2 - k2).^0.5;
[found] = find(imag(gam)>0);  gam(found) = -(gam(found)); clear found
[found] = find(imag(nu)>0);    nu(found) = -(nu(found));  clear found

% fase vertical
egammaN(:) = exp(-1i * gam * (z - m(e).z));
enuN(:) = exp(-1i * nu *  (z - m(e).z));
if (e ~= N+1)
    egammaP(:) = exp(1i * gam * (z - m(e+1).z));
    enuP(:) = exp(1i * nu * (z - m(e+1).z));
else
    egammaP(:) = 0;
    enuP(:) = 0;
end
% coeficientes de las ondas
if (e == N+1) % semiespacio de abajo
    %horizontal
    coeffsPSV(1:2,1,:) = Bpsv(4*N+1:4*N+2,1,:);
    coeffsPSV(3:4,1,:) = 0;
    coeffs_SH(1,1,:) = B_sh(2*N+1,1,:);
    coeffs_SH(2,1,:) = 0;
    %vertical
    coeffsPSV(1:2,2,:) = Bpsv(4*N+1:4*N+2,2,:);
    coeffsPSV(3:4,2,:) = 0;
else% estrato
    coeffsPSV(1:4,1,:) = Bpsv(4*(e-1)+1:4*(e-1)+4,1,:);
    coeffsPSV(1:4,2,:) = Bpsv(4*(e-1)+1:4*(e-1)+4,2,:);
    coeffs_SH(1:2,1,:) = B_sh(2*(e+1)+1:2*(e+1)+2,1,:);
end

% calcular elementos mecánicos sin fase horizontal en todos los k
mec.h.psv.u = zeros(3,f.nmax+1);
mec.h.psv.s = zeros(6,f.nmax+1);
mec.h.sh.u = zeros(2,f.nmax+1);
mec.h.sh.s = zeros(5,f.nmax+1);
mec.v.psv.u = zeros(3,f.nmax+1);
mec.v.psv.s = zeros(6,f.nmax+1);
for ik = 1:f.nmax+1
    % Fza Horizontal
    % uz ur ut
    mec.h.psv.u(:,ik) =  ...
        mtrcs{1,1}.sD0(:,1,ik,e).* (egammaN(ik)*coeffsPSV(1,1,ik))+...
        mtrcs{1,1}.sD0(:,2,ik,e).* (enuN(ik)*   coeffsPSV(2,1,ik))+...
        mtrcs{1,1}.sD0(:,3,ik,e).* (egammaP(ik)*coeffsPSV(3,1,ik))+...
        mtrcs{1,1}.sD0(:,4,ik,e).* (enuP(ik)*   coeffsPSV(4,1,ik));
    
    % szz szr szt srr srt stt
    mec.h.psv.s(:,ik) =  ...
        mtrcs{1,1}.sS0(:,1,ik,e).* (egammaN(ik)*coeffsPSV(1,1,ik))+...
        mtrcs{1,1}.sS0(:,2,ik,e).* (enuN(ik)*   coeffsPSV(2,1,ik))+...
        mtrcs{1,1}.sS0(:,3,ik,e).* (egammaP(ik)*coeffsPSV(3,1,ik))+...
        mtrcs{1,1}.sS0(:,4,ik,e).* (enuP(ik)*   coeffsPSV(4,1,ik));
    %SH:
    % ur ut
    mec.h.sh.u(:,ik) = ...
        mtrcs{1,1}.sD0sh(:,1,ik,e).* (enuN(ik)*coeffs_SH(1,1,ik))+...
        mtrcs{1,1}.sD0sh(:,2,ik,e).* (enuP(ik)*coeffs_SH(2,1,ik));
    
    % szr stt srr srt stt
    mec.h.sh.s(:,ik) = ...
        mtrcs{1,1}.sS0sh(:,1,ik,e).* (enuN(ik)*coeffs_SH(1,1,ik))+...
        mtrcs{1,1}.sS0sh(:,2,ik,e).* (enuP(ik)*coeffs_SH(2,1,ik));
    
    % Fza Vertical
    % uz ur ut
    mec.v.psv.u(:,ik) =  ...
        mtrcs{1,1}.sD0(:,1,ik,e).* (egammaN(ik)*coeffsPSV(1,2,ik))+...
        mtrcs{1,1}.sD0(:,2,ik,e).* (enuN(ik)*   coeffsPSV(2,2,ik))+...
        mtrcs{1,1}.sD0(:,3,ik,e).* (egammaP(ik)*coeffsPSV(3,2,ik))+...
        mtrcs{1,1}.sD0(:,4,ik,e).* (enuP(ik)*   coeffsPSV(4,2,ik));
    
    % szz szr szt srr srt stt
    mec.v.psv.s(:,ik) =  ...
        mtrcs{1,1}.sS0(:,1,ik,e).* (egammaN(ik)*coeffsPSV(1,2,ik))+...
        mtrcs{1,1}.sS0(:,2,ik,e).* (enuN(ik)*   coeffsPSV(2,2,ik))+...
        mtrcs{1,1}.sS0(:,3,ik,e).* (egammaP(ik)*coeffsPSV(3,2,ik))+...
        mtrcs{1,1}.sS0(:,4,ik,e).* (enuP(ik)*   coeffsPSV(4,2,ik));
end
