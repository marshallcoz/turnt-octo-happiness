%% funcion de Green 3D semiespacio estratificado
function [G,T,FK] = GijTij_HSestr_dwn(m,f,ops,p_x,pXi)
%m = m_vars;f=f_vars;
%% FK
FK = zeros(1,f.nmax);
%% se resuelve en cilíndricas y se transforma a cartesianas
G = zeros(3);
T = zeros(3);
[r,th,thp,z] = cilindricas(p_x,pXi);
k = 0:f.dk:f.dk*f.nmax; k(1) = 0.01*f.dk;
k2 = k.^2;
k3 = k.^3;
kr = k*r;
ome = f.come;


% J0kr= besselj(0,kr);
% J1kr= besselj(1,kr);
% ega = exp(-1i*gam*abs(z));
% enu = exp(-1i*nu*abs(z));
% om2be2 = ome^2/bet^2;


%% matriz de coeficientes polarizaciones P-SV
Mpsv = zeros(4*ops.N+2,4*ops.N+2);
% subMatD = zeros(2,4); subMatS = zeros(2,4);
% subMatD0= zeros(2,4); subMatS0= zeros(2,4);
% diagMat = zeros(4,4);

iR= 0;iC= 0;
for e = 1:ops.N+1
    alf = m(e).alfa;
    bet = m(e).beta;
    mu = m(e).amu;
    lam = m(e).lam;
% %        ! Se debe cumplir que la parte imaginaria del número de onda 
% %        ! vertical debe ser menor que cero. La parte imaginaria contri-
% %        ! buye ondas planas inhomogéneas con decaimiento expo-
% %        ! nencial a medida que z crece.
    gam =((ome/alf)^2 - k2).^0.5;
    nu = ((ome/bet)^2 - k2).^0.5;
    [found] = find(imag(gam)>0);  gam(found) = -(gam(found)); clear found
    [found] = find(imag(nu)>0);  nu(found) = -(nu(found)); clear found
    % figure;plot(real(gam),'r');hold on;plot(imag(gam),'b')
    % figure;plot(real(nu),'r');hold on;plot(imag(nu),'b')
    nu2 = nu.^2;
    ga2 = gam.^2;
    k2_ga = k2./gam;
    dKga = 2*k.*gam;
    dKnu = 2*k.*nu;
    xi = k.^2 - nu.^2;
    subMatD0 = [[-gam -k  gam  -k ];
                [-k   nu  -k   -nu]] * 1i;
    
    subMatS0 = [[xi    -dKnu xi   dKnu];
                [-dKga -xi   dKga -xi ]] * mu;
    
            ! la profundidad z de la frontera superior del estrato
!         z_i = Z(e)   ! e=1  ->  z = z0 = 0
!         z_f = Z(e+1) ! e=1  ->  z = z1 
        do bord = 0,1
          if ((e .eq. 0) .and. (bord .eq. 0)) cycle
          if (e+bord > N+1) exit 
          ! si 1+0;1+1;2+0;[2+1] > 2
        if (bord .eq. 0) then !---->---->---->---->---->---->
          if (e /= N+1) then !(radiation condition lower HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e)))
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))
          else
            ega = Z0
            enu = Z0
          end if
          diagMat = RESHAPE((/ UR, Z0, Z0, Z0, & 
                               Z0, UR, Z0, Z0, & 
                               Z0, Z0, ega, Z0, & 
                               Z0, Z0, Z0, enu /), &
                           (/ 4,4 /))
        else !bord .eq. 1 -----------------------------------
          if (e /= 0) then !(radiation condition upper HS)
            ega = exp(-UI * gamma * (Z(e+1)-Z(e))) 
            enu = exp(-UI * nu * (Z(e+1)-Z(e)))       
          else
            ega = Z0
            enu = Z0
          end if !<----<----<----<----<----<----<----<----<--
          diagMat = RESHAPE((/ ega, Z0, Z0, Z0, & 
                               Z0, enu, Z0, Z0, & 
                               Z0, Z0, UR, Z0, & 
                               Z0, Z0, Z0, UR /), &
                           (/ 4,4 /))
        end 
          ! desplazamientos
          subMatD = matmul(subMatD0,diagMat)
          ! esfuerzos
          subMatS = matmul(subMatS0,diagMat)
        !ensamble de la macro columna i
          !evaluadas en el borde SUPERIOR del layer i
          if (bord == 0 .AND. e /= 0 ) then 
           if (e /= N+1) then !(radiation condition)
            if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+4 ) = -subMatD
            end 
             this_A( iR+1 : iR+2 , iC+1 : iC+4 ) = -subMatS
           else
           if (e/=i) then !(only stress bound.cond. in the surface
             this_A( iR-1 : iR   , iC+1 : iC+2 ) = -subMatD(:,1:2)
           end 
             this_A( iR+1 : iR+2 , iC+1 : iC+2 ) = -subMatS(:,1:2)
            ! exit
           end 
          end 
          
          !evaluadas en el borde INFERIOR del layer i
          if (bord == 1 .AND. e /= N+1 ) then ! cond de radiación en HS
           if (e /= 0) then !(radiation condition upward)
            ! ambas ondas
            this_A( iR+3 : iR+4 , iC+1 : iC+4 ) = subMatD
            this_A( iR+5 : iR+6 , iC+1 : iC+4 ) = subMatS
           else
            ! solo onda hacia arriba
            this_A( iR+3 : iR+4 , iC+3 : iC+4 ) = subMatD(:,3:4)
            this_A( iR+5 : iR+6 , iC+3 : iC+4 ) = subMatS(:,3:4)
           end 
          end 
        end  !bord loop del borde i superior o nferior
          iR= iR+4 
          iC= iC+4
          end ! !{e} loop de las macro columnas para cada estrato
          
