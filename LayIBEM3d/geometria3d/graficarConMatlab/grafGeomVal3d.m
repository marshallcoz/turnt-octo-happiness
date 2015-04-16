%% graficar geometría
cd /Users/marshall/Documents/DOC/unit5/bugfree-sansa/tridilibem/graficarConMatlab
%open 'geom.txt';
G = importdata('geom.txt');
n = size(G,1);
figure;hold
for i=1:n
   plotCircle3D(G(i,1:3),G(i,4:6),G(i,7));
end
view(113,-14)
axis equal
disp(['radio maximo= ' num2str(max(G(:,7)))])
disp(['radio promedio= ' num2str(mean(G(:,7)))])
%% graficar observadores
cd /Users/marshall/Documents/DOC/unit5/bugfree-sansa/tridilibem
dat = importdata('VALL3D.txt');
n = str2double(dat.textdata{8,1});
px = zeros(n,3);
for i = 10:9+n
    px(i-9,1) =  str2double(dat.textdata{i,1});
    px(i-9,2) =  str2double(dat.textdata{i,2});
    px(i-9,3) =  str2double(dat.textdata{i,3});
end
scatter3(px(:,1),px(:,2),px(:,3),30*ones(n,1),'filled','g');
%% graficar en frecuencia
cd /Users/marshall/Documents/DOC/unit5/bugfree-sansa/tridilibem
dat = importdata('OndaP1.txt');
nf = 10; % la frecuencia nf 0,1,2
resini = 103;
resfin = n;
renglon = 2 + 4*nf;
fila = dat(renglon,:)';
ux = fila(1:2:end) + 1i*fila(2:2:end);
fila = dat(renglon+1,:)';
uy = fila(1:2:end) + 1i*fila(2:2:end);
fila = dat(renglon+2,:)';
uz = fila(1:2:end) + 1i*fila(2:2:end);

quiver3(px(resini:resfin,1),px(resini:resfin,2),px(resini:resfin,3),...
    real(ux(resini:resfin)),real(uy(resini:resfin)),real(uz(resini:resfin)),'r');
quiver3(px(resini:resfin,1),px(resini:resfin,2),px(resini:resfin,3),...
    imag(ux(resini:resfin)),imag(uy(resini:resfin)),imag(uz(resini:resfin)),'b');
        