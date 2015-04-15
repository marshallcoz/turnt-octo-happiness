function plotphi(j,Bou,phiVec,dir)
% plot the problem geometry
hfig = figure('Name','Force density');
maxphi = max(abs(phiVec(:,dir)));
subplot(1,2,1);hold on
phi = phiVec(1:size(phiVec(:,dir),1)/2)/maxphi/2;
l_m=1:3;
for i = 1:Bou.nBou
    plotCircle3D(Bou.pt{i}.center, ...
        Bou.pt{i}.normal, ...
        Bou.pt{i}.radio,...
        phi(l_m));
    l_m = l_m+3;
end
view(45,36)
maxx = max(get(gca,'xlim'));
maxy = max(get(gca,'ylim'));
maxz = max(get(gca,'zlim'));
minx = min(get(gca,'xlim'));
miny = min(get(gca,'ylim'));
minz = min(get(gca,'zlim'));

axis equal
xlim([minx maxx])
ylim([miny maxy])
zlim([minz maxz])
xlabel('x')
ylabel('y')
title('{\phi} E')
%
set(gca,'Box','off')
subplot(1,2,2);hold on
phi = phiVec(size(phiVec(:,dir),1)/2+1:size(phiVec(:,dir),1))/maxphi/2;
l_m=1:3;
for i = 1:Bou.nBou
    plotCircle3D(Bou.pt{i}.center, ...
        Bou.pt{i}.normal, ...
        Bou.pt{i}.radio,...
        phi(l_m));
    l_m = l_m+3;
end

view(45,36)
maxx = max(get(gca,'xlim'));
maxy = max(get(gca,'ylim'));
maxz = max(get(gca,'zlim'));
minx = min(get(gca,'xlim'));
miny = min(get(gca,'ylim'));
minz = min(get(gca,'zlim'));

axis equal
xlim([minx maxx])
ylim([miny maxy])
zlim([minz maxz])
xlabel('x')
ylabel('y')
title('{\phi} R')
%
set(gca,'Box','off')
cd out
cd phi
hgsave(hfig,['FigBoundary_frec' num2str(j) '_fza' num2str(dir) '.fig']);
cd ..
cd ..
close(hfig)

