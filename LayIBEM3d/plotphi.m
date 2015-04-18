function plotphi(j,Bou,phiVec,dir)
% plot the problem geometry
hfig = figure('Name',['Force density given F' num2str(dir)]);
maxphi = max(abs(phiVec(:,dir)));
subplot(1,2,1);hold on
tam = size(phiVec(:,1),1)/2;
phi = phiVec(1:tam,dir)/maxphi/2;
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
mTextBox = uicontrol('style','text');
set(mTextBox,'String',...
    ['max= ' num2str(max(phiVec(1:tam,dir)))]);
set(mTextBox,'Units','characters')
mTextBoxPosition = get(mTextBox,'Position');
mTextBoxPosition(1) = 10;
mTextBoxPosition(2) = 4;
mTextBoxPosition(3) = 25;
mTextBoxPosition(4) = 1;
set(mTextBox,'Position',mTextBoxPosition);
mTextBox = uicontrol('style','text');
set(mTextBox,'String',...
    ['min= ' num2str(min(phiVec(1:tam,dir)))]);
set(mTextBox,'Units','characters')
mTextBoxPosition = get(mTextBox,'Position');
mTextBoxPosition(1) = 10;
mTextBoxPosition(2) = 3;
mTextBoxPosition(3) = 25;
mTextBoxPosition(4) = 1;
set(mTextBox,'Position',mTextBoxPosition);
%
set(gca,'Box','off')
subplot(1,2,2);hold on
phi = phiVec(tam+1:2*tam,dir)/maxphi/2;
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
mTextBox = uicontrol('style','text');
set(mTextBox,'String',...
    ['max= ' num2str(max(phiVec(tam+1:2*tam,dir)))]);
set(mTextBox,'Units','characters')
mTextBoxPosition = get(mTextBox,'Position');
mTextBoxPosition(1) = 55;
mTextBoxPosition(2) = 4;
mTextBoxPosition(3) = 25;
mTextBoxPosition(4) = 1;
set(mTextBox,'Position',mTextBoxPosition);
mTextBox = uicontrol('style','text');
set(mTextBox,'String',...
    ['min= ' num2str(min(phiVec(tam+1:2*tam,dir)))]);
set(mTextBox,'Units','characters')
mTextBoxPosition = get(mTextBox,'Position');
mTextBoxPosition(1) = 55;
mTextBoxPosition(2) = 3;
mTextBoxPosition(3) = 25;
mTextBoxPosition(4) = 1;
set(mTextBox,'Position',mTextBoxPosition);
%
set(gca,'Box','off')
cd out
cd phi
hgsave(hfig,['FigBoundary_frec' num2str(j) '_fza' num2str(dir) '.fig']);
cd ..
cd ..
close(hfig)

