function [xli,yli,zli] = plotGeom(j,Bou,res,receptor,p0,plotReceiv)
% plot the problem geometry
hfig = figure('Name','Geometry');hold on
set(gca,'ZDir','reverse');
for i = 1:Bou.nBou
   plotPatch(Bou.pt{i}.vert,0.9);
   plotCircle3D(Bou.pt{i}.center, ...
                Bou.pt{i}.normal, ...
                Bou.pt{i}.radio,1);
end
maxx = -1000;
minx = 1000;
maxy = -1000;
miny = 1000;
maxz = -1000;
minz = 1000;
if (plotReceiv)
% ****************** RECEPTORES **************************
for iPx = 1:res.nrecep %for each station
    [p_x] = receptor{iPx};%pick_receptor(iPx,res);
    text(p_x.center(1),p_x.center(2),p_x.center(3),'{\color{blue}V}')
    maxx = max(maxx,p_x.center(1));
    maxy = max(maxy,p_x.center(2));
    maxz = max(maxz,p_x.center(3));
    minx = min(minx,p_x.center(1));
    miny = min(miny,p_x.center(2));
    minz = min(minz,p_x.center(3));
end
% ***************** FUENTE PUNTUAL **********************
p_x = p0;
%text(p_x.center(1),p_x.center(2),p_x.center(3),'{\color{red}X}')
plot3(p_x.center(1),p_x.center(2),p_x.center(3),'r.','MarkerSize',25)
    maxx = max(maxx,p0.center(1));
    maxy = max(maxy,p0.center(2));
    maxz = max(maxz,p0.center(3));
    minx = min(minx,p0.center(1));
    miny = min(miny,p0.center(2));
    minz = min(minz,p0.center(3));
end
view(45,36)
maxx = max(max(get(gca,'xlim')),maxx);
maxy = max(max(get(gca,'ylim')),maxy);
maxz = max(max(get(gca,'zlim')),maxz);
minx = min(min(get(gca,'xlim')),minx);
miny = min(min(get(gca,'ylim')),miny);
minz = min(min(get(gca,'zlim')),minz);
light('Position',[maxx*2 miny*2 maxz*5],'Style','infinite');
axis equal
xli = [minx maxx];
yli = [miny maxy];
zli = [minz maxz];
% ****************** SUPERFICIE LIBRE *************
hplane = fill3([minx minx maxx maxx],...
               [maxy miny miny maxy],...
               [0 0 0 0],'g');
set(hplane,'edgecolor','none')
alpha(hplane,0.3);
% lineas de centro
hline = line([minx maxx],[0 0]);
set(hline,'linewidth',0.1);
set(hline,'LineStyle','-');
set(hline,'Color','black');
hline = line([0 0],[miny maxy]);
set(hline,'linewidth',0.1);
set(hline,'LineStyle','-');
set(hline,'Color','black');
%
hline = line([minx maxx],[0 0],[maxz maxz]);
set(hline,'linewidth',0.1);
set(hline,'LineStyle','-');
set(hline,'Color','black');
hline = line([0 0],[miny maxy],[maxz maxz]);
set(hline,'linewidth',0.1);
set(hline,'LineStyle','-');
set(hline,'Color','black');

xlim(xli)
ylim(yli)
zlim(zli)
xlabel('x');ylabel('y');zlabel('z')
%
set(gca,'Box','off')
cd out
set(gcf,'Visible', 'on');
if (plotReceiv)
nombrefigura = ['FigReceivers' num2str(j)];
else
nombrefigura = ['FigBoundary' num2str(j)];
end
saveas(hfig,[nombrefigura '.jpg']);
%print -depsc2 profile.eps;
print('-depsc2',[nombrefigura '.eps'])
hgsave(hfig,[nombrefigura '.fig']);
cd ..
close(hfig)