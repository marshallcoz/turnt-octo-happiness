function plotGeom(j,Bou,res,p0)
% plot the problem geometry
hfig = figure('Name','Geometry');hold on
for i = 1:Bou.nBou
   plotCircle3D(Bou.pt{i}.center, ...
                Bou.pt{i}.normal, ...
                Bou.pt{i}.radio,2);
end
maxx = -1000;
minx = 1000;
maxy = -1000;
miny = 1000;
maxz = -1000;
minz = 1000;
for iPx = 1:res.nrecep %for each station
    [p_x] = pick_receptor(iPx,res);
    text(p_x.center(1),p_x.center(2),p_x.center(3),'{\color{blue}V}')
    maxx = max(maxx,p_x.center(1));
    maxy = max(maxy,p_x.center(2));
    maxz = max(maxz,p_x.center(3));
    minx = min(minx,p_x.center(1));
    miny = min(miny,p_x.center(2));
    minz = min(minz,p_x.center(3));
    
end
p_x = p0;
text(p_x.center(1),p_x.center(2),p_x.center(3),'{\color{red}X}')
    maxx = max(maxx,p0.center(1));
    maxy = max(maxy,p0.center(2));
    maxz = max(maxz,p0.center(3));
    minx = min(minx,p0.center(1));
    miny = min(miny,p0.center(2));
    minz = min(minz,p0.center(3));

view(45,36)
maxx = max(max(get(gca,'xlim')),maxx);
maxy = max(max(get(gca,'ylim')),maxy);
maxz = max(max(get(gca,'zlim')),maxz);
minx = min(min(get(gca,'xlim')),minx);
miny = min(min(get(gca,'ylim')),miny);
minz = min(min(get(gca,'zlim')),minz);
light('Position',[maxx*2 miny*2 maxz*5],'Style','infinite');
axis equal
xlim([minx maxx])
ylim([miny maxy])
zlim([minz maxz])
%
set(gca,'Box','off')
cd out
hgsave(hfig,['FigBoundary' num2str(j) '.fig']);
cd ..
close(hfig)