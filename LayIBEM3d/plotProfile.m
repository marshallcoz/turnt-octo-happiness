function plotProfile(m,N)
col = zeros(N+1,1);
maxcol = 0.4; % el más obscuro
mincol = 0.6; % el más claro
for i=1:N+1
    col(i) = m(i).beta0;
end
maxVel = max(col); minVel = min(col);
if (maxVel == minVel)
mvc = 0;    
else
mvc = (maxcol-mincol)/(maxVel-minVel);
end
for i=1:N+1
    col(i) = mincol + col(i)*mvc;
end

hfig = figure('Name','Layered Media');
set(gcf,'Visible', 'off');
set(gcf,'PaperPositionMode','auto')
subplot(1,2,1)
set(gca,'YDir','reverse');
hold on
for i=1:N
    rectangle('Position',[0,m(i).z,1,(m(i+1).z-m(i).z)],'FaceColor',col(i).*[1 1 1])
end
i=N+1;
rectangle('Position',[0,m(i).z,1,(1.1*m(i).z-m(i).z)],'FaceColor',col(i).*[1 1 1])
ylim([m(1).z 1.1*m(N+1).z])
set(gca,'XTick',[])
ylabel('Depth [meters]')
subplot(1,2,2)
set(gca,'YDir','reverse');
hold on
for i=1:N
    line([m(i).beta0 m(i).beta0], [m(i).z m(i+1).z],'Color','r')
    line([m(i).beta0 m(i+1).beta0], [m(i+1).z m(i+1).z],'Color','r')
    line([m(i).alfa0 m(i).alfa0], [m(i).z m(i+1).z],'Color','b')
    line([m(i).alfa0 m(i+1).alfa0], [m(i+1).z m(i+1).z],'Color','b')
end
line([m(i+1).beta0 m(i+1).beta0], [m(i+1).z 1.1*m(i+1).z],'Color','r')
line([m(i+1).alfa0 m(i+1).alfa0], [m(i+1).z 1.1*m(i+1).z],'Color','b')
xlabel('Wave velocity [m/s]')
set(gca,'YTick',[])
ylim([m(1).z 1.1*m(N+1).z])

p = plot(0,0,'r',0,0,'b');
legend(p,'beta', 'alpha');
chdir out
outfigsize = get(hfig, 'Position');
outfigsize(3) = outfigsize(3)/2;
set(hfig, 'Position', outfigsize);
set(gcf,'Visible', 'on');
saveas(hfig,'profile.jpg');
print -depsc2 profile.eps;
hgsave(hfig,'profile.fig');
close(hfig);
chdir ..