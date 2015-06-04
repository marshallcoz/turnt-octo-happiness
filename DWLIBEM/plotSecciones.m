%% graficar Secciones del tunel
cd '/Users/marshall/Documents/DOC/coco/turnt-octo-happiness/DWLIBEM/Workdir/outs/traces1'
geom = load('Geom1.txt','-ascii');
nlines = size(geom,1);
dat = load('Secciones1.txt','-ascii');

nsecciones = 2;
npersec = 6;
dt = 0.0014796;
ntmax = 10;

coord = dat(:,1:2);
nptstime = (size(dat,2)-2)/2;
if (ntmax > nptstime); ntmax = nptstime; end
stt = dat(:,3:2+nptstime);
srt = dat(:,3+nptstime:end);
%%
j = 1:npersec;
for i=1:nsecciones
    for t = 1:ntmax
    hfig = figure('Name',['Seccion ' num2str(i) ' t' num2str(t*dt)]);
    set(gcf,'PaperPositionMode','auto')
    set(hfig, 'Position', [0 0 800 420])
    set(gcf,'Visible', 'off')
    subplot(2,2,[1 3])
    hold on
    for l=1:nlines
        line(geom(l,[1 3]),geom(l,[2 4]))
    end
    line(coord(j,1),coord(j,2),'Color',[1 0 0])
    str1 = num2str(1);
    text(coord(j(1),1),coord(j(1),2),str1)
    str1 = num2str(npersec);
    text(coord(j(npersec),1),coord(j(npersec),2),str1)
    
    set(gca,'YDir','reverse');
    axis square
    subplot(2,2,2)
    barh(stt(j,t));
    title('s_{\theta \theta}');
    subplot(2,2,4)
    barh(srt(j,t));
    title('s_{r \theta}');
    %saveas(hfig,['Secc' num2str(i) '_t' num2str(t) '.jpg']);
    str1 = ['Secc' num2str(i) '_t' num2str(t) '.eps'];
    print('-depsc2', str1);
    close('all','hidden')
    end
    j = j + npersec;
end
