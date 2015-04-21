function [res] = initreceptores(res,f_vars)
hfig = figure('Name','ReceiversPlot');hold on
       
xran = res.box(1,1):res.Del:res.box(1,2);
yran = res.box(2,1):res.Del:res.box(2,2);
zran = res.box(3,1):res.Del:res.box(3,2);

res.nx = size(xran,2);
res.ny = size(yran,2);
res.nz = size(zran,2);

res.nrecep = res.nx * res.ny * res.nz;
disp(['Num de receptores: ' num2str(res.nrecep)])
res.receptor = cell(res.nrecep,1);
cont = 1;

polZ0c = importdata(res.BouFileZ0,' ',2);
xv = polZ0c.data(:,1);
yv = polZ0c.data(:,2);

for i = 1:res.nx;
    for j = 1:res.ny;
       for k = 1:res.nz;
           res.receptor{cont}.center = [xran(i) yran(j) zran(k)];
           res.receptor{cont}.p = [i,j,k]; %indice de la coordenada
           res.receptor{cont}.normal(1:3) = res.norm;
           res.receptor{cont}.region = 1;%'E';
           inout = inpolygon(xran(i),yran(j),xv,yv);
           if (inout == 1)
               res.receptor{cont}.region = 2;
               auxtx = ['{\color{blue}' num2str(cont) '}'];
           else
               auxtx = ['{\color{red}' num2str(cont) '}'];
           end%'R';end
               text(xran(i),yran(j),zran(k),auxtx)
           
           res.receptor{cont}.greenG = zeros(3,3,f_vars.NFREC); 
           res.receptor{cont}.greenT = zeros(3,3,f_vars.NFREC); 
           res.receptor{cont}.sism = zeros(3,3,f_vars.ntiempo);
           
           cont = cont + 1;
       end
    end
end

res.fotogramas = zeros(3,f_vars.ntiempo,res.nx,res.ny,res.nz);

view(45,36)
axis equal
xlim([min(xran) max(xran)])
ylim([min(yran) max(yran)])
maxz = max(zran);
minz = min(zran);
maxz = max(max(get(gca,'zlim')),maxz);
minz = min(min(get(gca,'zlim')),minz);
zlim([minz maxz])
set(gca,'Box','off')
cd out
hgsave(hfig,'ReceiversPlot.fig');
cd ..
close(hfig)