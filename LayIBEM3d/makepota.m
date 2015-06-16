function [potaOUT] = makepota(receptor,Bou)
% tabla que relaciona las profundidades de cada punto
% cada fila : la misma profundidad
% primera columna : número de sólo receptores
% segunda columna : número de 
% tercera columna : profundidad
% cuarta  columna : estrato
% y los apuntadores a receptor y Bou

pota = zeros(size(receptor,1)+Bou.nBou,4+size(receptor,1)+Bou.nBou);
% asumimos que el primer receptor está en E (region 1)
pota(1,1) = 1;
pota(1,2) = 0;
pota(1,3) = receptor{1,1}.center(3);
pota(1,4) = receptor{1,1}.layer;
pota(1,5) = 1;
nDeps = 1;

for i = 2:size(receptor,1)
    if (receptor{i}.region == 2) % receptor en R
        continue
    end
    nuevo = true;
    for ii = 1:nDeps
        if (abs(receptor{i,1}.center(3) - ...
             pota(ii,3)) <= 1/100)
           pota(ii,1) = pota(ii,1) + 1;
           pota(ii,4+pota(ii,1)) = i;
           nuevo = false;
           break;
        end
    end
    if (nuevo)
        nDeps = nDeps + 1;
        pota(nDeps,1) = 1;
        pota(nDeps,2) = 0;
        pota(nDeps,3) = receptor{i,1}.center(3);
        pota(nDeps,4) = receptor{i,1}.layer;
        pota(nDeps,5) = i;
    end
end
for i = 1:Bou.nBou
    if (Bou.pt{i,1}.tipoFrontera == 2)
        continue
    end
    nuevo = true;
    for ii = 1:nDeps
        if (abs(Bou.pt{i,1}.center(3) - ...
             pota(ii,3)) <= 1/100)
           pota(ii,2) = pota(ii,2) + 1;
           pota(ii,4+pota(ii,2)) = -i;
           nuevo = false;
           break;
        end
    end
    if (nuevo)
        nDeps = nDeps + 1;
        pota(nDeps,1) = 0;
        pota(nDeps,2) = 1;
        pota(nDeps,3) = Bou.pt{i,1}.center(3);
        pota(nDeps,4) = Bou.pt{i,1}.layer;
        pota(nDeps,5) = -i;
    end
end
%%
maxnum = pota(1,1) + pota(1,2);
for i=2:nDeps
    maxnum = max(maxnum,pota(i,1) + pota(i,2));
end
maxnum = maxnum +4;
potaOUT = pota(1:nDeps,1:maxnum);