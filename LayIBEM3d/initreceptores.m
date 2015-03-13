function [res] = initreceptores(res)
       
xran = res.box(1,1):res.Del:res.box(1,2);
yran = res.box(2,1):res.Del:res.box(2,2);
zran = res.box(3,1):res.Del:res.box(3,2);

res.nx = size(xran,2);
res.ny = size(yran,2);
res.nz = size(zran,2);

%[ops.X,ops.Y,ops.Z] = meshgrid(xran,yran,zran);
res.nrecep = res.nx * res.ny * res.nz;
res.receptor = cell(res.nrecep,1);
cont = 1;
for i = 1:res.nx;
    for j = 1:res.ny;
       for k = 1:res.nz;
           res.receptor{cont}.center = [xran(i) yran(j) zran(k)];
           res.receptor{cont}.p = [i,j,k];
           res.receptor{cont}.normal(1:3) = res.norm;
           cont = cont + 1;
       end
    end
end
