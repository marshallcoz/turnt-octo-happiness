function [Bou] = initBoundary
% leer archivo de frontera y guardar en variable
% 
%data = importdata('Bou.txt',' ',2);
% arreglo = cell(j);
% arreglo{1,1} = '0.21sphe.txt'; arreglo{1,2} = 0.21;
% arreglo{2,1} = '0.11sphe.txt'; arreglo{2,2} = 0.11;
% arreglo{3,1} = '0.05sphe.txt'; arreglo{3,2} = 0.05;
data = importdata('0.21sphe.txt',' ',2);
nt = data.textdata(2);
Bou.nBou = str2double(nt{1,1});
% 7 columnas, centro(3), normal(3), radio
Bou.pt = cell(Bou.nBou,1);
for i = 1:Bou.nBou
  Bou.pt{i}.center(1:3) = data.data(i,1:3);
  Bou.pt{i}.normal(1:3) = data.data(i,4:6);
  Bou.pt{i}.radio = data.data(i,7);
end