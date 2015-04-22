function [Bou] = initBoundary(res)
data = importdata(res.BouFile,' ',2);
nt = data.textdata(2);
Bou.nBou = str2double(nt{1,1});
% 7 columnas, centro(3), normal(3), radio
Bou.pt = cell(Bou.nBou,1);
iren = 1;
for i = 1:Bou.nBou
  %Bou.pt{i}.center(1:3) = data.data(iren,1:3);
  Bou.pt{i}.normal(1:3) = data.data(iren,4:6);
  Bou.pt{i}.radio = data.data(iren,7);
  iren = iren+1;
  nvert = data.data(iren,1);
  Bou.pt{i}.vert = zeros(nvert,3);
  for vert = 1:nvert
      iren = iren+1;
      Bou.pt{i}.vert(vert,1:3) = data.data(iren,1:3);
  end
  Bou.pt{i}.center(1) = mean(Bou.pt{i}.vert(:,1));
  Bou.pt{i}.center(2) = mean(Bou.pt{i}.vert(:,2));
  Bou.pt{i}.center(3) = mean(Bou.pt{i}.vert(:,3));
  iren = iren+1;
end

% % now the green functions from each segment to each station
% for i = 1:Bou.nBou
%     Bou.pt{i}.gG = zeros(res.nrecep,3,3); %just valid each frequency
%     Bou.pt{i}.gT = zeros(res.nrecep,3,3); 
% end
