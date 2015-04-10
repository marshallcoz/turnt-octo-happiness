function [res] = pack_receptor(J,i,res,G,T)
res.receptor{i,1}.greenG(:,:,J) = G;
res.receptor{i,1}.greenT(:,:,J) = T;