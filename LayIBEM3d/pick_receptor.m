function [p_x] = pick_receptor(i,res)
p_x.center(1:3) = res.receptor{i}.center;
p_x.p(1:3) = res.receptor{i}.p;
p_x.normal(1:3) =res.receptor{i}.normal;
