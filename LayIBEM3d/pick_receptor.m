function [p_x] = pick_receptor(i,res,f_vars)
p_x.center(1:3) = res.receptor{i}.center;
p_x.p(1:3) = res.receptor{i}.p;
p_x.normal(1:3) =res.receptor{i}.normal;

p_x.greenG = zeros(3,3,f_vars.NFREC); 
p_x.greenT = zeros(3,3,f_vars.NFREC); 
p_x.sism = zeros(3,3,f_vars.ntiempo);
