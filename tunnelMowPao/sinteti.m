%% load
clear
cd '/Users/marshall/Documents/DOC/unit5/tunel_MowPau'
outA
figure; hold on
plot(real(s_tt_1),'r')
plot(imag(s_tt_1),'b')
tp = 0.3; 
ts = 1.5; 
t0 = -1.5;
%% sinteti
f.ntiempo = 1024;
f.DFREC = 0.05;
f.NFREC = 200;
f.dt = 1 / (f.ntiempo * f.DFREC);
[rick] = ricker(f,ts,tp,t0);

% rick
s = zeros(1,f_vars.ntiempo);
s(1,1:f_vars.NFREC) =  receptor{iPx}.greenG(1:3,dirFza,1:f_vars.NFREC);
s(1,f_vars.ntiempo-((2:f_vars.NFREC)-2)) = conj(s(1:3,2:f_vars.NFREC));

u = s_tt_1.*rick;

