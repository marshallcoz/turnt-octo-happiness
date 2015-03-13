%% test
cd '/Users/marshall/Documents/DOC/unit6/LayIBEM3d'
p_x.center(1:3) =[1 1 0];
pXi.center(1:3) =[0 0 0];
r = distancia(p_x,pXi);
gamma = (p_x.center -pXi.center)./r;
%%
J = 2;
DFREC = 0.5;
omei = - 1*PI/TW
%%
figure;hold;
quiver3(0,0,0,0,0,1,0.5)
quiver3(1,1,1,0,0,1,0.5)
%%
%figure
sis = s(1,:);
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,real(sis),'r');hold;
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,imag(sis),'b');
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,abs(sis),'k')
%%
figure
plot((0:f_vars.ntiempo-1).*f_vars.DFREC,s,'k')

