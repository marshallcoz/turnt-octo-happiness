function [eXi,atInterf] = el_estrato_es(m,N,zXi)
eXi = 1;
atInterf = false;
for e=2:N+1
    if (zXi >= m(e).z)
        eXi = e;
        if (zXi == m(e).z),atInterf = true;end
        return
    end
end
