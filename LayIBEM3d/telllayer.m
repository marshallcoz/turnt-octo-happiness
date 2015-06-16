function [l,isonit] = telllayer(m,z,N)
l = 1;
isonit=false;
for i = 2:N+1
    if (z < m(i).z)
        break
    else
        l = i;
    end
end
if (abs(z - m(l).z) < 1/100)
    isonit = true;
end