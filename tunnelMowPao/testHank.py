import scipy as sp
from scipy import special
z = 4+5j
s = "z = " + str(z)
print s
# sp.info(special.hankel1)
k = special.jn(-1,z)
s = " J_-1 = " + str(k)
print s
k = special.yv(-1,z)
s = " Y_-1 = "+ str(k)
print s
k = special.hankel1(-1,z)
s = " H^1_-1(z) = "+ str(k)
print s
k = special.hankel2(-1,z)
s = " H^2_-1(z) = "+ str(k)
print s
print " "
k = special.jn(0,z)
s = " J_0 = " + str(k)
print s
k = special.yv(0,z)
s = " Y_0 = " + str(k)
print s
k = special.hankel1(0,z)
s = " H^1_0(z) = " + str(k)
print s
k = special.hankel2(0,z)
s = " H^2_0(z) = "+ str(k)
print s
print " "
k = special.jn(1,z)
s = " J_1 = " + str(k)
print s
k = special.yv(1,z)
s = " Y_1 = " + str(k)
print s
k = special.hankel1(1,z)
s = " H^1_1(z) = " + str(k)
print s
k = special.hankel2(1,z)
s = " H^2_1(z) = "+ str(k)
print s
print " "
k = special.jn(2,z)
s = " J_2 = " + str(k)
print s
k = special.yv(2,z)
s = " Y_2 = " + str(k)
print s
k = special.hankel1(2,z)
s = " H^1_2(z) = " + str(k)
print s
k = special.hankel2(2,z)
s = " H^2_2(z) = "+ str(k)
print s
print " "
k = special.jn(3,z)
s = " J_3 = " + str(k)
print s
k = special.yv(3,z)
s = " Y_3 = " + str(k)
print s
k = special.hankel1(3,z)
s = " H^1_3(z) = " + str(k)
print s
k = special.hankel2(3,z)
s = " H^2_3(z) = "+ str(k)
print s
print " "
k = special.jn(49,z)
s = " J_49 = " + str(k)
print s
k = special.yv(49,z)
s = " Y_49 = " + str(k)
print s
k = special.hankel1(49,z)
s = " H^1_49(z) = " + str(k)
print s
k = special.hankel2(49,z)
s = " H^2_49(z) = "+ str(k)
print s
print " "
k = special.jn(50,z)
s = " J_50 = " + str(k)
print s
k = special.yv(50,z)
s = " Y_50 = " + str(k)
print s
k = special.hankel1(50,z)
s = " H^1_50(z) = " + str(k)
print s
k = special.hankel2(50,z)
s = " H^2_50(z) = "+ str(k)
print s
print " "

