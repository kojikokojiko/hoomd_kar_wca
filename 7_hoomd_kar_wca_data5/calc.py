import math


n=0.8 
a_0=1.0


y=(math.pi/4.0)*n*a_0*a_0


Y=(1-7/16*y)/(1-y)**2


f=1+2*y*Y
dfdy=2*(1+y/8)/(1-y)**3


cs=math.sqrt(f+f**2+y*dfdy)
print(cs)
print(y)
print(Y)