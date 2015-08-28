import numpy as np
from scipy.stats.mstats import mode
import time

x, y, z = np.loadtxt('1051_caustic.txt', usecols=(0,1,2), unpack=True)
r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
xf = x; yf = y; zf = z; rf = r

cosphif = zf/rf; sinphif = np.sqrt(1.0-cosphif**2)
costhf = xf/rf/sinphif; sinthf = yf/rf/sinphif

I = 6; J = 6; K = 2

N = I + J + K - 3

t1 = time.time()

M = np.zeros((I*J*K,len(xf)))

for n in range(len(xf)):
    m = 0 
    for k in range(K):
        for j in range(J):
            for i in range(I):
                M[m][n] = np.power(sinphif[n],i+j)*np.power(cosphif[n],k)*np.power(sinthf[n],j)*np.power(costhf[n],i)
                m += 1

print M

rN1 = np.power(rf,N+1)
c = np.dot(rf,np.linalg.pinv(M))
t2 = time.time()
print "matrix solution execution time =", t2-t1
print c

"""
cosphi = z/r; sinphi = np.sqrt(1.0-cosphi**2)
costh = x/r/sinphi; sinth = y/r/sinphi

m = 0 
rho = np.zeros_like(r)
for k in range(K):
    for j in range(J):
        for i in range(I):
            rho += c[m]*np.power(sinphi,i+j)*np.power(cosphi,k)*np.power(sinth,j)*np.power(costh,i)
            m += 1

chi2 = np.sum(np.square(r-rho))/(len(r)-1-I*J*K)
#chi2 = np.sum(np.abs(rf-rho))/(len(rf)-1-I*J*K)

print "chi2 =",chi2

#c = np.dot(rN1[np.newaxis],np.dot(M.T,np.linalg.inv(np.dot(M,M.T))))
print "c=", c
pi = np.pi
cos = np.cos
sin = np.sin

nphi = 90; ntheta = 180
phig = np.linspace(0.,pi,nphi)
thetag = np.linspace(0.,2.*pi,ntheta)

phi, theta = np.meshgrid(phig,thetag)
    
rho = np.zeros_like(phi)

sinphi = sin(phi); cosphi = cos(phi); sinth = sin(theta); costh = cos(theta)

m = 0 
for k in range(K):
    for j in range(J):
        for i in range(I):
            rho += c[m]*np.power(sinphi,i+j)*np.power(cosphi,k)*np.power(sinth,j)*np.power(costh,i)
            m += 1

xs = rho * sin(phi) * cos(theta)
ys = rho * sin(phi) * sin(theta)
zs = rho * cos(phi)
"""
