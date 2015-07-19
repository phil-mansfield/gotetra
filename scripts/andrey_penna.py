import numpy as np

x, y, z = np.loadtxt("/home/mansfield/code/go/src/github.com/phil-mansfield/gotetra/scripts/1092_caustic.txt", usecols=(0,1,2), unpack=True)

#x = x - np.mean(x); y = y - np.mean(y); z = z - np.mean(z)

r = np.sqrt(np.square(x) + np.square(y) + np.square(z))

rmax = r.max()
ind = (r<1.5)
x = x[ind]; y = y[ind]; z = z[ind]; r = r[ind]

I = 5; J = 5; K = 2

N = I + J + K - 3

M = np.zeros((I*J*K,len(x)))

for n in range(len(x)):
    m = 0
    for k in range(K):
        for j in range(J):
            for i in range(I):
                M[m][n] = (np.power(r[n],N-i-j-k) *
                           np.power(x[n],i) *
                           np.power(y[n],j) *
                           np.power(z[n],k))
                m += 1

rN1 = np.power(r,N+1)
c = np.dot(rN1,np.linalg.pinv(M))
#c = np.dot(rN1[np.newaxis],np.dot(M.T,np.linalg.inv(np.dot(M,M.T))))
print c
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
            rho += (c[m]*np.power(sinphi,i+j) *
                    np.power(cosphi,k) *
                    np.power(sinth,j) *
                    np.power(costh,i))
            m += 1

xs = rho * sin(phi) * cos(theta)
ys = rho * sin(phi) * sin(theta)
zs = rho * cos(phi)

import mayavi.mlab as mylab
fig=mylab.figure(1,bgcolor=(0,0,0))
mylab.clf()
ps = 1.*np.ones_like(x)
mylab.points3d(x,y,z,(ps), scale_factor=0.025,colormap='Spectral',opacity=1.)
mylab.mesh(xs, ys, zs,  colormap="bone",opacity=0.5)
mylab.outline();
mylab.xlabel('X'); mylab.ylabel('Y'); mylab.zlabel('Z')
