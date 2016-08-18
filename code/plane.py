import math

L = 9.43E8 #UV bol lumgalex
alpha = 0.53 # kpc scale radius from noeske
R = 0.05
I = L/(4.0* math.pi * R * 1E3 * R * 1E3)
print 'R {}'.format(R)
print 'log I {}'.format(math.log10(I))
print 'log R {}'.format(math.log10(R))
print 'log L {}'.format(math.log10(L))
sigma = 10.0**((math.log10(I) + 1.21*math.log10(R)-0.55)/1.6)
print 'log sigma {}'.format(math.log10(sigma))
print 'sigma {}'.format(sigma)


obs = True

if obs:
    import numpy as np

    data = np.loadtxt('scale.dat')
    log_sigma = data[:,0]
    log_r = data[:,1]
    log_lum = data[:,2]
    ii = (log_sigma < np.log10(90.0)) & (log_sigma > np.log10(30.0)) 
    ii = ii & (log_lum> np.log10(0.33*L))
    print 10**log_sigma[ii]
    print 10**log_r[ii]
    print 10**log_lum[ii]
    #print log_sigma
