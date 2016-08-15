from numpy import *
from scipy.interpolate import interp1d
import emcee
import os

#Constants
c = 299792.458

#Functions

def model_to_data_interp(x_m, y_m, x_d):

    f = interp1d(x_m, y_m, kind="linear", fill_value=0.0, bounds_error=False)    
    x_m_new = x_d.copy()
    y_m_new = f(x_m_new)

    return x_m_new, y_m_new


def normalization(x_d, y_d, x_m, y_m):

    obs_int = trapz(y_d, x_d)
    teo_int = trapz(y_m, x_m)
    y_m = y_m * obs_int/teo_int

    return y_m

def model(logtau, vmax, theta, logT, x_d, y_d):

    filename = "model_lt{}_vm{}_th{}.dat".format(logtau, vmax, theta)
    os.system( ('./analytic_solution.x {} {} {} > {}'.format(10.0**logtau, vmax, theta, filename)))

    try:
        mod = genfromtxt(filename)

    except ValueError:
        mod = genfromtxt(filename, skip_footer=3)

    x_m_i = mod[:,0]
    y_m_i = mod[:,1]

    v_th = 12.85*sqrt((10.0**logT)/(1.0E4))
    x_m_i = x_m_i*v_th 

    x_m, y_m = model_to_data_interp(x_m_i, y_m_i, x_d)

    y_m = normalization(x_d, y_d, x_m, y_m)

    os.system(('rm {}'.format(filename)))

    return x_m, y_m



#emcee functions

def lnprior(param):

    logtau, vmax, theta, logT = param

    if 6.0 < logtau < 9.0 and 200.0 < vmax < 600.0 and 0 < theta < 90 and 4 < logT < 4.5 :
        return 0.0

    return -inf


def lnlike(param, x_d, y_d, y_sigma):

    logtau, vmax, theta, logT = param

    x_m, y_m = model(logtau, vmax, theta, logT, x_d, y_d)

    chi_squared = 0.5*sum(((y_d-y_m)/y_sigma)**2)

    print('chi squared {}'.format(chi_squared))
    return -chi_squared


def lnprob(param, x_d, y_d, y_sigma):

    lp = lnprior(param)
    if not isfinite(lp):
        return -inf

    return lp + lnlike(param, x_d, y_d, y_sigma)



#Observed data

tol = loadtxt('../data/obs/t1214-cleanf0_normalized.txt')

x_data = tol[:,0]
y_data = tol[:,1]
y_data_sigma = tol[:,2]


#First guess

logtau_0 = 7.0
vmax_0 = 300
theta_0 = 45.0
logT_0 = 4.3

#x_0, y_0 = model(logtau_0, vmax_0, theta_0, logT_0, x_data, y_data)

first_guess = [logtau_0, vmax_0, theta_0, logT_0]



#Running emcee

ndim = 4
nwalkers = 24
nsteps = 500

pos0 = [first_guess+ array([0.2,20.0,5.0,0.01])*random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_data, y_data, y_data_sigma), threads=24)

print("Running MCMC...")
sampler.run_mcmc(pos0, nsteps)# rstate0=random.get_state())
print("Done!")

#f = open("chain.dat", "w")
#f.close()

#for result in sampler.sample(pos0, iterations=nsteps, storechain=False):
#    position = result[0]
#    f = open("chain.dat", "a")
#    print position
#    for k in range(position.shape[0]):
#        f.write("{} {}\n".format(k, " ".join(str(v) for v in position[k])))
#    f.close()


# Saving results

samples_fc = sampler.flatchain
savetxt('sampler_flatchain.dat', samples_fc, delimiter=',')

print("Mean acceptance fraction: {0:.3f}".format(mean(sampler.acceptance_fraction)))

#Links:
# https://github.com/dfm/emcee/blob/master/examples/line.py
# http://dan.iel.fm/emcee/current/user/line/
