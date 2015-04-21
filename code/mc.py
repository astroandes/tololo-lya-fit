from numpy import *
from scipy.interpolate import interp1d
import emcee
import os

#Constants

z = 0.029836
c = 299792.458

#Functions

def model_to_data_interp(x_m, y_m, x_d):

    f = interp1d(x_m, y_m, kind="linear")

    try:
        x_m_new = x_d
        y_m_new = f(x_d)

        return x_m_new, y_m_new

    except ValueError:

        n = len(x_d)

        x_m_new = x_d
        y_m_new = zeros(n)

        for i in range(n):

            x_d_i = x_d[i]

            if(x_d_i >= x_m[0] and x_d_i <= x_m[-1]):
                y_m_new[i] = f(x_d_i)
            else:
                y_m_new[i] = 0.0

        return x_m_new, y_m_new


def normalization(y_d, y_m):

    obs_int = sum(y_d)
    teo_int = sum(y_m)
    y_m = y_m*obs_int/teo_int

    return y_m

def model(logtau, vmax, theta, logT, voff, x_d, y_d):

    os.system( ('./analytic_solution.x %f %f %f > "./model_lt%f_vm%f_th%f.dat"')%(10.0**logtau, vmax, theta, logtau, vmax, theta) )

    try:
	       mod = genfromtxt( ('./model_lt%f_vm%f_th%f.dat')%(logtau, vmax, theta) )

    except ValueError:
	       mod = genfromtxt( ('./model_lt%f_vm%f_th%f.dat')%(logtau, vmax, theta), skip_footer=3)

    x_m_i = mod[:,0]
    y_m_i = mod[:,1]

    v_th = 12.85*sqrt((10.0**logT)/(10.0**4.0))
    x_m_i = x_m_i*v_th + voff

    x_m, y_m = model_to_data_interp(x_m_i, y_m_i, x_d)

    y_m = normalization(y_d, y_m)

    os.system( ('rm model_lt%f_vm%f_th%f.dat')%(logtau, vmax, theta) )

    return x_m, y_m


#Observed data

lars = loadtxt('./LARS02.txt')

x_data = lars[:-17,0] - c*z
y_data = lars[:-17,1]

#emcee functions

def lnprior(param):

    logtau, vmax, theta, logT, voff = param

    if 4 < logtau < 6 and 0 < vmax < 100 and 0 < theta < 90 and 4 < logT < 4.5 and -20 < voff < 20:
        return 0.0

    return -inf


def lnlike(param, x_d, y_d):

    logtau, vmax, theta, logT, voff = param

    x_m, y_m = model(logtau, vmax, theta, logT, voff, x_d, y_d)

    chi_squared = (1.0/2.0)*sum((y_d-y_m)**2)

    return -chi_squared


def lnprob(param, x_d, y_d):

    lp = lnprior(param)
    if not isfinite(lp):
        return -inf

    return lp + lnlike(param, x_d, y_d)


#First guess

logtau_0 = 5
vmax_0 = 50
theta_0 = 45
logT_0 = 4
voff_0 = -11

x_0, y_0 = model(logtau_0, vmax_0, theta_0, logT_0, voff_0, x_data, y_data)

first_guess = [logtau_0, vmax_0, theta_0, logT_0, voff_0]


#Running emcee

ndim = 5
nwalkers = 10
nsteps = 2000

pos = [first_guess+ 1e-3*random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x_data, y_data), threads=8)

print("Running MCMC...")
sampler.run_mcmc(pos, nsteps, rstate0=random.get_state())
print("Done!")


# Saving results

samples_fc = sampler.flatchain

savetxt('sampler_flatchain.dat', samples_fc, delimiter=',')

print("Mean acceptance fraction: {0:.3f}".format(mean(sampler.acceptance_fraction)))

#Links:
# https://github.com/dfm/emcee/blob/master/examples/line.py
# http://dan.iel.fm/emcee/current/user/line/
