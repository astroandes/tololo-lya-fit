from numpy import *
from matplotlib.pyplot import *
from scipy.interpolate import interp1d
import triangle
import os

#Info
nwalkers = 10
nsteps = 2000
nparameters = 5

z = 0
c = 299792.458

#Data
samples = loadtxt('./sampler_flatchain.dat', delimiter=',')

logtaus = samples[:,0]
vmaxs = samples[:,1]
thetas = samples[:,2]
logTs = samples[:,3]
voffs = samples[:,4]

#Histograms
labels=[r"$\log{\tau}$", "$v_{max}$", r"$\theta$", "$\log{T}$", "$\Delta v_{off}$"]

#logtau
fig = figure()

n,b,p = hist(logtaus, bins=100, histtype='step')
logtau_b = b[argmax(n)]
vlines(logtau_b, 0, 400)
xlabel(labels[0])
ylabel('$\mathrm{N}$')

savefig('./histograms_logtau.png')

#vmax
fig = figure()

n,b,p = hist(vmaxs, bins=100, histtype='step')
vmax_b = b[argmax(n)]
vlines(vmax_b, 0, 700)
xlabel(labels[1])
ylabel('$\mathrm{N}$')

savefig('./histograms_vmax.png')

#theta
fig = figure()

n,b,p = hist(thetas, bins=100, histtype='step')
theta_b = b[argmax(n)]
vlines(theta_b, 0, 1200)
xlabel(labels[2])
ylabel('$\mathrm{N}$')

savefig('./histograms_theta.png')

#logT
fig = figure()

n,b,p = hist(logTs, bins=100, histtype='step')
logT_b = b[argmax(n)]
vlines(logT_b, 0, 450)
xlabel(labels[3])
ylabel('$\mathrm{N}$')

savefig('./histograms_logT.png')

#voff
fig = figure()

n,b,p = hist(voffs, bins=100, histtype='step')
voff_b = b[argmax(n)]
vlines(voff_b, 0, 800)
xlabel(labels[4])
ylabel('$\mathrm{N}$')

savefig('./histograms_voff.png')

#Best parameters
os.system('echo "'+str(logtau_b)+','+str(vmax_b)+','+str(theta_b)+','+str(logT_b)+','+str(voff_b)+'" > best_parameters.dat')
