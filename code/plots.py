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
tol = loadtxt('./tol.txt')
samples = loadtxt('./sampler_flatchain.dat', delimiter=',')

x_data = tol[:,0]
y_data = tol[:,1]

#Best parameters

f_best = open('./best_parameters.dat')
best = f_best.readline().split(',')
f_best.close()

logtau_b = float(best[0])
vmax_b = float(best[1])
theta_b = float(best[2])
logT_b = float(best[3])
voff_b = float(best[4].split(r'\n')[0])

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

    mod = genfromtxt( ('./model_lt%f_vm%f_th%f.dat')%(logtau, vmax, theta) )

    x_m_i = mod[:,0]
    y_m_i = mod[:,1]

    v_th = 12.85*sqrt((10.0**logT)/(10.0**4.0))
    x_m_i = x_m_i*v_th + voff

    x_m, y_m = model_to_data_interp(x_m_i, y_m_i, x_d)

    y_m = normalization(y_d, y_m)

    return x_m, y_m


#Plot fit and original
fig = figure(figsize=(12,6))

x_model, y_model = model(logtau_b, vmax_b, theta_b, logT_b, voff_b, x_data, y_data)

scatter(x_data, y_data, c='c', alpha=0.6, label='$\mathrm{TOL 1214-277}$')
plot(x_data, y_data, c='c', alpha=0.6)
plot(x_model, y_model, c='b', alpha=0.6, label='$\mathrm{Rotation\ model}$', lw=2)

xlabel('$\mathrm{Velocity (km/s)}$', fontsize=20)
ylabel('$\mathrm{Intensity (Arbitrary Units)}$', fontsize=20)
legend(loc='best', fontsize=20)

savefig('./mcmc.png')


# Walkers plot
clf()

fig, axes = subplots(5, 1, sharex=True, figsize=(12, 13))

axes[0].plot(samples[:,0].reshape(nsteps,nwalkers), color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(logtau_b, color="#888888", lw=2)
axes[0].set_ylabel(r"$\log{\tau}$")

axes[1].plot(samples[:,1].reshape(nsteps,nwalkers), color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(vmax_b, color="#888888", lw=2)
axes[1].set_ylabel("$v_{max}$")

axes[2].plot(samples[:,2].reshape(nsteps,nwalkers), color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(theta_b, color="#888888", lw=2)
axes[2].set_ylabel(r"$\theta$")

axes[3].plot(samples[:,3].reshape(nsteps,nwalkers), color="k", alpha=0.4)
axes[3].yaxis.set_major_locator(MaxNLocator(5))
axes[3].axhline(logT_b, color="#888888", lw=2)
axes[3].set_ylabel("$\log{T}$")

axes[4].plot(samples[:,4].reshape(nsteps,nwalkers), color="k", alpha=0.4)
axes[4].yaxis.set_major_locator(MaxNLocator(5))
axes[4].axhline(voff_b, color="#888888", lw=2)
axes[4].set_ylabel("$\Delta v_{off}$")
axes[4].set_xlabel("step number")

fig.tight_layout(h_pad=0.0)

savefig('./walkers.png')


# Triangle plot

fig = triangle.corner(samples, labels=[r"$\log{\tau}$", "$v_{max}$", r"$\theta$", "$\log{T}$", "$\Delta v_{off}$"],
                        truths=[logtau_b, vmax_b, theta_b, logT_b, voff_b])

savefig('./parameters.png')
