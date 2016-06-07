import numpy as np
#http://www.brera.inaf.it/utenti/marcella/masse.html
#estimate stellar mass from non detection in 2mass
m_l = 1.0 # mass to light ratio
d_pc = 100E6 # distance in parsec
Mag_sun = 3.41 # Sun magnitude in this filter
mag_gal  = 15.0
log_m_stellar = np.log10(m_l) + 2.0 * np.log10(d_pc) - 2.0 + 0.4 * Mag_sun - 0.4 * mag_gal
print log_m_stellar
