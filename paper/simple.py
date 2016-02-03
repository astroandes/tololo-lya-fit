tau = 1E7
v_rot = 300 # km/s
sigma_0 = 5.898E-14 #cm-2
n_values = [1,1E-3]
kpc_in_cm = 3.08E21

for n in n_values:
    D = 2.0 * tau / (sigma_0 * n) # in cm
    M_D = 2.16E5 * (v_rot)**2 * (D/kpc_in_cm)
    print("n= {} cm^-3 D={} kpc".format(n,D/kpc_in_cm, M_D/1E10)) 


D_kpc = [0.11, 0.34] # rango in kpc
for D in D_kpc:
    M_D = 2.16E5 * (v_rot)**2 * (D)
    print("M = {} 1E9 Msun".format(M_D/1E9)) 
