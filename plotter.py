import numpy as np
import matplotlib.pyplot as plt
import os, time

#make programs
#os.system('make')
os.system('make -f make_2')


#lufile = open('lu_err.dat','w')


#lufile.write('n & $\epsilon_i$ \\\ \hline \n')

n = 300
rhomax = (40,10,8,5)
omegas = (0.01, 0.5, 1, 5)

for i in range(len(omegas)):
    
    run = ' '.join(('./2_c.x', str(omegas[i]), str(rhomax[i])))
    os.system(run)
    #Fetch simulation
    wavefunction = np.loadtxt('ground.dat')
    probability = np.zeros(n+2)
    probability[1:-1] = wavefunction**2
    rho = np.linspace(0,rhomax[i], n+2)
    
    
    

    #Plot data against analytical solution

    
    plt.plot(rho,probability,'-')
    plt.hold(1)
                   

plt.legend(map(' '.join,zip((r'$\omega_r$ = '+ s for s in map(str,omegas)),
                            (r', $\rho_{max}$ = '+ s for s in map(str,rhomax)))))
plt.xlabel(r'$\rho$')
plt.ylabel(r'$|u(\rho)|^2$')
plt.savefig('probability_densities.png')


    #lufile.write(' '.join((str(steps),'&',str(eps),'\\\ \hline \n')))

