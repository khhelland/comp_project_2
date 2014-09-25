import numpy as np
import matplotlib.pyplot as plt
import os, time

#make programs
os.system('make')
os.system('make -f make_2')


#lufile = open('lu_err.dat','w')


#lufile.write('n & $\epsilon_i$ \\\ \hline \n')

n = 300
rhomax = 10
omegas = (0.01, 0.5, 1, 5)
for omega in omegas:
    
    run = ' '.join(('./2_c.x', str(omega)))
    os.system(run)
    #Fetch simulation
    wavefunction = np.loadtxt('ground.dat')
    probability = np.zeros(n+2)
    probability[1:-1] = wavefunction**2
    rho = np.linspace(0,rhomax, n+2)
    
    
    

    #Plot data against analytical solution

    
    plt.plot(rho,probability,'-')
    plt.hold(1)
                   

#plt.legend((str(omegas))
plt.xlabel(r'$\rho$')
plt.ylabel(r'$|u(\rho)|^2$')
plt.show()


    #lufile.write(' '.join((str(steps),'&',str(eps),'\\\ \hline \n')))

