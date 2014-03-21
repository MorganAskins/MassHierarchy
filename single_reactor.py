# Reactor Antineutrino energy spectrum
import math
import pylab as plt
import scipy.stats as stats
import numpy as np
import time

# normalizations
core_235U = 0.58
core_239Pu = 0.3
core_238U = 0.07
core_241Pu = 0.05
distance = 13000                 # meters

# constants
mass_electron = 0.511
mass_neutron = 939.56
mass_proton = 938.27
delta = mass_neutron - mass_proton

event_rate = 25                 # per day

def main():
    k = 0.035
    xmin = 0.001
    xmax = 10
    bins = 1000
    bin_width = (xmax - xmin) / bins

    x = []
    y = []
    z = []
    q = []
    r = []
    t = []
    ys = [0]*bins
    zs = [0]*bins
    for i in range(bins):
        xvalue = bin_width * i + xmin
        x.append(xvalue)
        t.append( cross_section(xvalue) * reactor(xvalue) )
        y.append( cross_section(xvalue) * reactor(xvalue)  * survival(xvalue, 'normal') )
        z.append( cross_section(xvalue) * reactor(xvalue)  * survival(xvalue, 'inverted') )
        for j in range(bins):
            smear_value = j * bin_width + xmin
            ys[j] += y[i] * bin_width * gaus(smear_value, xvalue, k*xvalue**0.5)
            zs[j] += z[i] * bin_width * gaus(smear_value, xvalue, k*xvalue**0.5)        

    for i in range(bins):
        if z[i] > 0:
            q.append( ys[i] / t[i] )
            r.append( zs[i] / t[i] )
            #q.append( ys[i] / zs[i] )
        else:
            q.append( 1 )
            r.append( 1 )

    
    normal = open("normal_spectrum.dat",'w')
    inverted = open("inverted_spectrum.dat",'w')
    no_oscillations = open("no_oscillation_spectrum.dat",'w')

        
    plt.plot(x,q,x,r)
    plt.plot(x,ys,'b',x,zs,'r',x,t,'g')
    plt.show()    

def gaus(x,mu,sigma):
    return (2 * math.pi * sigma**2 )**(-0.5) * math.e **( -(x-mu)**2 / (2 * sigma**2) )

def cross_section( energy ): # MeV
    if energy > delta + mass_electron:
        E_electron = energy + mass_proton - mass_neutron
        P_electron = ( E_electron**2 - mass_electron**2 )**0.5
        f = 1.6857                  # Statistical / Coulomb factor
        T = 881.5                    # Neutron mean lifetime
        return 2 * math.pi**2 * P_electron * E_electron / (mass_electron**5 * f * T )
    return 0

def reactor( energy ):
    return ( core_235U * math.e**(0.870 - 0.160 * energy - 0.091*energy**2) +
             core_239Pu * math.e**(0.896 - 0.239 * energy - 0.981*energy**2) +
             core_238U * math.e**(0.976 - 0.162 * energy - 0.0790*energy**2) +
             core_241Pu * math.e**(0.793 - 0.080 * energy - 0.1085*energy**2) )

def survival( energy, hierarchy ):
    s212 = 0.857
    s223 = 0.95
    s213 = 0.098
    c13 = (1 - s213)**0.5
    c12 = (1 - s212)**0.5
    
    m21 = 7.5E-5

    # if hierarchy == 'normal':
    #     m32 = 2.32E-3
    #     m31 = m21 + m32
    # else:
    #     m32 = 2.32E-3
    #     m31 = m32 - m21

    if hierarchy == 'normal':
        m32 = 2.32E-3
        m31 = m21 + m32
    else:
        m31 = 2.32E-3
        m32 = m21 + m31

    return ( 1 - c13**4 * s212 * math.sin( 1.27 * m21 * distance / energy )**2
             - s213 * ( c12**2 * math.sin( 1.27 * m31 * distance / energy )**2 + 
                        s212 * math.sin( 1.27 * m32 * distance / energy )**2 ) )

if __name__ == "__main__":
    main()
