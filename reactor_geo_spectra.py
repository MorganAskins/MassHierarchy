# Steve Dye gives spectrum in events / 10 keV / 10^32 protons / year
# Read in files and convert to events / 10 keV / kton / year

import numpy as np
import math
import pylab as plt

def main():
    
    density = 0.8               # g / cm^3
    # Mineral oil is ~ CH2
    atomic_weight = (12.0107 + 2 *1.008) / 6.022E23 # g / molecule
    # A molecule as 8 protons
    volume = (1000)*(100**3)    # cm^3
    protons = density * volume / atomic_weight * 8
    conversion = protons / 1E32 / 365.25

    print("Target Protons of WatchMan: ", protons)
    print("Conversion: ", conversion)

    uranium = np.loadtxt('imb_uran.dat')
    thorium = np.loadtxt('imb_thor.dat')
    reactor = np.loadtxt('imb_reactor_bkg.dat')

    mod_uranium = open('imb_uran_norm.dat','w')
    mod_thorium = open('imb_thor_norm.dat','w')
    mod_reactor = open('imb_reactor_norm.dat','w')

    for i in range(len(uranium)):
        mod_uranium.write( str( uranium[i] * conversion ) + '\n' )
    for i in range(len(thorium)):
        mod_thorium.write( str( thorium[i] * conversion ) + '\n' )
    for i in range(len(reactor)):
        mod_reactor.write( str( reactor[i] * conversion ) + '\n' )

if __name__ == "__main__":
    main()

