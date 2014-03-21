# Reactor Antineutrino energy spectrum
# Build up 10000 inverted / normal spectrums with smeared oscillation parameters
# Binning is events / 10 keV / day
# Operations are as follows:
# Read in the Li9, geonu, reactor, neutron backgrounds (All normalized)
# Output will be of the format:
# (Inverted/Normal ID) [ Histogrammed Data Array ] // One line for each experiment
# First two lines output the ideal (Inf Stat) experiments
# Following this are 20,000 experiments (same length) finite statistics

import math
import pylab as lab
import scipy.stats as stats
import numpy as np
import time
import random
import multiprocessing as mp
import sys

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

event_rate = 25                 # total per day
k = 0.035                   # Detector Energy Resolution
bins = 1000
xmin = 0.000001
xmax = 10                   # MeV
bin_width = (xmax - xmin) / bins

m_atm = 2.32E-3
m_atm_lower = 0.08E-3
m_atm_upper = 0.12E-3


ktony = int( sys.argv[1] )
run_length = 365 * ktony

def main():

    runs = 10
    m23_errors = random_asym_normal(m_atm, m_atm_lower, m_atm_upper)
    backgrounds = build_backgrounds()

    normal_file = open('normal_events_' + sys.argv[1] + 'ky.txt', 'w')
    inverted_file = open('inverted_events_' + sys.argv[1] + 'ky.txt', 'w')

    inf_normal = combine( backgrounds, build_spectrum( 'normal', 'no', m23_errors ) )
    inf_inverted = combine( backgrounds, build_spectrum( 'inverted', 'no', m23_errors ) )
    
    write_to_file( normal_file, inf_normal )
    write_to_file( normal_file, inf_inverted )
    write_to_file( inverted_file, inf_normal )
    write_to_file( inverted_file, inf_inverted )

    processes = 4

    norm_stat = []

    for run in range(runs):
        pool = mp.Pool( processes )
        normal_neutrino_pool = pool.starmap( full_build, ( ('normal', m23_errors, backgrounds) for i in range(processes) ) )
        inverted_neutrino_pool = pool.starmap( full_build, ( ('inverted', m23_errors, backgrounds) for i in range(processes) ) )
        pool.close()
        pool.join()
        for i in range( processes ):
            write_to_file( normal_file, normal_neutrino_pool[i] )
            #norm_stat.append( normal_neutrino_pool[i] )
            write_to_file( inverted_file, inverted_neutrino_pool[i] )
        if run % 10 == 0:
            print( "Event: ", run )        
        
    # lab.ion()
    # for i in range( len(norm_stat) ):
    #     lab.plot( norm_stat[i] )
    #     lab.draw()
    #     lab.pause(0.5)
    

##########################
## FUNCTION DEFINITIONS ##
##########################

def full_build( hierarchy, m23_errors, backgrounds):
    neutrinos = build_spectrum( hierarchy, 'yes', m23_errors )
    data_spectrum = combine( backgrounds, neutrinos )
    total_events = sum( data_spectrum ) * run_length    
    #real_spectrum = sample_spectrum( data_spectrum, int(total_events) )
    real_spectrum = sample_spectrum_simple( data_spectrum, int(total_events))
    return real_spectrum

def sample_spectrum_simple( spectrum, samples ):
    binned_spectrum = [0]*bins
    # renormalize
    # for i in range(len(spectrum)):
    #     spectrum[i] *= run_length
    
    for i in range(len( spectrum )):
        lamb = spectrum[i]
        binned_spectrum[i] = int(random.normalvariate( lamb, lamb**0.5 ))
        if binned_spectrum[i] < 0:
            binned_spectrum[i] = 0
    
    return binned_spectrum
    

def sample_spectrum( spectrum, samples ):
    binned_spectrum = [0]*bins
    top = 0
    for i in range( len(spectrum) ):
        if spectrum[i] > top:
            top = spectrum[i]

    for i in range(samples):
        value = random_from_spectrum( spectrum, bins, top )
        binned_spectrum[value] += 1
        
    return binned_spectrum

def random_from_spectrum( spectrum, bin_max, top ):
    random_x = int( random.uniform(0, bin_max) )
    random_y = random.uniform(0, top)
    if spectrum[random_x] > random_y:
        return random_x
    else:
        return random_from_spectrum( spectrum, bin_max, top )

def write_to_file( out_file, input_list ):
    for i in range( len(input_list) ):
        out_file.write( str( input_list[i] ) + ' ' )
    out_file.write('\n')

def build_backgrounds():
    # Fill with backgrounds
    # All energies should be the prompt signal and not the neutrino energy
    # For now I will assume no nuclear recoil, and thus E_nu = E_e + m_e + m_neutron - m_proton
    Li9_bkg = np.loadtxt(open('spectra/Li9_Spectrum.dat')) # positron energy !!
    reactor_bkg = np.loadtxt(open('spectra/imb_reactor_norm.dat'))
    thorium_bkg = np.loadtxt(open('spectra/imb_thor_norm.dat'))
    uranium_bkg = np.loadtxt(open('spectra/imb_uran_norm.dat'))
    ## MARC: ## neutron_bkg = np.loadtxt(open('spectra/neutrons.dat'))
    # Shift the spectra from neutrino to positron and add to the neutrino events
    shifted_reactor = []
    shifted_thorium = []
    shifted_uranium = []

    for i in range(bins):
        this_bin = int( i + (delta + mass_electron) * 100 )
        if this_bin < bins:
            shifted_reactor.append( reactor_bkg[this_bin] )
            shifted_thorium.append( thorium_bkg[this_bin] )
            shifted_uranium.append( uranium_bkg[this_bin] )
        else:
            shifted_reactor.append( 0 )
            shifted_thorium.append( 0 )
            shifted_uranium.append( 0 )
    ## Combine into one spectrum
    backgrounds = []
    for i in range(bins):
        backgrounds.append( Li9_bkg[i] + shifted_reactor[i] + shifted_thorium[i] + shifted_uranium[i] ) # + neutron_bkg[i] )

    return backgrounds
            
def combine( backgrounds, neutrinos ):
    spectra = [backgrounds, neutrinos]
    total_spectra = len(spectra)
    combination = [0]*bins
    for i in range(bins):
        xvalue = bin_width * i + xmin
        sigma = k*xvalue**0.5
        amplitude = 0
        for j in range(total_spectra):
            amplitude += spectra[j][i]*run_length
        if (xvalue - 3*sigma - xmin ) / bin_width > 0:
            bin_min = int( (xvalue - 3*sigma - xmin ) / bin_width )
        else:
            bin_min = 0
        if (xvalue + 3*sigma - xmin ) / bin_width < bins:
            bin_max = int( (xvalue + 3*sigma - xmin ) / bin_width )
        else:
            bin_max = 0
        for q in range(bin_min, bin_max):
            smear_value = q * bin_width + xmin
            combination[q] += amplitude * bin_width * gaus(smear_value, xvalue, sigma )

    return combination


def build_spectrum( hierarchy, wiggle, m23_errors ):
    spectrum = []
    inf_inverted = []

    # Fill with neutrino events
    for i in range(bins):
        xvalue = bin_width * i + xmin
        spectrum.append( cross_section(xvalue) * reactor(xvalue)  * survival(xvalue, hierarchy, wiggle, m23_errors) )

    # Normalize
    norm = sum(spectrum)
    for i in range(bins):
        spectrum[i] *= event_rate / norm

    shifted_spectrum = []
    # Shift Spectrum
    for i in range(bins):
        this_bin = int( i + (delta + mass_electron) * 100 )
        if this_bin < bins:
            shifted_spectrum.append( spectrum[this_bin] )
        else:
            shifted_spectrum.append(0)

    return shifted_spectrum

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

def survival( energy, hierarchy, wiggle, m23_errors ):

    theta_12 = ( 33.57 ) * ( math.pi / 180 )
    theta_12_error = ( 0.76 ) * ( math.pi / 180 )

    theta_13 = (8.71) * (math.pi / 180 )
    theta_13_error = ( 0.38 ) * (math.pi / 180)

    m_sol = 7.45E-5
    m_sol_error = 0.18E-5

    if wiggle == 'yes':
        # Smear the angles with a gaussian
        theta_12 = random.normalvariate( theta_12, theta_12_error )
        theta_13 = random.normalvariate( theta_13, theta_13_error )
        m_sol = random.normalvariate( m_sol, m_sol_error )

    if hierarchy == 'normal':
        m31 = 2.42E-3
        m31_error = 0.06E-3
        if wiggle == 'yes':
            m31 = random.normalvariate( m31, m31_error )
        m32 = m31 - m_sol

    if hierarchy == 'inverted':
        m32 = 2.42E-3
        m32_error = 0.06E-3
        if wiggle == 'yes':
            m32 = random.normalvariate( m32, m32_error )
        m31 = m32 - m_sol

    m21 = m_sol
    s212 = math.sin( 2*theta_12 )**2
    s213 = math.sin( 2*theta_13 )**2
    c13 = math.cos( theta_13 )
    c12 = math.cos( theta_12 )
    s12 = math.sin( theta_12 )

    return ( 1 - c13**4 * s212 * math.sin( 1.27 * m21 * distance / energy )**2
             - s213 * ( c12**2 * math.sin( 1.27 * m31 * distance / energy )**2 + 
                        s12**2 * math.sin( 1.27 * m32 * distance / energy )**2 ) )

class random_asym_normal:
    bins = 1000
    xmin = 0
    xmax = 0
    bin_width = 0
    left_max = 0
    right_max = 0
    dist = []
    x_content = []
    top = 1
    def __init__(self, mu, sigma_left, sigma_right):
        xmax = mu + 5 * sigma_right
        xmin = mu - 5 * sigma_left
        bin_width = (xmax - xmin) / bins
        left_max = self.gaus( mu, mu, sigma_left )
        right_max = self.gaus(mu, mu, sigma_right)
        for i in range(bins):
            xvalue = xmin + bin_width * i
            self.x_content.append(xvalue)
            if xvalue < mu:
                self.dist.append( self.gaus(xvalue, mu, sigma_left ) / left_max )
            else:
                self.dist.append( self.gaus(xvalue, mu, sigma_right ) / right_max )

        norm = sum(self.dist)
        for i in range(bins):
            self.dist[i] /= norm
        self.top = 1/norm

    def gaus(self,x,mu,sigma):
        return (2 * math.pi * sigma**2 )**(-0.5) * math.e **( -(x-mu)**2 / (2 * sigma**2) )

    def random(self):
        size = len( self.x_content )
        random_x = int(random.uniform(0, size))
        random_y = random.uniform(0, self.top)
    
        if self.dist[random_x] > random_y:
            return self.x_content[random_x]
        else:
            return self.random()


if __name__ == "__main__":
    main()
