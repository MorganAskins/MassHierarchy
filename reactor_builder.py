import pylab as lab
import numpy as np
import math
from matplotlib import animation

bins = 10000
xmin = .1
xmax = 10
bin_width = (xmax - xmin) / bins
# normalizations
core_235U = 0.58
core_239Pu = 0.3
core_238U = 0.07
core_241Pu = 0.05
distance = 50000                # meters

# constants
mass_electron = 0.511
mass_neutron = 939.56
mass_proton = 938.27
delta = mass_neutron - mass_proton
event_rate = 25                 # total per day
k = 0.035                   # Detector Energy Resolution



def main():

    m_atm = 2.42E-3
    m_atm_lower = 0.06E-3
    m_atm_upper = 0.06E-3

    d23 = m_atm - 3 * m_atm_lower
    m23_errors = random_asym_normal(m_atm, m_atm_lower, m_atm_upper)

    # lab.ion()
    fig = lab.figure(figsize=(16.0, 9.0))

    count = 0

    while d23 < m_atm + 3 * m_atm_upper:
        count += 1
        normal = build_spectrum( 'normal', 'no', m23_errors, m_atm )
        inverted = build_spectrum( 'inverted', 'no', m23_errors, d23 )
#        inverted = build_spectrum( 'inverted', 'no', m23_errors, d23 )
        diff = []
        ratio = []
        for i in range(len(normal)):
            diff.append( (normal[i] - inverted[i]) )
            if inverted[i] > 0:
                ratio.append( (normal[i] / inverted[i]) * 0.1 )
            else:
                ratio.append(0.1)
        d23 += 0.001E-3
        # print( d23 )
        lab.clf()
        lab.xlabel('Energy(10keV)')
        lab.title('Normal vs Inverted hierarchy at ' + str(int(distance/1000)) + 'km with changing $\Delta m_{32}^2$')
        lab.axis([0,1000,-0.02,0.14])
        lab.plot(ratio, label='Ratio')
        lab.plot(diff, label='Difference')
        lab.plot(normal, label='Normal')
        lab.plot(inverted, label='Inverted')
        lab.legend(fontsize='xx-large',bbox_to_anchor = (0.8, 0.6), fancybox=True, title=('$\Delta m_{32}^2 = $' + '%.6f'%(d23) ))
        fig.savefig('movies/fig_' + ('%04d' % count) + '.png')
        
        # lab.pause(0.001)
    
    # hmm = []
    # for i in range(len(normal)):
    #     if inverted[i] > 0:
    #         hmm.append(normal[i]/inverted[i])
    #     else:
    #         hmm.append(1)

    # lab.plot(normal)
    # lab.plot(inverted)
    # lab.plot(hmm)

def build_spectrum_LoverE( hierarchy, wiggle, m23_errors, d23):
    # Overloads normal build_spectrum with the bool LoverE
    spectrum = []
    inf_inverted = []
    
    # Convert to L/E
    new_xmax = distance / (xmin)
    new_xmin = distance / (xmax)
    new_bin_width = (new_xmax - new_xmin)/bins
    # Fill with neutrino events
    for i in range(bins):
        xvalue = distance / (new_bin_width * i + xmin + delta + mass_electron)
        spectrum.append( cross_section(xvalue) * reactor(xvalue)  * survival(xvalue, hierarchy, wiggle, m23_errors, d23))

    # Normalize
    norm = sum(spectrum)
    for i in range(bins):
        spectrum[i] *= event_rate / norm

    return spectrum

def build_spectrum( hierarchy, wiggle, m23_errors, d23 ):
    spectrum = []
    inf_inverted = []

    # Fill with neutrino events
    for i in range(bins):
        xvalue = bin_width * i + xmin + delta + mass_electron
        spectrum.append( cross_section(xvalue) * reactor(xvalue)  * survival(xvalue, hierarchy, wiggle, m23_errors, d23))

    # Normalize
    norm = sum(spectrum)
    for i in range(bins):
        spectrum[i] *= event_rate / norm

    return spectrum

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

def survival( energy, hierarchy, wiggle, m23_errors, d23 ):

    theta_12 = ( 33.57 ) * ( math.pi / 180 )
    theta_12_error = ( 0.76 ) * ( math.pi / 180 )

    theta_13 = (8.71) * (math.pi / 180 )
    theta_13_error = ( 0.38 ) * (math.pi / 180)

    m_sol = 7.45E-5
    m_sol_error = 0.18E-5

    # if hierarchy == 'normal':
    #     m31 = 2.42E-3
    #     m31_error = 0.06E-3
    #     m32 = m31 - m_sol
    # if hierarchy == 'inverted':
    #     m32 = 2.42E-3
    #     m32_error = 0.06E-3
    #     m31 = m32 - m_sol

    # if hierarchy == 'normal':
    #     m32 = 2.37E-3
    #     m32 = d23
    #     m32_error = 0.08E-3
    #     m31 = m32 - m_sol
    # if hierarchy == 'inverted':
    #     m32 = 2.43E-3
    #     m32_error = 0.10E-3
    #     m31 = m32 + m_sol

    if wiggle == 'yes':
        # Smear the angles with a gaussian
        theta_12 = random.normalvariate( theta_12, theta_12_error )
        theta_13 = random.normalvariate( theta_13, theta_13_error )
        m_sol = random.normalvariate( m_sol, m_sol_error )
        m_atm_adjusted = random.normalvariate( m32, m32_error )

    m21 = m_sol
    if hierarchy == 'normal':
        m31 = d23
        m32 = m31 - m21
    else:
        m32 = d23
        m31 = m32 - m21

    #changing it up a bit
    # if hierarchy == 'normal':
    #     m31 = 2.55E-3
    #     m32 = m31 - m21
    # else:
    #     m31 = 2.43E-3
    #     m32 = m31 - m21

    s212 = math.sin( 2*theta_12 )**2
    s213 = math.sin( 2*theta_13 )**2
    c13 = math.cos( theta_13 )
    c12 = math.cos( theta_12 )
    s12 = math.sin( theta_12 )

    # if hierarchy == 'normal':
    #     m32 = d23
    #     m31 = m21 + m32
    # else:
    #     m31 = d23
    #     m32 = m21 + m31

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
