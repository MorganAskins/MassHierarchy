# Li9 decay spectrum
# For each decay branch, must draw beta spectrum convolved with width
import numpy as np
import math
import pylab as plt
import scipy.special

simple = False

def main():
    # Some constants of nature
    m_electron = 0.511          # MeV
    Z = 3                       # protons in Lithium
    A = 9                       # Li9
    alpha = 1/137.035           # Fine structure
    S = (1 - (alpha*Z)**2)**0.5
    rho = 1.25E-15 * A**(1/3) / 6.582E-22
    
    # List The decay branching ratios and decays widths
    Q_beta = 13.6063
    # Four branches:
    Q_branch = [2.4294, 2.780, 7.940, 11.283]
    ratio_branch = [0.34,0.1,0.015,0.04]
    # normalize branching ratios:
    total_branch = 0
    for i in range(4):
        total_branch += ratio_branch[i]
    for i in range(4):
        ratio_branch[i] = ratio_branch[i] / total_branch
    
    width = [0.00077, 1.08, 1, 0.575]

    file = open('Li9_Spectrum.dat','w')

    xmin = 0.001
    xmax = 12                   # MeV
    bins = 1200
    bin_width = (xmax - xmin) / bins
    x = []
    branch = [[],[],[],[]]
    conv_gaus = [[0]*bins,[0]*bins,[0]*bins,[0]*bins]
    for i in range(0,bins):
        xvalue = i * bin_width + xmin
        x.append(xvalue)
        for j in range(0,4):
            if simple:
                branch[j].append( simple_beta( xvalue, Q_beta - Q_branch[j], m_electron ) )
            else:
                branch[j].append( beta( P( E( xvalue, m_electron), m_electron ), 
                                        E( xvalue, m_electron ), (Q_beta - Q_branch[j]), xvalue, A, Z, S, rho ) )
            for k in range(bins):
                conv_gaus[j][k] += branch[j][i] * gaus( k * bin_width + xmin, xvalue, width[j]/2. )
            #conv_gaus[j].append( gaus(xvalue, 0, ratio_branch[j] ) )
    
    # Normalize
    rate = 6.24                 # per kton / day
    for i in range(4):
        norm = sum(conv_gaus[i])
        norm2 = sum(branch[i])
        for j in range(len(conv_gaus[i])):
            conv_gaus[i][j] *= rate * ratio_branch[i] / norm
            branch[i][j] *= rate * ratio_branch[i] / norm2

    spectrum = []
    for i in range(len(conv_gaus[0])):
        spectrum.append( conv_gaus[0][i] + conv_gaus[1][i] + conv_gaus[2][i] + conv_gaus[3][i] )

    for i in range(4):
        plt.plot(x,conv_gaus[i])

    for i in range( len(spectrum) ):
        file.write( str(spectrum[i]) + '\n' )
    print("Events per kton per day: ", sum(spectrum) )

    plt.plot(x,spectrum)
    plt.show()

def simple_beta(x, Q, m):
    if Q - x > 0:
        return (x**2 + 2*x*m)**(0.5)*(Q-x)**2*(x+m)
    return 0

def E(T,m):
    return T + m

def P(e,m):
    return (e**2 - m**2)**0.5

def eta(Z,e,p):
    alpha = 1 / 137.035
    return alpha * Z * e * p

def fermi(Z,S,p,rho,eta):
    return ( 2 * (1+S) )/(math.gamma(1+2*S)**2)*(2*p*rho)**(2*S-2)*math.e**(math.pi * eta)*abs(scipy.special.gamma(complex(S,eta)))**2

def beta(p,e,Q,T,A,Z,S,rho):
    if Q - T > 0:
        return p**2 * (Q - T)**2 * fermi(Z,S,p,rho, eta(Z,e,p) )
    return 0

def spectrum(x,ratio,width,energy):
    return ratio * gaus(x,energy,width/2)

def gaus(x,mu,sigma):
    return (2*math.pi*sigma**2)**(-0.5)*math.e**(-(x-mu)**2/(2*sigma**2))

if __name__ == "__main__":
    main()
