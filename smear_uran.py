# Smear the uranium decay data
import random
import math
import pylab as plt
import numpy as np

def main():
    k = 0.06
    uran = np.loadtxt('imb_uran.dat')
    thor = np.loadtxt('imb_thor.dat')
    reactor = np.loadtxt('imb_reactor_bkg.dat')
    a = uran + thor + reactor
    size = len(a)
    x = []                      # x values
    y = [0]*size

    for i in range(0,size):
        x.append(i)
        for j in range(0,size):
            y[j] += a[i] * gaus(j,i,k*i**(0.5))

    print("signal: ", sum(a), " new: ", sum(y))
    #plt.plot(x,uran,'r',x,thor,'g',x,reactor,'b')
    plt.plot(x,a,'r',x,y,'g')
    plt.show()

def gaus(x, mu, sigma):
    if sigma == 0:
        return 0
    return (2*math.pi*sigma**2)**(-0.5)*math.e**(-(x-mu)**2/(2*sigma**2))


if __name__ == "__main__":
    main()
