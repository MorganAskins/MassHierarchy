# Smear a function with Root(E) resolution
import random
import math
import pylab as plt
import numpy as np

def main():
    k = 1.8
    bins = 2000
    # define the sample histogram bounded from (0,10)
    x = []
    y = []
    z = [0]*bins
    m = []
    xmax = 100
    xmin = 1
    bin_width = (xmax - xmin) / bins
    for i in range(0,bins):
        location = i * bin_width + xmin
        x.append(location)
        value = signal(location)
        y.append( value )
        # convolution stuff
        for j in range(0,bins):
            subloc = j * bin_width + xmin
            z[j] += value * bin_width * gaus(subloc,location,k*(location)**(0.5))

    print("signal: ", sum(y), " new: ", sum(z))
    #plt.plot(x,y)
    plt.plot(x,z,'r', x,y,'g')
    plt.show()


def gaus(x,mu,sigma):
    return (math.pi*sigma**2)**(-0.5)*math.e**(-(x-mu)**2/sigma**2)

def signal(x):
    return math.sin(x)**2

def congaus(x,sig):
    return (math.pi * x**2)**(-0.5)*math.e**(-(x-sig)**2/x)

# def signal(x):
#     if x > 1.9 and x < 2.1:
#         return 1
#     if x > 4.9 and x < 5.1:
#         return 1
#     if x > 8.9 and x < 9.1:
#         return 1
#     return 0

if __name__ == "__main__":
    main()
