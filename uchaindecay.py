# Simple script to plot a decay chain in python
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import time
import random
import math

def main():
    # ignore fast decays
    u238 = 100000
    u234 = 0
    th230 = 0
    ra226 = 0
    po210 = 0

    Pb = 0
    u238_hl = 1/10000.0
    u234_hl = 1/2455.0
    th230_hl = 1/753.800
    ra226_hl = 1/16.020
    po210_hl = 1/1.38376

    x = [1,2,3,4,5]
    y = []

    plt.ion()

    time = 0

    while True:
        time += 1
        u238, u234 = decay( u238, u234, u238_hl )
        u234, th230 = decay( u234, th230, u234_hl )
        th230, ra226 = decay( th230, ra226, th230_hl )
        ra226, po210 = decay( ra226, po210, ra226_hl )
        po210, Pb = decay( po210, Pb, po210_hl )
        #time.sleep(0.001)           # seconds?
        u238 = 100000
        if time%1000 == 0 :
            y.clear()
            y.append( u238 )
            y.append( u234 )
            y.append( th230 )
            y.append( ra226 )
            y.append( po210 )
            plt.clf()
            plt.plot(x,y, 'ro')
            plt.draw()
            print( "238U: ", u238, " 234U: ", u234, " 230Th: ", th230, " 226Ra ", ra226, " 210Po ", po210)
    

def decay(mother, daughter, mother_life):
    # mother will decay and give to daughter
    change = mother_life * mother
    mother -= change
    daughter += change
    return mother, daughter

if __name__ == "__main__":
    main()
