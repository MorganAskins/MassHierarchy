import pylab as lab
import numpy as np
import math
from matplotlib import animation
import scipy.fftpack

import reactor_builder as builder

#distance = 25000

def main():
    m_atm = 2.42E-3
    m_atm_lower = 0.06E-3
    m_atm_upper = 0.06E-3

    m23_errors = builder.random_asym_normal(m_atm, m_atm_lower, m_atm_upper)
    distance = builder.distance
    builder.distance = 50000

    d23 = m_atm
    normal = builder.build_spectrum_LoverE('normal', 'no', m23_errors, m_atm)
    inverted = builder.build_spectrum_LoverE('inverted', 'no', m23_errors, m_atm)
    diff = []
    ratio = []
    for inv, nor in zip(inverted, normal):
        diff.append( nor - inv )
        if inv > 0:
            ratio.append( nor * 0.17 / inv )
        else:
            ratio.append(0.17)

    normal_E = builder.build_spectrum('normal', 'no', m23_errors, m_atm)
    inverted_E = builder.build_spectrum('inverted', 'no', m23_errors, m_atm)
            
    lab.subplot(221)
    lab.xlabel('L/E '+str(builder.distance/1000)+' m/MeV')
    lab.title('Normal vs Inverted hierarchy at ' + str(int(distance/1000)) + 'km with changing $\Delta m_{32}^2$')
    lab.axis([0,1000,-0.02,0.2])
    lab.plot(ratio, label='Ratio')
    lab.plot(diff, label='Difference')
    lab.plot(normal, label='Normal')
    lab.plot(inverted, label='Inverted')
    #lab.legend(fontsize='xx-large',bbox_to_anchor = (0.8, 0.6), fancybox=True, title=('$\Delta m_{32}^2 = $' + '%.6f'%(d23) ))
    #lab.show()

    # Now perform a fourier spectrum on the data
    bins = len(normal)
    bin_width = 1
    fourier_depth = 1000000
    myfft_normal = abs(scipy.fft(normal, fourier_depth))**2
    myfft_inverted = abs(scipy.fft(inverted, fourier_depth))**2
    freqs = scipy.fftpack.fftfreq(myfft_normal.size, bin_width)
    
    # w=[]
    # for i in range(len(myfft)):
    #     w.append(i*1/(bins*bin_width))

    difference = np.array(myfft_normal) - np.array(myfft_inverted)
    print(len(difference), len(freqs))
    
    lab.subplot(222)
    lab.plot(freqs, difference)
    #lab.plot(normal_E, label='Normal')
    #lab.plot(inverted_E, label='Inverted')
    lab.subplot(223)
    lab.plot(freqs, myfft_normal)
    lab.plot(freqs, myfft_inverted)
    lab.xlim(-0.1, 0.1)
    lab.show()


if __name__ == "__main__":
    main()
