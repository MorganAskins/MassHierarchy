# Studying the effects of a changing dm32
import math
import pylab as lab

distance = 13000

def main():
    x = []
    normal = []
    inverted = []
    diff = []
    bins = 1000
    xmin = 0.001
    xmax = 10
    bin_width = (xmax - xmin)/bins

    for i in range(bins):
        xvalue = xmin + i * bin_width
        x.append(xvalue)
        normal.append( survival( xvalue, 'normal' ) )
        inverted.append( survival( xvalue, 'inverted' ) )
        diff.append( normal[i] - inverted[i] )

    lab.plot(x,normal)
    lab.plot(x,inverted)
    lab.show()

def survival( energy, hierarchy ):
    theta_12 = ( 33.57 ) * ( math.pi / 180 )
    theta_12_error = ( 0.76 ) * ( math.pi / 180 )

    theta_13 = ( 8.71 ) * ( math.pi / 180 )
    theta_13_error = ( 0.38 ) * ( math.pi / 180 )

    m_sol = 7.45E-5
    m_sol_error = 0.18E-5

    if hierarchy == 'normal':
        m31 = 2.42E-3
        m31_error = 0.06E-3
        m32 = m31 - m_sol

    if hierarchy == 'inverted':
        m32 = 2.42E-3
        m32_error = 0.06E-3
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

if __name__ == '__main__':
    main()
