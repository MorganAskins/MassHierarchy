# Generate the geoneutrino spectrum from Nikolai's thesis
import math
import matplotlib.pyplot as plt
import scipy.special as sp

# some constants
me = 0.511
hbar = 196.3269718
alpha = 1/137.0388

def main():
    # Th/U ratio is ~ 4.2
    thuratio = 0.5
    bins = 1000
    x = []
    y1 = []
    y2 = []
    max_x = 6
    for i in range(0, bins):
        x.append(i / bins * max_x)
        enu = i/bins * max_x
        y1.append( ( Dist(91, 2.194, enu, 10000, 234) ))  # 234Pa
                     #+ Dist(83, 3.279, enu, 15, 214 )   # 214Bi
                     #+ Dist(81, 5.482, enu, 15, 210 )   # 210Tl
                     #+ Dist(89, 2.124, enu, 10000 * thuratio, 228) # 228Ac
                     #+ Dist(83, 2.252, enu, 10000 * thuratio, 212) ))# 212Bi
                     #+ Dist(81, 4.998, enu, 10000 * thuratio, 208) ) # 208Tl
                   #* cross( enu ) )
        y2.append( Dist(91, 2.194, enu, 10000, 234)   # 234Pa
                   + Dist(83, 3.279, enu, 15, 214 )   # 214Bi
                   + Dist(81, 5.482, enu, 15, 210 )   # 210Tl
                   + Dist(89, 2.124, enu, 10000 * thuratio, 228) # 228Ac
                   + Dist(83, 2.252, enu, 10000 * thuratio, 212) )# 212Bi
                   #+ Dist(81, 4.998, enu, 10000 * thuratio, 208) ) # 208Tl

        #print("x: ", x[i], " y: ", y[i])
    plt.yscale('log')
    plt.plot(x,y1)
    #plt.plot(x,y2)
    plt.show()

def Pe( eo, enu ):              # eo is endpoint energy of beta decay
    return ( (eo - enu) * (eo- enu + 2 * me) )**0.5

def We( eo, enu ):
    return eo - enu + me

def Y( z, eo, enu ):
    return alpha * z * We(eo,enu) / Pe(eo,enu)

def gamma(z):
    return ( 1 - (alpha * z)**2 )**0.5

def Weo(enu):
    if enu > me + 1.3:
        return enu - 1.3
    return 0

def Peo(enu):
    a = ( Weo(enu)**2 - me**2 )**0.5
    return a

def cross(enu):
    return 0.0952 * Weo(enu) * Peo(enu)

def Fermi(z, eo, enu, L, A):
    value = L * 4 * (2 * Pe(eo, enu) * 1.2 * A**(1/3) / hbar )**( 2 * gamma(z) - 2 ) * math.e ** (math.pi * Y(z, eo, enu) )* abs( sp.gamma(complex(gamma(z), Y(z,eo,enu) ) ) )**2 / abs( sp.gamma(2*gamma(z) + 1) )**2
    return value

def Heaviside(x):
    if x==0: 
        return 0.5
    return 0 if x < 0 else 1

def Dist(z,eo,enu,L,A):
    if enu < eo:
        return Pe(eo, enu)*We(eo, enu)*(enu)**2*Fermi(z,eo,enu,L,A)
    return 0

if __name__ == "__main__":
    main()
