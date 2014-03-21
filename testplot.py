import pylab as lab
import sys

def main():
    data = open( sys.argv[1], 'r' )
    nline = data.readline()
    iline = data.readline()

    norm = []
    inv = []
    ratio = []

    for i in range( len(nline.split())):
        norm.append( float( nline.split()[i] ) )
        inv.append( float( iline.split()[i] ) )
        if inv[i] > 0:
            ratio.append( norm[i] / inv[i] )
        else:
            ratio.append(1)

    # Normalize:
    normal = sum(norm) / 300
    ratio = []
    for i in range( len(norm) ):
        norm[i] /= normal
        inv[i] /= normal
        if inv[i] > 0:
            ratio.append( norm[i] / inv[i] )
        else:
            ratio.append(1)

    lab.plot(norm, 'b', label='Normal Hierarchy')
    lab.plot(inv, 'g', label='Inverted Hierarchy')
    lab.plot(ratio, 'r', label='Ratio Normal / Inverted' )
    
    lab.xlabel('Energy (10 keV) ')
    lab.ylabel('Events/10 keV')
    lab.title('Infinite Statistics at 13 km')
    lab.legend()
    #lab.plot(ratio, 'r')
    lab.show()

if __name__ == "__main__":
    main()
