import numpy as np
import pylab as lab
import sys

def main():

    # normal = open( 'normal_events.txt', 'r' )
    print( sys.argv )
    data_file = open( str(sys.argv[1]) )
    
    # First two lines are normal, then inverted infinite statistics, no wiggle

    first_line = data_file.readline()
    second_line = data_file.readline()
    normal_true = []
    inverted_true = []

    for i in range( len(first_line.split()) ):
        normal_true.append( float( first_line.split()[i] ) )
        inverted_true.append( float( second_line.split()[i] ) )

    chi_squared_normal = []
    chi_squared_inverted = []
    test = []

    where_am_i = 0

    for line in data_file:
        data = []
        for i in range( len(line.split()) ):
            data.append( float( line.split()[i] ) )
        chi_squared_normal.append( chi_squared( data, normal_true ) )
        chi_squared_inverted.append( chi_squared( data, inverted_true ) )
        where_am_i += 1
        if where_am_i % 100 == 0:
            print("event: ", where_am_i)
    

    bins = 100

    print( "Average Chi_Squared_Normal: ", sum(chi_squared_normal)/len(chi_squared_normal) )
    print( "Average Chi_Squared_Inverted: ", sum(chi_squared_inverted)/len(chi_squared_inverted) )

    lab.hist( chi_squared_normal, bins, histtype='step', label='Chi^2 Normal Hierarchy')
    lab.hist( chi_squared_inverted, bins, histtype='step', label='Chi^2 Inverted Hierarchy')
    lab.xlabel('Chi^2 Value')
    lab.legend()
    lab.title('10 kton*years Normal Hierarchy Data')
    #lab.plot( normal_true )
    #lab.plot( inverted_true )
    lab.show()

def chi_squared( one, two ):
    sum_var = 0
    for i in range( len(one) ):
        if two[i] != 0:
            sum_var += (one[i] - two[i])**2 / two[i]

    return sum_var

if __name__ == "__main__":
    main()
