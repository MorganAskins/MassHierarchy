import math
import numpy
import pylab as lab

theta_12 = ( 33.57 ) * ( math.pi / 180 )
dm12 = 7.45E-5
dm32_error = 6E-5

distance = 13000

def main():

    # e1 = []
    # e2 = []
    # size = 0.1E-5
    # for i in range(600):
    #     value = i * size + dm32 - 4*dm32_error
    #     e1.append( equation1(value) )
    #     e2.append( equation2(value) )
    
    x1=[]
    x2=[]
    y=[]
    size = 0.1E-5
    dm32 = 2.42E-3 - dm12
    lab.ion()
    for i in range(1,100):
        xvalue = i / 10
        x1.clear()
        x2.clear()
        for j in range(100):
            xxvalue = j * size + dm32 - 4*dm32_error
            x1.append(equation1( xxvalue, dm32 ,xvalue ))
            x2.append(equation2( xxvalue, dm32, xvalue ))
        lab.clf()
        lab.plot(x1)
        lab.plot(x2)
        lab.pause(1)
        
        
        #y.append(intersect(dm32,xvalue))

    lab.plot(x,y)
    lab.show()
    # normal = []
    # inverted = []
    # size = 0.1E-5
    # top = int( 5 * dm32_error / size )
    # print(top)
    # for i in range(top):
    #     value = i * size + dm32 - 5*dm32_error
    #     normal.append(value)
    #     inverted.append( intersect(value) )
        
    # lab.rc('font', family='serif')
    # #lab.rc('text', usetex=True)
    # lab.xlabel('$\Delta m_{32}^2$ Inverted Hierarchy')
    # lab.ylabel('$\Delta m_{32}^2$ Normal Hierarchy')
    # lab.title('Normal / Inverted Hierarchy Degeneracy from $\Delta m_{32}^2$')
    # lab.plot(normal,inverted)
    # lab.show()

def equation1( delta, dm32, energy ):
    return math.sin( (dm32 + dm12)*1.27*distance/energy )**2 - math.sin( (delta - dm12)*1.27*distance/energy )**2

def equation2( delta, dm32, energy ):
    return math.tan(theta_12)**(-2) * ( math.sin(delta*1.27*distance/energy)**2 - math.sin(dm32*1.27*distance/energy)**2 )


def intersect(dm32,energy):
    step_size = 0.001E-3
    begin = dm32 - 4 * dm32_error
    end = dm32 + 4 * dm32_error
    position = begin
    mini = 1
    mini_spot = 0

    while position < end:
        diff = equation1( position, dm32, energy ) - equation2( position, dm32, energy )
        if math.fabs(diff) < mini:
            mini = math.fabs(diff)
            mini_spot = position
        position += step_size

    return mini_spot

if __name__ == '__main__':
    main()
