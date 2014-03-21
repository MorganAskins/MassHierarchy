import matplotlib.pyplot as plt
import math

gamma = 3.77
xi = 2.5
alpha = 0.268
epsilon = 618

def main():
    bins = 1000
    x = []
    y = []
    for i in range(1,bins):
        value = i / bins * 7
        x.append(value)
        y.append( ratio(value) )
    plt.plot(x,y)
    plt.show()



def ratio(distance):
    return gamma * Energy(distance) * math.e**(distance/xi) / ( epsilon * (math.e**(distance/xi) -1 ))

def Energy(distance):
    return alpha * distance

if __name__ == "__main__":
    main()
