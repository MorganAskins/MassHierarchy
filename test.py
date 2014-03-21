import multiprocessing as mp

def main():

    pool = mp.Pool(processes = 4)
    mymap = pool.map(f, range(10) )
    print ( mymap )
    
    

def f(x,y):
    r = 0
    for k in range(1, 5020000):
        r += x ** ( 1 / k**1.5 ) * y
    return r

if __name__ == '__main__':
    main()
