import numpy as np
import sys

if __name__ == "__main__":

    size = int(sys.argv[1])

    for x in range(size):
        if sys.argv[2] == 'sin':
            print(str(x) + ' ' + str(np.sin(x/10.)) + ' 0')
        elif sys.argv[2] == 'cos':
            print(str(x) + ' ' + str(np.cos(x/10.)) + ' 0')
        elif sys.argv[2] == 'line':
            print(str(x) + ' ' +  str(x) + ' 0')
        elif sys.argv[2] == 'rect':
            if x < size/2:
                print(str(x) + ' 1 0')
            else:
                print(str(x) + ' 0 0')
        else:
            raise Exception('Bad input argument:' + str(sys.argv))