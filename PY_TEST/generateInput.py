import sys
import random

# Serve per creare i file di input da dopo dare in pasto a FFT.c 
# generando (int) size numeri immaginari presi a random in un range da (int) lowRange e (int) highRange
# Utilizzo: python generateInput1.py size lowRange highRange
filenameToCreate = ""
size = 0
lowRange = 0
highRange = 0
print("Arguments: " + str(sys.argv) + "\n")
if (len(sys.argv) == 4):
    size = int(sys.argv[1])
    filenameToCreate = str(size) + "random"
    lowRange = int(sys.argv[2])
    highRange = int(sys.argv[3])
else:
    print("Error, bad arguments")
    sys.exit()

fileToCreate = open(filenameToCreate, "w+")
for i in range(size):
    realInt = random.randint(lowRange, highRange)
    realFloat = 0 #random.randint(0, 99)
    imgInt = random.randint(lowRange, highRange)
    imgFloat = 0 #random.randint(0, 99)
    real = str(realInt) + "." + str(realFloat)
    img = str(imgInt) + "." + str(imgFloat)
    strToWrite = str(i) + " " + real + " " + img + "\n"
    fileToCreate.write(strToWrite)



fileToCreate.close()
