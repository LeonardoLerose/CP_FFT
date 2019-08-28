import sys
import random

# Serve per creare i file di input da dopo dare in pasto a FFT.c 
# generando un vettore di dimensione numReps*patternSize da due pattern, Re e Im, da ripetere
# Utilizzo: python generateInput1.py numReps patternRe patternIm
# Esempio: generateinput2.py 4 1-2 3-4   Crea: [(1,3),(2,4),(1,3),(2,4),(1,3),(2,4),(1,3),(2,4)]

filenameToCreate = ""
numReps = 0
lowRange = 0
highRange = 0
print("Arguments: " + str(sys.argv) + "\n")
if (len(sys.argv) == 4):
    numReps = int(sys.argv[1])
    patternRe = str(sys.argv[2])
    patternIm = str(sys.argv[3])
    patternReInt = [int(n) for n in patternRe.split("-")] 
    patternImgInt = [int(n) for n in patternIm.split("-")] 
    filenameToCreate = "pattern" + str(len(patternReInt)) + "_" + str(numReps)

else:
    print("Error, bad arguments")
    sys.exit()

fileToCreate = open(filenameToCreate, "w+")
sample = []
for i, j in zip(patternReInt, patternImgInt):
    sample.append(str(i) + " " + str(j))

print(sample)

for i in range(numReps*len(sample)):
    strToWrite = str(i) + " " + sample[i%len(sample)]+ "\n"
    print(strToWrite)
    fileToCreate.write(strToWrite)



fileToCreate.close()

