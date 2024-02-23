# Model of non-trivial SQUID Circuits 

# CPR = Ic ^ (-10*n) & phi ^ (-2*n)

# Authors: Cliff Sun, Harshvardhan, Professor Bezryadin

import numpy as np
import matplotlib.pyplot as plt

arrayOfJunctions = [0, 0.001, 0.999, 1] # Where the junctions are located; must be even amount of junctions (includes 0 and 1)

criticalCurrents = np.ones(int(len(arrayOfJunctions) / 2)) # all the critical currents for each junction (Critical Currents are assumed to be 1 by default)

# prints out the elements in the junction in a better format

index = 0
while (index < (len(arrayOfJunctions) - 1)):
    if (index == len(arrayOfJunctions) - 2):
        print(str(arrayOfJunctions[index]) + " - " + str(arrayOfJunctions[index + 1]), end = " ")
    else:
        print(str(arrayOfJunctions[index]) + " - " + str(arrayOfJunctions[index + 1]) + ",", end = " ")
    index += 2

# Parameters for Current -> Magnetic Field, Junction Locations, Critical Currents, Initial Phase Difference

# B is the integer value of flux quanta present in the SQUID

def current(B, arrJ, arrC, y, numOfSegments): # y is initial phase difference of the whole circuit, B is the magnetic field, arrJ is the location of junctions, arrC is critical current associated with each junction

    curr = 0 # summation of all currents in the entire junction

    limit = int(len(arrJ) / 2) # number of junctions in the SQUID

    for n in range(limit):
        sizeOfSegment = float((arrJ[2 * n + 1] - arrJ[2 * n]) / numOfSegments)
        for i in range(numOfSegments):
            curr += float(arrC[n]/(10**n)) * np.sin((y + (2 * np.pi * B) * (arrJ[2 * n] + i * sizeOfSegment))/(2**n)) * (1/numOfSegments)

    # phase difference evolves according to 2 * pi * B

    # curr += (critical current element in array)(sin(y + (2 * pi * B) * length)

    return curr

def maxCurrent(B, arrayJ, arrayC, numOfSegments): # Spits out the maximum current by varying the gauge invariant phase of the left end (free parameter) gamma
    Y=np.linspace(0, 2*np.pi, 150)
    dummyArray=[]
    for gamma in Y:
        dummyArray.append(current(B, arrayJ, arrayC, gamma, numOfSegments))
    return max(dummyArray)

MagField = np.linspace(0, 10, 5000) # an array of Magnetic Fields ranging from 0 to 100 with 5000 total elements

# The 2 lines below is where the useful section of the code is for modeling a SQUID

IMaxPoint = []

numOfSegments = 40

for B in MagField:
    IMaxPoint.append(maxCurrent(B, arrayOfJunctions, criticalCurrents, numOfSegments) / int(len(arrayOfJunctions) / 2)) # This integer represents the number of segments you want to cut each junction up into (the higher the number, the better the approximation)

# The 2 lines above is where the useful section of the code is for modeling a SQUID

plt.plot(MagField, IMaxPoint, 'g-')
plt.xlabel('Flux')
plt.ylabel('Critical Current')  
plt.grid()             
plt.show()