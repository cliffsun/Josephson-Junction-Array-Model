import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, e

def EvenArrayOfJunctions(sigma, numOfJunctions, width, arrayJ = []): # Generates an array of junctions given some standard deviation & mean width
    junctionCenter = 1/(numOfJunctions - 1)
    arrOfWidthsDiv2 = np.zeros(numOfJunctions)
    arrOfWidths = np.zeros(numOfJunctions)
    arrOfJunctions =  np.zeros(numOfJunctions * 2)
    delta = 4 * (numOfJunctions - 1) / ((numOfJunctions) * (numOfJunctions + 2)) * sigma**2
    sum = 0
    if sigma == 0 and len(arrayJ) != 0:
        return arrayJ
    for i in range(1, numOfJunctions//2 + 1):
        sum += i * delta
        arrOfWidthsDiv2[i - 1] = np.sqrt(i * delta)
    for j in range(0, numOfJunctions, 2): 
        arrOfWidths[j] = width + arrOfWidthsDiv2[j//2]
        arrOfWidths[j + 1] = np.abs(width - arrOfWidthsDiv2[j//2])
    if (len(arrayJ) != 0):
        arrOfJunctions = np.array(arrayJ)
    arrOfJunctions[1] = arrOfWidths[0]
    for k in range(2, arrOfJunctions.size - 2, 2):
        arrOfJunctions[k] = arrayJ[k] - (arrOfWidths[k//2]/2 - width) if (len(arrayJ) != 0) else k//2 * junctionCenter - arrOfWidths[k//2]/2
        arrOfJunctions[k + 1] = arrayJ[k + 1] + (arrOfWidths[k//2]/2 - width) if (len(arrayJ) != 0) else k//2 * junctionCenter + arrOfWidths[k//2]/2
    arrOfJunctions[-1] = 1
    arrOfJunctions[-2] = arrayJ[-2] - (arrOfWidths[-1] - width) if len(arrayJ) != 0 else 1 - arrOfWidths[-1]
    return arrOfJunctions

def OddArrayOfJunctions(sigma, numOfJunctions, width, arrayJ = []): # Generates an array of junctions given some standard deviation & mean width
    junctionCenter = 1/(numOfJunctions - 1)
    arrOfWidthsDiv2 = np.zeros(numOfJunctions)
    arrOfWidths = np.zeros(numOfJunctions)
    arrOfJunctions = np.zeros(numOfJunctions * 2)
    delta = 4 / (numOfJunctions + 1) * sigma**2
    sum = 0
    if sigma == 0 and len(arrayJ) != 0:
        return arrayJ
    for i in range(1, (numOfJunctions - 1)//2 + 1):
        sum += i * delta
        arrOfWidthsDiv2[i - 1] = np.sqrt(i * delta)
    arrOfWidths[numOfJunctions - 1] = width
    for j in range(0, numOfJunctions - 1, 2): 
        arrOfWidths[j] = width + arrOfWidthsDiv2[j//2]
        arrOfWidths[j + 1] = np.abs(width - arrOfWidthsDiv2[j//2])
    if (len(arrayJ) != 0):
        arrOfJunctions = np.array(arrayJ)
    arrOfJunctions[1] = arrOfWidths[0]
    for k in range(2, arrOfJunctions.size - 2, 2):
        arrOfJunctions[k] = arrayJ[k] - (arrOfWidths[k//2]/2 - width) if (len(arrayJ) != 0) else k//2 * junctionCenter - arrOfWidths[k//2]/2
        arrOfJunctions[k + 1] = arrayJ[k + 1] + (arrOfWidths[k//2]/2 - width) if (len(arrayJ) != 0) else k//2 * junctionCenter + arrOfWidths[k//2]/2
    arrOfJunctions[-2] = arrayJ[-2] - (arrOfWidths[-1] - width) if len(arrayJ) != 0 else 1 - arrOfWidths[-1]
    arrOfJunctions[-1] = 1
    return arrOfJunctions

def ArrayOfJunctions(sigma, numOfJunctions, width, arrayJ):
    if (numOfJunctions % 2 == 0):
        return EvenArrayOfJunctions(sigma, numOfJunctions, width, arrayJ)
    else:
        return OddArrayOfJunctions(sigma, numOfJunctions, width, arrayJ)

def stateOfArray(arrJ):
    state = []
    state.append(len(arrJ)//2)
    mean = 0
    for i in range(len(arrJ)//2):
        mean += (arrJ[2*i + 1] - arrJ[2*i])
    state.append(mean/(len(arrJ)//2))
    return state

def meanOfArray(arrJ):
    percentage = 0
    for i in range(len(arrJ) // 2):
        percentage += (arrJ[2 * i + 1] - arrJ[2 * i])
    return  1 - percentage

def checkArray(arrJ):
    arr = arrJ.copy()
    if not np.array_equal(arr, sorted(arr)):
        if arr[1] > arr[2]:
            arr[1] = arr[2]
        for i in range(len(arr)//2 - 1):
            if arr[2 * i + 1] > arr[2 * i + 2]:
                arr[2*i + 1] = arr[2*i + 2]
            if arr[2 * i] > arr[2 * i + 1]:
                arr[2 * i] = arr[2 * i - 1]
    return arr

# Parameters for Current -> Magnetic Field, Junction Locations, Critical Currents, Initial Phase Difference

# B is the integer value of flux quanta present in the SQUID

def current(B, array, y): # This will give us the value of current flowing across the two superconductors depending on gamma (denoted as y), Magnetic Field B, and the location/widths of the junctions
    l=len(array)
    n=1+l/2
    i=0
    f=0
    x=y
    while(i<n):
        if(i==0): #Treating the first junction separately
            f+=np.sin(50*np.pi*B*array[0])*np.sin(x+50*np.pi*B*array[0])/(50*np.pi*B*array[0])
            # x+=2*50*np.pi*B*array[0]
        elif(i>0 and i<n-1):
            diff=array[2*i]-array[2*i-1] #The interior junctions
            f+=np.sin(50*np.pi*B*diff)*np.sin(x+2*50*np.pi*B*array[2*i-1]+50*np.pi*B*diff)/(50*np.pi*B*diff)
            # x+=2*50*np.pi*B*diff
        else:
            diff=1-array[2*i-1] #Last Junction
            # print(array[2*i-1])
            f+=np.sin(50*np.pi*B*diff)*np.sin(x+2*50*np.pi*B*array[2*i-1]+50*np.pi*B*diff)/(50*np.pi*B*diff)
            # x+=2*50*np.pi*B*diff
        i+=1
    return f

def maxCurrent(B, arrayJ): # Spits out the maximum current by varying the gauge invariant phase of the left end (free parameter) gamma
    Y=np.linspace(0, 2*np.pi, 150)
    dummyArray=[]
    for gamma in Y:
        dummyArray.append(current(B, arrayJ, gamma))
    return max(dummyArray)

FluxField = np.linspace(0.001, 0.2, 1000) # an array of Magnetic Fields ranging from 0 to 100 with 5000 total elements

#----------------------------------

arrayOfJunctions = []

# arrayOfJunctions = [0, 0.00001, 0.99999, 1]

# arrayOfJunctions = [0, 0.01, 0.03, 0.13, 0.14, 0.142, 0.16, 0.36, 0.38, 0.68, 0.681, 0.688, 0.7, 0.75, 0.77, 0.83, 0.84, 0.841, 0.89, 1]

number_of_junctions = 5

#----------------------------------

IMaxPointSigma1 = []

IMaxPointSigma2 = []

IMaxPointSigma3 = []

junctionNumber = stateOfArray(arrayOfJunctions)[0] if (len(arrayOfJunctions) != 0) else number_of_junctions

meanWidth = stateOfArray(arrayOfJunctions)[1] if (len(arrayOfJunctions) != 0) else 0.0001


#----------------------------------

Sigma1 = 0

arraySigma1 = checkArray(ArrayOfJunctions(Sigma1, junctionNumber, meanWidth, arrayOfJunctions))

percentageSigma1 = meanOfArray(arraySigma1)

#----------------------------------

Sigma2 = 0.005

arraySigma2 = checkArray(ArrayOfJunctions(Sigma2, junctionNumber, meanWidth, arrayOfJunctions))

percentageSigma2 = meanOfArray(arraySigma2)

#----------------------------------

Sigma3 = 0.01

arraySigma3 = checkArray(ArrayOfJunctions(Sigma3, junctionNumber, meanWidth, arrayOfJunctions))

percentageSigma3 = meanOfArray(arraySigma3)

#----------------------------------


for f in FluxField:
    IMaxPointSigma1.append(maxCurrent(f, arraySigma1[1:-1])) # This integer represents the number of segments you want to cut each junction up into (the higher the number, the better the approximation)
print(arraySigma1)
for f in FluxField:
    IMaxPointSigma2.append(maxCurrent(f, arraySigma2[1:-1]))
print(arraySigma2)
for f in FluxField:
    IMaxPointSigma3.append(maxCurrent(f, arraySigma3[1:-1]))
print(arraySigma3)

# Normalizes Flux_Field to be for flux quanta

FluxField = FluxField * 50

plt.figure(300)
plt.plot(FluxField, IMaxPointSigma1, 'black' , label="Sigma = " + str(Sigma1) + ": " + str(round(percentageSigma1 * 100, 2)) + "%")
plt.plot(FluxField, IMaxPointSigma2, '#ADD8E6', label="Sigma = " + str(Sigma2) + ": " + str(round(percentageSigma2 * 100, 2)) + "%")
plt.plot(FluxField, IMaxPointSigma3, 'green', label="Sigma = " + str(Sigma3) + ": " + str(round(percentageSigma3 * 100, 2)) + "%")
plt.legend(loc="upper right", frameon=True)
plt.xlabel('Magnetic Field')
plt.ylabel('Critical Current with ' + str(junctionNumber) + ' junctions')  
plt.grid()             
plt.show()