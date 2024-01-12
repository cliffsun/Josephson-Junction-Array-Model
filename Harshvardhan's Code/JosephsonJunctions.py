#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 15:28:46 2023

@author: harshvardhanmantry
"""

import numpy as np
import matplotlib.pyplot as plt

# n = int(input('Enter number of current flowing regions: '))
# randomList=np.random.random(2*(n-1)) # Number of points betweeen 0 and 1 will always be even
# randomList.sort()
# print(randomList)
# testList=[0.3, 0.5]
#test1List=[0.03205079, 0.0405224,  0.05092537, 0.06884022, 0.08417173, 0.1011885,
 # 0.11006377, 0.1144234,  0.13990064, 0.19354071, 0.20555725, 0.21358621,
# 0.21632374, 0.22668973, 0.25120565, 0.29798007, 0.3280096, 0.35133005,
 # 0.45859968, 0.52191115, 0.52657748, 0.54832122, 0.57996021, 0.59682881,
#  0.60895792, 0.66822352, 0.74195716, 0.75274705, 0.75416494, 0.77706997,
#  0.79543251, 0.83889362, 0.9228953,  0.95009882, 0.95421896, 0.95938356,
#  0.95951322, 0.98242441]
# test2=[0.01, 0.03, 0.13, 0.14, 0.142, 0.16, 0.36, 0.38, 0.68, 0.681, 0.688, 0.7, 0.75, 0.77, 0.83, 0.84, 0.841, 0.89]
test2 = [0.5, 0.5]

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



MagField=np.linspace(0.001, 0.4, 5000)


def maxCurrent(B, array):
    """ Spits out the maximum current by varying the gauge invariant phase of the left end (free parameter) gamma"""
    Y=np.linspace(0, 2*np.pi, 50)
    dummyArray=[]
    for gamma in Y:
        dummyArray.append(current(B, array, gamma))
    return max(dummyArray)

#Creating the max supercurrent array

IMax=[]
for B in MagField:
    IMax.append(maxCurrent(B, test2))



plt.figure(dpi=300)
plt.plot(MagField, IMax, 'g-')
plt.xlabel('Magnetic field (T)')          
plt.ylabel('I_max')  
plt.grid()             
plt.show()





# m= float(input("Enter the slope of the linear function: "))
# c=float(input("Enter the Y - intercept of the linear function: "))
# def linearFn(x, m, c):
#     return m*x+c


# def currentLinear(B, y, m, c): # This will give us a a general expression for current depending on gamma
#     array=np.linspace(0, 1, 101)
#     n=len(array)
#     dx=array[1]
#     i=0
#     f=0
#     gamma=y
#     while(i<n-1):
#        f+=linearFn(array[i], m, c)*np.sin(50*np.pi*B*dx)*np.sin(gamma+2*50*np.pi*B*array[i]+50*np.pi*B*dx)/B
#        i+=1
#     return f
   

# MagField=np.linspace(0.001, 5, 1000)


# def maxCurrent(B, m, c):
#     Y=np.linspace(0, 2*np.pi, 50)
#     dummyArray=[]
#     for gamma in Y:
#         dummyArray.append(currentLinear(B, gamma, m, c))
#     return max(dummyArray)

# IMax=[]
# for B in MagField:
#     IMax.append(maxCurrent(B, m, c))

# plt.plot(MagField, IMax, 'g-')
# plt.xlabel('B')          
# plt.ylabel('I_max')  
# plt.grid()             
# plt.show()

# plt.plot(np.linspace(0, 1, 11), linearFn(np.linspace(0, 1, 11), m, c), "k-")
# plt.xlabel('x')          
# plt.ylabel('J(x)')  
# plt.grid()             
# plt.show()

    
                




