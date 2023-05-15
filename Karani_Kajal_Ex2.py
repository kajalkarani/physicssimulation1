# -*- coding: utf-8 -*-
"""
@author: Kajal

"""



import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#variables

Cd = 1
p0 = 1200 #g/m^3
A = 0.5 #m^3
g = 9.81 #m/s^2
dT = 1 #s  
vi = 0 #m/s
ti = 0 #s
m = 100000 #g
h = 7640 #m

#function definitions

def AnalyticalPrediction(v0,y0,t0): 
    velocity1, position1, time1 = [],[],[] 
    v=0.0
    t=0.0
    y=y0
    k = (Cd*p0*A)/2 #drag coefficient
    while (y>0):
        v = -1*math.sqrt(m*g/k)*math.tanh((math.sqrt(k*g/m)*t))
        y = y0-((m/(2*k)))*math.log(math.cosh(math.sqrt(k*g/m)*t)**2)
        #append elements in lists
        velocity1.append(abs(v))
        position1.append(y)
        time1.append(t)
        t +=dT
    return velocity1, position1, time1


def Euler(v,y,t):
    k = (Cd*p0*A)/2 #drag coefficient 
    velocity, position, time = [],[],[]
    while (y>0):
        vnext = v-dT*(g+(k/m)*abs(v)*(v))
        ynext = y+dT*v
        tnext = t+dT
        #following variables
        v = vnext
        y = ynext
        t = tnext
        #append elements in lists
        velocity.append(abs(v))
        position.append(y)
        time.append(t)
    return velocity, position, time


def modifiedEuler(v,y,t):
    velocity2, position2, time2 = [],[],[]
    k = (Cd*p0*A)/2 #drag coefficient
    while (y>0):
        vmidpoint = v-(dT/2)*(g-(v**2)*(k/m))
        vnext = v-dT*(g-(vmidpoint**2)*(k/m))
        ynext = y+vmidpoint*dT
        tnext = t+dT
        #append elements in lists
        velocity2.append(abs(v))
        position2.append(y)
        time2.append(t)
        #following variables
        v = vnext
        y = ynext
        t = tnext
    return velocity2, position2, time2


def modifiedEulerDensity(v,y,t):
    velocity3, position3, time3 = [],[],[]
    while (y>0):
        p = p0*math.exp(-y/h) #density
        k = (Cd*p*A)/2 #drag coefficient
        vmidpoint = v-(dT/2)*(g-(v**2)*(k/m))
        vnext = v-dT*(g-(vmidpoint**2)*(k/m))
        ynext = y+vmidpoint*dT
        tnext = t+dT
        #append elements in lists
        velocity3.append(abs(v))
        position3.append(y)
        time3.append(t)
        #following variables
        v = vnext
        y = ynext
        t = tnext
    return velocity3, position3, time3


def Variations(v,y,t,A,Cd): #equation from lecture 4 (range kutta)
    velocity4, position4, time4 = [],[],[]
    while (y>0):
        p = p0*math.exp(-y/h) #density
        k = (Cd*p*A)/2 #drag coefficient
        vmidpoint1 = v-(dT/2)*(g-(v**2)*(k/m))
        vmidpoint2 = v-(dT/2)*(g-(vmidpoint1**2)*(k/m))
        vnext = v-(dT/2)*(g-(vmidpoint2**2)*(k/m))
        ynext = y+(dT/6)*(v-(dT/2)*(g-(v**2)*(k/m))+2*vmidpoint1+2*vmidpoint2+vnext)
        tnext = t+dT 
        #append elements in lists
        velocity4.append(abs(v))
        position4.append(y)
        time4.append(t)
        #following variables
        v = vnext
        y = ynext
        t = tnext
    return velocity4,position4,time4


#list definitions
    
velocity, position, time = Euler(vi,1000,ti)
velocity1, position1, time1 = AnalyticalPrediction(vi,1000,ti)
velocity2, position2, time2 = modifiedEuler(vi,1000,ti)
velocity3, position3, time3 = modifiedEulerDensity(vi,39045,ti)
velocity4, position4, time4 = Variations(vi,39045,ti,A,Cd)

#main program 

MyInput = '0'

#menu

while MyInput != 'q':
    MyInput = input('Enter a choice, "a", "b", "c", "d", "e": ')
    print('You entered the choice ', MyInput)

#part a
    
    if MyInput == 'a':
        
        yi = 1000 #initial height 
        #velocity against time
        plt.plot(time1,velocity1) 
        plt.title('predicted velocity as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        plt.show()
        #displacement against time
        plt.plot(time1,position1)
        plt.title('predicted position as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('position (m)')
        plt.show()

#part b
        
    elif MyInput == 'b': 
        
        yi = 1000 #initial height 
        #velocity against time
        plt.plot(time,velocity, 'b-') #Euler
        plt.plot(time1,velocity1,'r-') #prediction
        plt.title('velocity as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        #legend 
        blue_patch = mpatches.Patch(color='blue', label='Eulers method')
        red_patch = mpatches.Patch(color='red', label='analytical prediction')
        plt.legend(handles=[red_patch, blue_patch])
        plt.show()
        #displacement against time
        plt.plot(time,position,'b-') #Euler
        plt.plot(time1,position1,'r-') #prediction
        plt.title('position as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('position(m)')
        #legend
        blue_patch = mpatches.Patch(color='blue', label='Eulers method')
        red_patch = mpatches.Patch(color='red', label='analytical prediction')
        plt.legend(handles=[red_patch, blue_patch])
        plt.show()
        
    elif MyInput == 'c': 
        
        yi = 1000 #initial height 
        #velocity against time
        plt.plot(time,velocity,'b-') #Euler
        plt.plot(time1,velocity1,'r-') #prediction
        plt.plot(time2,velocity2,'g-') #modifiedEuler
        plt.title('velocity as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        #legend
        blue_patch = mpatches.Patch(color='blue', label='Eulers method')
        red_patch = mpatches.Patch(color='red', label='Analytical prediction')
        green_patch = mpatches.Patch(color='green', label='Modified Eulers method')
        plt.legend(handles=[blue_patch, red_patch, green_patch])
        plt.show()
        #displacement against time
        plt.plot(time,position,'b-') #Euler
        plt.plot(time1,position1,'r-') #prediction
        plt.plot(time2,position2,'g-') #modifiedEuler
        plt.title('position as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('position (m)')
        #legend
        blue_patch = mpatches.Patch(color='blue', label='Eulers method')
        red_patch = mpatches.Patch(color='red', label='Analytical prediction')
        green_patch = mpatches.Patch(color='green', label='Modified Eulers method')
        plt.legend(handles=[blue_patch, red_patch, green_patch])
        plt.show()
        
        #for the report, factor k/m is increased and the step size dT is varied after this part for analysis
        #code for creating legends taken from https://matplotlib.org/3.1.1/tutorials/intermediate/legend_guide.html

#part c

    elif MyInput == 'd': 
        
        yi = 39045 #initial height 
        #velocity against time
        plt.plot(time3,velocity3) 
        plt.title('velocity as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        plt.show()
        #displacement against time
        plt.plot(time3,position3) 
        plt.title('position as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('position (m)')
        plt.show()
        
#part d

    elif MyInput == 'e':
        
        yi = 39045 #initial height 
        V,Y,T = Variations(vi,yi,ti,0.3,1)
        V1,Y1,T1 = Variations(vi,yi,ti,0.4,1.1)
        V2,Y2,T2 = Variations(vi,yi,ti,0.5,1.2)
        V3,Y3,T3 = Variations(vi,yi,ti,0.6,1.3)
        V4,Y4,T4 = Variations(vi,yi,ti,0.7,1.4)
        #velocity against time for a range of values
        plt.plot(T,V) 
        plt.plot(T1,V1)
        plt.plot(T2,V2)
        plt.plot(T3,V3)
        plt.plot(T4,V4)
        plt.title('velocity as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('velocity (m/s)')
        plt.show()
        #position against time for a range of values
        plt.plot(T,Y) 
        plt.plot(T1,Y1)
        plt.plot(T2,Y2)
        plt.plot(T3,Y3)
        plt.plot(T4,Y4)
        plt.title('position as a function of time')
        plt.xlabel('time (s)')
        plt.ylabel('position (m)')
        plt.show()
        












    
