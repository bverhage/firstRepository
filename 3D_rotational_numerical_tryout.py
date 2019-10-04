# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 16:58:36 2019

@author: billy
"""

# -*- coding: utf-8 -*-
"""
This is a numeric differential equation solver
"""
print('------ begin of code ------')
## preperation
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

## world variables

## list of the active forces

##functions
## functions for numeric diff.eq solving
def EF(Time,w,dt):
    ### Euler forward numerical integration method
    ### w_(n+1)=w_n+dt*f(tn,wn)
    ans=w[:,-1:]+dt*f(Time[-1],w[:,-1:]);
    return(ans)

def TZ(Time,w,dt):
    ### trapezodial numerical integration method
    ### w_(n+1)=w_n+dt*(f(t_n,w_n)+f(tn+dt,w*_(n+1)))/2
    ### with w*_(n+1)=w_n+dt*f(tn,wn)
    
    ans= w[:,-1:]+dt/2*(f(Time[-1],w[:,-1:])+f(Time[-1]+dt,w[:,-1:]+dt*f(Time[-1],w[:,-1:])));
    
    return(ans)
    
## diff.eq  
    
    
    ## the differential eq      
def f(t,w):
    #creating a place matrix
    o=np.delete(w,3,axis=0)
    o=np.delete(o,3,axis=0)
    o=np.delete(o,3,axis=0)
    #creating a speed matrix
    o_dot=np.delete(w,0,axis=0)
    o_dot=np.delete(o_dot,0,axis=0)
    o_dot=np.delete(o_dot,0,axis=0)
    
    

    
    #the homogeneus part of de Differential eqation
    #A=np.array([[0,1,0,0,0,0],
    #            [0,0,0,0,0,0],
    #            [0,0,0,1,0,0],
    #            [0,0,0,0,0,0],
    #            [0,0,0,0,0,1],
    #            [0,0,0,0,0,0]])
    
    A=np.block([
                [np.zeros((3,3)),np.identity(3)],
                [np.zeros((3,3)),np.zeros((3,3))]
                ])
    
    
    #the non homogenus part
    ans=A.dot(w)+Internal(o_dot)+External()
    return(ans)      




## Forces
def Internal(o_dot):
    #interial matrix 
    I=np.array([[2,0,0],
                [0,2,0],
                [0,0,1]])
    
    #linear part of the rotional acceleration
    Internal=-1*np.linalg.inv(I).dot(np.transpose(np.cross(np.transpose(I.dot(o_dot)),np.transpose(o_dot))))

    ans=np.block([
                [np.zeros((3,1))],
                [Internal]
                ])
    return(ans)
    
def External():
    #air resitance
    I=np.array([[2,0,0],
                [0,2,0],
                [0,0,1]])
    
    M=np.array([[0],
                [0],
                [0]])
    
    ans=np.block([
                [np.zeros((3,1))],
                [np.linalg.inv(I).dot(M)]
                ])
    return(ans)


## variables
t0=0
Time=[t0]
tE=10
dt=0.5
w0=np.array([[0],
             [0],
             [0],
             [5],
             [10],
             [0]])
w=w0

## The excecution
while(Time[-1]<tE):
    #solution
    wn=EF(Time,w,dt)
    w=np.append(w,wn,axis=1)

    #time step
    Time.append(Time[-1]+dt)
    

    
    
## plots
    #creating a place matrix
    o=np.delete(w,3,axis=0)
    o=np.delete(o,3,axis=0)
    o=np.delete(o,3,axis=0)
    #creating a speed matrix
    o_dot=np.delete(w,0,axis=0)
    o_dot=np.delete(o_dot,0,axis=0)
    o_dot=np.delete(o_dot,0,axis=0)

#for i in range(len(Time)):
#    Forces.append(FSUM(w[:,i:]))


fig1 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax.scatter(o[0,:], o[1,:], o[2,:], c='r', marker='o')
fig2 = plt.figure()
ax = fig2.add_subplot(111, projection='3d')
ax.scatter(Time, o[1,:], o[2,:], c='r', marker='o')
fig2 = plt.figure()
plt.plot(Time,o_dot[0,:],'-',color='red')
plt.plot(Time,o_dot[1,:],'-',color='blue')
plt.plot(Time,o_dot[2,:],'-',color='yello')
#plt.subplot(312)
#plt.plot(Time,v[2,:],'-')
#plt.subplot(313)
#plt.plot(Time[0:-1],Forces,'-')
#plt.show
print('------ end of code ------')    


